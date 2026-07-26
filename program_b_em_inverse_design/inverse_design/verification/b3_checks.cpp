#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include "b3_checks.hpp"
#include "program_b/inverse_design/band_objective.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

// ---------------------------------------------------------------------------
// Gate B3-PHYSICAL-ADMISSIBILITY.
//
// SCOPE, stated up front and repeated in the report: a passing result is a
// passive, causal parametric constitutive model that approximates the target
// response over a stated frequency band. It is NOT a fabricated medium or an
// experimental device. There is no unit-cell geometry, no hardware control model
// and no Maxwell propagation anywhere in B3.
// ---------------------------------------------------------------------------

namespace program_b {
namespace checks {

using namespace program_b::inverse_design;
using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace {

using CD = Complex<double>;

std::string sci(double v) {
    std::ostringstream o;
    o.setf(std::ios::scientific);
    o.precision(4);
    o << v;
    return o.str();
}

struct Setup {
    LaboratoryMetric<double> g_lab;
    ObserverFrame<double> frame;
    ObjectiveWeights<double> weights;
    BandSpec band;
};

Setup default_setup() {
    Setup s;
    s.g_lab = LaboratoryMetric<double>::minkowski();
    double V[4] = {1.0, 0.0, 0.0, 0.0};
    s.frame = build_observer_frame(
        s.g_lab, make_laboratory_observer(s.g_lab, V, "lab_rest"));
    // The band SPANS the resonances at 3.0, 3.6 and 4.4. That is not cosmetic:
    // probing far below every resonance makes the model NON-IDENTIFIABLE -- see
    // b3.identifiability.requires_probing_the_resonance -- because F and w0 can
    // then be traded against each other at nearly constant response.
    s.band.omega_min = 2.0;
    s.band.omega_max = 5.0;
    s.band.samples = 31;
    return s;
}

// A damped, passive, single-pole reference medium well away from resonance.
DispersiveModelSpec one_pole_spec(ParameterMode mode) {
    DispersiveModelSpec spec;
    spec.eps_poles = 1;
    spec.mu_poles = 1;
    spec.drude = false;
    spec.magnetoelectric = false;
    spec.mode = mode;
    return spec;
}

DispersiveModelSpec two_pole_spec(ParameterMode mode) {
    DispersiveModelSpec spec = one_pole_spec(mode);
    spec.eps_poles = 2;
    return spec;
}

DispersiveModelSpec magnetoelectric_spec(ParameterMode mode) {
    DispersiveModelSpec spec = one_pole_spec(mode);
    spec.magnetoelectric = true;
    return spec;
}

void fill_physical(double p[kPhysicalSlots]) {
    for (int i = 0; i < kPhysicalSlots; ++i) p[i] = 0.0;
}

// eps_inf, one eps pole, mu_inf, one mu pole. Resonances placed ABOVE the band
// so the medium is low-loss inside it.
void one_pole_values(double p[kPhysicalSlots]) {
    fill_physical(p);
    p[kSlotEpsInf] = 2.10;
    p[kSlotEpsPole0F] = 1.60; p[kSlotEpsPole0W0] = 3.00; p[kSlotEpsPole0Gamma] = 0.18;
    p[kSlotMuInf] = 1.15;
    p[kSlotMuPole0F] = 0.45; p[kSlotMuPole0W0] = 3.60; p[kSlotMuPole0Gamma] = 0.22;
}

void two_pole_values(double p[kPhysicalSlots]) {
    one_pole_values(p);
    p[kSlotEpsPole1F] = 0.90; p[kSlotEpsPole1W0] = 4.40; p[kSlotEpsPole1Gamma] = 0.30;
}

void magnetoelectric_values(double p[kPhysicalSlots]) {
    one_pole_values(p);
    p[kSlotZetaAmp] = 0.06; p[kSlotZetaW0] = 3.20; p[kSlotZetaGamma] = 0.25;
}

template <typename S>
PhysicalDesignVariables<S> make_design(const DispersiveModelSpec& spec,
                                       const double p[kPhysicalSlots]) {
    return design_from_physical<S>(spec, p);
}

double diff4x4(const Rank4<double>& a, const Rank4<double>& b) {
    double w = 0.0;
    for (int i = 0; i < 256; ++i) w = std::fmax(w, std::fabs(a.t[i] - b.t[i]));
    return w;
}

// Deterministic perturbation of a physical parameter set.
void perturb_physical(const double in[kPhysicalSlots], double out[kPhysicalSlots],
                      double amp, int seed) {
    unsigned int s = static_cast<unsigned int>(seed) * 2654435761u + 7u;
    for (int i = 0; i < kPhysicalSlots; ++i) {
        s = s * 1664525u + 1013904223u;
        const double u = static_cast<double>((s >> 8) & 0xFFFFu) / 65535.0;
        out[i] = in[i] * (1.0 + amp * (2.0 * u - 1.0));
    }
}

}  // namespace

void run_b3_checks(VerificationReport& report, const std::string& out_dir,
                  B3Metrics& metrics) {
    const Setup S = default_setup();
    const FrequencyConventions conv;
    metrics.fourier_convention = conv.fourier_kernel;
    metrics.model_family = "Lorentz + optional Drude, isotropic; mu parameterized "
                           "directly, mu^-1 derived; zeta = z(w) K, reciprocal";

    // =======================================================================
    // 1. Complex arithmetic.
    // =======================================================================
    {
        Stopwatch sw;
        const CD a(2.0, -3.0), b(-1.5, 0.75);
        const CD sum = a + b, prod = a * b, quo = a / b;

        // Hand values.
        const double pr = 2.0 * -1.5 - (-3.0) * 0.75;      // -3 + 2.25 = -0.75
        const double pi = 2.0 * 0.75 + (-3.0) * (-1.5);    // 1.5 + 4.5 = 6.0
        const double d = 1.5 * 1.5 + 0.75 * 0.75;          // 2.8125
        const double qr = (2.0 * -1.5 + (-3.0) * 0.75) / d;
        const double qi = ((-3.0) * (-1.5) - 2.0 * 0.75) / d;

        double worst = 0.0;
        worst = std::fmax(worst, std::fabs(sum.re - 0.5));
        worst = std::fmax(worst, std::fabs(sum.im + 2.25));
        worst = std::fmax(worst, std::fabs(prod.re - pr));
        worst = std::fmax(worst, std::fabs(prod.im - pi));
        worst = std::fmax(worst, std::fabs(quo.re - qr));
        worst = std::fmax(worst, std::fabs(quo.im - qi));

        // (a/b)*b == a, and a*conj(a) == |a|^2 real.
        const CD back = quo * b;
        worst = std::fmax(worst, std::fabs(back.re - a.re));
        worst = std::fmax(worst, std::fabs(back.im - a.im));
        const CD aa = a * conj(a);
        worst = std::fmax(worst, std::fabs(aa.re - norm_squared(a)));
        worst = std::fmax(worst, std::fabs(aa.im));

        report.check_residual(
            "b3.complex.arithmetic",
            "addition, multiplication and division against hand-computed values; "
            "(a/b)*b = a; a conj(a) is real and equals |a|^2",
            VerificationLevel::A, "n/a", "a = 2-3i, b = -1.5+0.75i", worst, 1e-15,
            sw.seconds(),
            "written from the definitions rather than std::complex, whose "
            "division and abs branch on magnitudes -- which a hypercomplex AD "
            "scalar cannot tolerate");
    }

    {
        Stopwatch sw;
        CMat3<double> A;
        A(0, 0) = CD(2.0, 0.3);  A(0, 1) = CD(0.4, -0.2); A(0, 2) = CD(-0.1, 0.5);
        A(1, 0) = CD(0.4, -0.2); A(1, 1) = CD(1.7, 0.25); A(1, 2) = CD(0.2, 0.1);
        A(2, 0) = CD(-0.1, 0.5); A(2, 1) = CD(0.2, 0.1);  A(2, 2) = CD(2.4, 0.15);

        // A is symmetric but NOT Hermitian, which is the whole point: the two
        // residuals must disagree.
        const double sym = transpose_residual(A);
        const double herm = hermiticity_residual(A);

        CD det;
        const CMat3<double> Ai = invert(A, det);
        const CMat3<double> I = A * Ai;
        double worst_inv = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                worst_inv = std::fmax(worst_inv,
                                      std::fabs(I(i, j).re - (i == j ? 1.0 : 0.0)));
                worst_inv = std::fmax(worst_inv, std::fabs(I(i, j).im));
            }

        report.check_predicate(
            "b3.complex.transpose_and_adjoint_are_distinct",
            "for a complex symmetric matrix the plain transpose residual vanishes "
            "while the Hermitian adjoint residual does not",
            VerificationLevel::A, sym < 1e-15 && herm > 1e-2,
            "|A - A^T| = 0 and |A - A^dagger| > 0",
            "transpose residual " + sci(sym) + ", adjoint residual " + sci(herm),
            sw.seconds(),
            "Lorentz reciprocity is the UNCONJUGATED transpose and passivity is "
            "the conjugate one. Using the adjoint for reciprocity is a physics "
            "error that only appears once the medium is lossy -- which is where "
            "B3 begins");

        report.check_residual(
            "b3.complex.matrix_inverse_roundtrip", "A A^-1 = I for a complex 3x3",
            VerificationLevel::A, "n/a", "adjugate form, branch-free", worst_inv,
            1e-14, sw.seconds());
    }

    {
        Stopwatch sw;
        // Real limit: with zero imaginary parts the complex path must reproduce
        // the real one bit for bit.
        DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        one_pole_values(p);
        p[kSlotEpsPole0Gamma] = 0.0;   // lossless
        p[kSlotMuPole0Gamma] = 0.0;
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        double worst = 0.0;
        for (double w : S.band.frequencies()) {
            const ResponseAtFrequency<double> r = evaluate_response(v, w);
            worst = std::fmax(worst, std::fabs(r.epsilon.im));
            worst = std::fmax(worst, std::fabs(r.mu.im));
            worst = std::fmax(worst, std::fabs(r.mu_inverse.im));
            // and the real part must equal the closed form
            const double den = 3.0 * 3.0 - w * w;
            worst = std::fmax(worst, std::fabs(r.epsilon.re - (2.10 + 1.60 / den)));
        }
        report.check_residual(
            "b3.complex.real_limit_recovery",
            "with every damping set to zero the response is exactly real and "
            "matches the undamped closed form",
            VerificationLevel::B, "n/a", "21 frequencies, gamma = 0", worst, 1e-14,
            sw.seconds(),
            "the imaginary parts are exactly zero, not small: a lossless model "
            "must not leak loss through the arithmetic");
    }

    // =======================================================================
    // 2. Conventions.
    // =======================================================================
    {
        Stopwatch sw;
        // The sign that pins the convention. Under e^{-iwt} a Lorentz term with
        // positive F, gamma has Im eps > 0 (passive). Under e^{+iwt} the same
        // algebraic form would give Im eps < 0, i.e. gain -- so the model forms
        // the milestone specifies are consistent with e^{-iwt} and with nothing
        // else.
        const double F = 1.6, w0 = 3.0, gam = 0.18;
        double smallest_im = 1e300;
        for (double w : S.band.frequencies()) {
            const CD e = lorentz_term<double>(F, w0, gam, w);
            smallest_im = std::fmin(smallest_im, e.im);
        }
        const CD d = drude_term<double>(1.2, 0.1, 1.0);

        std::ostringstream got;
        got << "kernel " << conv.fourier_kernel << "; min Im(Lorentz) over band = "
            << sci(smallest_im) << "; Im(Drude) at w=1 = " << sci(d.im);

        report.check_predicate(
            "b3.convention.sign_of_loss_is_pinned",
            "under the recorded e^{-i omega t} kernel both specified model forms "
            "have STRICTLY POSITIVE imaginary part for positive parameters, which "
            "is passive; under the opposite kernel both would be gain",
            VerificationLevel::B, smallest_im > 0.0 && d.im > 0.0,
            "Im > 0 for the Lorentz term across the band and for the Drude term",
            got.str(), sw.seconds(),
            "this is the check that fixes the Fourier convention rather than "
            "assuming it. Every passivity sign downstream depends on it");

        report.check_predicate(
            "b3.convention.recorded",
            "the Fourier kernel, time derivative, constitutive relation, "
            "dissipated power, Poynting vector, passive sign, causal half-plane "
            "and reality condition are all recorded",
            VerificationLevel::A,
            std::string(conv.fourier_kernel) == "exp(-i*omega*t)" &&
                std::string(conv.passive_sign).size() > 0 &&
                std::string(conv.poynting).size() > 0 &&
                std::string(conv.causal_half_plane).size() > 0,
            "all conventions recorded",
            std::string(conv.fourier_kernel) + " | " + conv.dissipated_power,
            sw.seconds());
    }

    // =======================================================================
    // 3. Dispersive models.
    // =======================================================================
    {
        Stopwatch sw;
        const DispersiveModelSpec spec = two_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        two_pole_values(p);
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        double worst = 0.0;
        for (double w : {0.7, 1.0, 1.5}) {
            const ResponseAtFrequency<double> r = evaluate_response(v, w);
            // Hand-computed two-pole Lorentz sum.
            auto term = [&](double F, double w0, double g) {
                const double dr = w0 * w0 - w * w, di = -g * w;
                const double m = dr * dr + di * di;
                return CD(F * dr / m, -F * di / m);
            };
            const CD e0 = term(1.60, 3.00, 0.18);
            const CD e1 = term(0.90, 4.40, 0.30);
            worst = std::fmax(worst, std::fabs(r.epsilon.re - (2.10 + e0.re + e1.re)));
            worst = std::fmax(worst, std::fabs(r.epsilon.im - (e0.im + e1.im)));
        }
        report.check_residual(
            "b3.dispersion.lorentz_closed_form",
            "the two-pole Lorentz permittivity matches a hand-written closed form",
            VerificationLevel::B, "n/a", "3 frequencies", worst, 1e-14, sw.seconds());

        // mu is parameterized, mu^-1 derived. mu^-1 * mu = 1 exactly, and the
        // RESONANCE lives in mu.
        double worst_inv = 0.0, mu_peak = 0.0, mu_inv_peak = 0.0;
        for (double w = 0.5; w <= 5.0; w += 0.01) {
            const ResponseAtFrequency<double> r = evaluate_response(v, w);
            const CD one = r.mu * r.mu_inverse;
            worst_inv = std::fmax(worst_inv, std::fabs(one.re - 1.0));
            worst_inv = std::fmax(worst_inv, std::fabs(one.im));
            if (std::fabs(r.mu.im) > mu_peak) mu_peak = std::fabs(r.mu.im);
            if (std::fabs(r.mu_inverse.im) > mu_inv_peak)
                mu_inv_peak = std::fabs(r.mu_inverse.im);
        }
        report.check_residual(
            "b3.dispersion.mu_inverse_is_derived_from_mu",
            "[mu^-1](w) mu(w) = 1 at every sampled frequency",
            VerificationLevel::A, "n/a", "w from 0.5 to 5", worst_inv, 1e-13,
            sw.seconds(),
            "mu is the quantity PARAMETERIZED and mu^-1 is derived by inversion. "
            "Parameterizing mu^-1 as a sum of Lorentz terms would place poles of "
            "mu^-1 -- zeros of mu -- at the resonances, which is a different and "
            "generally non-passive medium");

        // Reality condition eps(-w) = conj(eps(w)).
        double worst_reality = 0.0;
        for (double w : {0.4, 1.0, 2.5, 3.0}) {
            const ResponseAtFrequency<double> a = evaluate_response(v, w);
            const ResponseAtFrequency<double> b = evaluate_response(v, -w);
            worst_reality = std::fmax(worst_reality, std::fabs(a.epsilon.re - b.epsilon.re));
            worst_reality = std::fmax(worst_reality, std::fabs(a.epsilon.im + b.epsilon.im));
        }
        report.check_residual(
            "b3.dispersion.reality_condition",
            "eps(-w) = conj(eps(w)), the condition for a real time-domain response",
            VerificationLevel::A, "n/a", "4 frequencies", worst_reality, 1e-14,
            sw.seconds());
    }

    {
        Stopwatch sw;
        const DispersiveModelSpec spec = magnetoelectric_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        magnetoelectric_values(p);
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        double worst = 0.0, biggest = 0.0;
        for (double w : S.band.frequencies()) {
            const ComplexBlocks<double> b = blocks_at_frequency(v, w);
            worst = std::fmax(worst, cmat3_max_diff(b.zeta_adjoint, transpose(b.zeta)));
            biggest = std::fmax(biggest, cmat3_max_abs(b.zeta));
        }
        report.check_residual(
            "b3.dispersion.magnetoelectric_is_reciprocal_by_construction",
            "zeta_adjoint = zeta^T at every frequency, exactly",
            VerificationLevel::A, "n/a",
            "21 frequencies; max|zeta| = " + sci(biggest), worst, 0.0, sw.seconds(),
            "EXACT zero: zeta = z(w) K with K a fixed real antisymmetric shape, so "
            "the transpose relation is structural. Independent arbitrary complex "
            "entries would break reciprocity by construction and are not used");
    }

    // =======================================================================
    // 5. Causality -- structural.
    // =======================================================================
    {
        Stopwatch sw;
        bool all_lower = true;
        std::ostringstream poles;
        poles.setf(std::ios::fixed);
        poles.precision(4);
        struct PC { const char* name; double w0, g; };
        const PC cases[4] = {{"eps pole 0", 3.00, 0.18},
                             {"eps pole 1", 4.40, 0.30},
                             {"mu pole 0", 3.60, 0.22},
                             {"overdamped", 0.20, 3.00}};
        for (const PC& c : cases) {
            const PolePair pp = lorentz_poles(c.w0, c.g);
            all_lower = all_lower && pp.both_in_lower_half_plane;
            poles << c.name << ": " << pp.real_part_a << (pp.imag_part_a < 0 ? "" : "+")
                  << pp.imag_part_a << "i, " << pp.real_part_b
                  << (pp.imag_part_b < 0 ? "" : "+") << pp.imag_part_b << "i"
                  << (pp.overdamped ? " (overdamped)" : "") << "; ";
        }
        metrics.pole_locations = poles.str();

        report.check_predicate(
            "b3.causality.poles_in_lower_half_plane",
            "every pole of every Lorentz term sits in the LOWER half omega-plane, "
            "including the overdamped case where both roots are purely imaginary",
            VerificationLevel::B, all_lower,
            "Im(pole) < 0 for all 8 roots", poles.str(), sw.seconds(),
            "under e^{-i omega t} causality requires analyticity in the UPPER half "
            "plane. This is enforced by the parameterization -- gamma is strictly "
            "positive by construction -- and not argued from sampled data");

        // The transform that makes it structural.
        double smallest = 1e300;
        for (double u : {-1e6, -100.0, -1.0, 0.0, 1.0, 100.0})
            smallest = std::fmin(smallest, positive_transform<double>(u, 1e-6));
        report.check_predicate(
            "b3.causality.positive_transform_is_strictly_positive",
            "the damping transform is strictly positive for every input, including "
            "extreme negative values, so no optimizer step can produce an "
            "acausal model",
            VerificationLevel::A, smallest > 0.0, "positive for all inputs",
            "smallest value over u in {-1e6 .. 1e2} = " + sci(smallest),
            sw.seconds(),
            "this is HARD feasibility for causality: not a penalty that can be "
            "traded away, but a parameterization the infeasible region is not in");

        // Round trip through the transform.
        double worst_rt = 0.0;
        for (double value : {1e-4, 0.05, 0.18, 1.0, 25.0}) {
            const double u = positive_transform_inverse(value, 1e-6);
            worst_rt = std::fmax(worst_rt,
                                 std::fabs(positive_transform<double>(u, 1e-6) - value));
        }
        report.check_residual(
            "b3.causality.transform_roundtrip",
            "positive_transform(positive_transform_inverse(x)) = x, so a known "
            "physical model can be expressed exactly as a design vector",
            VerificationLevel::A, "n/a", "5 values", worst_rt, 1e-12, sw.seconds());
    }

    {
        Stopwatch sw;
        // Kramers-Kronig, DIAGNOSTIC ONLY.
        const DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        one_pole_values(p);
        p[kSlotMuPole0F] = 0.0;
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        // Reported as a CONVERGENCE TREND rather than a single tolerance. A
        // fixed threshold here would be a statement about the quadrature, not
        // about causality; what actually distinguishes a truncation error from a
        // physics violation is that widening the integration band drives it
        // down. So that is what is measured.
        struct KK { double lo, hi; int n; double residual; };
        KK trials[3] = {{1e-2, 1e2, 20001, 0.0},
                        {1e-3, 1e3, 40001, 0.0},
                        {1e-4, 1e4, 80001, 0.0}};
        for (KK& k : trials) {
            double w = 0.0;
            for (double probe : {0.8, 1.2, 1.6})
                w = std::fmax(w, kramers_kronig_residual(v, probe, k.lo, k.hi,
                                                         k.n, 2.10));
            k.residual = w;
        }
        const bool decreasing = trials[1].residual < trials[0].residual &&
                                trials[2].residual < trials[1].residual;

        std::ostringstream kkgot;
        kkgot.setf(std::ios::scientific);
        kkgot.precision(3);
        for (const KK& k : trials)
            kkgot << "[" << k.lo << "," << k.hi << "] -> " << k.residual << "  ";

        report.check_predicate(
            "b3.causality.kramers_kronig_diagnostic",
            "the Kramers-Kronig residual DECREASES monotonically as the "
            "integration band widens, and is small on the widest band -- the "
            "signature of quadrature truncation rather than a causality "
            "violation",
            VerificationLevel::D, decreasing && trials[2].residual < 5e-3,
            "monotone decrease with band width, and below 5e-3 on the widest",
            kkgot.str(), sw.seconds(),
            "A DIAGNOSTIC, NOT the causality proof. A fixed tolerance here would "
            "be measuring the principal-value quadrature, not the physics -- "
            "which is why the convergence trend is asserted instead. Causality "
            "is established structurally by pole placement, and that case is "
            "b3.causality.poles_in_lower_half_plane");
    }

    // =======================================================================
    // 4. Passivity.
    // =======================================================================
    {
        Stopwatch sw;
        const DispersiveModelSpec spec = magnetoelectric_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        magnetoelectric_values(p);
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);
        const std::vector<double> ws = S.band.frequencies();
        const PassivityReport pr = assess_passivity(v, ws);

        report.check_residual(
            "b3.passivity.loss_matrix_is_hermitian",
            "L = (C - C^dagger)/(2i) is Hermitian at every band frequency",
            VerificationLevel::A, "n/a", "21 frequencies",
            pr.worst_hermiticity_residual, 1e-13, sw.seconds(),
            "Hermitian by construction, so this catches an error in the "
            "assembly rather than in the physics");

        report.check_residual(
            "b3.passivity.reciprocity_in_EH_representation",
            "xi_EH = zeta_EH^T, the reciprocity relation DERIVED for this "
            "project's (E,B) sign convention",
            VerificationLevel::B, "n/a", "21 frequencies",
            pr.worst_reciprocity_residual, 1e-13, sw.seconds(),
            "many references write xi = -zeta^T; that sign traces to the opposite "
            "convention for zeta in D = eps E - zeta B. This project's sign was "
            "derived at B1C and cross-checked against Minkowski's moving-medium "
            "relations, so the relation is verified here rather than imported");
    }

    {
        Stopwatch sw;
        // Vacuum: lossless, so the smallest loss eigenvalue is exactly zero.
        DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        spec.eps_poles = 0;
        spec.mu_poles = 0;
        double p[kPhysicalSlots];
        fill_physical(p);
        p[kSlotEpsInf] = 1.0;
        p[kSlotMuInf] = 1.0;
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);
        const PassivityReport pr = assess_passivity(v, S.band.frequencies());

        report.check_predicate(
            "b3.passivity.vacuum_is_lossless",
            "the nondispersive vacuum limit has loss eigenvalues exactly zero: "
            "passive, and on the boundary rather than inside it",
            VerificationLevel::B,
            std::fabs(pr.min_loss_eigenvalue) < 1e-14 && pr.passive_across_band,
            "min eigenvalue of L is 0",
            "min loss eigenvalue " + sci(pr.min_loss_eigenvalue), sw.seconds());
    }

    {
        Stopwatch sw;
        // A damped passive Lorentz oscillator.
        const DispersiveModelSpec spec = two_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        two_pole_values(p);
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);
        const PassivityReport pr = assess_passivity(v, S.band.frequencies());
        metrics.min_loss_eigenvalue = pr.min_loss_eigenvalue;
        metrics.max_passivity_violation = pr.max_violation;
        metrics.worst_passivity_frequency = pr.worst_frequency;

        report.check_predicate(
            "b3.passivity.damped_lorentz_is_passive",
            "a damped two-pole Lorentz medium is passive across the whole band: "
            "the smallest eigenvalue of the loss matrix is non-negative",
            VerificationLevel::B, pr.passive_across_band && pr.min_loss_eigenvalue >= 0.0,
            "min eigenvalue of L >= 0 at every frequency",
            "min loss eigenvalue " + sci(pr.min_loss_eigenvalue) + " at w = " +
                sci(pr.worst_frequency),
            sw.seconds());
    }

    {
        Stopwatch sw;
        // An intentionally ACTIVE oscillator, which must fail.
        DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        one_pole_values(p);
        p[kSlotEpsPole0Gamma] = -0.18;   // negative damping: gain
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);
        const PassivityReport pr = assess_passivity(v, S.band.frequencies());
        const PolePair pp = lorentz_poles(3.00, -0.18);

        report.check_predicate(
            "b3.passivity.active_oscillator_is_rejected",
            "a Lorentz term with NEGATIVE damping is detected as active -- the "
            "loss matrix acquires a negative eigenvalue -- and its poles move "
            "into the upper half plane, so it is acausal as well",
            VerificationLevel::B,
            !pr.passive_across_band && pr.max_violation > 1e-6 &&
                !pp.both_in_lower_half_plane,
            "passivity violated and poles in the upper half plane",
            "min loss eigenvalue " + sci(pr.min_loss_eigenvalue) + ", violation " +
                sci(pr.max_violation) + ", pole Im = " + sci(pp.imag_part_a),
            sw.seconds(),
            "the test must be able to FAIL a medium, or a passing verdict on a "
            "passive one means nothing");
    }

    {
        Stopwatch sw;
        // The case the milestone warns about: Im eps > 0 and Im mu > 0
        // everywhere, yet the medium is NOT passive, because the magnetoelectric
        // coupling drives an off-diagonal block of L.
        DispersiveModelSpec spec = magnetoelectric_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        one_pole_values(p);
        p[kSlotEpsPole0Gamma] = 1e-4;   // almost lossless
        p[kSlotMuPole0Gamma] = 1e-4;
        p[kSlotZetaAmp] = 6.0;          // strong, and genuinely complex
        p[kSlotZetaW0] = 1.10;
        p[kSlotZetaGamma] = 0.40;
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        const std::vector<double> ws = S.band.frequencies();
        double smallest_im_eps = 1e300, smallest_im_mu = 1e300;
        for (double w : ws) {
            const ResponseAtFrequency<double> r = evaluate_response(v, w);
            smallest_im_eps = std::fmin(smallest_im_eps, r.epsilon.im);
            smallest_im_mu = std::fmin(smallest_im_mu, r.mu.im);
        }
        const PassivityReport pr = assess_passivity(v, ws);

        std::ostringstream got;
        got << "min Im eps = " << sci(smallest_im_eps) << ", min Im mu = "
            << sci(smallest_im_mu) << ", yet min eigenvalue of L = "
            << sci(pr.min_loss_eigenvalue) << " at w = " << sci(pr.worst_frequency);

        report.check_predicate(
            "b3.passivity.diagonal_imaginary_parts_are_not_sufficient",
            "a medium with STRICTLY POSITIVE Im eps and Im mu at every frequency "
            "is nevertheless NOT passive when the magnetoelectric coupling is "
            "strong: the full loss matrix has a negative eigenvalue",
            VerificationLevel::B,
            smallest_im_eps > 0.0 && smallest_im_mu > 0.0 &&
                !pr.passive_across_band,
            "positive diagonal imaginary parts, negative loss eigenvalue",
            got.str(), sw.seconds(),
            "this is why 'all imaginary parts are positive' must NOT be equated "
            "with passivity for a coupled bianisotropic medium. The off-diagonal "
            "block L_12 = (xi_EH - zeta_EH^dagger)/(2i) carries the coupling, and "
            "with near-lossless diagonal blocks it alone sets the spectrum");
    }

    // =======================================================================
    // 6. The band objective.
    // =======================================================================
    {
        Stopwatch sw;
        // (A11) must be real-linear, which is what licenses evaluating Z(omega)
        // as A11(Re blocks) + i A11(Im blocks) rather than writing a second
        // complex reconstruction.
        ConstitutiveDesignVariables<double> a, b;
        a.family = DesignFamily::UnrestrictedReciprocal;
        b.family = DesignFamily::UnrestrictedReciprocal;
        for (int k = 0; k < 21; ++k) {
            a.q[k] = 0.3 + 0.11 * k;
            b.q[k] = -0.7 + 0.05 * k * k;
        }
        const double al = 1.7, be = -0.9;
        ConstitutiveDesignVariables<double> c = a;
        for (int k = 0; k < 21; ++k) c.q[k] = al * a.q[k] + be * b.q[k];

        const Rank4<double> Za = reconstruct_design_tensor(a, S.g_lab, S.frame);
        const Rank4<double> Zb = reconstruct_design_tensor(b, S.g_lab, S.frame);
        const Rank4<double> Zc = reconstruct_design_tensor(c, S.g_lab, S.frame);
        double worst = 0.0;
        for (int i = 0; i < 256; ++i)
            worst = std::fmax(worst, std::fabs(Zc.t[i] - (al * Za.t[i] + be * Zb.t[i])));

        report.check_residual(
            "b3.band.a11_is_real_linear",
            "(A11) is homogeneous linear in the blocks with real coefficients",
            VerificationLevel::A, "n/a", "alpha = 1.7, beta = -0.9", worst, 1e-12,
            sw.seconds(),
            "this is what licenses Z(w) = A11(Re blocks) + i A11(Im blocks), so "
            "the verified real reconstruction is reused rather than a second "
            "complex implementation being written");
    }

    {
        Stopwatch sw;
        const DispersiveModelSpec spec = two_pole_spec(ParameterMode::BoundedTransform);
        double p[kPhysicalSlots];
        two_pole_values(p);
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);
        const BandTarget<double> t =
            band_target_from_model(v, S.g_lab, S.frame, S.band, "two_pole");

        BandConstraints<double> c;
        const BandEvaluation<double> e =
            evaluate_band_design(v, t, c, S.frame, S.weights);

        report.check_predicate(
            "b3.band.zero_at_exact_target_and_reported_separately",
            "J_band vanishes at the exact target, and the centre, mean and worst "
            "frequency mismatches are reported as three distinct numbers",
            VerificationLevel::A,
            std::fabs(real_part(e.J_band)) < 1e-28 && e.centre_mismatch < 1e-28 &&
                e.mean_mismatch < 1e-28 && e.worst_mismatch < 1e-28,
            "all mismatch measures zero at the target",
            "centre " + sci(e.centre_mismatch) + ", mean " + sci(e.mean_mismatch) +
                ", worst " + sci(e.worst_mismatch) + " at w = " +
                sci(e.worst_mismatch_frequency),
            sw.seconds(),
            "the contributions are never collapsed into one number: a model that "
            "is excellent at band centre and poor at the edges is a different "
            "object from one that is uniformly mediocre");
    }

    // =======================================================================
    // 8. Gradients.
    // =======================================================================
    {
        Stopwatch sw;
        // Analytic derivatives of a single Lorentz pole:
        //   eps = eps_inf + F/D,  D = (w0^2 - w^2) - i gamma w
        //   d eps / dF     = 1/D
        //   d eps / dw0    = -2 F w0 / D^2
        //   d eps / dgamma = i F w / D^2
        const double F = 1.6, w0 = 3.0, gam = 0.18, w = 1.3;
        const CD D(w0 * w0 - w * w, -gam * w);
        const CD invD = reciprocal(D);
        const CD invD2 = invD * invD;
        const CD dF = invD;
        const CD dw0 = invD2 * (-2.0 * F * w0);
        const CD dgam = CD(0.0, 1.0) * invD2 * (F * w);

        auto eval = [&](double f, double ww0, double g) {
            return lorentz_term<double>(f, ww0, g, w);
        };

        double worst = 0.0;
        for (double h : {1e-4, 1e-5, 1e-6}) {
            const CD gF = (eval(F + h, w0, gam) - eval(F - h, w0, gam)) *
                          (1.0 / (2 * h));
            const CD g0 = (eval(F, w0 + h, gam) - eval(F, w0 - h, gam)) *
                          (1.0 / (2 * h));
            const CD gg = (eval(F, w0, gam + h) - eval(F, w0, gam - h)) *
                          (1.0 / (2 * h));
            worst = std::fmax(worst, std::fabs(gF.re - dF.re));
            worst = std::fmax(worst, std::fabs(gF.im - dF.im));
            worst = std::fmax(worst, std::fabs(g0.re - dw0.re));
            worst = std::fmax(worst, std::fabs(g0.im - dw0.im));
            worst = std::fmax(worst, std::fabs(gg.re - dgam.re));
            worst = std::fmax(worst, std::fabs(gg.im - dgam.im));
        }
        metrics.worst_gradient_analytic = worst;
        report.check_residual(
            "b3.gradient.analytic_one_lorentz_pole",
            "d eps/dF = 1/D, d eps/dw0 = -2 F w0/D^2, d eps/dgamma = i F w/D^2 "
            "match centred differences",
            VerificationLevel::B, "n/a", "3 step sizes, w = 1.3", worst, 1e-6,
            sw.seconds(),
            "hand-derived from the closed form, for both real and imaginary parts");
    }

    {
        Stopwatch sw;
        // COMPLEX-STEP differentiation, instantiating the model at
        // Scalar = Complex<double>: the outer imaginary unit is the perturbation.
        // No subtractive cancellation, so it is exact for arbitrarily small h.
        const DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        one_pole_values(p);
        const PhysicalDesignVariables<double> vd = make_design<double>(spec, p);

        const double w = 1.3, h = 1e-30;
        double worst = 0.0;
        const std::vector<int> slots = spec.active_slots();
        for (std::size_t k = 0; k < slots.size(); ++k) {
            PhysicalDesignVariables<CD> vc;
            vc.spec = spec;
            for (std::size_t j = 0; j < slots.size(); ++j)
                vc.q[j] = CD(real_part(vd.q[j]), 0.0);
            vc.q[k] = CD(real_part(vd.q[k]), h);

            const ResponseAtFrequency<CD> rc = evaluate_response(vc, w);
            const double d_re = rc.epsilon.re.im / h;
            const double d_im = rc.epsilon.im.im / h;

            // Against centred differences of the plain-double model.
            PhysicalDesignVariables<double> hi = vd, lo = vd;
            const double fd_h = 1e-6;
            hi.q[k] += fd_h;
            lo.q[k] -= fd_h;
            const ResponseAtFrequency<double> rh = evaluate_response(hi, w);
            const ResponseAtFrequency<double> rl = evaluate_response(lo, w);
            worst = std::fmax(worst,
                              std::fabs(d_re - (rh.epsilon.re - rl.epsilon.re) / (2 * fd_h)));
            worst = std::fmax(worst,
                              std::fabs(d_im - (rh.epsilon.im - rl.epsilon.im) / (2 * fd_h)));
        }
        metrics.worst_gradient_complex_step = worst;
        report.check_residual(
            "b3.gradient.complex_step",
            "complex-step derivatives, from instantiating the model at "
            "Scalar = Complex<double> with a step of 1e-30, agree with centred "
            "differences",
            VerificationLevel::D, "n/a",
            std::to_string(slots.size()) + " parameters, h = 1e-30 against "
            "centred h = 1e-6",
            worst, 1e-6, sw.seconds(),
            "the agreement is limited by the CENTRED DIFFERENCE, not by the "
            "complex step: the latter has no subtractive cancellation and is "
            "exact to machine precision. The scalar seam is what makes this "
            "possible without a second implementation");
    }

    // Batched HyPAD gradient of the band objective.
    {
        Stopwatch sw;
        const DispersiveModelSpec spec = two_pole_spec(ParameterMode::BoundedTransform);
        double p[kPhysicalSlots], pp[kPhysicalSlots];
        two_pole_values(p);
        perturb_physical(p, pp, 0.15, 3);
        const PhysicalDesignVariables<double> truth = make_design<double>(spec, p);
        const PhysicalDesignVariables<double> x = make_design<double>(spec, pp);
        const BandTarget<double> t =
            band_target_from_model(truth, S.g_lab, S.frame, S.band, "two_pole");
        BandConstraints<double> c;

        const int n = x.size();
        double fd[kPhysicalSlots] = {};
        for (int k = 0; k < n; ++k) {
            PhysicalDesignVariables<double> hi = x, lo = x;
            const double h = 1e-6;
            hi.q[k] += h;
            lo.q[k] -= h;
            fd[k] = (real_part(evaluate_band_design(hi, t, c, S.frame, S.weights).J_total) -
                     real_part(evaluate_band_design(lo, t, c, S.frame, S.weights).J_total)) /
                    (2 * h);
        }
        metrics.worst_gradient_fd = 0.0;

#ifdef WARP_HYPAD_HAVE_OTILIB
        {
            namespace sens = warp_hypad::sensitivity;
            using sens::OTIScalar;
            metrics.oti_available = true;

            const int batch = sens::kOtiDirections;
            double g[kPhysicalSlots] = {};
            int batches = 0;
            std::ostringstream active_report;

            for (int start = 0; start < n; start += batch) {
                PhysicalDesignVariables<OTIScalar> vo;
                vo.spec = spec;
                for (int j = 0; j < n; ++j) vo.q[j] = OTIScalar(real_part(x.q[j]));
                int count = 0;
                active_report << "[";
                for (int b = 0; b < batch && start + b < n; ++b) {
                    vo.q[start + b] = sens::seed_parameter(real_part(x.q[start + b]), b);
                    active_report << (b ? "," : "") << (start + b);
                    ++count;
                }
                active_report << "]";

                Mat4<OTIScalar> lab_o;
                for (int a2 = 0; a2 < 4; ++a2)
                    for (int b2 = 0; b2 < 4; ++b2)
                        lab_o(a2, b2) = OTIScalar(S.g_lab.value(a2, b2));
                const LaboratoryMetric<OTIScalar> lab(lab_o);
                OTIScalar Vo[4];
                for (int a2 = 0; a2 < 4; ++a2) Vo[a2] = OTIScalar(S.frame.V[a2]);
                const ObserverFrame<OTIScalar> fo = build_observer_frame(
                    lab, make_laboratory_observer(lab, Vo, "lab_rest"));

                BandTarget<OTIScalar> to;
                to.g_lab = lab;
                to.origin = t.origin;
                to.omega = t.omega;
                for (std::size_t q2 = 0; q2 < t.Z.size(); ++q2) {
                    ComplexRank4<OTIScalar> z;
                    for (int i = 0; i < 256; ++i) {
                        z.re.t[i] = OTIScalar(t.Z[q2].re.t[i]);
                        z.im.t[i] = OTIScalar(t.Z[q2].im.t[i]);
                    }
                    to.Z.push_back(z);
                }
                BandConstraints<OTIScalar> co;
                co.enforce_passivity_penalty = c.enforce_passivity_penalty;
                co.passivity_weight = c.passivity_weight;
                ObjectiveWeights<OTIScalar> wo;
                wo.norm = S.weights.norm;
                wo.delta_Z = S.weights.delta_Z;

                const BandEvaluation<OTIScalar> e =
                    evaluate_band_design(vo, to, co, fo, wo);
                for (int b = 0; b < count; ++b)
                    g[start + b] = sens::extract_gradient(e.J_total, b);
                ++batches;
            }
            metrics.oti_batches = batches;

            double worst = 0.0;
            for (int k = 0; k < n; ++k) worst = std::fmax(worst, std::fabs(g[k] - fd[k]));
            metrics.worst_gradient_oti = worst;

            report.check_predicate(
                "b3.gradient.oti_is_batched",
                "the band objective's gradient is assembled from batches of at "
                "most kOtiDirections = 3 active directions",
                VerificationLevel::A, batches == (n + batch - 1) / batch,
                std::to_string((n + batch - 1) / batch) + " batches",
                std::to_string(batches) + " batches over " + std::to_string(n) +
                    " parameters, active " + active_report.str(),
                sw.seconds(),
                "the whole chain q -> C(w) -> Z(w) -> J_band is differentiated "
                "through the same scalar-generic objective, not a second "
                "implementation");

            report.check_residual(
                "b3.gradient.oti_matches_finite_difference",
                "the batch-assembled HyPAD gradient of J_band matches centred "
                "differences",
                VerificationLevel::D, "n/a",
                std::to_string(n) + " parameters", worst, 1e-6, sw.seconds());
        }
#else
        {
            report.check_skipped(
                "b3.gradient.oti_is_batched",
                "batched HyPAD gradient of the band objective",
                "the OTIlib backend is disabled; the analytic, complex-step and "
                "finite-difference gradient cases still ran, and the optimizer "
                "runs on finite differences",
                sw.seconds());
            report.check_skipped(
                "b3.gradient.oti_matches_finite_difference",
                "HyPAD band gradient against centred differences",
                "the OTIlib backend is disabled", sw.seconds());
        }
#endif
    }

    // =======================================================================
    // 9. The verification problems.
    // =======================================================================
    std::ofstream hist(out_dir + "/b3_objective_histories.csv", std::ios::binary);
    if (hist) hist << "problem,iteration,objective,gradient_norm,step,backtracks\n";
    auto record = [&](const std::string& name, const BandOptimizerResult& r) {
        if (hist) {
            hist.setf(std::ios::scientific);
            hist.precision(10);
            for (const IterationRecord& it : r.history)
                hist << name << ',' << it.iteration << ',' << it.objective << ','
                     << it.gradient_norm << ',' << it.step_length << ','
                     << it.backtracks << '\n';
        }
    };

    BandOptimizerSettings opt;
    opt.max_iterations = 1500;

    // Problem 5: exact recovery of a manufactured one-pole target.
    {
        Stopwatch sw;
        const DispersiveModelSpec spec = one_pole_spec(ParameterMode::BoundedTransform);
        double p[kPhysicalSlots], pp[kPhysicalSlots];
        one_pole_values(p);
        perturb_physical(p, pp, 0.20, 11);
        const PhysicalDesignVariables<double> truth = make_design<double>(spec, p);
        const PhysicalDesignVariables<double> start = make_design<double>(spec, pp);
        const BandTarget<double> t =
            band_target_from_model(truth, S.g_lab, S.frame, S.band, "one_pole");
        BandConstraints<double> c;

        const BandOptimizerResult r =
            optimize_band_design<double>(start, t, c, S.frame, S.weights, opt);
        record("one_pole_recovery", r);
        metrics.one_pole_recovery_error = r.best_objective;
        metrics.termination_reason = to_string(r.reason);

        double worst_param = 0.0;
        double pr_[kPhysicalSlots], pt_[kPhysicalSlots];
        physical_values(r.best, pr_);
        physical_values(truth, pt_);
        for (int slot : spec.active_slots())
            worst_param = std::fmax(worst_param, std::fabs(pr_[slot] - pt_[slot]));
        metrics.worst_parameter_recovery_error = worst_param;

        report.check_predicate(
            "b3.recovery.one_pole",
            "a manufactured single-pole target is recovered from a 20 percent "
            "perturbation, in physical parameters and not merely in objective",
            VerificationLevel::B,
            r.best_objective < 1e-16 && worst_param < 1e-4,
            "J_band below 1e-16 and every physical parameter within 1e-4",
            "J: " + sci(r.initial_objective) + " -> " + sci(r.best_objective) +
                ", max|dp| = " + sci(worst_param) + ", " + to_string(r.reason),
            sw.seconds());
    }

    // Problem 6: two-pole recovery.
    {
        Stopwatch sw;
        const DispersiveModelSpec spec = two_pole_spec(ParameterMode::BoundedTransform);
        double p[kPhysicalSlots], pp[kPhysicalSlots];
        two_pole_values(p);
        perturb_physical(p, pp, 0.10, 23);
        const PhysicalDesignVariables<double> truth = make_design<double>(spec, p);
        const PhysicalDesignVariables<double> start = make_design<double>(spec, pp);
        const BandTarget<double> t =
            band_target_from_model(truth, S.g_lab, S.frame, S.band, "two_pole");
        BandConstraints<double> c;

        const BandOptimizerResult r =
            optimize_band_design<double>(start, t, c, S.frame, S.weights, opt);
        record("two_pole_recovery", r);
        metrics.two_pole_recovery_error = r.best_objective;

        const BandEvaluation<double> e =
            evaluate_band_design(r.best, t, c, S.frame, S.weights);
        metrics.centre_mismatch = e.centre_mismatch;
        metrics.mean_mismatch = e.mean_mismatch;
        metrics.max_mismatch = e.worst_mismatch;
        metrics.worst_mismatch_frequency = e.worst_mismatch_frequency;

        report.check_predicate(
            "b3.recovery.two_pole",
            "a two-pole target is recovered to a band mismatch below 1e-12",
            VerificationLevel::B, r.best_objective < 1e-14,
            "J_band below 1e-14",
            "J: " + sci(r.initial_objective) + " -> " + sci(r.best_objective) +
                "; centre " + sci(e.centre_mismatch) + ", mean " +
                sci(e.mean_mismatch) + ", worst " + sci(e.worst_mismatch) +
                " at w = " + sci(e.worst_mismatch_frequency) + ", " +
                to_string(r.reason),
            sw.seconds());
    }

    // Problem 7: a one-pole model against a two-pole target.
    {
        Stopwatch sw;
        double p2[kPhysicalSlots];
        two_pole_values(p2);
        const DispersiveModelSpec spec2 = two_pole_spec(ParameterMode::BoundedTransform);
        const PhysicalDesignVariables<double> truth = make_design<double>(spec2, p2);
        const BandTarget<double> t =
            band_target_from_model(truth, S.g_lab, S.frame, S.band, "two_pole");

        const DispersiveModelSpec spec1 = one_pole_spec(ParameterMode::BoundedTransform);
        double p1[kPhysicalSlots];
        one_pole_values(p1);
        const PhysicalDesignVariables<double> start = make_design<double>(spec1, p1);
        BandConstraints<double> c;

        const BandOptimizerResult r =
            optimize_band_design<double>(start, t, c, S.frame, S.weights, opt);
        record("restricted_one_pole", r);
        metrics.restricted_residual = r.best_objective;

        report.check_predicate(
            "b3.restricted_one_pole_approximation",
            "a one-pole model fitted to a two-pole target leaves a STRICTLY "
            "NONZERO residual, and still improves substantially on its start",
            VerificationLevel::B,
            r.best_objective > 1e-9 && r.best_objective < r.initial_objective,
            "nonzero residual, and a genuine improvement",
            "J: " + sci(r.initial_objective) + " -> " + sci(r.best_objective) +
                ", " + to_string(r.reason),
            sw.seconds(),
            "the residual is a property of the model family, not an optimizer "
            "failure -- the same statement B2 made about restricted design "
            "families, now in the frequency domain");
    }

    // Problem 8: hard-bound versus soft-penalty.
    {
        Stopwatch sw;
        const DispersiveModelSpec hard_spec = two_pole_spec(ParameterMode::BoundedTransform);
        const DispersiveModelSpec soft_spec = two_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        two_pole_values(p);
        const PhysicalDesignVariables<double> truth =
            make_design<double>(hard_spec, p);
        const BandTarget<double> t =
            band_target_from_model(truth, S.g_lab, S.frame, S.band, "two_pole");

        // A bound the target genuinely violates: eps_inf must not exceed 1.6,
        // but the target needs 2.10.
        BandConstraints<double> c;
        c.bound_parameters = true;
        for (int i = 0; i < kPhysicalSlots; ++i) {
            c.slot_lo[i] = -1e30;
            c.slot_hi[i] = 1e30;
        }
        c.slot_lo[kSlotEpsInf] = 1.0;
        c.slot_hi[kSlotEpsInf] = 1.6;
        c.bound_weight = 1e3;

        double pstart[kPhysicalSlots];
        two_pole_values(pstart);
        pstart[kSlotEpsInf] = 1.3;

        BandOptimizerSettings hard = opt;
        hard.mode = FeasibilityMode::HardFeasible;
        const BandOptimizerResult rh = optimize_band_design<double>(
            make_design<double>(hard_spec, pstart), t, c, S.frame, S.weights, hard);
        record("hard_feasible", rh);

        BandOptimizerSettings soft = opt;
        soft.mode = FeasibilityMode::SoftPenalty;
        const BandOptimizerResult rs = optimize_band_design<double>(
            make_design<double>(soft_spec, pstart), t, c, S.frame, S.weights, soft);
        record("soft_penalized", rs);

        metrics.hard_bound_violation = rh.max_bound_violation;
        metrics.soft_bound_violation = rs.max_bound_violation;

        std::ostringstream got;
        got << to_string(rh.mode) << ": violation " << sci(rh.max_bound_violation)
            << ", J = " << sci(rh.best_objective) << "; " << to_string(rs.mode)
            << ": violation " << sci(rs.max_bound_violation) << ", J = "
            << sci(rs.best_objective);

        report.check_predicate(
            "b3.hard_feasible_versus_soft_penalty",
            "the hard-feasible run satisfies the bound EXACTLY -- zero violation "
            "-- and is labelled feasible, while the soft-penalised run is "
            "reported under a different name whatever the size of its violation",
            VerificationLevel::B,
            rh.max_bound_violation <= 0.0 && rh.hard_feasible && !rs.hard_feasible,
            "hard run exactly feasible; soft run never labelled feasible",
            got.str(), sw.seconds(),
            "B2 admitted a small violation because a smooth hinge trades "
            "feasibility for differentiability. Here feasibility comes from a "
            "bounded parameter transform plus projection, so it is exact, and a "
            "soft result is NEVER called feasible merely because its violation "
            "is small");
    }

    // Problem 9: bandwidth tradeoff.
    {
        Stopwatch sw;
        double p2[kPhysicalSlots];
        two_pole_values(p2);
        const DispersiveModelSpec spec2 = two_pole_spec(ParameterMode::BoundedTransform);
        const PhysicalDesignVariables<double> truth = make_design<double>(spec2, p2);
        const DispersiveModelSpec spec1 = one_pole_spec(ParameterMode::BoundedTransform);
        double p1[kPhysicalSlots];
        one_pole_values(p1);

        std::ofstream bw(out_dir + "/b3_bandwidth_tradeoff.csv", std::ios::binary);
        if (bw) bw << "half_width,J_band,centre_mismatch,mean_mismatch,worst_mismatch\n";

        double narrow = 0.0, wide = 0.0;
        bool monotone = true;
        double previous = -1.0;
        for (double half : {0.15, 0.5, 1.0, 1.6}) {
            BandSpec band;
            band.omega_min = 3.4 - half;   // centred between the eps resonances
            band.omega_max = 3.4 + half;
            band.samples = 25;
            const BandTarget<double> t =
                band_target_from_model(truth, S.g_lab, S.frame, band, "two_pole");
            BandConstraints<double> c;
            const BandOptimizerResult r = optimize_band_design<double>(
                make_design<double>(spec1, p1), t, c, S.frame, S.weights, opt);
            const BandEvaluation<double> e =
                evaluate_band_design(r.best, t, c, S.frame, S.weights);
            if (bw) {
                bw.setf(std::ios::scientific);
                bw.precision(10);
                bw << half << ',' << r.best_objective << ',' << e.centre_mismatch
                   << ',' << e.mean_mismatch << ',' << e.worst_mismatch << '\n';
            }
            if (half == 0.15) narrow = r.best_objective;
            if (half == 1.6) wide = r.best_objective;
            if (previous >= 0.0 && r.best_objective < previous) monotone = false;
            previous = r.best_objective;
        }

        report.check_predicate(
            "b3.bandwidth_tradeoff",
            "a fixed one-pole model fits a two-pole target progressively worse as "
            "the band widens, and the degradation is monotone",
            VerificationLevel::B, monotone && wide > narrow && narrow >= 0.0,
            "monotone increase in residual with bandwidth",
            "J at half-width 0.15 = " + sci(narrow) + ", at 1.6 = " + sci(wide),
            sw.seconds(),
            "the tradeoff a fixed model family imposes, measured rather than "
            "asserted");
    }

    // Problem 10: reciprocal magnetoelectric model.
    {
        Stopwatch sw;
        const DispersiveModelSpec spec = magnetoelectric_spec(ParameterMode::BoundedTransform);
        double p[kPhysicalSlots], pp[kPhysicalSlots];
        magnetoelectric_values(p);
        perturb_physical(p, pp, 0.12, 31);
        const PhysicalDesignVariables<double> truth = make_design<double>(spec, p);
        const BandTarget<double> t =
            band_target_from_model(truth, S.g_lab, S.frame, S.band, "magnetoelectric");
        BandConstraints<double> c;
        const BandOptimizerResult r = optimize_band_design<double>(
            make_design<double>(spec, pp), t, c, S.frame, S.weights, opt);
        record("magnetoelectric", r);

        double worst_recip = 0.0;
        for (double w : S.band.frequencies()) {
            const ComplexBlocks<double> b = blocks_at_frequency(r.best, w);
            worst_recip = std::fmax(worst_recip,
                                    cmat3_max_diff(b.zeta_adjoint, transpose(b.zeta)));
        }
        const PassivityReport pr = assess_passivity(r.best, S.band.frequencies());

        report.check_predicate(
            "b3.recovery.reciprocal_magnetoelectric",
            "a dispersive magnetoelectric target is recovered while reciprocity "
            "holds exactly at every frequency throughout",
            VerificationLevel::B,
            r.best_objective < 1e-14 && worst_recip == 0.0,
            "recovered, with zeta_adjoint = zeta^T exactly",
            "J -> " + sci(r.best_objective) + ", reciprocity residual " +
                sci(worst_recip) + ", min loss eigenvalue " +
                sci(pr.min_loss_eigenvalue),
            sw.seconds());
    }

    // Problem 11: coordinate-change consistency at each frequency.
    {
        Stopwatch sw;
        Mat4<double> L = Mat4<double>::identity();
        L(0, 0) = 1.15; L(0, 1) = 0.22;
        L(1, 0) = 0.18; L(1, 1) = 1.08;
        L(2, 1) = 0.11; L(2, 2) = 1.04; L(2, 3) = 0.07;
        L(3, 0) = 0.04; L(3, 2) = 0.13; L(3, 3) = 0.97;
        const warp_hypad::InverseResult<double> Li = warp_hypad::inverse(L);

        Mat4<double> gp;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                double acc = 0.0;
                for (int c2 = 0; c2 < 4; ++c2)
                    for (int d = 0; d < 4; ++d)
                        acc += Li.inverse(c2, a) * Li.inverse(d, b) * S.g_lab.value(c2, d);
                gp(a, b) = acc;
            }
        const LaboratoryMetric<double> g_prime(gp);
        double Vp[4];
        for (int a = 0; a < 4; ++a) {
            double acc = 0.0;
            for (int b = 0; b < 4; ++b) acc += L(a, b) * S.frame.V[b];
            Vp[a] = acc;
        }
        // The frame is PUSHED FORWARD, not rebuilt. Rebuilding it by
        // Gram-Schmidt in the new coordinates would pick different spatial
        // directions, and the anisotropic magnetoelectric shape K singles those
        // out -- see the companion case below, which measures exactly that.
        const ObserverFrame<double> frame_p =
            push_frame_forward(S.frame, L, g_prime);
        const ObserverFrame<double> frame_rebuilt = build_observer_frame(
            g_prime, make_laboratory_observer(g_prime, Vp, "lab_rest_primed"));

        const DispersiveModelSpec spec = magnetoelectric_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        magnetoelectric_values(p);
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        double worst = 0.0;
        for (double w : S.band.frequencies()) {
            const ComplexRank4<double> Z =
                reconstruct_Z_at_frequency(v, S.g_lab, S.frame, w);
            const ComplexRank4<double> Zp =
                reconstruct_Z_at_frequency(v, g_prime, frame_p, w);
            // Push the unprimed tensor forward and compare, part by part.
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int c2 = 0; c2 < 4; ++c2)
                        for (int d = 0; d < 4; ++d) {
                            double ar = 0.0, ai = 0.0;
                            for (int e = 0; e < 4; ++e)
                                for (int f = 0; f < 4; ++f)
                                    for (int g = 0; g < 4; ++g)
                                        for (int h = 0; h < 4; ++h) {
                                            const double lam = L(a, e) * L(b, f) *
                                                               L(c2, g) * L(d, h);
                                            ar += lam * Z.re(e, f, g, h);
                                            ai += lam * Z.im(e, f, g, h);
                                        }
                            worst = std::fmax(worst, std::fabs(Zp.re(a, b, c2, d) - ar));
                            worst = std::fmax(worst, std::fabs(Zp.im(a, b, c2, d) - ai));
                        }
        }
        metrics.worst_coordinate_change = worst;

        // The companion measurement: the same coefficients on an INDEPENDENTLY
        // rebuilt frame give a different tensor, because the magnetoelectric
        // shape is orientation dependent.
        double rebuilt_gap = 0.0;
        {
            const double w = S.band.centre();
            const ComplexRank4<double> Za =
                reconstruct_Z_at_frequency(v, g_prime, frame_p, w);
            const ComplexRank4<double> Zb =
                reconstruct_Z_at_frequency(v, g_prime, frame_rebuilt, w);
            for (int i = 0; i < 256; ++i) {
                rebuilt_gap = std::fmax(rebuilt_gap, std::fabs(Za.re.t[i] - Zb.re.t[i]));
                rebuilt_gap = std::fmax(rebuilt_gap, std::fabs(Za.im.t[i] - Zb.im.t[i]));
            }
        }

        report.check_predicate(
            "b3.frame_orientation_matters_for_anisotropic_shapes",
            "the same coefficients on an independently rebuilt triad give a "
            "DIFFERENT tensor, because the magnetoelectric shape K singles out "
            "particular spatial directions",
            VerificationLevel::A, rebuilt_gap > 1e-3,
            "a measurable difference between pushed-forward and rebuilt frames",
            "max|Z(pushed) - Z(rebuilt)| = " + sci(rebuilt_gap), sw.seconds(),
            "this is physics, not a defect. The isotropic parts are triad "
            "independent -- sum_i e_i^a e_i^b = h^{ab} for any triad -- but an "
            "anisotropic shape is not, so transforming the setup means "
            "transforming the FRAME, which is what the covariance case above "
            "does");

        report.check_residual(
            "b3.coordinate_change_consistency_per_frequency",
            "the same physical model built in transformed coordinates gives the "
            "pushed-forward Z at EVERY band frequency, real and imaginary parts "
            "alike",
            VerificationLevel::A, "n/a", "21 frequencies, non-orthogonal L",
            worst, 1e-11, sw.seconds(),
            "the frequency dependence lives entirely in the blocks; the geometry "
            "carries none, so covariance holds frequency by frequency");
    }

    // Problem 12: the real limit reduces to B2.
    {
        Stopwatch sw;
        DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        spec.eps_poles = 0;
        spec.mu_poles = 0;
        double p[kPhysicalSlots];
        fill_physical(p);
        p[kSlotEpsInf] = 2.75;
        p[kSlotMuInf] = 1.0 / 0.4;    // so mu^-1 = 0.4
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        // The same medium as a B2 design.
        ConstitutiveDesignVariables<double> b2;
        b2.family = DesignFamily::Isotropic;
        b2.q[0] = 2.75;
        b2.q[1] = 0.4;
        const Rank4<double> Z_b2 = reconstruct_design_tensor(b2, S.g_lab, S.frame);

        double worst = 0.0, worst_im = 0.0;
        for (double w : S.band.frequencies()) {
            const ComplexRank4<double> Z =
                reconstruct_Z_at_frequency(v, S.g_lab, S.frame, w);
            worst = std::fmax(worst, diff4x4(Z.re, Z_b2));
            for (int i = 0; i < 256; ++i)
                worst_im = std::fmax(worst_im, std::fabs(Z.im.t[i]));
        }
        metrics.real_limit_residual = worst;

        report.check_predicate(
            "b3.real_limit_reduces_to_b2",
            "a nondispersive lossless model reproduces the B2 isotropic design's "
            "Z exactly at every frequency, with identically zero imaginary part",
            VerificationLevel::B, worst < 1e-13 && worst_im == 0.0,
            "B3 real limit equals the B2 tensor, imaginary part exactly zero",
            "max|Z_B3 - Z_B2| = " + sci(worst) + ", max|Im Z| = " + sci(worst_im),
            sw.seconds(),
            "B3 must contain B2, or the frequency-domain layer has changed the "
            "physics rather than extended it");
    }

    // Vacuum nondispersive limit, and a lossless Lorentz away from resonance.
    {
        Stopwatch sw;
        DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        spec.eps_poles = 0;
        spec.mu_poles = 0;
        double p[kPhysicalSlots];
        fill_physical(p);
        p[kSlotEpsInf] = 1.0;
        p[kSlotMuInf] = 1.0;
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);

        const Rank4<double> Zvac =
            covariant_constitutive(EffectiveMetric<double>(S.g_lab.value), S.g_lab).Z;
        double worst = 0.0;
        for (double w : S.band.frequencies()) {
            const ComplexRank4<double> Z =
                reconstruct_Z_at_frequency(v, S.g_lab, S.frame, w);
            worst = std::fmax(worst, diff4x4(Z.re, Zvac));
        }
        report.check_residual(
            "b3.limit.vacuum_nondispersive",
            "eps = mu = 1 with no poles reproduces the vacuum Z at every "
            "frequency",
            VerificationLevel::B, "n/a", "21 frequencies", worst, 1e-14,
            sw.seconds());
    }

    {
        Stopwatch sw;
        DispersiveModelSpec spec = one_pole_spec(ParameterMode::Direct);
        double p[kPhysicalSlots];
        one_pole_values(p);
        p[kSlotEpsPole0Gamma] = 0.0;
        p[kSlotMuPole0Gamma] = 0.0;
        const PhysicalDesignVariables<double> v = make_design<double>(spec, p);
        BandSpec low;
        low.omega_min = 0.6; low.omega_max = 1.6; low.samples = 21;
        const PassivityReport pr = assess_passivity(v, low.frequencies());

        double worst_im = 0.0;
        for (double w : low.frequencies()) {
            const ResponseAtFrequency<double> r = evaluate_response(v, w);
            worst_im = std::fmax(worst_im, std::fabs(r.epsilon.im));
        }
        report.check_predicate(
            "b3.limit.lossless_lorentz_away_from_resonance",
            "an undamped Lorentz medium evaluated well below its resonance is "
            "exactly lossless and sits on the passivity boundary",
            VerificationLevel::B,
            worst_im == 0.0 && std::fabs(pr.min_loss_eigenvalue) < 1e-12,
            "Im eps exactly zero and min loss eigenvalue zero",
            "max|Im eps| = " + sci(worst_im) + ", min loss eigenvalue " +
                sci(pr.min_loss_eigenvalue),
            sw.seconds(),
            "band 0.6 to 1.6 against resonances at 3.0 and 3.6, so the model is "
            "far from its poles and the response is real");
    }

    // Determinism.
    {
        Stopwatch sw;
        const DispersiveModelSpec spec = one_pole_spec(ParameterMode::BoundedTransform);
        double p[kPhysicalSlots], pp[kPhysicalSlots];
        one_pole_values(p);
        perturb_physical(p, pp, 0.2, 11);
        const PhysicalDesignVariables<double> truth = make_design<double>(spec, p);
        const PhysicalDesignVariables<double> start = make_design<double>(spec, pp);
        const BandTarget<double> t =
            band_target_from_model(truth, S.g_lab, S.frame, S.band, "one_pole");
        BandConstraints<double> c;
        const BandOptimizerResult a =
            optimize_band_design<double>(start, t, c, S.frame, S.weights, opt);
        const BandOptimizerResult b =
            optimize_band_design<double>(start, t, c, S.frame, S.weights, opt);

        bool identical = a.iterations == b.iterations &&
                         a.best_objective == b.best_objective &&
                         a.history.size() == b.history.size();
        for (int k = 0; k < a.best.size() && identical; ++k)
            if (a.best.q[k] != b.best.q[k]) identical = false;

        report.check_predicate(
            "b3.deterministic_repeatability",
            "two band-optimizer runs from the same initialization are "
            "bit-identical",
            VerificationLevel::A, identical, "bit-identical runs",
            identical ? "identical over " + std::to_string(a.history.size()) +
                            " iterations"
                      : "runs diverged",
            sw.seconds());
    }

    // Passivity sweep export.
    {
        std::ofstream os(out_dir + "/b3_passivity_sweep.csv", std::ios::binary);
        if (os) {
            os << "omega,re_eps,im_eps,re_mu,im_mu,re_zeta,im_zeta,"
                  "min_loss_eigenvalue\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            const DispersiveModelSpec spec = magnetoelectric_spec(ParameterMode::Direct);
            double p[kPhysicalSlots];
            magnetoelectric_values(p);
            const PhysicalDesignVariables<double> v = make_design<double>(spec, p);
            for (double w = 0.2; w <= 5.0; w += 0.02) {
                const ResponseAtFrequency<double> r = evaluate_response(v, w);
                os << w << ',' << r.epsilon.re << ',' << r.epsilon.im << ','
                   << r.mu.re << ',' << r.mu.im << ',' << r.zeta.re << ','
                   << r.zeta.im << ',' << passivity_margin_at(v, w) << '\n';
            }
        }
    }
}

}  // namespace checks
}  // namespace program_b

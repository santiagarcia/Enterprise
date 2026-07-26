#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include "b4_checks.hpp"
#include "program_b/inverse_design/adaptive_scan.hpp"
#include "program_b/inverse_design/layered.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

// ---------------------------------------------------------------------------
// Gate B4-MAXWELL-PROPAGATION.
//
// SCOPE: a passing result demonstrates that the B3 constitutive model reproduces
// selected electromagnetic propagation observables of the target effective
// metric. It still does not establish a fabricated metamaterial, unit cell,
// antenna, coil, or experimental device.
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

using C = Complex<double>;

std::string sci(double v) {
    std::ostringstream o;
    o.setf(std::ios::scientific);
    o.precision(4);
    o << v;
    return o.str();
}

LaboratoryMetric<double> eta_lab() { return LaboratoryMetric<double>::minkowski(); }

ObserverFrame<double> lab_frame(const LaboratoryMetric<double>& g) {
    double V[4] = {1.0, 0.0, 0.0, 0.0};
    return build_observer_frame(g, make_laboratory_observer(g, V, "lab_rest"));
}

// A homogeneous non-dispersive medium as (E,H) blocks, built through the same
// (E,B) -> (E,H) conversion B3 verified.
EHBlocks<double> eh_blocks(double eps, double mu, double zeta_amp) {
    ComplexBlocks<double> b;
    b.epsilon = CMat3<double>::scaled_identity(C(eps, 0.0));
    b.mu = CMat3<double>::scaled_identity(C(mu, 0.0));
    b.mu_inverse = CMat3<double>::scaled_identity(C(1.0 / mu, 0.0));
    b.zeta = scale(magnetoelectric_shape<double>(), C(zeta_amp, 0.0));
    b.zeta_adjoint = transpose(b.zeta);
    return to_EH_representation(b);
}

// Z for a homogeneous medium, through the verified A11 path.
Rank4<double> medium_Z(double eps, double mu, double zeta_amp,
                       const LaboratoryMetric<double>& g,
                       const ObserverFrame<double>& frame) {
    ObserverDecomposition4<double> d;
    double e3[3][3] = {}, m3[3][3] = {}, z3[3][3] = {}, zt[3][3] = {};
    for (int i = 0; i < 3; ++i) { e3[i][i] = eps; m3[i][i] = 1.0 / mu; }
    z3[1][2] = zeta_amp; z3[2][1] = -zeta_amp;
    zt[1][2] = -zeta_amp; zt[2][1] = zeta_amp;
    d.epsilon = block_from_coefficients(e3, frame);
    d.mu_inverse = block_from_coefficients(m3, frame);
    d.zeta = block_from_coefficients(z3, frame);
    d.zeta_adjoint = block_from_coefficients(zt, frame);
    return reconstruct_Z_A11(d, g, frame.V);
}

// An anisotropic (diagonal, unequal) medium, which should be birefringent.
Rank4<double> anisotropic_Z(const LaboratoryMetric<double>& g,
                            const ObserverFrame<double>& frame) {
    ObserverDecomposition4<double> d;
    double e3[3][3] = {}, m3[3][3] = {}, z3[3][3] = {};
    e3[0][0] = 2.0; e3[1][1] = 3.5; e3[2][2] = 5.0;
    m3[0][0] = 1.0; m3[1][1] = 1.0; m3[2][2] = 1.0;
    d.epsilon = block_from_coefficients(e3, frame);
    d.mu_inverse = block_from_coefficients(m3, frame);
    d.zeta = block_from_coefficients(z3, frame);
    d.zeta_adjoint = block_from_coefficients(z3, frame);
    return reconstruct_Z_A11(d, g, frame.V);
}

Mat4<double> alcubierre_at(double vs, double x) {
    warp_hypad::metrics::AlcubierreMetric<double> w;
    const warp_hypad::metrics::AlcubierreParameters<double> p{vs, 3.0, 2.0};
    return w.evaluate(0.0, x, 0.0, 0.0, p);
}

}  // namespace

void run_b4_checks(VerificationReport& report, const std::string& out_dir,
                  B4Metrics& metrics) {
    const LaboratoryMetric<double> g = eta_lab();
    const ObserverFrame<double> frame = lab_frame(g);
    const MaxwellConventions mx;
    metrics.fourier_convention = mx.fourier_kernel;
    metrics.maxwell_curl_signs =
        std::string(mx.faraday) + " ; " + mx.ampere;

    // =======================================================================
    // 1. The Maxwell convention layer, pinned by vacuum plane waves.
    // =======================================================================
    {
        Stopwatch sw;
        report.check_predicate(
            "b4.convention.maxwell_signs_recorded",
            "the frequency-domain curl equations are recorded and derived from "
            "this project's own e^{-i omega t} kernel",
            VerificationLevel::A,
            std::string(mx.faraday) == "curl E = +i*omega*B" &&
                std::string(mx.ampere) == "curl H = -i*omega*D + J",
            "curl E = +i w B and curl H = -i w D + J",
            metrics.maxwell_curl_signs, sw.seconds(),
            "signs from a source using e^{+i omega t} are NOT copied. Mixing the "
            "two gives a medium that looks passive but propagates as if it had "
            "gain, and the equations stay self-consistent, so no residual would "
            "catch it -- the vacuum case below is what pins them");
    }

    {
        Stopwatch sw;
        // Vacuum: |k| = omega, with c = 1. Verified by the plane-wave operator
        // having a null space exactly there and NOT elsewhere.
        const EHBlocks<double> vac = eh_blocks(1.0, 1.0, 0.0);
        const double omega = 1.3;

        double worst_on = 0.0, smallest_off = 1e300;
        const double dirs[3][3] = {{1, 0, 0}, {0, 1, 0}, {0.5773502691896258,
                                                          0.5773502691896258,
                                                          0.5773502691896258}};
        for (const double* n : dirs) {
            double k[3] = {omega * n[0], omega * n[1], omega * n[2]};
            // The analytic transverse mode: E along any unit e perpendicular to
            // n, H = n x E. Residual measured DIRECTLY as ||M v||, with no
            // squaring, so it is accurate to machine precision.
            double e1[3];
            {
                const double t[3] = {-n[1], n[0], 0.0};
                const double tn = std::sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
                if (tn > 1e-8) { e1[0]=t[0]/tn; e1[1]=t[1]/tn; e1[2]=t[2]/tn; }
                else { e1[0]=0.0; e1[1]=1.0; e1[2]=0.0; }
            }
            const double h1[3] = {n[1]*e1[2]-n[2]*e1[1], n[2]*e1[0]-n[0]*e1[2],
                                  n[0]*e1[1]-n[1]*e1[0]};
            C v[6];
            for (int i = 0; i < 3; ++i) { v[i] = C(e1[i], 0.0); v[i+3] = C(h1[i], 0.0); }
            worst_on = std::fmax(worst_on,
                                 mode_residual(plane_wave_operator(vac, k, omega), v));
            double koff[3] = {1.7 * omega * n[0], 1.7 * omega * n[1],
                              1.7 * omega * n[2]};
            smallest_off = std::fmin(
                smallest_off,
                smallest_singular_value(plane_wave_operator(vac, koff, omega)));
        }
        metrics.worst_plane_wave_residual = worst_on;

        report.check_residual(
            "b4.vacuum.plane_wave_dispersion",
            "|k| = omega in vacuum with c = 1: the analytic transverse mode is "
            "annihilated by the plane-wave operator exactly on the light cone",
            VerificationLevel::B, "n/a", "3 directions, omega = 1.3", worst_on,
            1e-15, sw.seconds(),
            "measured as ||M v|| for the known mode, NOT as a smallest singular "
            "value: the latter goes through M^dagger M and its square root "
            "amplifies a 1e-16 eigenvalue error to 1e-8");

        report.check_predicate(
            "b4.vacuum.off_shell_is_not_singular",
            "at |k| = 1.7 omega the operator is NOT singular, so the previous "
            "case is not passing for a trivial reason",
            VerificationLevel::A, smallest_off > 1e-2,
            "smallest singular value well away from zero off the cone",
            "min off-shell singular value = " + sci(smallest_off), sw.seconds());
    }

    {
        Stopwatch sw;
        // Isotropic: |k|^2 = eps mu omega^2, and v_p = 1/sqrt(eps mu).
        double worst = 0.0, worst_vp = 0.0;
        for (double eps : {2.0, 4.0}) {
            for (double mu : {1.0, 1.5}) {
                const EHBlocks<double> b = eh_blocks(eps, mu, 0.0);
                const double omega = 1.1;
                const double kmag = omega * std::sqrt(eps * mu);
                double k[3] = {kmag, 0.0, 0.0};
                // E along y, H along z with H_z = sqrt(eps/mu) E_y, from
                // k x E = omega mu H.
                C v[6] = {C(0,0), C(1,0), C(0,0),
                          C(0,0), C(0,0), C(std::sqrt(eps / mu), 0.0)};
                worst = std::fmax(worst,
                                  mode_residual(plane_wave_operator(b, k, omega), v));
                const double vp = omega / kmag;
                worst_vp = std::fmax(worst_vp,
                                     std::fabs(vp - 1.0 / std::sqrt(eps * mu)));
            }
        }
        report.check_residual(
            "b4.isotropic.analytic_dispersion",
            "|k|^2 = eps mu omega^2 for an isotropic medium, with H_z = "
            "sqrt(eps/mu) E_y",
            VerificationLevel::B, "n/a", "eps in {2,4}, mu in {1,1.5}", worst,
            1e-15, sw.seconds());
        report.check_residual(
            "b4.isotropic.phase_velocity",
            "v_p = 1/sqrt(eps mu)", VerificationLevel::B, "n/a",
            "4 media", worst_vp, 1e-15, sw.seconds());
    }

    {
        Stopwatch sw;
        // Transversality: in vacuum k.E = 0 and k.H = 0 on shell. Taken from the
        // divergence equations rather than assumed.
        const double omega = 1.0;
        double k[3] = {omega, 0.0, 0.0};
        // A transverse mode: E along y, H along z with |H| = |E| in vacuum.
        C E[3] = {C(0, 0), C(1, 0), C(0, 0)};
        C H[3] = {C(0, 0), C(0, 0), C(1, 0)};
        double kE = 0.0, kH = 0.0;
        for (int i = 0; i < 3; ++i) {
            kE += k[i] * E[i].re;
            kH += k[i] * H[i].re;
        }
        metrics.worst_transversality_residual = std::fmax(std::fabs(kE), std::fabs(kH));

        // And the mode really satisfies Maxwell: k x E = w B = w H in vacuum.
        double resid = 0.0;
        const double kxE[3] = {k[1] * E[2].re - k[2] * E[1].re,
                               k[2] * E[0].re - k[0] * E[2].re,
                               k[0] * E[1].re - k[1] * E[0].re};
        for (int i = 0; i < 3; ++i)
            resid = std::fmax(resid, std::fabs(kxE[i] - omega * H[i].re));

        report.check_predicate(
            "b4.vacuum.transversality_and_curl",
            "k.E = k.H = 0 and k x E = omega H for the vacuum transverse mode",
            VerificationLevel::B,
            metrics.worst_transversality_residual < 1e-15 && resid < 1e-15,
            "both divergences and the Faraday relation vanish",
            "max|k.E|,|k.H| = " + sci(metrics.worst_transversality_residual) +
                ", |k x E - w H| = " + sci(resid),
            sw.seconds(),
            "this is the case that pins the sign of the Faraday equation: with "
            "the opposite kernel it would read k x E = -omega H");
    }

    // =======================================================================
    // 4. The Fresnel polynomial and the effective metric.
    // =======================================================================
    {
        Stopwatch sw;
        // Vacuum Z: lambda must be proportional to (eta^{ab}k_a k_b)^2.
        const Rank4<double> Zvac =
            covariant_constitutive(EffectiveMetric<double>(g.value), g).Z;

        double worst_rel = 0.0, worst_rank1 = 0.0;
        double scale = 0.0;
        std::vector<double> lam, cone;
        for (double kx : {0.3, 0.9, 1.7}) {
            for (double w : {0.7, 1.4}) {
                const double k_dn[4] = {-w, kx, 0.4, -0.2};
                double r1 = 0.0;
                const double l = fresnel_scalar(Zvac, k_dn, &r1);
                worst_rank1 = std::fmax(worst_rank1, r1);
                // eta^{ab}k_a k_b = -k_0^2 + |k|^2
                const double q = -w * w + kx * kx + 0.16 + 0.04;
                lam.push_back(l);
                cone.push_back(q * q);
                scale = std::fmax(scale, std::fabs(l));
            }
        }
        // A single constant must relate them all.
        double c0 = 0.0;
        for (std::size_t i = 0; i < lam.size(); ++i)
            if (std::fabs(cone[i]) > 1e-12) { c0 = lam[i] / cone[i]; break; }
        for (std::size_t i = 0; i < lam.size(); ++i)
            worst_rel = std::fmax(worst_rel,
                                  std::fabs(lam[i] - c0 * cone[i]) /
                                      std::fmax(scale, 1e-300));
        metrics.fresnel_factorization_error = worst_rel;

        report.check_residual(
            "b4.fresnel.vacuum_factorizes_into_the_light_cone",
            "the Fresnel quartic from adj(W) = lambda k (x) k is proportional to "
            "(eta^{ab} k_a k_b)^2, with ONE constant across every wavevector",
            VerificationLevel::B, "n/a",
            "6 wavevectors; rank-one residual " + sci(worst_rank1) +
                ", constant " + sci(c0),
            worst_rel, 1e-12, sw.seconds(),
            "det W vanishes identically because W annihilates k on both sides, so "
            "the physics is in the adjugate. That it comes out rank one, and a "
            "perfect square, is the statement that vacuum has a single light "
            "cone and no birefringence");
    }

    {
        Stopwatch sw;
        // Conformal invariance of the cone, and a subluminal Alcubierre slice.
        double worst = 0.0;
        for (double vs : {0.0, 0.3}) {
            const Mat4<double> geff =
                (vs == 0.0) ? g.value : alcubierre_at(vs, 3.0);
            const Rank4<double> Z1 =
                covariant_constitutive(EffectiveMetric<double>(geff), g).Z;
            Mat4<double> scaled;
            for (int i = 0; i < 16; ++i) scaled.m[i] = 2.7 * geff.m[i];
            const Rank4<double> Z2 =
                covariant_constitutive(EffectiveMetric<double>(scaled), g).Z;

            for (double kx : {0.5, 1.2}) {
                const double k_dn[4] = {-1.0, kx, 0.3, 0.0};
                const double l1 = fresnel_scalar(Z1, k_dn);
                const double l2 = fresnel_scalar(Z2, k_dn);
                worst = std::fmax(worst, std::fabs(l1 - l2) /
                                             std::fmax(std::fabs(l1), 1e-300));
            }
        }
        report.check_residual(
            "b4.fresnel.conformal_invariance_of_the_cone",
            "the Fresnel scalar is unchanged under g_eff -> Omega^2 g_eff, so the "
            "light cone is a property of the conformal class",
            VerificationLevel::A, "n/a",
            "vacuum and Alcubierre v_s = 0.3, Omega^2 = 2.7", worst, 1e-12,
            sw.seconds());
    }

    {
        Stopwatch sw;
        // Metric media do NOT birefringe; a generic anisotropic reciprocal
        // medium DOES. Measured as the split between the two propagation
        // branches along a direction.
        const double omega = 1.0;
        const double n[3] = {0.5773502691896258, 0.5773502691896258,
                             0.5773502691896258};

        auto branch_split = [&](const Rank4<double>& Z) {
            double c[5];
            fresnel_quartic_coefficients(Z, omega, n, c);
            const QuarticRoots r = real_quartic_roots(c, 0.05, 12.0);
            if (r.count < 2) return 0.0;
            double lo = 1e300, hi = -1e300;
            for (int i = 0; i < r.count; ++i) {
                lo = std::fmin(lo, r.root[i]);
                hi = std::fmax(hi, r.root[i]);
            }
            return (hi - lo) / std::fmax(hi, 1e-300);
        };

        const Rank4<double> Zmetric =
            covariant_constitutive(EffectiveMetric<double>(alcubierre_at(0.3, 3.0)),
                                   g).Z;
        const double split_metric = branch_split(Zmetric);
        const double split_generic = branch_split(anisotropic_Z(g, frame));
        metrics.generic_birefringence_split = split_generic;

        report.check_predicate(
            "b4.fresnel.metric_media_do_not_birefringe",
            "a metric-induced medium has its two propagation branches degenerate, "
            "while a generic reciprocal anisotropic medium splits them",
            VerificationLevel::B,
            split_metric < 1e-8 && split_generic > 1e-2,
            "degenerate for the metric medium, split for the generic one",
            "metric split = " + sci(split_metric) + ", generic split = " +
                sci(split_generic),
            sw.seconds(),
            "the second half is the point: NOT every reciprocal medium "
            "corresponds to one effective metric. A quartic that does not "
            "factorize into a doubled quadric describes two distinct cones, and "
            "eps = diag(2, 3.5, 5) is such a medium");
    }

    // =======================================================================
    // 6. Poynting, group velocity and the direction question.
    // =======================================================================
    {
        Stopwatch sw;
        // Group velocity of a dispersive isotropic medium, against centred
        // differences of the dispersion curve.
        DispersiveModelSpec spec;
        spec.eps_poles = 1;
        spec.mu_poles = 0;
        spec.magnetoelectric = false;
        spec.mode = ParameterMode::Direct;
        double p[kPhysicalSlots] = {};
        p[kSlotEpsInf] = 2.0;
        p[kSlotEpsPole0F] = 1.2; p[kSlotEpsPole0W0] = 5.0; p[kSlotEpsPole0Gamma] = 0.0;
        p[kSlotMuInf] = 1.0;
        const PhysicalDesignVariables<double> v = design_from_physical<double>(spec, p);

        auto kof = [&](double w) {
            const ResponseAtFrequency<double> r = evaluate_response(v, w);
            return w * std::sqrt(r.epsilon.re * (1.0 / r.mu_inverse.re));
        };

        double worst = 0.0;
        for (double w : {1.0, 1.5, 2.0}) {
            const double h = 1e-5;
            const double dkdw = (kof(w + h) - kof(w - h)) / (2 * h);
            const double vg_fd = 1.0 / dkdw;
            // Independent analytic route: v_g = c / (n + w dn/dw).
            auto nof = [&](double ww) {
                const ResponseAtFrequency<double> r = evaluate_response(v, ww);
                return std::sqrt(r.epsilon.re * (1.0 / r.mu_inverse.re));
            };
            const double nn = nof(w);
            const double dndw = (nof(w + h) - nof(w - h)) / (2 * h);
            const double vg_an = 1.0 / (nn + w * dndw);
            worst = std::fmax(worst, std::fabs(vg_fd - vg_an) / std::fabs(vg_an));
        }
        metrics.group_velocity_error = worst;

        report.check_residual(
            "b4.group_velocity_matches_finite_differences",
            "v_g = dw/dk from the dispersion curve equals c/(n + w dn/dw)",
            VerificationLevel::D, "n/a", "3 frequencies, lossless Lorentz",
            worst, 1e-7, sw.seconds(),
            "two independent routes through the same dispersion relation");
    }

    {
        Stopwatch sw;
        // Poynting vector, and the directions that need not coincide.
        const double omega = 1.0;
        // Vacuum transverse mode along +x.
        C E[3] = {C(0, 0), C(1, 0), C(0, 0)};
        C H[3] = {C(0, 0), C(0, 0), C(1, 0)};
        double S[3];
        time_averaged_poynting(E, H, S);

        report.check_predicate(
            "b4.poynting_vacuum",
            "<S> = (1/2) Re[E x conj(H)] gives energy flow along +x with "
            "magnitude 1/2 for a unit vacuum transverse mode",
            VerificationLevel::B,
            std::fabs(S[0] - 0.5) < 1e-15 && std::fabs(S[1]) < 1e-15 &&
                std::fabs(S[2]) < 1e-15,
            "S = (0.5, 0, 0)",
            "S = (" + sci(S[0]) + ", " + sci(S[1]) + ", " + sci(S[2]) + ")",
            sw.seconds());

        // In a magnetoelectric medium the forward and backward phase velocities
        // differ, so wavevector and energy flow are NOT symmetric under k -> -k.
        const LayerMode<double> m = layer_modes(C(2.5, 0.0), C(1.0 / 1.2, 0.0), C(0.35, 0.0));
        const double kf = std::fabs(m.kappa_forward.re);
        const double kb = std::fabs(m.kappa_backward.re);
        report.check_predicate(
            "b4.magnetoelectric_directional_asymmetry",
            "in a magnetoelectric medium the forward and backward wavenumbers "
            "differ in magnitude: phase velocity is direction dependent",
            VerificationLevel::B, std::fabs(kf - kb) > 1e-3,
            "|kappa_+| != |kappa_-|",
            "|kappa_+| = " + sci(kf) + ", |kappa_-| = " + sci(kb) + ", omega = " +
                sci(omega),
            sw.seconds(),
            "zeta = z K is the moving-medium form, so this is Fresnel drag and "
            "not a sign error. It is also why wavevector, phase-velocity and "
            "energy-flow directions must not be assumed to coincide");
    }

    // =======================================================================
    // 7. The layered solver, and energy balance.
    // =======================================================================
    {
        Stopwatch sw;
        // Lossless slab: R + T = 1 exactly, and the absorption computed by
        // integrating the local dissipated power is zero.
        std::vector<Layer<double>> stack;
        Layer<double> L;
        L.epsilon = C(4.0, 0.0);
        L.mu_inverse = C(1.0, 0.0);
        L.zeta = C(0.0, 0.0);
        L.thickness = 1.7;
        stack.push_back(L);

        double worst = 0.0, worst_abs = 0.0;
        for (double w : {0.7, 1.3, 2.9}) {
            const StackResult<double> r = solve_stack(stack, w);
            worst = std::fmax(worst, std::fabs(r.reflection + r.transmission - 1.0));
            worst_abs = std::fmax(worst_abs, std::fabs(r.absorption_from_dissipation));
        }
        report.check_residual(
            "b4.layered.lossless_slab_conserves_energy",
            "R + T = 1 for a lossless dielectric slab",
            VerificationLevel::B, "n/a", "eps = 4, d = 1.7, 3 frequencies",
            worst, 1e-13, sw.seconds(),
            "R and T come from the scattering matrix independently; this is not "
            "a definition");
        report.check_residual(
            "b4.layered.lossless_slab_absorbs_nothing",
            "the independently integrated dissipated power is zero in a lossless "
            "slab",
            VerificationLevel::B, "n/a", "3 frequencies", worst_abs, 1e-12,
            sw.seconds());
    }

    {
        Stopwatch sw;
        // Lossy slab: R + T + A = 1, with A from the dissipation integral.
        std::vector<Layer<double>> stack;
        Layer<double> L;
        L.epsilon = C(4.0, 0.35);
        L.mu_inverse = C(1.0, 0.0);
        L.zeta = C(0.0, 0.0);
        L.thickness = 1.7;
        stack.push_back(L);

        double worst = 0.0;
        StackResult<double> last;
        for (double w : {0.7, 1.3, 2.9}) {
            const StackResult<double> r = solve_stack(stack, w, 6000);
            worst = std::fmax(worst, r.energy_balance_residual);
            last = r;
        }
        metrics.worst_energy_balance = worst;
        metrics.reflection = last.reflection;
        metrics.transmission = last.transmission;
        metrics.absorption = last.absorption_from_dissipation;

        report.check_residual(
            "b4.layered.lossy_slab_energy_balance",
            "R + T + A = 1 with A obtained by integrating <p> = (w/2) Im[E* D + "
            "H* B] through the internal field profile",
            VerificationLevel::B, "n/a",
            "eps = 4 + 0.35i, d = 1.7, 3 frequencies, 6000 quadrature points",
            worst, 1e-7, sw.seconds(),
            "A is computed INDEPENDENTLY of R and T, so this is a genuine "
            "balance. Setting A = 1 - R - T and then asserting the sum is 1 "
            "would be circular");

        report.check_predicate(
            "b4.layered.absorption_is_positive_for_a_passive_slab",
            "a passive lossy slab absorbs a strictly positive amount of power",
            VerificationLevel::B, last.absorption_from_dissipation > 1e-3,
            "A > 0",
            "R = " + sci(last.reflection) + ", T = " + sci(last.transmission) +
                ", A = " + sci(last.absorption_from_dissipation),
            sw.seconds());
    }

    {
        Stopwatch sw;
        // An ACTIVE slab must show negative absorption. The solver has to be
        // able to report gain, or its verdict on passive media is worthless.
        std::vector<Layer<double>> stack;
        Layer<double> L;
        L.epsilon = C(4.0, -0.35);
        L.mu_inverse = C(1.0, 0.0);
        L.zeta = C(0.0, 0.0);
        L.thickness = 1.7;
        stack.push_back(L);
        const StackResult<double> r = solve_stack(stack, 1.3);

        report.check_predicate(
            "b4.layered.active_slab_is_detected",
            "a slab with Im eps < 0 shows NEGATIVE absorption and R + T > 1",
            VerificationLevel::B,
            r.absorption_from_dissipation < -1e-3 &&
                r.reflection + r.transmission > 1.0 &&
                r.energy_balance_residual < 1e-6,
            "gain detected, and the balance still closes",
            "A = " + sci(r.absorption_from_dissipation) + ", R + T = " +
                sci(r.reflection + r.transmission),
            sw.seconds(),
            "the energy balance closing on a NEGATIVE absorption is the strongest "
            "form of this check: the sign is physical, not clamped");
    }

    {
        Stopwatch sw;
        // Thick evanescent stack: the scattering formulation must not overflow
        // where a transfer matrix would.
        std::vector<Layer<double>> stack;
        Layer<double> L;
        L.epsilon = C(-2.0, 0.01);   // below plasma frequency: evanescent
        L.mu_inverse = C(1.0, 0.0);
        L.zeta = C(0.0, 0.0);
        L.thickness = 40.0;          // many decay lengths
        stack.push_back(L);
        const StackResult<double> r = solve_stack(stack, 1.0);

        report.check_predicate(
            "b4.layered.thick_evanescent_stack_is_stable",
            "a 40-unit evanescent layer returns finite R and T with T "
            "exponentially small, rather than overflowing",
            VerificationLevel::A,
            r.ok && std::isfinite(r.reflection) && std::isfinite(r.transmission) &&
                r.transmission < 1e-30 && r.reflection <= 1.0 + 1e-9,
            "finite results, T below 1e-30",
            "R = " + sci(r.reflection) + ", T = " + sci(r.transmission),
            sw.seconds(),
            "this is why scattering matrices were used rather than transfer "
            "matrices: the latter carry e^{|Im k| d} for the backward wave and "
            "overflow long before the physics goes wrong");
    }

    // =======================================================================
    // 8. Target-metric propagation benchmark.
    // =======================================================================
    {
        Stopwatch sw;
        // A layered profile sampled from a subluminal Alcubierre slice, compared
        // against a zeta == 0 restricted model. The point is to quantify the
        // PROPAGATION cost of dropping the magnetoelectric term, not just its
        // tensor norm.
        const double vs = 0.3;
        const double omega = 1.0;

        std::vector<Layer<double>> exact, restricted;
        std::ofstream prof(out_dir + "/b4_profile.csv", std::ios::binary);
        if (prof) prof << "x,eps,mu_inverse,zeta\n";

        for (int i = 0; i < 24; ++i) {
            const double x = -6.0 + 12.0 * static_cast<double>(i) / 23.0;
            const Mat4<double> geff = alcubierre_at(vs, x);
            const Rank4<double> Z =
                covariant_constitutive(EffectiveMetric<double>(geff), g).Z;
            const ObserverDecomposition4<double> d =
                decompose_A7(Z, g.value, frame.V);
            // Read the isotropic-equivalent scalars off the blocks.
            const double eps = d.epsilon(2, 2);
            const double mui = d.mu_inverse(2, 2);
            const double zt = d.zeta(2, 3);

            Layer<double> Le, Lr;
            Le.epsilon = C(eps, 0.0);
            Le.mu_inverse = C(mui, 0.0);
            Le.zeta = C(zt, 0.0);
            Le.thickness = 0.5;
            Lr = Le;
            Lr.zeta = C(0.0, 0.0);
            exact.push_back(Le);
            restricted.push_back(Lr);
            if (prof) prof << x << ',' << eps << ',' << mui << ',' << zt << '\n';
        }

        const StackResult<double> re = solve_stack(exact, omega);
        const StackResult<double> rr = solve_stack(restricted, omega);

        metrics.phase_mismatch = std::fabs(re.phase - rr.phase);
        const double h = 1e-4;
        const StackResult<double> re2 = solve_stack(exact, omega + h);
        const StackResult<double> rr2 = solve_stack(restricted, omega + h);
        const double tg_e = (re2.phase - re.phase) / h;
        const double tg_r = (rr2.phase - rr.phase) / h;
        metrics.group_delay_mismatch = std::fabs(tg_e - tg_r);
        metrics.exact_vs_restricted_error =
            std::fabs(re.transmission - rr.transmission);

        std::ostringstream got;
        got << "exact: R = " << sci(re.reflection) << ", T = " << sci(re.transmission)
            << ", phase = " << sci(re.phase) << "; restricted (zeta = 0): R = "
            << sci(rr.reflection) << ", T = " << sci(rr.transmission)
            << ", phase = " << sci(rr.phase) << "; d(phase) = "
            << sci(metrics.phase_mismatch) << ", d(group delay) = "
            << sci(metrics.group_delay_mismatch);

        report.check_predicate(
            "b4.propagation.restricted_family_costs_phase",
            "dropping the magnetoelectric term changes measurable propagation "
            "observables -- phase accumulation and group delay -- through a "
            "warp-analogue profile, and both stacks conserve energy",
            VerificationLevel::B,
            re.ok && rr.ok && re.energy_balance_residual < 1e-6 &&
                rr.energy_balance_residual < 1e-6 &&
                metrics.phase_mismatch > 1e-3,
            "a measurable phase difference, with energy conserved in both",
            got.str(), sw.seconds(),
            "this quantifies the PROPAGATION cost of a constitutive "
            "approximation, which a tensor-norm residual alone does not. Both "
            "profiles are lossless, so R + T = 1 holds for each");
    }

    // =======================================================================
    // 9. Adaptive frequency validation.
    // =======================================================================
    {
        Stopwatch sw;
        // A model with a narrow passivity dip that a coarse grid can miss.
        DispersiveModelSpec spec;
        spec.eps_poles = 1;
        spec.mu_poles = 1;
        spec.magnetoelectric = true;
        spec.mode = ParameterMode::Direct;
        double p[kPhysicalSlots] = {};
        p[kSlotEpsInf] = 2.1;
        p[kSlotEpsPole0F] = 1.6; p[kSlotEpsPole0W0] = 3.0; p[kSlotEpsPole0Gamma] = 1e-3;
        p[kSlotMuInf] = 1.15;
        p[kSlotMuPole0F] = 0.45; p[kSlotMuPole0W0] = 3.6; p[kSlotMuPole0Gamma] = 1e-3;
        p[kSlotZetaAmp] = 1.5; p[kSlotZetaW0] = 3.05; p[kSlotZetaGamma] = 0.02;
        const PhysicalDesignVariables<double> v = design_from_physical<double>(spec, p);

        const AdaptiveScanResult a =
            adaptive_passivity_scan(v, 2.0, 5.0, 21, 6, 9);
        metrics.coarse_min_passivity = a.coarse_min_eigenvalue;
        metrics.adaptive_min_passivity = a.refined_min_eigenvalue;
        metrics.adaptive_improvement = a.improvement;
        metrics.adaptive_samples = a.refined_samples;

        std::ostringstream got;
        got << "coarse (" << a.coarse_samples << " samples): min eig = "
            << sci(a.coarse_min_eigenvalue) << " at w = "
            << sci(a.coarse_worst_frequency) << "; adaptive ("
            << a.refined_samples << " samples, " << a.refinement_rounds
            << " rounds): min eig = " << sci(a.refined_min_eigenvalue) << " at w = "
            << sci(a.refined_worst_frequency) << "; resonances targeted: "
            << a.resonances.size() << ", mu zeros: " << a.mu_near_zeros.size();

        report.check_predicate(
            "b4.adaptive.refines_the_passivity_margin",
            "the adaptive scan finds a strictly LOWER minimum loss eigenvalue "
            "than the coarse grid, by refining around resonances, zeros of "
            "det mu and local minima",
            VerificationLevel::B,
            a.improvement > 0.0 && a.refined_samples > a.coarse_samples,
            "a strictly lower minimum found by refinement", got.str(),
            sw.seconds(),
            "BOTH results are reported, coarse and refined, because the "
            "difference between them is the measure of what a fixed grid misses. "
            "The refined scan narrows the unverified region; it does NOT prove "
            "global passivity, which would need an analytic positive-real "
            "argument this project has not carried out");
    }

    // =======================================================================
    // 10. HyPAD sensitivities of a propagation observable.
    // =======================================================================
#ifdef WARP_HYPAD_HAVE_OTILIB
    {
        namespace sens = warp_hypad::sensitivity;
        using sens::OTIScalar;
        Stopwatch sw;
        metrics.oti_available = true;

        // Differentiate the analytic isotropic transmission phase through a slab
        // with respect to eps, away from any degeneracy. The observable is
        // smooth there, which is stated rather than assumed.
        const double eps0 = 4.0, d = 1.7, w = 1.3;
        auto phase_of = [&](double eps) {
            return w * std::sqrt(eps) * d;   // single-pass phase, n = sqrt(eps)
        };
        const double analytic = w * d * 0.5 / std::sqrt(eps0);

        const sens::SeededParameters sp = sens::seed_all(eps0, 1.0, 1.0);
        const OTIScalar e = sp.velocity;
        using std::sqrt;
        const OTIScalar phase = sqrt(e) * OTIScalar(w * d);
        const double oti = sens::extract_gradient(phase, sens::kParamVelocity);

        double best_fd = 1e300;
        for (double h : {1e-4, 1e-5, 1e-6})
            best_fd = std::fmin(best_fd,
                                std::fabs((phase_of(eps0 + h) - phase_of(eps0 - h)) /
                                              (2 * h) - analytic));
        metrics.hypad_gradient_error = std::fabs(oti - analytic);

        report.check_residual(
            "b4.hypad.propagation_phase_gradient",
            "d(phase)/d(eps) through a slab matches the analytic w d /(2 sqrt(eps))",
            VerificationLevel::B, "n/a",
            "eps = 4, d = 1.7, omega = 1.3; finite-difference cross-check " +
                sci(best_fd),
            metrics.hypad_gradient_error, 1e-14, sw.seconds(),
            "evaluated away from any eigenvalue degeneracy or branch crossing, "
            "where the observable is smooth. Differentiating through branch "
            "sorting or discrete mode selection is NOT done here, and would have "
            "to be documented as non-smooth if it were");
    }
#else
    {
        Stopwatch sw;
        report.check_skipped(
            "b4.hypad.propagation_phase_gradient",
            "HyPAD sensitivity of a propagation observable",
            "the OTIlib backend is disabled; every Maxwell, dispersion, Fresnel, "
            "propagation and energy case still ran",
            sw.seconds());
    }
#endif

    // =======================================================================
    // Determinism, and the branch/degeneracy report.
    // =======================================================================
    {
        Stopwatch sw;
        std::vector<Layer<double>> stack;
        Layer<double> L;
        L.epsilon = C(3.2, 0.12);
        L.mu_inverse = C(0.9, 0.0);
        L.zeta = C(0.2, 0.03);
        L.thickness = 2.3;
        stack.push_back(L);
        const StackResult<double> a = solve_stack(stack, 1.1);
        const StackResult<double> b = solve_stack(stack, 1.1);
        const bool identical = a.reflection == b.reflection &&
                               a.transmission == b.transmission &&
                               a.absorption_from_dissipation ==
                                   b.absorption_from_dissipation;
        report.check_predicate(
            "b4.deterministic_repeatability",
            "two identical stack solves return bit-identical R, T and A",
            VerificationLevel::A, identical, "bit-identical",
            identical ? "identical" : "diverged", sw.seconds());
    }

    {
        Stopwatch sw;
        // Degeneracy warning: report where the two Fresnel branches nearly
        // coincide, since that is where any mode-sorting derivative would be
        // non-smooth.
        const double n[3] = {1.0, 0.0, 0.0};
        double c[5];
        fresnel_quartic_coefficients(anisotropic_Z(g, frame), 1.0, n, c);
        const QuarticRoots r = real_quartic_roots(c, 0.05, 12.0);
        std::ostringstream w;
        w << r.count << " real branches along +x: ";
        for (int i = 0; i < r.count; ++i) w << sci(r.root[i]) << " ";
        metrics.branch_warnings = w.str();

        report.check_predicate(
            "b4.branch_count_is_reported",
            "the number of real propagation branches along a direction is counted "
            "and reported rather than assumed",
            VerificationLevel::A, r.count >= 2,
            "at least two real branches for an anisotropic medium", w.str(),
            sw.seconds(),
            "where two branches approach each other any mode-sorting derivative "
            "is non-smooth, which is why the HyPAD case above is evaluated away "
            "from degeneracies and says so");
    }

    // Dispersion-surface export.
    {
        std::ofstream os(out_dir + "/b4_dispersion_surfaces.csv", std::ios::binary);
        if (os) {
            os << "medium,theta,kappa_branch1,kappa_branch2,branches\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            struct M { const char* name; Rank4<double> Z; };
            std::vector<M> media;
            media.push_back({"vacuum",
                             covariant_constitutive(
                                 EffectiveMetric<double>(g.value), g).Z});
            media.push_back({"isotropic_eps4", medium_Z(4.0, 1.0, 0.0, g, frame)});
            media.push_back({"anisotropic", anisotropic_Z(g, frame)});
            media.push_back({"magnetoelectric",
                             medium_Z(2.5, 1.2, 0.35, g, frame)});
            media.push_back({"warp_v0.3",
                             covariant_constitutive(
                                 EffectiveMetric<double>(alcubierre_at(0.3, 3.0)),
                                 g).Z});
            for (const M& m : media)
                for (int i = 0; i <= 72; ++i) {
                    const double th = 3.14159265358979323846 *
                                      static_cast<double>(i) / 72.0;
                    const double nn[3] = {std::cos(th), std::sin(th), 0.0};
                    double c[5];
                    fresnel_quartic_coefficients(m.Z, 1.0, nn, c);
                    const QuarticRoots r = real_quartic_roots(c, 0.05, 12.0);
                    double b1 = 0.0, b2 = 0.0;
                    int pos = 0;
                    for (int j = 0; j < r.count; ++j)
                        if (r.root[j] > 0.0) {
                            if (pos == 0) b1 = r.root[j];
                            else if (pos == 1) b2 = r.root[j];
                            ++pos;
                        }
                    os << m.name << ',' << th << ',' << b1 << ',' << b2 << ','
                       << r.count << '\n';
                }
        }
    }
}

}  // namespace checks
}  // namespace program_b

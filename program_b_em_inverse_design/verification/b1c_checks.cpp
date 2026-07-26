#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "b1c_checks.hpp"
#include "program_b/effective_metric/bel_extraction.hpp"
#include "program_b/effective_metric/covariant_constitutive.hpp"
#include "program_b/effective_metric/field_reconstruction.hpp"
#include "program_b/effective_metric/metric_medium_map.hpp"
#include "program_b/effective_metric/metric_medium_status.hpp"
#include "program_b/effective_metric/z_reconstruction.hpp"
#include "warp_hypad/geometry/signature.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

// ---------------------------------------------------------------------------
// The AGGREGATE B1C gate.
//
// It does NOT include legacy_unverified.hpp. Every constitutive block here comes
// from decompose_A7(); every reconstruction from reconstruct_Z_A11(); every
// observer from the typed laboratory/effective-geometry constructors. See
// docs/b1c_aggregate_gate_audit.md for the case-by-case record of what was
// replaced and why.
//
// SCOPE. The sub-gates B1C-S23 and B1C-S456 are separate binaries and separate
// CTest tests, and this gate DEPENDS on them rather than repeating them. What
// lives here is only what is cross-cutting:
//
//   * the construction of Z from a pair of metrics -- vacuum, conformal
//     invariance, tensoriality. The sub-gates take Z as given.
//   * the warp analogue's blocks against closed forms, over the profile.
//   * the ergoregion correction, and what does and does not diverge there.
//   * the two velocity-sensitivity modes.
//   * objective invariance.
//   * the five-way status classification.
//   * HyPAD channels through Z.
// ---------------------------------------------------------------------------

namespace program_b {
namespace checks {

using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;
using warp_hypad::metrics::AlcubierreMetric;
using warp_hypad::metrics::AlcubierreParameters;

namespace {

const double kR = 3.0, kSigma = 2.0;

LaboratoryMetric<double> eta_lab() {
    return LaboratoryMetric<double>::minkowski();
}
Mat4<double> minkowski() { return Mat4<double>::minkowski(); }

Mat4<double> alcubierre_at(double vs, double x, double t = 0.0) {
    AlcubierreMetric<double> w;
    w.time_treatment = warp_hypad::metrics::TimeTreatment::InertialTimeDependent;
    const AlcubierreParameters<double> p{vs, kR, kSigma};
    return w.evaluate(t, x, 0.0, 0.0, p);
}

double shape(double r) {
    return (std::tanh(kSigma * (r + kR)) - std::tanh(kSigma * (r - kR))) /
           (2 * std::tanh(kSigma * kR));
}

double diff4(const Mat4<double>& a, const Mat4<double>& b) {
    double w = 0.0;
    for (int i = 0; i < 16; ++i) w = std::fmax(w, std::fabs(a.m[i] - b.m[i]));
    return w;
}
double max_abs4(const Mat4<double>& a) {
    double w = 0.0;
    for (int i = 0; i < 16; ++i) w = std::fmax(w, std::fabs(a.m[i]));
    return w;
}
double diff4x4(const Rank4<double>& a, const Rank4<double>& b) {
    double w = 0.0;
    for (int i = 0; i < 256; ++i) w = std::fmax(w, std::fabs(a.t[i] - b.t[i]));
    return w;
}
double max_abs4x4(const Rank4<double>& a) {
    double w = 0.0;
    for (int i = 0; i < 256; ++i) w = std::fmax(w, std::fabs(a.t[i]));
    return w;
}

// Eigenvalue diagnostics of the SPATIAL 3x3 part of a block, valid for
// g_lab = eta with V = (1,0,0,0), where the spatial slice is indices 1..3.
struct Spectrum {
    double min_abs = 0.0;
    double max_abs = 0.0;
    double determinant = 0.0;
    double condition = 0.0;
    double max_component = 0.0;
};

Spectrum spectrum_of(const Mat4<double>& block) {
    double m[3][3];
    Spectrum s;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            m[i][j] = block(i + 1, j + 1);
            s.max_component = std::fmax(s.max_component, std::fabs(m[i][j]));
        }
    double eig[3];
    symmetric_eigenvalues_3x3(m, eig);
    s.min_abs = std::fmin(std::fabs(eig[0]),
                          std::fmin(std::fabs(eig[1]), std::fabs(eig[2])));
    s.max_abs = std::fmax(std::fabs(eig[0]),
                          std::fmax(std::fabs(eig[1]), std::fabs(eig[2])));
    s.determinant = eig[0] * eig[1] * eig[2];
    s.condition = (s.min_abs > 0.0) ? s.max_abs / s.min_abs
                                    : std::numeric_limits<double>::infinity();
    return s;
}

std::string sci(double v) {
    std::ostringstream o;
    o.setf(std::ios::scientific);
    o.precision(4);
    o << v;
    return o.str();
}

}  // namespace

void run_b1c_checks(VerificationReport& report, const std::string& out_dir,
                    VelocityMetrics& velocity) {
    const LaboratoryMetric<double> lab = eta_lab();
    const Mat4<double> eta = minkowski();
    double V_rest[4] = {1.0, 0.0, 0.0, 0.0};
    const Observer<double> lab_rest =
        make_laboratory_observer(lab, V_rest, "lab_rest");

    // =======================================================================
    // 1. Covariant Z for Minkowski. Pins every factor in the definition.
    // =======================================================================
    {
        Stopwatch sw;
        const CovariantConstitutive<double> cc =
            covariant_constitutive(EffectiveMetric<double>(eta), lab);

        Mat4<double> F;
        const double e[3] = {0.3, -0.7, 0.5}, b[3] = {0.2, 0.9, -0.4};
        F(0, 1) = e[0]; F(1, 0) = -e[0];
        F(0, 2) = e[1]; F(2, 0) = -e[1];
        F(0, 3) = e[2]; F(3, 0) = -e[2];
        F(2, 3) = b[0]; F(3, 2) = -b[0];
        F(3, 1) = b[1]; F(1, 3) = -b[1];
        F(1, 2) = b[2]; F(2, 1) = -b[2];

        const double worst =
            cc.ok ? diff4(excitation_direct(cc.Z, F), raise_two_form(F, lab)) : 1.0;

        report.check_residual(
            "b1c.covariant_minkowski_is_vacuum",
            "H^{ab} = Z^{abcd} F_{cd} reproduces H = F exactly for Minkowski, "
            "on a general antisymmetric F",
            VerificationLevel::A, "n/a", "Heaviside-Lorentz, c = 1", worst,
            1e-15, sw.seconds(),
            "pins the 1/2 in the definition of Z and the absence of a second "
            "1/2 in the constitutive relation");
    }

    // =======================================================================
    // 2. Conformal invariance of Z.
    // =======================================================================
    {
        Stopwatch sw;
        const Mat4<double> g = alcubierre_at(0.3, 2.5);
        const CovariantConstitutive<double> base =
            covariant_constitutive(EffectiveMetric<double>(g), lab);

        double worst = base.ok ? 0.0 : 1.0;
        for (double om2 : {0.25, 0.5, 2.0, 7.0}) {
            Mat4<double> gs;
            for (int i = 0; i < 16; ++i) gs.m[i] = om2 * g.m[i];
            const CovariantConstitutive<double> cs =
                covariant_constitutive(EffectiveMetric<double>(gs), lab);
            if (!cs.ok) { worst = 1.0; continue; }
            worst = std::fmax(worst, diff4x4(cs.Z, base.Z));
        }
        report.check_residual(
            "b1c.covariant_conformal_invariance",
            "Z^{abcd} is unchanged under g_eff -> Omega^2 g_eff",
            VerificationLevel::A, "n/a", "Omega^2 = 0.25, 0.5, 2, 7", worst,
            1e-14, sw.seconds(),
            "density factor Omega^4 against Omega^-4 from the two inverse "
            "metrics");
    }

    // =======================================================================
    // 3. Tensorial transformation under a coordinate change.
    // =======================================================================
    {
        Stopwatch sw;
        Mat4<double> L = Mat4<double>::identity();
        const double c = std::cos(0.6), s = std::sin(0.6);
        L(1, 1) = c;  L(1, 2) = -s;
        L(2, 1) = s;  L(2, 2) = c;
        L(3, 3) = 1.2;   // anisotropic rescale, so this is not merely a rotation

        const warp_hypad::InverseResult<double> Li = warp_hypad::inverse(L);
        const Mat4<double> g = alcubierre_at(0.3, 2.5);

        auto push = [&](const Mat4<double>& m) {
            Mat4<double> o;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    double acc = 0.0;
                    for (int p = 0; p < 4; ++p)
                        for (int q = 0; q < 4; ++q)
                            acc += Li.inverse(p, a) * Li.inverse(q, b) * m(p, q);
                    o(a, b) = acc;
                }
            return o;
        };

        const CovariantConstitutive<double> Zc =
            covariant_constitutive(EffectiveMetric<double>(g), lab);
        const CovariantConstitutive<double> Zp = covariant_constitutive(
            EffectiveMetric<double>(push(g)), LaboratoryMetric<double>(push(eta)));

        double worst = (Zc.ok && Zp.ok && Li.ok) ? 0.0 : 1.0;
        if (Zc.ok && Zp.ok && Li.ok) {
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int cc2 = 0; cc2 < 4; ++cc2)
                        for (int d = 0; d < 4; ++d) {
                            double acc = 0.0;
                            for (int p = 0; p < 4; ++p)
                                for (int q = 0; q < 4; ++q)
                                    for (int r = 0; r < 4; ++r)
                                        for (int t2 = 0; t2 < 4; ++t2)
                                            acc += L(a, p) * L(b, q) * L(cc2, r) *
                                                   L(d, t2) * Zc.Z(p, q, r, t2);
                            worst = std::fmax(worst, std::fabs(Zp.Z(a, b, cc2, d) - acc));
                        }
        }
        report.check_residual(
            "b1c.covariant_transforms_tensorially",
            "Z'^{abcd} = L^a_p L^b_q L^c_r L^d_s Z^{pqrs} under a linear change "
            "of coordinates applied to both g_eff and g_lab",
            VerificationLevel::A, "n/a", "rotation by 0.6 rad plus a z rescale",
            worst, 1e-13, sw.seconds(),
            "a coordinate artefact must not be mistaken for invariant "
            "non-realizability, so the tensorial character is checked directly");
    }

    // =======================================================================
    // 4. The warp analogue's blocks against CLOSED FORMS.
    //
    // This is what replaces the two withdrawn Plebanski equality cases. For
    // g_lab = eta and the laboratory rest observer V = (1,0,0,0), the (A7)
    // blocks of the Alcubierre analogue have exact closed forms, derivable
    // because the metric's spatial block is delta_ij and det g = -1 identically:
    //
    //     eps^{ab}     = h^{ab} = diag(0,1,1,1)                     EXACTLY,
    //                    for every v_s and f, including v_s f > 1
    //     [mu^-1]^{ab} = diag(0, 1, -g_00, -g_00),  -g_00 = 1 - v^2 f^2
    //     zeta^{23}    = -g_01 = v f,  zeta^{32} = +g_01,  rest zero
    //
    // eps = h follows from eps^{ij} = -g_eff^{00} (gamma^{-1})^{ij} with
    // gamma_ij = delta_ij and g_eff^{00} = det(gamma)/det(g_eff) = -1.
    // =======================================================================
    {
        Stopwatch sw;
        double w_eps = 0.0, w_mu = 0.0, w_zeta = 0.0, w_recip = 0.0;
        const Mat4<double> h = projector_upper(eta, V_rest);

        for (double vs : {0.2, 0.4, 1.0, 2.0})
            for (double x : {0.0, 1.5, 2.5, 3.0, 3.385, 4.0, 6.0}) {
                const Mat4<double> g = alcubierre_at(vs, x);
                const Rank4<double> Z =
                    covariant_constitutive(EffectiveMetric<double>(g), lab).Z;
                const ObserverDecomposition4<double> d =
                    decompose_A7(Z, eta, lab_rest.V);

                Mat4<double> mu_ref, zeta_ref;
                mu_ref(1, 1) = 1.0;
                mu_ref(2, 2) = -g(0, 0);
                mu_ref(3, 3) = -g(0, 0);
                zeta_ref(2, 3) = -g(0, 1);
                zeta_ref(3, 2) = g(0, 1);

                w_eps = std::fmax(w_eps, diff4(d.epsilon, h));
                w_mu = std::fmax(w_mu, diff4(d.mu_inverse, mu_ref));
                w_zeta = std::fmax(w_zeta, diff4(d.zeta, zeta_ref));
                w_recip = std::fmax(w_recip, reciprocity_residual(d));
            }

        const std::string params = "v_s in {0.2,0.4,1,2}, 7 positions each";
        report.check_residual(
            "b1c.warp_epsilon_closed_form",
            "eps^{ab} = h^{ab} exactly for the Alcubierre analogue with "
            "g_lab = eta and the laboratory rest observer, at every v_s",
            VerificationLevel::B, "n/a", params, w_eps, 1e-13, sw.seconds(),
            "derived: eps^{ij} = -g^{00}(gamma^-1)^{ij}, and the Alcubierre "
            "spatial block is delta_ij with det g = -1, so g^{00} = -1. This is "
            "why the withdrawn cond(eps) divergence test could never fire");
        report.check_residual(
            "b1c.warp_mu_inverse_closed_form",
            "[mu^-1]^{ab} = diag(0, 1, -g_00, -g_00) with -g_00 = 1 - v^2 f^2",
            VerificationLevel::B, "n/a", params, w_mu, 1e-13, sw.seconds(),
            "the transverse eigenvalue is exactly -g_00 and therefore passes "
            "through ZERO on the ergosurface: mu itself diverges there. That is "
            "the quantity which actually blows up for this observer");
        report.check_residual(
            "b1c.warp_zeta_closed_form",
            "zeta^{23} = -g_01 = v f, zeta^{32} = +g_01, all other components "
            "zero",
            VerificationLevel::B, "n/a", params, w_zeta, 1e-13, sw.seconds(),
            "an exact signed closed form for the warp medium's magnetoelectric "
            "block. Replaces the withdrawn min(error_to_xi, error_to_zeta) "
            "comparison, which asserted no sign at all");
        report.check_residual(
            "b1c.warp_blocks_are_reciprocal",
            "zeta_adjoint = zeta^T for the warp analogue, as required by the "
            "pair symmetry Z^{abcd} = Z^{cdab}",
            VerificationLevel::A, "n/a", params, w_recip, 1e-14, sw.seconds());
    }

    // The handedness of the warp medium's magnetoelectric response, asserted in
    // ONE direction.
    {
        Stopwatch sw;
        const double vs = 0.4, x = 3.0;
        const Mat4<double> g = alcubierre_at(vs, x);
        const Rank4<double> Z =
            covariant_constitutive(EffectiveMetric<double>(g), lab).Z;
        const ObserverDecomposition4<double> d = decompose_A7(Z, eta, lab_rest.V);

        const double E3[3] = {0.0, 0.0, 0.0};
        const double B3[3] = {0.0, 0.0, 1.0};   // B = +z
        const Mat4<double> F = field_strength_from_3vectors<double>(E3, B3);
        const ObserverFields<double> f = decompose_fields(d, F, lab, lab_rest.V);

        const double expected = g(0, 1);   // D_y = -zeta^{23} B_z = g_01 = -v f
        std::ostringstream got;
        got.setf(std::ios::fixed);
        got.precision(9);
        got << "D = (" << f.D.up[1] << ", " << f.D.up[2] << ", " << f.D.up[3]
            << "), expected D_y = " << expected;

        report.check_predicate(
            "b1c.warp_magnetoelectric_handedness",
            "for the warp analogue with the bubble moving along +x, B = +B0 "
            "z-hat and E = 0, the induced displacement points along -y",
            VerificationLevel::B,
            f.D.up[2] < -1e-6 && std::fabs(f.D.up[1]) < 1e-13 &&
                std::fabs(f.D.up[3]) < 1e-13,
            "D_y strictly negative, D_x = D_z = 0", got.str(), sw.seconds(),
            "one direction is asserted; a +y result fails. Same handedness as "
            "the moving isotropic medium of B1C-S456, which is checked against "
            "Minkowski's moving-medium relations");
        report.check_value(
            "b1c.warp_magnetoelectric_magnitude",
            "D_y = g_01 = -v f B0 exactly", VerificationLevel::B, "n/a",
            "v_s = 0.4 at the wall, B0 = 1", expected, f.D.up[2], 1e-13,
            sw.seconds());
    }

    // The Plebanski blocks and the covariant blocks belong to DIFFERENT
    // observers, and this records the exact factor between them. It is what the
    // two withdrawn equality cases were really measuring.
    {
        Stopwatch sw;
        const double vs = 0.2;
        double w_pleb = 0.0, w_mu_identity = 0.0, biggest_gap = 0.0;

        for (double x : {0.0, 1.5, 2.5, 3.0, 4.0, 6.0}) {
            const Mat4<double> g = alcubierre_at(vs, x);
            const Rank4<double> Z =
                covariant_constitutive(EffectiveMetric<double>(g), lab).Z;
            const ObserverDecomposition4<double> d =
                decompose_A7(Z, eta, lab_rest.V);
            const ForwardResult<double> pl = metric_to_medium(g);
            if (!pl.ok) { w_pleb = 1.0; continue; }

            const double vf = vs * shape(std::fabs(x));
            const double factor = 1.0 / (1.0 - vf * vf);

            // Plebanski's own closed form, in ITS observer's normalisation.
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    const double want =
                        (i != j) ? 0.0 : (i == 0 ? 1.0 : factor);
                    w_pleb = std::fmax(w_pleb,
                                       std::fabs(pl.response.epsilon(i, j) - want));
                }

            // [mu^-1] from the laboratory observer IS the inverse of the
            // Plebanski mu, componentwise -- an identity special to this metric
            // family, not a statement that the two observers agree.
            const Mat3Inverse<double> mi = mat3_inverse(pl.response.mu);
            if (!mi.ok) { w_mu_identity = 1.0; continue; }
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    w_mu_identity = std::fmax(
                        w_mu_identity,
                        std::fabs(d.mu_inverse(i + 1, j + 1) - mi.inverse(i, j)));

            biggest_gap = std::fmax(biggest_gap, factor - 1.0);
        }

        report.check_residual(
            "b1c.plebanski_epsilon_closed_form",
            "the Plebanski eps is diag(1, 1/(1-v^2f^2), 1/(1-v^2f^2)) for the "
            "Alcubierre metric",
            VerificationLevel::B, "n/a", "v_s = 0.2, six positions", w_pleb,
            1e-13, sw.seconds(),
            "its own closed form, in its own g_00 normalisation. Compare with "
            "b1c.warp_epsilon_closed_form, which is the IDENTITY: the two "
            "differ by exactly 1/(1-v^2f^2) and cannot both be eps for one "
            "observer");
        report.check_predicate(
            "b1c.plebanski_differs_by_observer_normalisation",
            "the covariant eps for the laboratory observer and the Plebanski eps "
            "differ, by exactly the ratio of the two observers' normalisations",
            VerificationLevel::B, biggest_gap > 1e-3,
            "a measurable, quantified discrepancy",
            "max(eps_pleb/eps_lab - 1) = " + sci(biggest_gap) +
                " = 1/(1-v_s^2 f^2) - 1",
            sw.seconds(),
            "THE WITHDRAWN EQUALITY. The old case asserted these were equal to "
            "1e-12 and failed by 4.1667e-2 = 1/24 = 1/(1-0.04) - 1 at v_s=0.2, "
            "f=1. Plebanski's 3+1 tensors are normalised by g_00 of g_eff; the "
            "(A7) blocks are normalised by g_lab. They are different observers "
            "and must not be compared for equality");
        report.check_residual(
            "b1c.plebanski_mu_inverse_coincides",
            "[mu^-1] from the (A7) laboratory decomposition equals the 3x3 "
            "inverse of the Plebanski mu, componentwise",
            VerificationLevel::D, "n/a", "v_s = 0.2, six positions",
            w_mu_identity, 1e-12, sw.seconds(),
            "this coincidence is real and is retained as a cross-check between "
            "two independent routes, but it is NOT evidence that the two "
            "observers agree -- eps demonstrably does not. It is special to a "
            "metric whose spatial block is delta_ij");
    }

    // =======================================================================
    // 5. Observer change over the warp profile. The authoritative replacement
    //    for the withdrawn observer-equality comparison.
    // =======================================================================
    {
        Stopwatch sw;
        double U_raw[4] = {1.0, 0.0, 0.5, 0.0};
        const Observer<double> U =
            make_laboratory_observer(lab, U_raw, "boosted_beta_y_0.5");

        double worst = 0.0, smallest_block_gap = 1e300;
        int pairs = 0;
        for (double vs : {0.2, 0.4, 2.0})
            for (double x : {0.0, 1.5, 2.5, 3.0, 4.0, 6.0}) {
                const Mat4<double> g = alcubierre_at(vs, x);
                const Rank4<double> Z =
                    covariant_constitutive(EffectiveMetric<double>(g), lab).Z;

                const ObserverDecomposition4<double> bV =
                    decompose_A7(Z, eta, lab_rest.V);
                const ObserverDecomposition4<double> bU = decompose_A7(Z, eta, U.V);

                const Rank4<double> ZV = reconstruct_Z_A11(bV, lab, lab_rest.V);
                const Rank4<double> ZU = reconstruct_Z_A11(bU, lab, U.V);

                worst = std::fmax(worst, diff4x4(ZV, Z));
                worst = std::fmax(worst, diff4x4(ZU, Z));
                worst = std::fmax(worst, diff4x4(ZV, ZU));

                // Only count the gap where the medium is not vacuum; at x = 6
                // the bubble has died away and every observer sees vacuum.
                const double gap = std::fmax(diff4(bV.mu_inverse, bU.mu_inverse),
                                             diff4(bV.zeta, bU.zeta));
                if (max_abs4(bV.zeta) > 1e-3) {
                    smallest_block_gap = std::fmin(smallest_block_gap, gap);
                    ++pairs;
                }
            }

        report.check_predicate(
            "b1c.observer_change_blocks_differ_over_profile",
            "at every profile point where the medium is not vacuum, two "
            "admissible laboratory observers see different blocks",
            VerificationLevel::A, pairs > 0 && smallest_block_gap > 1e-3,
            "smallest block difference above 1e-3 over all such points",
            std::to_string(pairs) + " magnetoelectric points; smallest gap " +
                sci(smallest_block_gap),
            sw.seconds(),
            "without this the agreement below would say nothing about observer "
            "independence");
        report.check_residual(
            "b1c.observer_change_reconstruction_over_profile",
            "Zhat[V] = Z, Zhat[U] = Z and Zhat[V] = Zhat[U] across the warp "
            "profile: the blocks are observer dependent, the tensor is not",
            VerificationLevel::A, "n/a",
            "v_s in {0.2,0.4,2}, six positions each, V = lab rest and "
            "U = boosted beta_y = 0.5",
            worst, 1e-12, sw.seconds(),
            "this is what replaces the withdrawn comparison of one observer's "
            "blocks against another's. Observer-dependent quantities are not "
            "asserted equal; the tensor they rebuild is");
    }

    // =======================================================================
    // 6. THE CORRECTION. The coordinate-time observer of g_eff fails at
    //    g_00 >= 0 while the covariant tensor remains perfectly valid.
    // =======================================================================
    {
        Stopwatch sw;
        int inside = 0, cov_ok_inside = 0, obs_bad_inside = 0, lorentzian_inside = 0;

        for (int i = -320; i <= 320; ++i) {
            const double x = 0.025 * i;
            const Mat4<double> g = alcubierre_at(2.0, x);
            if (!(g(0, 0) >= 0.0)) continue;  // only the ergoregion
            ++inside;

            const CovariantConstitutive<double> cc =
                covariant_constitutive(EffectiveMetric<double>(g), lab);
            if (cc.ok) ++cov_ok_inside;
            if (warp_hypad::real_part(warp_hypad::det(g)) < 0.0) ++lorentzian_inside;

            // The observer built on the EFFECTIVE metric is the one that fails.
            const EffectiveGeometryObserver<double> obs =
                make_effective_geometry_observer(EffectiveMetric<double>(g), V_rest,
                                                 "effective_coordinate_time");
            if (!obs.admissible) ++obs_bad_inside;
        }

        std::ostringstream got;
        got << inside << " points with g_00 >= 0; covariant Z valid at "
            << cov_ok_inside << "; g_eff Lorentzian at " << lorentzian_inside
            << "; coordinate-time observer of g_eff inadmissible at "
            << obs_bad_inside;

        report.check_predicate(
            "b1c.ergoregion_covariant_valid_observer_fails",
            "inside v_s f > 1 the effective metric stays Lorentzian and "
            "Z^{abcd} stays valid at EVERY point, while the coordinate-time "
            "observer of g_eff is inadmissible at every point",
            VerificationLevel::A,
            inside > 0 && cov_ok_inside == inside &&
                lorentzian_inside == inside && obs_bad_inside == inside,
            "covariant valid everywhere, observer invalid everywhere",
            got.str(), sw.seconds(),
            "THIS IS THE CORRECTION. B1 reported 'no equivalent medium exists' "
            "here. What actually fails is one observer's decomposition. The "
            "surface is an ERGOSURFACE, not an invariant realizability wall");
    }

    // The metric is Lorentzian for every velocity tested, including v_s = 10.
    {
        Stopwatch sw;
        int bad = 0, total = 0;
        double worst_det = 0.0;
        for (double vs : {0.2, 1.0, 2.0, 5.0, 10.0})
            for (int i = -200; i <= 200; i += 5) {
                ++total;
                const Mat4<double> g = alcubierre_at(vs, 0.04 * i);
                const auto sig = warp_hypad::geometry::analyse_signature(g);
                if (!sig.lorentzian) ++bad;
                worst_det = std::fmax(worst_det,
                                      std::fabs(warp_hypad::real_part(warp_hypad::det(g)) + 1.0));
            }
        std::ostringstream got;
        got << bad << " of " << total << " points non-Lorentzian; worst "
            << "|det g + 1| = " << worst_det;
        report.check_predicate(
            "b1c.alcubierre_lorentzian_at_all_velocities",
            "the Alcubierre metric has signature (-,+,+,+) and det g = -1 at "
            "every point for v_s up to 10",
            VerificationLevel::A, bad == 0 && worst_det < 1e-12,
            "0 non-Lorentzian points", got.str(), sw.seconds(),
            "the (t,x) block has determinant -1 identically, so its two "
            "eigenvalues always have opposite signs however large v_s f is. "
            "This is why no invariant obstruction exists at v_s f = 1");
    }

    // =======================================================================
    // 7. Alternative observers inside the ergoregion.
    // =======================================================================
    {
        Stopwatch sw;
        int inside = 0, eulerian_ok = 0, lab_coord_ok = 0, lab_blocks_ok = 0;

        for (int i = -200; i <= 200; ++i) {
            const double x = 0.03 * i;
            const Mat4<double> g = alcubierre_at(2.0, x);
            if (!(g(0, 0) >= 0.0)) continue;
            ++inside;

            // Eulerian observer of the EFFECTIVE metric: normal to the slicing.
            if (eulerian_observer(g).admissible) ++eulerian_ok;
            // The laboratory observer is always fine: g_lab = eta.
            if (lab_rest.admissible) ++lab_coord_ok;

            // And the laboratory decomposition itself is finite here.
            const Rank4<double> Z =
                covariant_constitutive(EffectiveMetric<double>(g), lab).Z;
            const ObserverDecomposition4<double> d =
                decompose_A7(Z, eta, lab_rest.V);
            const bool finite = std::isfinite(max_abs4(d.epsilon)) &&
                                std::isfinite(max_abs4(d.mu_inverse)) &&
                                std::isfinite(max_abs4(d.zeta)) &&
                                max_abs4(d.epsilon) < 1e3;
            if (finite) ++lab_blocks_ok;
        }

        std::ostringstream got;
        got << inside << " ergoregion points: Eulerian observer of g_eff "
            << "admissible at " << eulerian_ok << ", laboratory observer at "
            << lab_coord_ok << ", laboratory (A7) blocks finite at "
            << lab_blocks_ok;

        report.check_predicate(
            "b1c.alternative_observers_in_ergoregion",
            "inside the ergoregion the Eulerian observer of g_eff and the "
            "laboratory observer both remain admissible where the g_eff "
            "coordinate-time observer does not, and the laboratory blocks stay "
            "finite",
            VerificationLevel::A,
            inside > 0 && eulerian_ok == inside && lab_coord_ok == inside &&
                lab_blocks_ok == inside,
            "both alternatives admissible and the blocks finite throughout",
            got.str(), sw.seconds(),
            "demonstrates directly that the failure is observer-dependent: "
            "changing observer removes it");
    }

    // Boosted laboratory observers give valid but DIFFERENT blocks. Migrated
    // from the legacy tetrad decomposition to decompose_A7.
    {
        Stopwatch sw;
        const Mat4<double> g = alcubierre_at(0.2, 3.0);
        const Rank4<double> Z =
            covariant_constitutive(EffectiveMetric<double>(g), lab).Z;
        const ObserverDecomposition4<double> ref =
            decompose_A7(Z, eta, lab_rest.V);

        int ok_count = 0;
        double biggest_difference = 0.0;
        for (double beta : {0.1, 0.3, 0.6}) {
            const Observer<double> ob = laboratory_boosted_observer(lab, beta);
            if (!ob.admissible) continue;
            const ObserverDecomposition4<double> d = decompose_A7(Z, eta, ob.V);
            if (!std::isfinite(max_abs4(d.mu_inverse))) continue;
            ++ok_count;
            biggest_difference =
                std::fmax(biggest_difference,
                          std::fmax(diff4(d.mu_inverse, ref.mu_inverse),
                                    diff4(d.zeta, ref.zeta)));
        }

        report.check_predicate(
            "b1c.boosted_observer_blocks_differ",
            "a boosted inertial laboratory observer gives a valid but DIFFERENT "
            "(A7) decomposition of the same Z",
            VerificationLevel::A, ok_count == 3 && biggest_difference > 1e-6,
            "3 valid decompositions, all differing from the rest-frame one",
            std::to_string(ok_count) + " valid; largest block difference " +
                sci(biggest_difference),
            sw.seconds(),
            "eps, mu^-1 and zeta are observer-dependent by construction. If "
            "they had come out identical, the decomposition would be ignoring "
            "V. Supersedes b1c.boosted_observer_decomposition, which measured "
            "the same property through the quarantined tetrad path");
    }

    // =======================================================================
    // 8. What does and does not diverge at the ergosurface.
    //
    // WITHDRAWN: the old case asserted that cond(eps) must diverge. It cannot.
    // For this observer eps^{ab} = h^{ab} identically, so cond(eps) = 1
    // everywhere -- a tensor may scale uniformly towards infinity with its
    // condition number pinned at 1, and here eps does not scale at all. The
    // measured value was 1 both near and far, which is the correct answer to
    // the question the test asked and not evidence about the physics.
    //
    // What DOES diverge, under the same observer and the same conventions, is
    // the PERMEABILITY: the transverse eigenvalue of [mu^-1] equals -g_00
    // exactly and therefore reaches zero on the ergosurface, so mu = (mu^-1)^-1
    // diverges. Every diagnostic the correction task lists is computed and
    // exported; divergence is asserted only for the quantities predicted to
    // diverge.
    // =======================================================================
    {
        Stopwatch sw;
        std::ofstream os(out_dir + "/ergosurface_conditioning.csv", std::ios::binary);
        if (os) {
            os << "x,g00,f,vf,covariant_valid,effective_observer_admissible,"
                  "eps_cond,eps_eig_min,eps_eig_max,eps_det,eps_max_component,"
                  "mu_inv_cond,mu_inv_eig_min,mu_inv_eig_max,mu_inv_det,"
                  "mu_inv_max_component,zeta_norm,Z_norm\n";
            os.setf(std::ios::scientific);
            os.precision(10);
        }

        // Diagnostics at the sampled point closest to the ergosurface, and far
        // outside it.
        Spectrum eps_near, eps_far, mu_near, mu_far;
        double zeta_near = 0.0, zeta_far = 0.0;
        double Z_near = 0.0, Z_far = 0.0;
        double g00_near = 1e300;
        double worst_mu_eig_vs_g00 = 0.0, worst_eps_vs_projector = 0.0;
        double Z_min = 1e300, Z_max = 0.0;
        const Mat4<double> h = projector_upper(eta, V_rest);

        for (int i = 0; i <= 600; ++i) {
            const double x = 2.0 + 0.005 * i;   // sweep outward through r = 3
            const Mat4<double> g = alcubierre_at(2.0, x);
            const double f = shape(std::fabs(x));
            const CovariantConstitutive<double> cc =
                covariant_constitutive(EffectiveMetric<double>(g), lab);
            const EffectiveGeometryObserver<double> eobs =
                make_effective_geometry_observer(EffectiveMetric<double>(g), V_rest,
                                                 "effective_coordinate_time");
            if (!cc.ok) continue;

            const ObserverDecomposition4<double> d =
                decompose_A7(cc.Z, eta, lab_rest.V);
            const Spectrum se = spectrum_of(d.epsilon);
            const Spectrum sm = spectrum_of(d.mu_inverse);
            const double zn = max_abs4(d.zeta);
            const double zz = max_abs4x4(cc.Z);

            worst_eps_vs_projector =
                std::fmax(worst_eps_vs_projector, diff4(d.epsilon, h));
            // [mu^-1] has eigenvalues {1, -g_00, -g_00}, so the SMALLEST
            // magnitude is min(1, |g_00|) -- not |g_00| itself. Deep inside the
            // ergoregion |g_00| = |1 - v^2f^2| exceeds 1 (it reaches 2.97 at
            // x = 2 for v_s = 2) and the unit longitudinal eigenvalue is then
            // the smaller one.
            worst_mu_eig_vs_g00 = std::fmax(
                worst_mu_eig_vs_g00,
                std::fabs(sm.min_abs - std::fmin(1.0, std::fabs(g(0, 0)))));
            Z_min = std::fmin(Z_min, zz);
            Z_max = std::fmax(Z_max, zz);

            if (os)
                os << x << ',' << g(0, 0) << ',' << f << ',' << 2.0 * f << ','
                   << 1 << ',' << (eobs.admissible ? 1 : 0) << ','
                   << se.condition << ',' << se.min_abs << ',' << se.max_abs << ','
                   << se.determinant << ',' << se.max_component << ','
                   << sm.condition << ',' << sm.min_abs << ',' << sm.max_abs << ','
                   << sm.determinant << ',' << sm.max_component << ',' << zn << ','
                   << zz << '\n';

            if (std::fabs(g(0, 0)) < g00_near) {
                g00_near = std::fabs(g(0, 0));
                eps_near = se;
                mu_near = sm;
                zeta_near = zn;
                Z_near = zz;
            }
            if (std::fabs(x - 5.0) < 0.0026) {
                eps_far = se;
                mu_far = sm;
                zeta_far = zn;
                Z_far = zz;
            }
        }

        report.check_residual(
            "b1c.ergosurface.epsilon_is_the_projector",
            "eps^{ab} = h^{ab} at every point of the sweep, so its condition "
            "number is identically 1 and cannot diverge",
            VerificationLevel::A, "n/a",
            "v_s = 2, x from 2 to 5 through the ergosurface",
            worst_eps_vs_projector, 1e-13, sw.seconds(),
            "THE WITHDRAWN ASSERTION, restated as what is actually true. The "
            "old case demanded cond(eps) > 10 x its far value and measured "
            "1 versus 1. A uniformly scaling tensor keeps condition number 1; "
            "eps here does not even scale");

        {
            std::ostringstream got;
            got << "cond(eps) = " << eps_near.condition << " nearest the "
                << "ergosurface (|g_00| = " << g00_near << ") versus "
                << eps_far.condition << " at x = 5";
            report.check_predicate(
                "b1c.ergosurface.epsilon_conditioning_does_not_diverge",
                "cond(eps) is 1 both at the ergosurface and far from it, "
                "recorded so the withdrawn expectation is not silently dropped",
                VerificationLevel::A,
                std::fabs(eps_near.condition - 1.0) < 1e-9 &&
                    std::fabs(eps_far.condition - 1.0) < 1e-9,
                "cond(eps) = 1 at both ends", got.str(), sw.seconds(),
                "this case exists to keep the withdrawal visible in the gate "
                "output rather than only in the audit document");
        }

        report.check_residual(
            "b1c.ergosurface.mu_inverse_min_eigenvalue_closed_form",
            "the smallest eigenvalue magnitude of [mu^-1] equals min(1, |g_00|) "
            "exactly along the sweep",
            VerificationLevel::B, "n/a",
            "v_s = 2, x from 2 to 5 through the ergosurface",
            worst_mu_eig_vs_g00, 1e-13, sw.seconds(),
            "the spectrum is {1, -g_00, -g_00}, so the transverse pair reaches "
            "zero precisely on the ergosurface and mu = (mu^-1)^-1 diverges "
            "there. It is min(1, |g_00|) and NOT |g_00|: deep inside the "
            "ergoregion |g_00| = |1 - v^2f^2| exceeds 1 -- 2.97 at x = 2 for "
            "v_s = 2 -- and the unit longitudinal eigenvalue is then the smaller "
            "one. This is the quantity actually predicted to diverge for this "
            "observer");

        {
            std::ostringstream got;
            got << "nearest the ergosurface (|g_00| = " << g00_near
                << "): min|eig mu^-1| = " << mu_near.min_abs << ", cond(mu^-1) = "
                << mu_near.condition << ", |det mu^-1| = "
                << std::fabs(mu_near.determinant) << "; at x = 5: "
                << mu_far.min_abs << ", " << mu_far.condition << ", "
                << std::fabs(mu_far.determinant);
            report.check_predicate(
                "b1c.ergosurface.permeability_diverges",
                "the permeability diverges at the ergosurface: min|eig mu^-1| "
                "falls below 1e-4, cond(mu^-1) exceeds 1e3 and |det mu^-1| "
                "falls below 1e-8, while all three are of order 1 far away",
                VerificationLevel::B,
                mu_near.min_abs < 1e-4 && mu_near.condition > 1e3 &&
                    std::fabs(mu_near.determinant) < 1e-8 &&
                    std::fabs(mu_far.min_abs - 1.0) < 1e-3 &&
                    std::fabs(mu_far.condition - 1.0) < 1e-3,
                "mu^-1 degenerate near the surface, unit far away", got.str(),
                sw.seconds(),
                "asserted for a quantity predicted to diverge, under the same "
                "observer and the same conventions as the decomposition. "
                "Replaces b1c.eigenvalue_conditioning_diverges_at_ergosurface");
        }

        {
            std::ostringstream got;
            got << "max|Z| = " << Z_near << " nearest the ergosurface versus "
                << Z_far << " at x = 5, difference "
                << std::fabs(Z_near - Z_far) << "; range over the whole sweep ["
                << Z_min << ", " << Z_max << "]; max|zeta| = " << zeta_near
                << " versus " << zeta_far;
            report.check_predicate(
                "b1c.ergosurface.covariant_tensor_stays_bounded",
                "max|Z^{abcd}| is bounded over the whole sweep and is unchanged "
                "to 1e-4 across the ergosurface itself, and the magnetoelectric "
                "norm stays finite, while the permeability diverges",
                VerificationLevel::A,
                Z_max < 2.0 && Z_min > 0.1 &&
                    std::fabs(Z_near - Z_far) < 1e-4 && zeta_near < 2.0 &&
                    zeta_far < 2.0,
                "max|Z| in (0.1, 2) throughout, and near == far to 1e-4; "
                "max|zeta| finite",
                got.str(), sw.seconds(),
                "the two halves of the correction in one measurement: the "
                "covariant tensor does not notice the ergosurface, one "
                "observer's material parameters do. max|Z| is bounded but NOT "
                "constant over the sweep -- it rises to about 1.49 deep inside "
                "the ergoregion, where g_eff^{11} = 1 - v^2f^2 reaches -2.98 -- "
                "so boundedness is what is asserted, not constancy");
        }
    }

    // =======================================================================
    // 9. Velocity sensitivity, rebuilt. See b1c_velocity_checks.cpp.
    // =======================================================================
    run_velocity_checks(report, velocity);

    // =======================================================================
    // 10. Invariant versus non-invariant objectives.
    // =======================================================================
    {
        namespace M = warp_hypad::metrics;
        Stopwatch sw;
        const M::AlcubierreParameters<double> p{0.35, kR, kSigma};
        M::AlcubierreMetric<double> w;
        w.time_treatment = M::TimeTreatment::InertialTimeDependent;

        auto invariant_objective = [&](double t) {
            double acc = 0.0;
            const double dx = 0.02;
            for (double x = -40.0; x <= 40.0; x += dx) {
                const double v = w.evaluate(t, x, 0.0, 0.0, p)(0, 1);
                acc += v * v * dx;
            }
            return acc;
        };
        auto sensor_objective = [&](double t) {
            const double v = w.evaluate(t, 0.0, 0.0, 0.0, p)(0, 1);
            return v * v;
        };

        const double i0 = invariant_objective(0.0), i2 = invariant_objective(12.0);
        const double s0 = sensor_objective(0.0), s2 = sensor_objective(12.0);
        const double inv_rel = std::fabs(i2 - i0) / std::fmax(std::fabs(i0), 1e-300);
        const double sen_rel = std::fabs(s2 - s0) / std::fmax(std::fabs(s0), 1e-300);

        report.check_residual(
            "b1c.invariant_objective_is_time_translation_invariant",
            "an objective with a translation-invariant weight over a domain with "
            "remote boundaries is unchanged between t = 0 and t = 12",
            VerificationLevel::B, "n/a", "domain |x| <= 40, bubble moves 4.2",
            inv_rel, 1e-6, sw.seconds(),
            "this is the condition under which optimizing at t = 0 is "
            "legitimate -- it is a property of the objective, not a general "
            "fact");
        report.check_predicate(
            "b1c.laboratory_objective_is_not_invariant",
            "an objective read by a FIXED laboratory sensor at x = 0 changes "
            "substantially between t = 0 and t = 12",
            VerificationLevel::B, sen_rel > 0.5,
            "sensor reading changes by more than 50%",
            "relative change " + std::to_string(sen_rel), sw.seconds(),
            "so 'optimize at t = 0' is NOT generally valid. Program B's "
            "measurement and hardware objectives default to fixed laboratory "
            "coordinates for exactly this reason");
    }

    // =======================================================================
    // 11. Corrected status map export, on the production decomposition.
    // =======================================================================
    {
        std::ofstream os(out_dir + "/corrected_status_map.csv", std::ios::binary);
        if (os) {
            os << "v_s,x,g00,lorentzian,covariant_valid,effective_observer_ok,"
                  "eulerian_ok,eps_max,mu_inv_max,zeta_max,status\n";
            for (int vi = 1; vi <= 30; ++vi) {
                const double vs = 0.1 * vi;
                for (int i = -160; i <= 160; i += 2) {
                    const double x = 0.05 * i;
                    const Mat4<double> g = alcubierre_at(vs, x);
                    const EffectiveMetric<double> g_eff(g);
                    const CovariantConstitutive<double> cc =
                        covariant_constitutive(g_eff, lab);
                    const EffectiveGeometryObserver<double> eo =
                        make_effective_geometry_observer(g_eff, V_rest,
                                                         "effective_coordinate_time");
                    const auto sig = warp_hypad::geometry::analyse_signature(g);
                    const MetricMediumReport rep =
                        analyse_metric_medium(g_eff, lab, lab_rest, eo);
                    os << vs << ',' << x << ',' << g(0, 0) << ','
                       << (sig.lorentzian ? 1 : 0) << ',' << (cc.ok ? 1 : 0) << ','
                       << (eo.admissible ? 1 : 0) << ','
                       << (eulerian_observer(g).admissible ? 1 : 0) << ','
                       << rep.epsilon_max_component << ','
                       << rep.mu_inverse_max_component << ','
                       << rep.zeta_max_component << ',' << to_string(rep.status)
                       << '\n';
                }
            }
        }
    }

    // The five statuses must be distinguishable, not collapsed.
    {
        Stopwatch sw;
        int distinct = 0;
        std::ostringstream got;

        auto status_of = [&](const Mat4<double>& g_eff_raw) {
            const EffectiveMetric<double> g_eff(g_eff_raw);
            const EffectiveGeometryObserver<double> eo =
                make_effective_geometry_observer(g_eff, V_rest,
                                                 "effective_coordinate_time");
            return analyse_metric_medium(g_eff, lab, lab_rest, eo).status;
        };

        if (status_of(alcubierre_at(0.2, 3.0)) ==
            MetricMediumStatus::CovariantMetricMediumValid) {
            ++distinct;
            got << "valid ";
        }
        if (status_of(alcubierre_at(2.0, 0.0)) ==
            MetricMediumStatus::SelectedObserverDecompositionInvalid) {
            ++distinct;
            got << "observer-invalid ";
        }
        if (status_of(Mat4<double>::identity()) ==
            MetricMediumStatus::MetricDegenerateOrNonLorentzian) {
            ++distinct;
            got << "metric-degenerate ";
        }

        report.check_predicate(
            "b1c.statuses_are_distinguishable",
            "the valid, observer-invalid and metric-degenerate statuses are each "
            "produced by an appropriate input and are not collapsed together",
            VerificationLevel::A, distinct == 3, "3 distinct statuses observed",
            std::to_string(distinct) + " of 3: " + got.str(), sw.seconds(),
            "the original B1 error was reporting all of these as 'no equivalent "
            "medium'. Now computed through analyse_metric_medium(), which takes "
            "the laboratory and effective-geometry observers as separately typed "
            "arguments so the two questions cannot be conflated");
    }

#ifdef WARP_HYPAD_HAVE_OTILIB
    // OTI through the covariant tensor.
    {
        namespace sens = warp_hypad::sensitivity;
        using sens::OTIScalar;
        Stopwatch sw;

        const double vs = 0.2, x = 3.0;
        const sens::SeededParameters sp = sens::seed_all(vs, kR, kSigma);
        const AlcubierreParameters<OTIScalar> p{sp.velocity, sp.radius, sp.sigma};
        AlcubierreMetric<OTIScalar> w;
        w.time_treatment = warp_hypad::metrics::TimeTreatment::InertialTimeDependent;

        Mat4<OTIScalar> eta_o;
        eta_o(0, 0) = OTIScalar(-1.0);
        for (int i = 1; i < 4; ++i) eta_o(i, i) = OTIScalar(1.0);

        const CovariantConstitutive<OTIScalar> cc = covariant_constitutive(
            EffectiveMetric<OTIScalar>(w.evaluate(0.0, x, 0.0, 0.0, p)),
            LaboratoryMetric<OTIScalar>(eta_o));

        auto Z0101_double = [&](double v, double rr, double ss) {
            AlcubierreMetric<double> wd;
            wd.time_treatment = warp_hypad::metrics::TimeTreatment::InertialTimeDependent;
            const AlcubierreParameters<double> q{v, rr, ss};
            return covariant_constitutive(
                       EffectiveMetric<double>(wd.evaluate(0.0, x, 0.0, 0.0, q)), lab)
                .Z(0, 1, 0, 1);
        };

        const double base[3] = {vs, kR, kSigma};
        const double steps[4] = {1e-3, 1e-4, 1e-5, 1e-6};
        double worst = cc.ok ? 0.0 : 1.0;
        for (int d = 0; d < 3 && cc.ok; ++d) {
            const double oti = sens::extract_gradient(cc.Z(0, 1, 0, 1), d);
            double best = 1e300;
            for (double hh : steps) {
                double lo[3] = {base[0], base[1], base[2]};
                double hi[3] = {base[0], base[1], base[2]};
                lo[d] -= hh; hi[d] += hh;
                const double fd = (Z0101_double(hi[0], hi[1], hi[2]) -
                                   Z0101_double(lo[0], lo[1], lo[2])) / (2 * hh);
                best = std::fmin(best, std::fabs(fd - oti) /
                                           std::fmax(std::fabs(oti), 1e-9));
            }
            worst = std::fmax(worst, best);
        }

        report.check_residual(
            "b1c.oti_through_covariant_tensor",
            "d Z^{0101} / d(v, R, sigma) matches centred finite differences",
            VerificationLevel::D, "n/a", "v=0.2, R=3, sigma=2 at the wall",
            worst, 1e-6, sw.seconds(),
            "the covariant tensor is templated, so parameter channels flow "
            "through it unchanged");
    }
#else
    {
        Stopwatch sw;
        report.check_skipped(
            "b1c.oti_through_covariant_tensor",
            "HyPAD derivatives through the covariant constitutive tensor",
            "the OTIlib backend is disabled", sw.seconds());
    }
#endif
}

}  // namespace checks
}  // namespace program_b

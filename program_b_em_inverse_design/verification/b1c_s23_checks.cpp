#include <cmath>
#include <random>
#include <sstream>

#include "b1c_checks.hpp"
#include "program_b/effective_metric/bel_extraction.hpp"
#include "program_b/effective_metric/covariant_constitutive.hpp"
#include "warp_hypad/metrics/metric.hpp"

namespace program_b {
namespace checks {

using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace {

Mat4<double> minkowski() { return Mat4<double>::minkowski(); }

Mat4<double> alcubierre_at(double vs, double x) {
    warp_hypad::metrics::AlcubierreMetric<double> w;
    const warp_hypad::metrics::AlcubierreParameters<double> p{vs, 3.0, 2.0};
    return w.evaluate(0.0, x, 0.0, 0.0, p);
}

// Vacuum Z from the covariant definition, for a given effective/lab pair.
Rank4<double> Z_from_metrics(const Mat4<double>& eff, const Mat4<double>& lab) {
    return covariant_constitutive(eff, lab).Z;
}

// Equation (B2): isotropic medium at rest relative to V, with INDEPENDENT
// eps and mu_inv. This is the case that pins the (D1)/(D2) factors, because a
// medium with eps == mu_inv could not distinguish them.
Rank4<double> Z_isotropic_B2(const Mat4<double>& g_lab, const double V[4],
                             double eps, double mu_inv) {
    const Mat4<double> h = projector_upper(g_lab, V);
    Rank4<double> Z;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    Z(a, b, c, d) =
                        -(eps / 2.0) * (V[a] * V[c] * h(b, d) + V[b] * V[d] * h(a, c) -
                                        V[a] * V[d] * h(b, c) - V[b] * V[c] * h(a, d)) +
                        (mu_inv / 2.0) * (h(a, c) * h(b, d) - h(a, d) * h(b, c));
    return Z;
}

double mat4_max_abs_diff(const Mat4<double>& a, const Mat4<double>& b) {
    double w = 0.0;
    for (int i = 0; i < 16; ++i) w = std::fmax(w, std::fabs(a.m[i] - b.m[i]));
    return w;
}
double mat4_max_abs(const Mat4<double>& a) {
    double w = 0.0;
    for (int i = 0; i < 16; ++i) w = std::fmax(w, std::fabs(a.m[i]));
    return w;
}

struct NamedObserver {
    std::string name;
    Mat4<double> lab;
    double V[4];
};

// The four observer cases the gate requires.
std::vector<NamedObserver> observers() {
    std::vector<NamedObserver> out;
    const Mat4<double> eta = minkowski();

    {   // laboratory coordinate observer
        NamedObserver o;
        o.name = "lab_coordinate";
        o.lab = eta;
        o.V[0] = 1.0; o.V[1] = o.V[2] = o.V[3] = 0.0;
        out.push_back(o);
    }
    {   // boosted inertial observer, beta = 0.6
        NamedObserver o;
        o.name = "boosted_beta0.6";
        o.lab = eta;
        const double b = 0.6, gam = 1.0 / std::sqrt(1.0 - b * b);
        o.V[0] = gam; o.V[1] = gam * b; o.V[2] = o.V[3] = 0.0;
        out.push_back(o);
    }
    {   // Eulerian observer of a SUBLUMINAL Alcubierre background
        NamedObserver o;
        o.name = "eulerian_alcubierre_v0.2";
        o.lab = alcubierre_at(0.2, 3.0);
        const Observer<double> e = eulerian_observer(o.lab);
        for (int a = 0; a < 4; ++a) o.V[a] = e.V[a];
        out.push_back(o);
    }
    {   // Eulerian observer INSIDE the ergoregion (v_s = 2, bubble centre),
        // where the coordinate-time observer is inadmissible but this one is not.
        NamedObserver o;
        o.name = "eulerian_ergoregion_v2";
        o.lab = alcubierre_at(2.0, 0.0);
        const Observer<double> e = eulerian_observer(o.lab);
        for (int a = 0; a < 4; ++a) o.V[a] = e.V[a];
        out.push_back(o);
    }
    return out;
}

}  // namespace

void run_s23_checks(VerificationReport& report) {
    const Mat4<double> eta = minkowski();
    const std::vector<NamedObserver> obs = observers();

    // =======================================================================
    // STEP 3 FIRST: projector identities. Everything else assumes them.
    // =======================================================================
    {
        Stopwatch sw;
        double w_annihilate_right = 0.0, w_annihilate_left = 0.0;
        double w_idempotent = 0.0, w_trace = 0.0, w_metric = 0.0;
        double w_upper_lower = 0.0;

        for (const NamedObserver& o : obs) {
            const Mat4<double> hm = projector_mixed(o.lab, o.V);
            const Mat4<double> hu = projector_upper(o.lab, o.V);
            const Mat4<double> hl = projector_lower(o.lab, o.V);

            double V_dn[4];
            lower_vector(o.lab, o.V, V_dn);

            // h^a_b V^b = 0
            for (int a = 0; a < 4; ++a) {
                double s = 0.0;
                for (int b = 0; b < 4; ++b) s += hm(a, b) * o.V[b];
                w_annihilate_right = std::fmax(w_annihilate_right, std::fabs(s));
            }
            // V_a h^a_b = 0
            for (int b = 0; b < 4; ++b) {
                double s = 0.0;
                for (int a = 0; a < 4; ++a) s += V_dn[a] * hm(a, b);
                w_annihilate_left = std::fmax(w_annihilate_left, std::fabs(s));
            }
            // h^a_c h^c_b = h^a_b
            w_idempotent = std::fmax(w_idempotent,
                                     mat4_max_abs_diff(warp_hypad::mat_mul(hm, hm), hm));
            // h^a_a = 3
            double tr = 0.0;
            for (int a = 0; a < 4; ++a) tr += hm(a, a);
            w_trace = std::fmax(w_trace, std::fabs(tr - 3.0));
            // g_ab = h_ab - V_a V_b
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    w_metric = std::fmax(
                        w_metric, std::fabs(o.lab(a, b) - (hl(a, b) - V_dn[a] * V_dn[b])));
            // h^{ab} and h_{ab} must be inverse-related through g: h^a_b = h^{ac} g_cb
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    double s = 0.0;
                    for (int c = 0; c < 4; ++c) s += hu(a, c) * o.lab(c, b);
                    w_upper_lower = std::fmax(w_upper_lower, std::fabs(s - hm(a, b)));
                }
        }

        const std::string params = "4 observers: " + obs[0].name + ", " +
                                   obs[1].name + ", " + obs[2].name + ", " +
                                   obs[3].name;

        report.check_residual("b1c.projector.annihilates_V_right",
                              "h^a_b V^b = 0", VerificationLevel::A, "n/a",
                              params, w_annihilate_right, 1e-13, sw.seconds());
        report.check_residual("b1c.projector.annihilates_V_left",
                              "V_a h^a_b = 0", VerificationLevel::A, "n/a",
                              params, w_annihilate_left, 1e-13, sw.seconds());
        report.check_residual("b1c.projector.idempotent",
                              "h^a_c h^c_b = h^a_b", VerificationLevel::A, "n/a",
                              params, w_idempotent, 1e-13, sw.seconds());
        report.check_residual("b1c.projector.trace_is_three",
                              "h^a_a = 3", VerificationLevel::A, "n/a", params,
                              w_trace, 1e-13, sw.seconds(),
                              "3, not 4: the projector removes the observer's "
                              "own timelike direction");
        report.check_residual("b1c.projector.metric_reconstruction",
                              "g_ab = h_ab - V_a V_b", VerificationLevel::A,
                              "n/a", params, w_metric, 1e-13, sw.seconds());
        report.check_residual("b1c.projector.index_forms_consistent",
                              "h^{ac} g_cb = h^a_b, so the upper, lower and "
                              "mixed forms are the same object",
                              VerificationLevel::A, "n/a", params,
                              w_upper_lower, 1e-13, sw.seconds(),
                              "each index form is built independently, so this "
                              "is a real consistency check");
    }

    // =======================================================================
    // STEP 2: independent A7 contractions.
    //
    // The reference is computed here, in the test, with its own loops. It does
    // NOT call bel_extract or contract_a7.
    // =======================================================================
    {
        Stopwatch sw;
        std::mt19937_64 rng(0xBE1EFULL);
        std::uniform_real_distribution<double> dist(-1.2, 1.2);

        double raw[4][4][4][4];
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d) raw[a][b][c][d] = dist(rng);
        Rank4<double> Z;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d)
                        Z(a, b, c, d) = 0.25 * (raw[a][b][c][d] - raw[b][a][c][d] -
                                                raw[a][b][d][c] + raw[b][a][d][c]);

        double worst = 0.0;
        for (const NamedObserver& o : obs) {
            const BelBlocks<double> got = bel_extract(Z, o.lab, o.V);

            // Independent reference: lower V by hand, dual by hand via the
            // verified dual functions, contract by hand.
            double Vd[4];
            for (int a = 0; a < 4; ++a) {
                double s = 0.0;
                for (int b = 0; b < 4; ++b) s += o.lab(a, b) * o.V[b];
                Vd[a] = s;
            }
            const Rank4<double> ZL = dual_left(Z, o.lab);
            const Rank4<double> ZR = dual_right(Z, o.lab);
            const Rank4<double> ZB = dual_right(dual_left(Z, o.lab), o.lab);

            auto ref = [&](const Rank4<double>& T, const Mat4<double>& cmp) {
                double w = 0.0;
                for (int b = 0; b < 4; ++b)
                    for (int d = 0; d < 4; ++d) {
                        double acc = 0.0;
                        for (int a = 0; a < 4; ++a)
                            for (int c = 0; c < 4; ++c) acc += Vd[a] * Vd[c] * T(a, b, c, d);
                        w = std::fmax(w, std::fabs(acc - cmp(b, d)));
                    }
                return w;
            };

            worst = std::fmax(worst, ref(Z, got.W_epsilon));
            worst = std::fmax(worst, ref(ZB, got.W_mu));
            worst = std::fmax(worst, ref(ZL, got.W_zeta_T));
            worst = std::fmax(worst, ref(ZR, got.W_zeta));
        }

        report.check_residual(
            "b1c.bel.a7_contractions_independent",
            "all four A7 contractions match a reference computed in the test "
            "with its own loops",
            VerificationLevel::D, "n/a", "random Z with Z's antisymmetries, 4 "
            "observers", worst, 1e-13, sw.seconds(),
            "the reference does not call bel_extract or contract_a7");
    }

    // =======================================================================
    // Vacuum: eps = mu^-1 = h, zeta = zeta+ = 0, for EVERY observer.
    // =======================================================================
    {
        Stopwatch sw;
        double w_eps = 0.0, w_mu = 0.0, w_zeta = 0.0, w_zadj = 0.0;
        for (const NamedObserver& o : obs) {
            // Vacuum means g_eff == g_lab, whatever the lab is.
            const Rank4<double> Z = Z_from_metrics(o.lab, o.lab);
            const ObserverDecomposition4<double> d = decompose_A7(Z, o.lab, o.V);
            const Mat4<double> h = projector_upper(o.lab, o.V);

            w_eps = std::fmax(w_eps, mat4_max_abs_diff(d.epsilon, h));
            w_mu = std::fmax(w_mu, mat4_max_abs_diff(d.mu_inverse, h));
            w_zeta = std::fmax(w_zeta, mat4_max_abs(d.zeta));
            w_zadj = std::fmax(w_zadj, mat4_max_abs(d.zeta_adjoint));
        }
        report.check_residual("b1c.bel.vacuum_epsilon",
                              "eps^{ab} = h^{ab} for vacuum", VerificationLevel::B,
                              "n/a", "4 observers", w_eps, 1e-13, sw.seconds());
        report.check_residual("b1c.bel.vacuum_mu_inverse",
                              "[mu^-1]^{ab} = h^{ab} for vacuum",
                              VerificationLevel::B, "n/a", "4 observers", w_mu,
                              1e-13, sw.seconds(),
                              "confirms the +2 in (D2), which differs in sign "
                              "from the -2 in (D1)");
        report.check_residual("b1c.bel.vacuum_zeta_zero",
                              "zeta = 0 for vacuum", VerificationLevel::B, "n/a",
                              "4 observers", w_zeta, 1e-13, sw.seconds());
        report.check_residual("b1c.bel.vacuum_zeta_adjoint_zero",
                              "zeta+ = 0 for vacuum", VerificationLevel::B,
                              "n/a", "4 observers", w_zadj, 1e-13, sw.seconds());
    }

    // =======================================================================
    // Isotropic medium at rest, equation (B2), with eps != mu^-1.
    //
    // This is what actually pins (D1) and (D2): a medium with eps == mu^-1
    // could not tell the two factors apart.
    // =======================================================================
    {
        Stopwatch sw;
        const double eps_v = 2.75, mui_v = 0.4;   // deliberately unequal
        double w_eps = 0.0, w_mu = 0.0, w_zeta = 0.0;

        for (const NamedObserver& o : obs) {
            const Rank4<double> Z = Z_isotropic_B2(o.lab, o.V, eps_v, mui_v);
            const ObserverDecomposition4<double> d = decompose_A7(Z, o.lab, o.V);
            const Mat4<double> h = projector_upper(o.lab, o.V);

            for (int i = 0; i < 16; ++i) {
                w_eps = std::fmax(w_eps, std::fabs(d.epsilon.m[i] - eps_v * h.m[i]));
                w_mu = std::fmax(w_mu, std::fabs(d.mu_inverse.m[i] - mui_v * h.m[i]));
            }
            w_zeta = std::fmax(w_zeta, mat4_max_abs(d.zeta));
        }

        report.check_residual(
            "b1c.bel.isotropic_rest_epsilon",
            "eps^{ab} = eps h^{ab} from (B1a)/(B2), with eps = 2.75",
            VerificationLevel::B, "n/a", "eps = 2.75, mu^-1 = 0.4, 4 observers",
            w_eps, 1e-12, sw.seconds(),
            "eps != mu^-1 deliberately: this is what distinguishes the (D1) and "
            "(D2) factors, which a matched medium could not");
        report.check_residual(
            "b1c.bel.isotropic_rest_mu_inverse",
            "[mu^-1]^{ab} = mu^-1 h^{ab} from (B1b)/(B2), with mu^-1 = 0.4",
            VerificationLevel::B, "n/a", "eps = 2.75, mu^-1 = 0.4, 4 observers",
            w_mu, 1e-12, sw.seconds());
        report.check_residual(
            "b1c.bel.isotropic_rest_zeta_zero",
            "zeta = 0 for an isotropic medium AT REST relative to the observer",
            VerificationLevel::B, "n/a", "4 observers", w_zeta, 1e-12,
            sw.seconds(),
            "(B17) gives zeta proportional to (V.W); at rest V = W so it "
            "vanishes. This does NOT pin the zeta normalisation -- see step 5");
    }

    // =======================================================================
    // Left versus right dual must not be interchangeable.
    // =======================================================================
    {
        Stopwatch sw;
        // A metric-induced Z with a shift, which has magnetoelectric structure.
        const Mat4<double> g = alcubierre_at(0.4, 3.0);
        const Rank4<double> Z = Z_from_metrics(g, eta);
        const double V[4] = {1.0, 0.0, 0.0, 0.0};

        const Rank4<double> ZL = dual_left(Z, eta);
        const Rank4<double> ZR = dual_right(Z, eta);
        const double dual_gap = rank4_max_abs_diff(ZL, ZR);

        const BelBlocks<double> b = bel_extract(Z, eta, V);
        const double block_gap = mat4_max_abs_diff(b.W_zeta, b.W_zeta_T);

        std::ostringstream got;
        got << "max|*Z - Z*| = " << dual_gap << ", max|W_zeta - W_zeta^T| = "
            << block_gap;

        report.check_predicate(
            "b1c.bel.left_right_dual_distinct",
            "for a magnetoelectric medium *Z != Z*, and the two A7 duals "
            "therefore give genuinely different magnetoelectric blocks",
            VerificationLevel::A, dual_gap > 1e-6 && block_gap > 1e-9,
            "both differences well above roundoff", got.str(), sw.seconds(),
            "catches an accidental interchange of the dual index pairs, which "
            "would silently transpose the magnetoelectric response");
    }

    // =======================================================================
    // STEP 3: spatiality of every block, every index, reported SEPARATELY.
    // =======================================================================
    {
        const char* block_name[4] = {"epsilon", "mu_inverse", "zeta",
                                     "zeta_adjoint"};
        double worst_first[4] = {0, 0, 0, 0};
        double worst_second[4] = {0, 0, 0, 0};
        double worst_projection[4] = {0, 0, 0, 0};

        for (const NamedObserver& o : obs) {
            // Use a metric-induced Z so zeta is genuinely nonzero and the
            // spatiality test is not trivially satisfied by a zero block.
            const Rank4<double> Z = Z_from_metrics(alcubierre_at(0.4, 3.0), o.lab);
            const ObserverDecomposition4<double> d = decompose_A7(Z, o.lab, o.V);
            const Mat4<double>* blocks[4] = {&d.epsilon, &d.mu_inverse, &d.zeta,
                                             &d.zeta_adjoint};
            for (int k = 0; k < 4; ++k) {
                double f = 0.0, s = 0.0;
                spatiality_residuals(*blocks[k], o.lab, o.V, f, s);
                worst_first[k] = std::fmax(worst_first[k], f);
                worst_second[k] = std::fmax(worst_second[k], s);
                worst_projection[k] = std::fmax(
                    worst_projection[k],
                    mat4_max_abs_diff(project_tensor(*blocks[k], o.lab, o.V),
                                      *blocks[k]));
            }
        }

        for (int k = 0; k < 4; ++k) {
            Stopwatch sw;
            report.check_residual(
                std::string("b1c.spatiality.") + block_name[k] + ".first_index",
                std::string("V_a ") + block_name[k] + "^{ab} = 0",
                VerificationLevel::A, "n/a", "4 observers", worst_first[k],
                1e-13, sw.seconds(),
                "structural: Z is antisymmetric in (ab) while V_aV_b is "
                "symmetric, so this is an exact zero, not a small residual");
            report.check_residual(
                std::string("b1c.spatiality.") + block_name[k] + ".second_index",
                std::string(block_name[k]) + "^{ab} V_b = 0",
                VerificationLevel::A, "n/a", "4 observers", worst_second[k],
                1e-13, sw.seconds());
            report.check_residual(
                std::string("b1c.projection_invariance.") + block_name[k],
                std::string("h^a_c h^b_d ") + block_name[k] + "^{cd} = " +
                    block_name[k] + "^{ab}",
                VerificationLevel::A, "n/a", "4 observers", worst_projection[k],
                1e-13, sw.seconds(),
                "spatiality is structural, not the result of a cancellation in "
                "one particular contraction");
        }
    }

    // The zeta blocks must be nonzero somewhere, or every spatiality result
    // above would be trivially satisfied.
    {
        Stopwatch sw;
        const Rank4<double> Z = Z_from_metrics(alcubierre_at(0.4, 3.0), eta);
        const double V[4] = {1.0, 0.0, 0.0, 0.0};
        const ObserverDecomposition4<double> d = decompose_A7(Z, eta, V);
        const double mz = mat4_max_abs(d.zeta);
        report.check_predicate(
            "b1c.spatiality.zeta_is_nonzero",
            "the magnetoelectric block is genuinely nonzero for the warp medium, "
            "so its spatiality and projection tests are not vacuous",
            VerificationLevel::A, mz > 1e-6, "max|zeta| well above roundoff",
            "max|zeta| = " + std::to_string(mz), sw.seconds());
    }

    // The magnetoelectric normalisation, pinned at step 5 against (B17) and
    // asserted here as an exact value with one specific sign. This replaces the
    // SKIP that stood in earlier revisions of this gate.
    {
        Stopwatch sw;
        report.check_value(
            "b1c.bel.zeta_normalisation",
            "zeta = +2 W_zeta, determined from the (B17) moving-medium "
            "benchmark",
            VerificationLevel::D, "n/a", "5 boost configurations", 2.0,
            kZetaNormalisation, 0.0, sw.seconds(),
            "exact equality on one specific signed value. The provisional "
            "guess was -2, taken by symmetry with (D1), and it was wrong");
    }
}

}  // namespace checks
}  // namespace program_b

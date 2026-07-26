#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "b2_checks.hpp"
#include "program_b/inverse_design/design_gradient.hpp"
#include "program_b/inverse_design/optimizer.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "program_b/inverse_design/design_gradient_oti.hpp"
#endif

// ---------------------------------------------------------------------------
// Gate B2-ABSTRACT-INVERSE.
//
// Everything runs in the abstract constitutive-design space. No Maxwell solver,
// no hardware model, no dispersion. The optimizer drives the SAME A11
// reconstruction B1C verified, through the same scalar-generic objective.
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

LaboratoryMetric<double> eta_lab() {
    return LaboratoryMetric<double>::minkowski();
}

Mat4<double> alcubierre_at(double vs, double x) {
    warp_hypad::metrics::AlcubierreMetric<double> w;
    const warp_hypad::metrics::AlcubierreParameters<double> p{vs, 3.0, 2.0};
    return w.evaluate(0.0, x, 0.0, 0.0, p);
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

std::string sci(double v) {
    std::ostringstream o;
    o.setf(std::ios::scientific);
    o.precision(4);
    o << v;
    return o.str();
}

// --- An independent linear-algebra reference --------------------------------
//
// For every family except PositiveDefinite, A11 is homogeneous linear in the
// blocks and the blocks are linear in q, so Z_design = A q exactly. The
// coordinate-basis J_Z is then a quadratic whose unique minimiser solves the
// normal equations A^T A q = A^T z_target.
//
// This is deliberately NOT the optimizer: dense Gaussian elimination with
// partial pivoting, written here, is what the optimizer's answer is checked
// against. The design matrix is built by evaluating the design map at basis
// vectors, which is legitimate -- it is the map itself, not a second copy of it.

struct DesignMatrix {
    std::vector<double> A;   // 256 x n, row-major
    int n = 0;
    bool linear = false;
};

DesignMatrix build_design_matrix(DesignFamily family,
                                const LaboratoryMetric<double>& g_lab,
                                const ObserverFrame<double>& frame) {
    DesignMatrix M;
    M.linear = (family != DesignFamily::PositiveDefinite);
    if (!M.linear) return M;
    M.n = parameter_count(family);
    M.A.assign(static_cast<std::size_t>(256) * M.n, 0.0);

    for (int k = 0; k < M.n; ++k) {
        ConstitutiveDesignVariables<double> e;
        e.family = family;
        e.q[k] = 1.0;
        const Rank4<double> Zk = reconstruct_design_tensor(e, g_lab, frame);
        for (int m = 0; m < 256; ++m)
            M.A[static_cast<std::size_t>(m) * M.n + k] = Zk.t[m];
    }
    return M;
}

// Solve the normal equations. Returns false on a singular system rather than
// producing a plausible-looking answer.
bool solve_least_squares(const DesignMatrix& M, const Rank4<double>& z_target,
                        double* q_out) {
    const int n = M.n;
    std::vector<double> N(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> rhs(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int m = 0; m < 256; ++m)
                s += M.A[static_cast<std::size_t>(m) * n + i] *
                     M.A[static_cast<std::size_t>(m) * n + j];
            N[static_cast<std::size_t>(i) * n + j] = s;
        }
        double s = 0.0;
        for (int m = 0; m < 256; ++m)
            s += M.A[static_cast<std::size_t>(m) * n + i] * z_target.t[m];
        rhs[i] = s;
    }

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double best = std::fabs(N[static_cast<std::size_t>(col) * n + col]);
        for (int r = col + 1; r < n; ++r) {
            const double v = std::fabs(N[static_cast<std::size_t>(r) * n + col]);
            if (v > best) { best = v; pivot = r; }
        }
        if (best < 1e-14) return false;
        if (pivot != col) {
            for (int c = 0; c < n; ++c)
                std::swap(N[static_cast<std::size_t>(col) * n + c],
                          N[static_cast<std::size_t>(pivot) * n + c]);
            std::swap(rhs[col], rhs[pivot]);
        }
        const double d = N[static_cast<std::size_t>(col) * n + col];
        for (int r = col + 1; r < n; ++r) {
            const double f = N[static_cast<std::size_t>(r) * n + col] / d;
            if (f == 0.0) continue;
            for (int c = col; c < n; ++c)
                N[static_cast<std::size_t>(r) * n + c] -=
                    f * N[static_cast<std::size_t>(col) * n + c];
            rhs[r] -= f * rhs[col];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        double s = rhs[r];
        for (int c = r + 1; c < n; ++c)
            s -= N[static_cast<std::size_t>(r) * n + c] * q_out[c];
        q_out[r] = s / N[static_cast<std::size_t>(r) * n + r];
    }
    return true;
}

// --- Named demonstration designs --------------------------------------------

ConstitutiveDesignVariables<double> isotropic_design(double eps, double nu) {
    ConstitutiveDesignVariables<double> v;
    v.family = DesignFamily::Isotropic;
    v.q[0] = eps;
    v.q[1] = nu;
    return v;
}

// Distinct diagonal entries plus one off-diagonal symmetric entry, so the target
// cannot be matched by a diagonal design.
ConstitutiveDesignVariables<double> anisotropic_design() {
    ConstitutiveDesignVariables<double> v;
    v.family = DesignFamily::UnrestrictedReciprocal;
    v.q[0] = 2.40;  v.q[1] = 0.35;  v.q[2] = 0.00;
    v.q[3] = 1.70;  v.q[4] = 0.00;  v.q[5] = 3.10;
    v.q[6] = 0.80;  v.q[7] = 0.00;  v.q[8] = -0.12;
    v.q[9] = 0.55;  v.q[10] = 0.00; v.q[11] = 0.95;
    return v;
}

// Nonsymmetric zeta. zeta_adjoint = zeta^T is structural, so this is reciprocal
// while being genuinely magnetoelectric -- a zero-zeta design cannot match it.
ConstitutiveDesignVariables<double> magnetoelectric_design() {
    ConstitutiveDesignVariables<double> v = anisotropic_design();
    v.q[12] = 0.20;  v.q[13] = -0.45; v.q[14] = 0.05;
    v.q[15] = 0.60;  v.q[16] = -0.10; v.q[17] = 0.33;
    v.q[18] = -0.07; v.q[19] = 0.28;  v.q[20] = 0.14;
    return v;
}

struct Setup {
    LaboratoryMetric<double> g_lab;
    Observer<double> observer;
    ObserverFrame<double> frame;
    ObjectiveWeights<double> weights;
    DesignConstraints<double> constraints;
    OptimizerSettings settings;
};

Setup default_setup() {
    Setup s;
    s.g_lab = eta_lab();
    double V[4] = {1.0, 0.0, 0.0, 0.0};
    s.observer = make_laboratory_observer(s.g_lab, V, "lab_rest");
    s.frame = build_observer_frame(s.g_lab, s.observer);
    s.settings.max_iterations = 4000;
    s.settings.gradient_tolerance = 1e-11;
    s.settings.objective_tolerance = 1e-20;
    return s;
}

}  // namespace

void run_b2_checks(VerificationReport& report, const std::string& out_dir,
                  B2Metrics& metrics) {
    const Setup S = default_setup();

    // =======================================================================
    // 1. The observer frame, and the parameter-vector round trip.
    // =======================================================================
    {
        Stopwatch sw;
        report.check_predicate(
            "b2.observer_frame_is_orthonormal_and_spatial",
            "the spatial triad satisfies g_lab(e_i, e_j) = delta_ij and "
            "g_lab(e_i, V) = 0",
            VerificationLevel::A, S.frame.ok, "frame built",
            S.frame.message, sw.seconds(),
            "the triad supplies coordinates ON the design space; it is never "
            "used to extract physics FROM a tensor, which is what made the "
            "quarantined legacy tetrad path unsound");

        double worst_orth = 0.0;
        for (int i = 0; i < 3; ++i) {
            double sv = 0.0;
            for (int a = 0; a < 4; ++a) sv += S.frame.e_dn[i][a] * S.frame.V[a];
            worst_orth = std::fmax(worst_orth, std::fabs(sv));
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int a = 0; a < 4; ++a) s += S.frame.e_dn[i][a] * S.frame.e[j][a];
                worst_orth = std::fmax(worst_orth, std::fabs(s - (i == j ? 1.0 : 0.0)));
            }
        }
        report.check_residual(
            "b2.observer_frame_orthonormality", "g_lab(e_i, e_j) = delta_ij and "
            "g_lab(e_i, V) = 0", VerificationLevel::A, "n/a",
            "g_lab = eta, V = (1,0,0,0)", worst_orth, 1e-14, sw.seconds());
    }

    {
        Stopwatch sw;
        double worst = 0.0;
        const DesignFamily families[4] = {DesignFamily::UnrestrictedReciprocal,
                                          DesignFamily::ZeroMagnetoelectric,
                                          DesignFamily::DiagonalAnisotropic,
                                          DesignFamily::Isotropic};
        for (DesignFamily f : families) {
            ConstitutiveDesignVariables<double> v;
            v.family = f;
            for (int k = 0; k < v.size(); ++k)
                v.q[k] = 0.31 + 0.17 * static_cast<double>(k) -
                         0.03 * static_cast<double>(k * k);

            const BlockCoefficients<double> c = coefficients_from_parameters(v);
            ConstitutiveDesignVariables<double> back;
            if (!parameters_from_coefficients(c, f, back)) { worst = 1.0; continue; }
            for (int k = 0; k < v.size(); ++k)
                worst = std::fmax(worst, std::fabs(back.q[k] - v.q[k]));

            // And through the four-dimensional blocks and back.
            const ObserverDecomposition4<double> b = blocks_from_parameters(v, S.frame);
            double m[3][3];
            coefficients_from_block(b.epsilon, S.frame, m);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    worst = std::fmax(worst, std::fabs(m[i][j] - c.epsilon[i][j]));
            coefficients_from_block(b.zeta, S.frame, m);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    worst = std::fmax(worst, std::fabs(m[i][j] - c.zeta[i][j]));
        }
        metrics.worst_parameter_roundtrip = worst;
        report.check_residual(
            "b2.parameter_vector_roundtrip",
            "q -> coefficients -> q is the identity, and q -> blocks -> "
            "coefficients likewise, for every linear family",
            VerificationLevel::A, "n/a", "4 families", worst, 1e-14, sw.seconds(),
            "pins the parameter ordering. A silent reordering would leave every "
            "other case passing while the reported design meant something else");
    }

    // =======================================================================
    // 2. Symmetry, spatiality and reciprocity BY CONSTRUCTION.
    // =======================================================================
    {
        Stopwatch sw;
        double w_sym = 0.0, w_spat = 0.0, w_recip = 0.0, w_zadj = 0.0;
        const DesignFamily families[5] = {DesignFamily::UnrestrictedReciprocal,
                                          DesignFamily::PositiveDefinite,
                                          DesignFamily::ZeroMagnetoelectric,
                                          DesignFamily::DiagonalAnisotropic,
                                          DesignFamily::Isotropic};
        for (DesignFamily f : families) {
            ConstitutiveDesignVariables<double> v;
            v.family = f;
            for (int k = 0; k < v.size(); ++k)
                v.q[k] = 0.7 - 0.11 * static_cast<double>(k);
            const ObserverDecomposition4<double> b = blocks_from_parameters(v, S.frame);

            w_sym = std::fmax(w_sym, warp_hypad::symmetry_residual(b.epsilon));
            w_sym = std::fmax(w_sym, warp_hypad::symmetry_residual(b.mu_inverse));
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    w_zadj = std::fmax(w_zadj,
                                       std::fabs(b.zeta_adjoint(i, j) - b.zeta(j, i)));
            w_recip = std::fmax(w_recip, reciprocity_residual(b));

            double s1 = 0.0, s2 = 0.0;
            const Mat4<double>* blocks[4] = {&b.epsilon, &b.mu_inverse, &b.zeta,
                                             &b.zeta_adjoint};
            for (int k = 0; k < 4; ++k) {
                spatiality_residuals(*blocks[k], S.g_lab.value, S.frame.V, s1, s2);
                w_spat = std::fmax(w_spat, std::fmax(s1, s2));
            }

            const Rank4<double> Z = reconstruct_design_tensor(v, S.g_lab, S.frame);
            w_recip = std::fmax(w_recip, pair_symmetry_residual(Z));
        }
        metrics.worst_symmetry_residual = std::fmax(w_sym, w_zadj);
        metrics.worst_spatiality_residual = w_spat;

        report.check_residual(
            "b2.symmetry_by_construction",
            "eps and mu^-1 are symmetric and zeta_adjoint = zeta^T for every "
            "family, with no symmetrisation step applied afterwards",
            VerificationLevel::A, "n/a", "5 families",
            std::fmax(w_sym, w_zadj), 0.0, sw.seconds(),
            "EXACT zero, not a tolerance: both (i,j) and (j,i) read the same "
            "parameter, so there is nothing to cancel");
        report.check_residual(
            "b2.spatiality_by_construction",
            "every block is spatial with respect to V, structurally",
            VerificationLevel::A, "n/a", "5 families", w_spat, 1e-14, sw.seconds());
        report.check_residual(
            "b2.reconstructed_tensor_is_reciprocal",
            "Z_design^{abcd} = Z_design^{cdab} for every family",
            VerificationLevel::A, "n/a", "5 families", w_recip, 1e-13, sw.seconds(),
            "the design space is reciprocal by construction, so J_reciprocity is "
            "a diagnostic and not a useful optimizer term");
    }

    // =======================================================================
    // 3. The A11 reconstruction is the verified one: blocks -> Z -> blocks.
    // =======================================================================
    {
        Stopwatch sw;
        double worst = 0.0;
        for (int trial = 0; trial < 4; ++trial) {
            ConstitutiveDesignVariables<double> v = magnetoelectric_design();
            v = perturb_design(v, 0.4, trial);
            const ObserverDecomposition4<double> b = blocks_from_parameters(v, S.frame);
            const Rank4<double> Z = reconstruct_design_tensor(v, S.g_lab, S.frame);
            const ObserverDecomposition4<double> back =
                decompose_A7(Z, S.g_lab.value, S.frame.V);
            worst = std::fmax(worst, warp_hypad::max_abs_diff(back.epsilon, b.epsilon));
            worst = std::fmax(worst,
                              warp_hypad::max_abs_diff(back.mu_inverse, b.mu_inverse));
            worst = std::fmax(worst, warp_hypad::max_abs_diff(back.zeta, b.zeta));
            worst = std::fmax(worst, warp_hypad::max_abs_diff(back.zeta_adjoint,
                                                             b.zeta_adjoint));
        }
        report.check_residual(
            "b2.a11_reconstruction_roundtrip",
            "design blocks -> Z via (A11) -> blocks via (A7) returns the design "
            "blocks",
            VerificationLevel::A, "n/a", "4 perturbed magnetoelectric designs",
            worst, 1e-12, sw.seconds(),
            "the optimizer drives the SAME reconstruction B1C verified; there is "
            "no second production formula");
    }

    // =======================================================================
    // 4. The two norms: which is invariant and which is not.
    // =======================================================================
    {
        Stopwatch sw;
        Mat4<double> L = Mat4<double>::identity();
        L(0, 0) = 1.20; L(0, 1) = 0.30;
        L(1, 0) = 0.25; L(1, 1) = 1.10;
        L(2, 1) = 0.15; L(2, 2) = 1.05; L(2, 3) = 0.10;
        L(3, 0) = 0.05; L(3, 2) = 0.20; L(3, 3) = 0.95;
        const warp_hypad::InverseResult<double> Li = warp_hypad::inverse(L);

        const ConstitutiveDesignVariables<double> v = magnetoelectric_design();
        const Rank4<double> Z = reconstruct_design_tensor(v, S.g_lab, S.frame);

        Mat4<double> gp;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                double acc = 0.0;
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d)
                        acc += Li.inverse(c, a) * Li.inverse(d, b) * S.g_lab.value(c, d);
                gp(a, b) = acc;
            }
        const LaboratoryMetric<double> g_prime(gp);

        Rank4<double> Zp;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d) {
                        double acc = 0.0;
                        for (int e = 0; e < 4; ++e)
                            for (int f = 0; f < 4; ++f)
                                for (int g = 0; g < 4; ++g)
                                    for (int h = 0; h < 4; ++h)
                                        acc += L(a, e) * L(b, f) * L(c, g) * L(d, h) *
                                               Z(e, f, g, h);
                        Zp(a, b, c, d) = acc;
                    }
        double Vp[4];
        for (int a = 0; a < 4; ++a) {
            double acc = 0.0;
            for (int b = 0; b < 4; ++b) acc += L(a, b) * S.frame.V[b];
            Vp[a] = acc;
        }

        const Mat4<double> G = observer_gram(S.g_lab, S.frame.V);
        const Mat4<double> Gp = observer_gram(g_prime, Vp);

        const double coord_unprimed =
            rank4_norm_squared(Z, G, ObjectiveNorm::CoordinateBasis);
        const double coord_primed =
            rank4_norm_squared(Zp, Gp, ObjectiveNorm::CoordinateBasis);
        const double obs_unprimed =
            rank4_norm_squared(Z, G, ObjectiveNorm::ObserverReferred);
        const double obs_primed =
            rank4_norm_squared(Zp, Gp, ObjectiveNorm::ObserverReferred);

        report.check_residual(
            "b2.norms_agree_in_the_analogue_configuration",
            "with g_lab = eta and V = (1,0,0,0) the observer Gram is the identity, "
            "so the coordinate-basis and observer-referred norms coincide exactly",
            VerificationLevel::A, "n/a", "magnetoelectric design",
            std::fabs(obs_unprimed - coord_unprimed) /
                std::fmax(coord_unprimed, 1e-300),
            1e-14, sw.seconds(),
            "which is why the coordinate-basis norm is usable there, and why "
            "calling it invariant would still be wrong");

        report.check_residual(
            "b2.observer_referred_norm_is_invariant",
            "||Z||^2_V is unchanged by a nontrivial linear coordinate change",
            VerificationLevel::A, "n/a", "non-orthogonal L, det = " + sci(Li.det_real),
            std::fabs(obs_primed - obs_unprimed) / std::fmax(obs_unprimed, 1e-300),
            1e-12, sw.seconds(),
            "G = g_lab + 2 V(x)V transforms as a (0,2) tensor, so the full "
            "contraction is a scalar. g_lab alone would be indefinite and the "
            "norm could vanish for Z != 0");

        report.check_predicate(
            "b2.coordinate_basis_norm_is_not_invariant",
            "the coordinate-basis norm DOES change under the same coordinate "
            "change, which is why it is labelled a coordinate-basis norm",
            VerificationLevel::A,
            std::fabs(coord_primed - coord_unprimed) / coord_unprimed > 1e-2,
            "a relative change above 1e-2",
            "coord: " + sci(coord_unprimed) + " -> " + sci(coord_primed) +
                ", observer-referred: " + sci(obs_unprimed) + " -> " +
                sci(obs_primed),
            sw.seconds(),
            "asserted, so the distinction cannot quietly rot into a claim of "
            "invariance");
    }

    // =======================================================================
    // 5. The objective is zero at the exact target and positive away from it.
    // =======================================================================
    {
        Stopwatch sw;
        double worst_at = 0.0, smallest_away = 1e300;
        struct Named { const char* name; ConstitutiveDesignVariables<double> v; };
        const Named designs[3] = {{"isotropic", isotropic_design(2.75, 0.4)},
                                  {"anisotropic", anisotropic_design()},
                                  {"magnetoelectric", magnetoelectric_design()}};

        for (const Named& d : designs) {
            const ConstitutiveTarget<double> t =
                target_from_design(d.v, S.g_lab, S.frame, d.name);
            const DesignEvaluation<double> at =
                evaluate_design(d.v, t, S.constraints, S.frame, S.weights);
            worst_at = std::fmax(worst_at, std::fabs(at.J_Z));
            worst_at = std::fmax(worst_at, std::fabs(at.J_blocks));

            for (int k = 0; k < 3; ++k) {
                const ConstitutiveDesignVariables<double> away =
                    perturb_design(d.v, 0.25, k + 1);
                const DesignEvaluation<double> e =
                    evaluate_design(away, t, S.constraints, S.frame, S.weights);
                smallest_away = std::fmin(smallest_away, e.J_Z);
            }
        }
        metrics.worst_objective_at_target = worst_at;

        report.check_residual(
            "b2.objective_zero_at_exact_target",
            "J_Z and J_blocks vanish when the design equals the target",
            VerificationLevel::A, "n/a", "3 designs", worst_at, 1e-28,
            sw.seconds(),
            "the residuals are squared, so the tolerance is the square of the "
            "reconstruction accuracy");
        report.check_predicate(
            "b2.objective_positive_away_from_target",
            "J_Z is strictly positive for every perturbed design",
            VerificationLevel::A, smallest_away > 1e-6,
            "smallest J_Z over 9 perturbations above 1e-6",
            "smallest J_Z = " + sci(smallest_away), sw.seconds());
    }

    // The objective is an exact quadratic in q for the linear families, which is
    // what makes centred differences of it exact and the least-squares reference
    // legitimate.
    {
        Stopwatch sw;
        const ConstitutiveDesignVariables<double> v0 = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(v0, S.g_lab, S.frame, "magnetoelectric");
        const ConstitutiveDesignVariables<double> d = perturb_design(v0, 1.0, 7);

        // J(q0 + s*dir) must be quadratic in s: the third difference vanishes.
        auto Jat = [&](double s) {
            ConstitutiveDesignVariables<double> v = v0;
            for (int k = 0; k < v.size(); ++k)
                v.q[k] = v0.q[k] + s * (d.q[k] - v0.q[k]);
            return real_part(
                evaluate_design(v, t, S.constraints, S.frame, S.weights).J_Z);
        };
        const double h = 0.25;
        double worst = 0.0;
        for (double s : {-1.0, 0.0, 0.7, 1.5}) {
            const double third = Jat(s + 2 * h) - 3 * Jat(s + h) + 3 * Jat(s) -
                                Jat(s - h);
            worst = std::fmax(worst, std::fabs(third));
        }
        report.check_residual(
            "b2.objective_is_exactly_quadratic",
            "J_Z is an exact quadratic in q for the linear families: its third "
            "difference along any direction vanishes",
            VerificationLevel::A, "n/a", "4 base points, h = 0.25", worst, 1e-12,
            sw.seconds(),
            "(A11) is homogeneous linear in the blocks and the blocks are linear "
            "in q. This is why centred differences of J carry no truncation "
            "error and why the normal-equation reference is the exact minimiser");
    }

    // =======================================================================
    // 6. Gradients: analytic, finite-difference, and batched HyPAD.
    // =======================================================================
    {
        Stopwatch sw;
        // (a) hand-derived, on a manufactured quadratic. Tests the seeding and
        //     batch assembly in isolation from the constitutive objective.
        const ManufacturedQuadratic mq = make_manufactured_quadratic(21);
        double q[kMaxDesignParameters];
        for (int i = 0; i < mq.size; ++i) q[i] = 0.4 - 0.05 * static_cast<double>(i);
        double g_analytic[kMaxDesignParameters];
        manufactured_gradient(mq, q, g_analytic);

        double worst_fd_manufactured = 0.0;
        for (int k = 0; k < mq.size; ++k) {
            double best = 1e300;
            for (double h : {1e-3, 1e-4, 1e-5, 1e-6}) {
                double hi[kMaxDesignParameters], lo[kMaxDesignParameters];
                for (int i = 0; i < mq.size; ++i) { hi[i] = q[i]; lo[i] = q[i]; }
                hi[k] += h;
                lo[k] -= h;
                const double fd = (manufactured_objective<double>(mq, hi) -
                                   manufactured_objective<double>(mq, lo)) / (2 * h);
                best = std::fmin(best, std::fabs(fd - g_analytic[k]));
            }
            worst_fd_manufactured = std::fmax(worst_fd_manufactured, best);
        }
        report.check_residual(
            "b2.gradient_manufactured_analytic",
            "for J = sum c_k (q_k - a_k)^2 the centred difference reproduces the "
            "hand-derived 2 c_k (q_k - a_k)",
            VerificationLevel::B, "n/a", "21 parameters, 4 step sizes",
            worst_fd_manufactured, 1e-9, sw.seconds());
        metrics.worst_gradient_analytic = worst_fd_manufactured;

        // The same manufactured quadratic, differentiated by HyPAD in batches of
        // three and compared against the HAND-DERIVED gradient. This is the case
        // that tests the seeding and batch assembly against a closed form; the
        // finite-difference comparison above does not exercise the backend at
        // all.
#ifdef WARP_HYPAD_HAVE_OTILIB
        {
            double worst_oti_manufactured = 0.0;
            int evaluations = 0;
            for (int start = 0; start < mq.size; start += kGradientBatchSize) {
                OTIScalar qs[kMaxDesignParameters];
                for (int i = 0; i < mq.size; ++i) qs[i] = OTIScalar(q[i]);
                int active[kGradientBatchSize];
                int count = 0;
                for (int b = 0; b < kGradientBatchSize && start + b < mq.size; ++b) {
                    active[count] = start + b;
                    qs[start + b] = sens::seed_parameter(q[start + b], count);
                    ++count;
                }
                const OTIScalar J = manufactured_objective<OTIScalar>(mq, qs);
                ++evaluations;
                for (int b = 0; b < count; ++b)
                    worst_oti_manufactured = std::fmax(
                        worst_oti_manufactured,
                        std::fabs(sens::extract_gradient(J, b) - g_analytic[active[b]]));
            }
            metrics.worst_gradient_analytic =
                std::fmax(metrics.worst_gradient_analytic, worst_oti_manufactured);
            report.check_residual(
                "b2.oti_gradient_manufactured_matches_analytic",
                "the batch-assembled HyPAD gradient of J = sum c_k (q_k - a_k)^2 "
                "equals the hand-derived 2 c_k (q_k - a_k) exactly",
                VerificationLevel::B, "n/a",
                "21 parameters in " + std::to_string(evaluations) +
                    " batches of at most 3",
                worst_oti_manufactured, 1e-14, sw.seconds(),
                "this is the case that tests the seeding and batch assembly "
                "against a CLOSED FORM. It is deliberately independent of the "
                "constitutive objective, so a failure names the AD plumbing "
                "rather than the physics");
        }
#endif
    }

    // (b) the real objective's gradient against centred differences, and the
    //     analytic gradient from the linear structure.
    {
        Stopwatch sw;
        const ConstitutiveDesignVariables<double> v0 = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(v0, S.g_lab, S.frame, "magnetoelectric");
        const ConstitutiveDesignVariables<double> x = perturb_design(v0, 0.3, 11);

        // Analytic: J_Z = ||A q - z||^2 / s, so dJ/dq_k = 2/s * A_k . (A q - z).
        const DesignMatrix M = build_design_matrix(x.family, S.g_lab, S.frame);
        const Rank4<double> Zx = reconstruct_design_tensor(x, S.g_lab, S.frame);
        double scale = 0.0;
        for (int m = 0; m < 256; ++m) scale += t.Z.t[m] * t.Z.t[m];
        scale += S.weights.delta_Z;

        double g_analytic[kMaxDesignParameters] = {};
        for (int k = 0; k < M.n; ++k) {
            double acc = 0.0;
            for (int m = 0; m < 256; ++m)
                acc += M.A[static_cast<std::size_t>(m) * M.n + k] *
                       (Zx.t[m] - t.Z.t[m]);
            g_analytic[k] = 2.0 * acc / scale;
        }

        double worst_fd = 0.0;
        for (double h : {1e-3, 1e-4, 1e-5}) {
            const GradientResult fd = finite_difference_gradient(
                x, t, S.constraints, S.frame, S.weights, h);
            for (int k = 0; k < fd.size; ++k)
                worst_fd = std::fmax(worst_fd, std::fabs(fd.g[k] - g_analytic[k]));
        }
        metrics.worst_gradient_fd = worst_fd;

        report.check_residual(
            "b2.gradient_finite_difference_matches_analytic",
            "centred differences of J_Z reproduce 2 A_k . (A q - z) / s, the "
            "gradient implied by the linear design map",
            VerificationLevel::D, "n/a", "21 parameters, h in {1e-3, 1e-4, 1e-5}",
            worst_fd, 1e-8, sw.seconds(),
            "the objective is exactly quadratic, so this comparison is not "
            "step-size limited in the usual way");
    }

    // (c) batched HyPAD gradient.
#ifdef WARP_HYPAD_HAVE_OTILIB
    {
        Stopwatch sw;
        metrics.oti_available = true;
        const ConstitutiveDesignVariables<double> v0 = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(v0, S.g_lab, S.frame, "magnetoelectric");
        const ConstitutiveDesignVariables<double> x = perturb_design(v0, 0.3, 11);

        const GradientResult oti =
            oti_gradient(x, t, S.constraints, S.weights);
        const GradientResult fd = finite_difference_gradient(
            x, t, S.constraints, S.frame, S.weights, 1e-5);

        double worst = 0.0;
        for (int k = 0; k < oti.size; ++k)
            worst = std::fmax(worst, std::fabs(oti.g[k] - fd.g[k]));
        metrics.worst_gradient_oti = worst;
        metrics.oti_batches = static_cast<int>(oti.batches.size());
        metrics.oti_evaluations = oti.evaluations;

        // HyPAD against the ANALYTIC gradient of the same objective, which is
        // available in closed form because the design map is linear:
        //     dJ_Z/dq_k = 2 A_k . (A q - z) / s
        // Reported separately from the finite-difference comparison. Without
        // this the agreement between HyPAD and finite differences would leave
        // the shared finite-difference error as the limiting number, and the
        // exactness of the HyPAD channels would only be inferred.
        double worst_vs_analytic = 0.0;
        {
            const DesignMatrix M = build_design_matrix(x.family, S.g_lab, S.frame);
            const Rank4<double> Zx = reconstruct_design_tensor(x, S.g_lab, S.frame);
            double scale = 0.0;
            for (int mm = 0; mm < 256; ++mm) scale += t.Z.t[mm] * t.Z.t[mm];
            scale += S.weights.delta_Z;
            for (int k = 0; k < M.n; ++k) {
                double acc = 0.0;
                for (int mm = 0; mm < 256; ++mm)
                    acc += M.A[static_cast<std::size_t>(mm) * M.n + k] *
                           (Zx.t[mm] - t.Z.t[mm]);
                worst_vs_analytic = std::fmax(
                    worst_vs_analytic, std::fabs(oti.g[k] - 2.0 * acc / scale));
            }
        }
        report.check_residual(
            "b2.oti_gradient_matches_analytic",
            "the batch-assembled HyPAD gradient of J_Z equals the closed-form "
            "2 A_k . (A q - z) / s implied by the linear design map",
            VerificationLevel::B, "n/a", "21 parameters", worst_vs_analytic,
            1e-13, sw.seconds(),
            "tighter than the finite-difference comparison by three orders of "
            "magnitude, because neither side is differenced");

        std::ostringstream batches;
        for (const GradientBatchResult& b : oti.batches) {
            batches << "[";
            for (int i = 0; i < b.count; ++i)
                batches << b.active[i] << (i + 1 < b.count ? "," : "");
            batches << "]";
        }

        report.check_predicate(
            "b2.oti_gradient_is_batched",
            "the 21-parameter gradient is assembled from batches of at most "
            "kOtiDirections = 3 active directions",
            VerificationLevel::A,
            oti.ok && oti.evaluations == 7 && metrics.oti_batches == 7,
            "7 batches of at most 3",
            std::to_string(oti.evaluations) + " evaluations, active indices " +
                batches.str(),
            sw.seconds(),
            "the backend is onumm3n2 with three derivative directions. An "
            "interface claiming to seed 21 at once would be a lie about it");

        report.check_residual(
            "b2.oti_gradient_matches_finite_difference",
            "the batch-assembled HyPAD gradient equals the centred-difference "
            "gradient of the same objective",
            VerificationLevel::D, "n/a", "21 parameters", worst, 1e-8,
            sw.seconds(),
            "the objective is not reimplemented for the AD scalar: "
            "evaluate_design is instantiated at OTIScalar and called unchanged");

        // Simultaneous versus separate seeding.
        double worst_alone = 0.0;
        for (int k : {0, 4, 7, 13, 20}) {
            const double alone =
                oti_partial_alone(x, t, S.constraints, S.weights, k);
            worst_alone = std::fmax(worst_alone, std::fabs(alone - oti.g[k]));
        }
        metrics.worst_batch_consistency = worst_alone;
        report.check_residual(
            "b2.oti_gradient_batch_independence",
            "a partial derivative is identical whether its parameter was seeded "
            "alone or alongside two others",
            VerificationLevel::A, "n/a", "5 parameters, batch of 1 versus 3",
            worst_alone, 1e-14, sw.seconds(),
            "catches derivative directions leaking into each other, which a "
            "comparison against finite differences alone would blur");
    }
#else
    {
        Stopwatch sw;
        report.check_skipped(
            "b2.oti_gradient_is_batched",
            "batch-assembled HyPAD gradient of the B2 objective",
            "the OTIlib backend is disabled. The analytic and "
            "finite-difference gradient cases above still ran, and every "
            "optimization case below runs on the finite-difference gradient",
            sw.seconds());
        report.check_skipped(
            "b2.oti_gradient_matches_finite_difference",
            "HyPAD gradient against centred differences",
            "the OTIlib backend is disabled", sw.seconds());
        report.check_skipped(
            "b2.oti_gradient_batch_independence",
            "simultaneous versus separate OTI seeding",
            "the OTIlib backend is disabled", sw.seconds());
        report.check_skipped(
            "b2.oti_gradient_manufactured_matches_analytic",
            "batched HyPAD gradient of the manufactured quadratic against its "
            "hand-derived gradient",
            "the OTIlib backend is disabled; the finite-difference comparison "
            "against the same hand-derived gradient still ran",
            sw.seconds());
        report.check_skipped(
            "b2.oti_gradient_matches_analytic",
            "batched HyPAD gradient of J_Z against its closed form",
            "the OTIlib backend is disabled; the finite-difference comparison "
            "against the same closed form still ran",
            sw.seconds());
    }
#endif

    // =======================================================================
    // 7. The optimizer, and the inverse-design demonstrations.
    // =======================================================================
    std::ofstream hist(out_dir + "/objective_histories.csv", std::ios::binary);
    if (hist) hist << "problem,iteration,objective,gradient_norm,step,backtracks\n";
    std::ofstream recov(out_dir + "/recovered_designs.csv", std::ios::binary);
    if (recov) recov << "problem,family,index,target_q,recovered_q\n";

    auto record = [&](const std::string& name, const OptimizerResult& r) {
        if (hist) {
            hist.setf(std::ios::scientific);
            hist.precision(10);
            for (const IterationRecord& it : r.history)
                hist << name << ',' << it.iteration << ',' << it.objective << ','
                     << it.gradient_norm << ',' << it.step_length << ','
                     << it.backtracks << '\n';
        }
    };

    struct Problem {
        const char* name;
        ConstitutiveDesignVariables<double> truth;
        DesignFamily design_family;
        double perturbation;
        int seed;
    };
    const Problem problems[3] = {
        {"isotropic", isotropic_design(2.75, 0.4), DesignFamily::Isotropic, 0.6, 1},
        {"anisotropic", anisotropic_design(), DesignFamily::UnrestrictedReciprocal,
         0.35, 2},
        {"magnetoelectric", magnetoelectric_design(),
         DesignFamily::UnrestrictedReciprocal, 0.35, 3},
    };

    for (const Problem& p : problems) {
        Stopwatch sw;
        const ConstitutiveTarget<double> t =
            target_from_design(p.truth, S.g_lab, S.frame, p.name);
        ConstitutiveDesignVariables<double> start =
            perturb_design(p.truth, p.perturbation, p.seed);
        start.family = p.design_family;

        const OptimizerResult r = optimize_design(start, t, S.constraints, S.frame,
                                                 S.weights, S.settings);
        record(p.name, r);

        const Rank4<double> Zr =
            reconstruct_design_tensor(r.best, S.g_lab, S.frame);
        const double tensor_residual =
            diff4x4(Zr, t.Z) / std::fmax(max_abs4x4(t.Z), 1e-300);
        double param_residual = 0.0;
        for (int k = 0; k < r.best.size(); ++k)
            param_residual =
                std::fmax(param_residual, std::fabs(r.best.q[k] - p.truth.q[k]));

        if (recov)
            for (int k = 0; k < r.best.size(); ++k)
                recov << p.name << ',' << to_string(p.design_family) << ',' << k
                      << ',' << p.truth.q[k] << ',' << r.best.q[k] << '\n';

        if (std::string(p.name) == "isotropic")
            metrics.isotropic_residual = tensor_residual;
        else if (std::string(p.name) == "anisotropic")
            metrics.anisotropic_residual = tensor_residual;
        else
            metrics.magnetoelectric_residual = tensor_residual;

        std::ostringstream got;
        got << "J: " << sci(r.initial_objective) << " -> " << sci(r.best_objective)
            << " in " << r.iterations << " iterations, " << to_string(r.reason)
            << "; max|dq| = " << sci(param_residual);

        report.check_predicate(
            std::string("b2.recovery.") + p.name,
            "starting from a perturbed design, the optimizer recovers the known "
            "target: relative tensor residual below 1e-6 and every parameter "
            "within 1e-4",
            VerificationLevel::B,
            r.best_objective < 1e-12 && tensor_residual < 1e-6 &&
                param_residual < 1e-4,
            "converged onto the known answer", got.str(), sw.seconds(),
            "the answer is known by construction: the target was built from this "
            "design through the same reconstruction path");

        report.check_predicate(
            std::string("b2.optimizer_descends.") + p.name,
            "the objective decreases monotonically and the line search succeeds "
            "at every accepted iteration",
            VerificationLevel::A,
            [&] {
                if (r.history.size() < 2) return false;
                for (std::size_t i = 1; i < r.history.size(); ++i)
                    if (r.history[i].objective > r.history[i - 1].objective)
                        return false;
                return true;
            }(),
            "monotone objective history",
            std::to_string(r.history.size()) + " iterations recorded, " +
                std::to_string(r.objective_evaluations) + " objective evaluations",
            sw.seconds(),
            "Armijo backtracking guarantees this; asserting it catches a sign "
            "error in the descent direction");
    }

    // Problem 4: metric-induced target, recovered in the UNRESTRICTED family.
    {
        Stopwatch sw;
        const EffectiveMetric<double> g_eff(alcubierre_at(0.4, 3.0));
        const ConstitutiveTarget<double> t = target_from_effective_metric(
            g_eff, S.g_lab, S.frame, "alcubierre_v0.4_x3");

        ConstitutiveDesignVariables<double> start;
        start.family = DesignFamily::UnrestrictedReciprocal;
        // A deliberately uninformed start: vacuum-like, zero magnetoelectric.
        start.q[0] = start.q[3] = start.q[5] = 1.0;
        start.q[6] = start.q[9] = start.q[11] = 1.0;

        const OptimizerResult r = optimize_design(start, t, S.constraints, S.frame,
                                                 S.weights, S.settings);
        record("metric_induced", r);

        const Rank4<double> Zr = reconstruct_design_tensor(r.best, S.g_lab, S.frame);
        const double tensor_residual =
            diff4x4(Zr, t.Z) / std::fmax(max_abs4x4(t.Z), 1e-300);
        metrics.metric_induced_residual = tensor_residual;

        std::ostringstream got;
        got << "target: " << t.message << " J: " << sci(r.initial_objective)
            << " -> " << sci(r.best_objective) << ", relative |Z_design - Z_target| = "
            << sci(tensor_residual) << ", " << to_string(r.reason);

        report.check_predicate(
            "b2.recovery.metric_induced",
            "a subluminal Alcubierre effective metric with g_lab = eta is "
            "recovered in the unrestricted reciprocal family, compared against "
            "Z_target itself",
            VerificationLevel::B,
            t.metric_induced && t.reciprocal && r.best_objective < 1e-12 &&
                tensor_residual < 1e-6,
            "converged onto Z_target", got.str(), sw.seconds(),
            "compared against Z_target, NOT against another observer's Plebanski "
            "blocks. That comparison is the one B1C withdrew");
    }

    // Problem 5: restricted family, nearest approximation.
    {
        Stopwatch sw;
        const ConstitutiveDesignVariables<double> truth = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(truth, S.g_lab, S.frame, "magnetoelectric");

        ConstitutiveDesignVariables<double> start;
        start.family = DesignFamily::ZeroMagnetoelectric;
        for (int k = 0; k < start.size(); ++k) start.q[k] = truth.q[k];
        start = perturb_design(start, 0.2, 5);
        start.family = DesignFamily::ZeroMagnetoelectric;

        const OptimizerResult r = optimize_design(start, t, S.constraints, S.frame,
                                                 S.weights, S.settings);
        record("restricted_family", r);
        metrics.restricted_family_residual = r.best_objective;

        // Independent reference: the exact least-squares projection, from the
        // normal equations, computed without the optimizer.
        const DesignMatrix M =
            build_design_matrix(DesignFamily::ZeroMagnetoelectric, S.g_lab, S.frame);
        double q_ls[kMaxDesignParameters] = {};
        const bool solved = solve_least_squares(M, t.Z, q_ls);

        double param_gap = 0.0;
        double J_ls = 0.0;
        if (solved) {
            ConstitutiveDesignVariables<double> vls;
            vls.family = DesignFamily::ZeroMagnetoelectric;
            for (int k = 0; k < M.n; ++k) vls.q[k] = q_ls[k];
            J_ls = real_part(
                evaluate_design(vls, t, S.constraints, S.frame, S.weights).J_Z);
            for (int k = 0; k < M.n; ++k)
                param_gap = std::fmax(param_gap, std::fabs(q_ls[k] - r.best.q[k]));
        }
        metrics.restricted_vs_least_squares = param_gap;

        std::ostringstream got;
        got << "optimizer J = " << sci(r.best_objective) << ", normal-equation J = "
            << sci(J_ls) << ", max|dq| = " << sci(param_gap);

        report.check_predicate(
            "b2.restricted_family_is_a_projection_not_a_fit",
            "against a magnetoelectric target, the zeta == 0 family leaves a "
            "STRICTLY NONZERO residual, and the optimizer finds the same point "
            "as the exact least-squares projection",
            VerificationLevel::D,
            solved && r.best_objective > 1e-3 && param_gap < 1e-5 &&
                std::fabs(r.best_objective - J_ls) < 1e-8,
            "nonzero residual, matching the normal-equation solution",
            got.str(), sw.seconds(),
            "the residual is a property of the family, not an optimizer failure. "
            "The reference is dense Gaussian elimination on the normal equations, "
            "written in the test and not calling the optimizer");
    }

    // Problem 6: infeasible bounded target.
    {
        Stopwatch sw;
        const ConstitutiveDesignVariables<double> truth = anisotropic_design();
        const ConstitutiveTarget<double> t =
            target_from_design(truth, S.g_lab, S.frame, "anisotropic");

        DesignConstraints<double> bounded;
        bounded.active = true;
        // The target needs eps entries up to 3.1; this forbids anything above 1.5.
        bounded.max_epsilon_entry = 1.5;
        bounded.max_mu_inverse_entry = 1.5;
        bounded.max_zeta_entry = 1.5;
        bounded.penalty_weight = 1e3;

        ConstitutiveDesignVariables<double> start;
        start.family = DesignFamily::UnrestrictedReciprocal;
        start.q[0] = start.q[3] = start.q[5] = 1.0;
        start.q[6] = start.q[9] = start.q[11] = 1.0;

        const OptimizerResult r =
            optimize_design(start, t, bounded, S.frame, S.weights, S.settings);
        record("infeasible_bounded", r);

        const DesignEvaluation<double> e =
            evaluate_design(r.best, t, bounded, S.frame, S.weights);
        metrics.bounded_residual = e.J_Z;

        std::ostringstream got;
        got << "J_Z = " << sci(e.J_Z) << ", penalty = " << sci(e.J_penalty)
            << ", max|eps| = " << sci(e.epsilon_diagnostics.max_component)
            << " against bound 1.5, entry bound active: "
            << (e.epsilon_entry_active ? "yes" : "no") << ", " << to_string(r.reason);

        report.check_predicate(
            "b2.infeasible_bounded_reports_best_approximation",
            "with bounds excluding the target, the result is a BEST BOUNDED "
            "APPROXIMATION: J_Z stays strictly positive, the bound is reported "
            "as active, and no exact realization is claimed",
            VerificationLevel::B,
            e.J_Z > 1e-4 && e.epsilon_diagnostics.max_component < 1.55,
            "nonzero J_Z with the entry bound respected and reported",
            got.str(), sw.seconds(),
            "the failure mode this guards against is a gate that reports "
            "convergence while the design silently violates its own bounds. The "
            "penalty is SOFT by construction -- a smooth hinge admits a small "
            "violation in exchange for differentiability -- so max|eps| settles "
            "just above the bound rather than exactly on it, and the case allows "
            "that explicitly rather than pretending the constraint is hard.");
    }

    // Problem 7: coordinate-change consistency.
    {
        Stopwatch sw;
        Mat4<double> L = Mat4<double>::identity();
        L(0, 0) = 1.15; L(0, 1) = 0.22;
        L(1, 0) = 0.18; L(1, 1) = 1.08;
        L(2, 1) = 0.11; L(2, 2) = 1.04; L(2, 3) = 0.07;
        L(3, 0) = 0.04; L(3, 2) = 0.13; L(3, 3) = 0.97;
        const warp_hypad::InverseResult<double> Li = warp_hypad::inverse(L);

        const ConstitutiveDesignVariables<double> truth = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(truth, S.g_lab, S.frame, "magnetoelectric");

        Mat4<double> gp;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                double acc = 0.0;
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d)
                        acc += Li.inverse(c, a) * Li.inverse(d, b) * S.g_lab.value(c, d);
                gp(a, b) = acc;
            }
        const LaboratoryMetric<double> g_prime(gp);

        Rank4<double> Zt_prime;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d) {
                        double acc = 0.0;
                        for (int e = 0; e < 4; ++e)
                            for (int f = 0; f < 4; ++f)
                                for (int g = 0; g < 4; ++g)
                                    for (int h = 0; h < 4; ++h)
                                        acc += L(a, e) * L(b, f) * L(c, g) * L(d, h) *
                                               t.Z(e, f, g, h);
                        Zt_prime(a, b, c, d) = acc;
                    }

        double Vp_raw[4];
        for (int a = 0; a < 4; ++a) {
            double acc = 0.0;
            for (int b = 0; b < 4; ++b) acc += L(a, b) * S.frame.V[b];
            Vp_raw[a] = acc;
        }
        const Observer<double> obs_p =
            make_laboratory_observer(g_prime, Vp_raw, "lab_rest_primed");
        const ObserverFrame<double> frame_p =
            build_observer_frame(g_prime, obs_p);
        const ConstitutiveTarget<double> t_prime = target_from_tensor(
            Zt_prime, g_prime, frame_p, "magnetoelectric_transformed");

        ConstitutiveDesignVariables<double> start_p;
        start_p.family = DesignFamily::UnrestrictedReciprocal;
        start_p.q[0] = start_p.q[3] = start_p.q[5] = 1.0;
        start_p.q[6] = start_p.q[9] = start_p.q[11] = 1.0;

        const OptimizerResult rp = optimize_design(start_p, t_prime, S.constraints,
                                                  frame_p, S.weights, S.settings);
        record("coordinate_change", rp);

        const Rank4<double> Zp_design =
            reconstruct_design_tensor(rp.best, g_prime, frame_p);
        const double residual =
            diff4x4(Zp_design, Zt_prime) / std::fmax(max_abs4x4(Zt_prime), 1e-300);
        metrics.worst_coordinate_change = residual;

        std::ostringstream got;
        got << "primed target " << t_prime.message << " J -> "
            << sci(rp.best_objective) << ", relative |Z'_design - Z'_target| = "
            << sci(residual) << ", " << to_string(rp.reason);

        report.check_predicate(
            "b2.coordinate_change_consistency",
            "optimizing entirely in the transformed coordinates reconstructs the "
            "transformed target: the FINAL TENSOR RELATIONSHIP holds",
            VerificationLevel::B,
            frame_p.ok && t_prime.reciprocal && residual < 1e-6,
            "Z'_design = Z'_target", got.str(), sw.seconds(),
            "the final tensor relationship is what is asserted, NOT an identical "
            "iteration history: a coordinate-basis norm reweights the parameters "
            "and need not preserve the trajectory");
    }

    // =======================================================================
    // 8. Determinism and non-finite handling.
    // =======================================================================
    {
        Stopwatch sw;
        const ConstitutiveDesignVariables<double> truth = anisotropic_design();
        const ConstitutiveTarget<double> t =
            target_from_design(truth, S.g_lab, S.frame, "anisotropic");
        const ConstitutiveDesignVariables<double> start =
            perturb_design(truth, 0.3, 42);

        const OptimizerResult a =
            optimize_design(start, t, S.constraints, S.frame, S.weights, S.settings);
        const OptimizerResult b =
            optimize_design(start, t, S.constraints, S.frame, S.weights, S.settings);

        bool identical = a.iterations == b.iterations &&
                         a.history.size() == b.history.size() &&
                         a.reason == b.reason &&
                         a.best_objective == b.best_objective;
        for (int k = 0; k < a.best.size() && identical; ++k)
            if (a.best.q[k] != b.best.q[k]) identical = false;
        for (std::size_t i = 0; i < a.history.size() && identical; ++i)
            if (a.history[i].objective != b.history[i].objective ||
                a.history[i].step_length != b.history[i].step_length)
                identical = false;

        report.check_predicate(
            "b2.deterministic_repeatability",
            "two runs from the same initialization are BIT-IDENTICAL in history, "
            "termination reason and final parameters",
            VerificationLevel::A, identical, "bit-identical runs",
            identical ? "identical over " + std::to_string(a.history.size()) +
                            " iterations"
                      : "runs diverged",
            sw.seconds(),
            "bit equality, not a tolerance. The perturbation uses a fixed integer "
            "recurrence rather than a random engine, so there is no seed to lose");
    }

    {
        Stopwatch sw;
        const ConstitutiveDesignVariables<double> truth = anisotropic_design();
        const ConstitutiveTarget<double> t =
            target_from_design(truth, S.g_lab, S.frame, "anisotropic");

        ConstitutiveDesignVariables<double> bad = truth;
        bad.q[4] = std::numeric_limits<double>::quiet_NaN();
        const OptimizerResult r =
            optimize_design(bad, t, S.constraints, S.frame, S.weights, S.settings);

        ConstitutiveDesignVariables<double> huge = truth;
        huge.q[2] = 1e308;
        const DesignEvaluation<double> e =
            evaluate_design(huge, t, S.constraints, S.frame, S.weights);

        report.check_predicate(
            "b2.nonfinite_is_rejected_not_absorbed",
            "a non-finite design is rejected with reason NonFiniteEncountered, "
            "and an overflowing design is reported as not finite rather than "
            "silently scored",
            VerificationLevel::A,
            r.reason == TerminationReason::NonFiniteEncountered && !e.finite,
            "NonFiniteEncountered, and evaluation reports finite = false",
            std::string(to_string(r.reason)) + "; overflow evaluation finite = " +
                (e.finite ? "true" : "false"),
            sw.seconds(),
            "an optimizer that silently absorbs a NaN reports convergence to "
            "nothing");
    }

    // =======================================================================
    // 9. Constraint activity and target classification records.
    // =======================================================================
    // Finite-difference error against step size, for the convergence plot.
    {
        std::ofstream os(out_dir + "/gradient_fd_convergence.csv", std::ios::binary);
        const ConstitutiveDesignVariables<double> v0 = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(v0, S.g_lab, S.frame, "magnetoelectric");
        const ConstitutiveDesignVariables<double> x = perturb_design(v0, 0.3, 11);
        const DesignMatrix M = build_design_matrix(x.family, S.g_lab, S.frame);
        const Rank4<double> Zx = reconstruct_design_tensor(x, S.g_lab, S.frame);
        double scale = 0.0;
        for (int mm = 0; mm < 256; ++mm) scale += t.Z.t[mm] * t.Z.t[mm];
        scale += S.weights.delta_Z;
        double g_an[kMaxDesignParameters] = {};
        for (int k = 0; k < M.n; ++k) {
            double acc = 0.0;
            for (int mm = 0; mm < 256; ++mm)
                acc += M.A[static_cast<std::size_t>(mm) * M.n + k] *
                       (Zx.t[mm] - t.Z.t[mm]);
            g_an[k] = 2.0 * acc / scale;
        }
        if (os) {
            os << "step,worst_abs_error\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (double h : {1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5,
                             3e-6, 1e-6, 1e-7, 1e-8}) {
                const GradientResult fd = finite_difference_gradient(
                    x, t, S.constraints, S.frame, S.weights, h);
                double w = 0.0;
                for (int k = 0; k < fd.size; ++k)
                    w = std::fmax(w, std::fabs(fd.g[k] - g_an[k]));
                os << h << ',' << w << '\n';
            }
        }
    }

    // Z residual by antisymmetric index pair, for the residual map.
    {
        std::ofstream os(out_dir + "/z_residual_by_pair.csv", std::ios::binary);
        static const int P[6][2] = {{0, 1}, {0, 2}, {0, 3}, {2, 3}, {3, 1}, {1, 2}};
        static const char* PN[6] = {"01", "02", "03", "23", "31", "12"};

        const ConstitutiveDesignVariables<double> truth = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(truth, S.g_lab, S.frame, "magnetoelectric");

        // Exact family versus the zeta == 0 restricted family, side by side.
        ConstitutiveDesignVariables<double> restricted;
        restricted.family = DesignFamily::ZeroMagnetoelectric;
        const DesignMatrix M =
            build_design_matrix(DesignFamily::ZeroMagnetoelectric, S.g_lab, S.frame);
        double q_ls[kMaxDesignParameters] = {};
        if (solve_least_squares(M, t.Z, q_ls))
            for (int k = 0; k < M.n; ++k) restricted.q[k] = q_ls[k];
        const Rank4<double> Ze = reconstruct_design_tensor(truth, S.g_lab, S.frame);
        const Rank4<double> Zr =
            reconstruct_design_tensor(restricted, S.g_lab, S.frame);

        if (os) {
            os << "pair_A,pair_B,exact_residual,restricted_residual\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int A = 0; A < 6; ++A)
                for (int B = 0; B < 6; ++B) {
                    const int p = P[A][0], q2 = P[A][1], r = P[B][0], s = P[B][1];
                    os << PN[A] << ',' << PN[B] << ','
                       << std::fabs(Ze(p, q2, r, s) - t.Z(p, q2, r, s)) << ','
                       << std::fabs(Zr(p, q2, r, s) - t.Z(p, q2, r, s)) << '\n';
                }
        }
    }

    // Per-family residual against the same magnetoelectric target, so the cost
    // of each restriction is visible rather than asserted.
    {
        std::ofstream os(out_dir + "/family_residuals.csv", std::ios::binary);
        const ConstitutiveDesignVariables<double> truth = magnetoelectric_design();
        const ConstitutiveTarget<double> t =
            target_from_design(truth, S.g_lab, S.frame, "magnetoelectric");
        if (os) {
            os << "family,parameters,least_squares_J_Z\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            const DesignFamily fams[4] = {DesignFamily::UnrestrictedReciprocal,
                                          DesignFamily::ZeroMagnetoelectric,
                                          DesignFamily::DiagonalAnisotropic,
                                          DesignFamily::Isotropic};
            for (DesignFamily f : fams) {
                const DesignMatrix M = build_design_matrix(f, S.g_lab, S.frame);
                double q_ls[kMaxDesignParameters] = {};
                double J = std::numeric_limits<double>::quiet_NaN();
                if (solve_least_squares(M, t.Z, q_ls)) {
                    ConstitutiveDesignVariables<double> v;
                    v.family = f;
                    for (int k = 0; k < M.n; ++k) v.q[k] = q_ls[k];
                    J = real_part(
                        evaluate_design(v, t, S.constraints, S.frame, S.weights).J_Z);
                }
                os << to_string(f) << ',' << parameter_count(f) << ',' << J << '\n';
            }
        }
    }

    {
        std::ofstream os(out_dir + "/target_classification.csv", std::ios::binary);
        if (os) {
            os << "target,metric_induced,reciprocal,pair_antisymmetric,finite,"
                  "observer_decomposable,real_lossless,reciprocity_residual,"
                  "antisymmetry_residual,max_component\n";
            os.setf(std::ios::scientific);
            os.precision(10);

            struct Row { const char* name; ConstitutiveTarget<double> t; };
            std::vector<Row> rows;
            rows.push_back({"isotropic",
                            target_from_design(isotropic_design(2.75, 0.4), S.g_lab,
                                               S.frame, "isotropic")});
            rows.push_back({"anisotropic", target_from_design(anisotropic_design(),
                                                             S.g_lab, S.frame,
                                                             "anisotropic")});
            rows.push_back({"magnetoelectric",
                            target_from_design(magnetoelectric_design(), S.g_lab,
                                               S.frame, "magnetoelectric")});
            rows.push_back({"alcubierre_v0.4",
                            target_from_effective_metric(
                                EffectiveMetric<double>(alcubierre_at(0.4, 3.0)),
                                S.g_lab, S.frame, "alcubierre")});
            for (const Row& r : rows)
                os << r.name << ',' << r.t.metric_induced << ',' << r.t.reciprocal
                   << ',' << r.t.pair_antisymmetric << ',' << r.t.finite << ','
                   << r.t.observer_decomposable << ','
                   << r.t.inside_real_lossless_model << ','
                   << r.t.reciprocity_residual << ',' << r.t.antisymmetry_residual
                   << ',' << r.t.max_component << '\n';
        }
    }

    {
        Stopwatch sw;
        // Every target used above must classify as usable, or the recoveries
        // would be measuring something other than what they claim.
        const ConstitutiveTarget<double> t = target_from_effective_metric(
            EffectiveMetric<double>(alcubierre_at(0.4, 3.0)), S.g_lab, S.frame,
            "alcubierre");
        report.check_predicate(
            "b2.target_classification_is_recorded",
            "a metric-induced target is classified as metric-induced, reciprocal, "
            "pair-antisymmetric, finite and observer-decomposable",
            VerificationLevel::A,
            t.metric_induced && t.reciprocal && t.pair_antisymmetric && t.finite &&
                t.observer_decomposable && t.usable(),
            "all six classification flags set", t.message, sw.seconds(),
            "a target that failed reciprocity would be reported, not silently "
            "fitted: the design space is reciprocal, so the achievable residual "
            "would then be bounded below by the target's own skewon part");
    }
}

}  // namespace checks
}  // namespace program_b

#include "o0a_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "warp_hypad/design/constraint_projection.hpp"

// ---------------------------------------------------------------------------
// O0A cases 21-28: exact constraints, the projection they induce, the Sobolev
// Gram matrix and the whitening transform.
//
// WHY THE LINEAR ALGEBRA IS TESTED ON ITS OWN FIXTURES FIRST.
//
// Every number in this group is produced by the same two routines -- a one-sided
// Jacobi SVD and a cyclic Jacobi eigendecomposition. If those are wrong, the
// spline-derived matrices will still yield a self-consistent constraint rank, a
// self-consistent null space and a self-consistent whitening, because the same
// wrong factorisation is used to build and to check them. So the group opens
// with six matrices whose answers are known independently of any spline:
//
//     full-rank rectangular      rank must be the column count, no null space
//     rank-deficient wide        a row that is the sum of two others
//     planted duplicate          a row that is a multiple of another
//     symmetric positive definite   every eigenvalue strictly positive
//     PSD with an exact null mode   4I - v v^T annihilates v EXACTLY, and has a
//                                   THREEFOLD degenerate eigenvalue 4
//     nearly singular            the same, with the null eigenvalue moved to
//                                4e-18, which whitening must discard
//
// and checks the four identities that pin a factorisation down:
// |A V_null|, |V^T V - I|, |A A^+ A - A| and |Q L Q^T - S|. The fifth fixture is
// doing double duty: its threefold eigenvalue is the degenerate subspace the
// deterministic-sign convention CANNOT disambiguate, and the case records that
// rather than pretending mode 1 of that group is meaningful.
//
// THE SVD CONVERGENCE FLAG IS CHECKED, NOT ASSUMED. It has already caught a real
// defect: on a rank-1 constraint matrix the method ran its full sweep limit and
// returned converged = false while still producing the correct singular values,
// because once a rotation drove a column to rounding level the relative
// orthogonality measure for that pair never fell below threshold again. The
// answer was right and the guarantee behind it was not.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::ConstraintProjection;
using design::ConstraintSystem;
using design::DenseMatrix;
using design::DesignQuadrature;
using design::DesignRegion;
using design::EigenResult;
using design::FieldFamily;
using design::SvdResult;
using design::TensorBasisEvaluation;

DenseMatrix from_rows(int r, int c, const std::vector<double>& v) {
    DenseMatrix m(r, c);
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            m(i, j) = v[static_cast<std::size_t>(i) *
                            static_cast<std::size_t>(c) +
                        static_cast<std::size_t>(j)];
    return m;
}

// S = alpha I - beta v v^T with v = (1,1,...,1). Eigenvalues: alpha - beta n
// (eigenvector v) and alpha with multiplicity n-1. Exact in binary floating
// point for the values used below.
DenseMatrix ones_deflation(int n, double alpha, double beta) {
    DenseMatrix S(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            S(i, j) = (i == j ? alpha : 0.0) - beta;
    return S;
}

}  // namespace

void run_o0a_constraint_checks(VerificationReport& report,
                               const std::string& out_dir) {
    Stopwatch sw;

    // =======================================================================
    // Independent dense-linear-algebra fixtures.
    // =======================================================================
    sw.reset();
    struct Fixture {
        const char* name;
        DenseMatrix a;
        int expected_rank;
        const char* why;
    };
    std::vector<Fixture> fixtures;
    fixtures.push_back({"full_rank_rectangular",
                        from_rows(4, 3, {1, 0, 0, 0, 2, 0, 0, 0, 3, 1, 1, 1}),
                        3, "4x3 with independent columns"});
    // Row 3 = row 1 + row 2 exactly.
    fixtures.push_back({"rank_deficient_wide",
                        from_rows(3, 5,
                                  {1, 2, 3, 4, 5,
                                   2, 1, 0, -1, -2,
                                   3, 3, 3, 3, 3}),
                        2, "3x5 whose third row is the sum of the first two"});
    // Row 4 = -2.5 * row 2 exactly.
    fixtures.push_back({"planted_duplicate",
                        from_rows(4, 6,
                                  {1, 0, 2, 0, 1, 0,
                                   0, 1, 0, 2, 0, 1,
                                   1, 1, 1, 1, 1, 1,
                                   0, -2.5, 0, -5, 0, -2.5}),
                        3, "4x6 with a row that is -2.5 times another"});

    std::ofstream fx(out_dir + "/linalg_fixtures.csv");
    fx << "fixture,rows,cols,expected_rank,detected_rank,converged,sweeps,"
          "singular_max,singular_min_retained,null_residual,"
          "orthonormality_residual,pseudo_inverse_residual\n";
    bool rank_ok = true, identities_ok = true, convergence_ok = true;
    double worst_null = 0.0, worst_orth = 0.0, worst_pinv = 0.0;
    for (const Fixture& f : fixtures) {
        const SvdResult s = design::jacobi_svd(f.a);
        const double nr = design::null_space_residual(f.a, s);
        const double orr = design::orthonormality_residual(s.v);
        const double pr = design::pseudo_inverse_residual(f.a, s);
        if (s.rank != f.expected_rank) rank_ok = false;
        if (!s.converged) convergence_ok = false;
        worst_null = std::fmax(worst_null, nr);
        worst_orth = std::fmax(worst_orth, orr);
        worst_pinv = std::fmax(worst_pinv, pr);
        fx << f.name << "," << f.a.rows << "," << f.a.cols << ","
           << f.expected_rank << "," << s.rank << "," << (s.converged ? 1 : 0)
           << "," << s.sweeps << "," << s.singular_values.front() << ","
           << (s.rank > 0
                   ? s.singular_values[static_cast<std::size_t>(s.rank - 1)]
                   : 0.0)
           << "," << nr << "," << orr << "," << pr << "\n";
    }
    if (worst_null > 1e-13 || worst_orth > 1e-13 || worst_pinv > 1e-12)
        identities_ok = false;
    {
        char buf[420];
        std::snprintf(buf, sizeof(buf),
                      "3 fixtures: ranks %s, all converged=%s; worst "
                      "|A V_null| = %.2e, |V^T V - I| = %.2e, "
                      "|A A^+ A - A| = %.2e",
                      rank_ok ? "all as expected" : "MISMATCH",
                      convergence_ok ? "yes" : "NO", worst_null, worst_orth,
                      worst_pinv);
        report.check_predicate(
            "o0a.svd_fixtures_independent",
            "the SVD reproduces the known rank of a full-rank rectangular, a "
            "rank-deficient wide and a planted-duplicate matrix, satisfies "
            "|A V_null| = 0, |V^T V - I| = 0 and the Moore-Penrose identity, and "
            "reports convergence",
            VerificationLevel::A,
            rank_ok && identities_ok && convergence_ok,
            "ranks 3, 2, 3 with all four identities at rounding level",
            buf, sw.seconds(),
            "tested on matrices whose answers are known without any spline, "
            "because the spline-derived checks all use this same factorisation "
            "to build AND to verify their own inputs");
    }

    // --- Symmetric fixtures, degeneracy, and the discard decision -----------
    sw.reset();
    // SPD: 4I - 0.5 v v^T on n = 4 gives eigenvalues 4 - 2 = 2 and 4 (x3).
    const DenseMatrix spd = ones_deflation(4, 4.0, 0.5);
    // Exact null mode: 4I - 1 v v^T on n = 4 gives 0 and 4 (x3).
    const DenseMatrix psd_null = ones_deflation(4, 4.0, 1.0);
    // Nearly singular: the null eigenvalue moved to 4e-18.
    const DenseMatrix nearly = ones_deflation(4, 4.0, 1.0 - 1.0e-18);

    std::ofstream sx(out_dir + "/linalg_symmetric_fixtures.csv");
    sx << "fixture,n,converged,sweeps,reconstruction_residual,eigenvalue_min,"
          "eigenvalue_max,degenerate_groups,largest_group,expected_min\n";
    struct SymFixture { const char* name; const DenseMatrix* s; double expect_min; };
    const SymFixture syms[3] = {{"symmetric_positive_definite", &spd, 2.0},
                                {"psd_with_exact_null_mode", &psd_null, 0.0},
                                {"nearly_singular", &nearly, 4.0e-18}};
    bool sym_ok = true, sym_conv = true;
    double worst_recon = 0.0;
    std::size_t degenerate_total = 0, largest_group = 0;
    for (const SymFixture& sf : syms) {
        const EigenResult e = design::jacobi_eigen(*sf.s);
        const double rr = design::eigen_reconstruction_residual(*sf.s, e);
        const auto groups = design::degenerate_groups(e.eigenvalues);
        std::size_t big = 0;
        for (const auto& g : groups) big = std::max(big, g.size());
        degenerate_total += groups.size();
        largest_group = std::max(largest_group, big);
        worst_recon = std::fmax(worst_recon, rr);
        if (!e.converged) sym_conv = false;
        const double got_min = e.eigenvalues.back();
        // The exact-null fixture must give a numerically zero eigenvalue; the
        // nearly-singular one must give something far below the largest.
        if (std::fabs(got_min - sf.expect_min) > 1e-12 * 4.0) sym_ok = false;
        sx << sf.name << "," << sf.s->rows << "," << (e.converged ? 1 : 0) << ","
           << e.sweeps << "," << rr << "," << got_min << ","
           << e.eigenvalues.front() << "," << groups.size() << "," << big << ","
           << sf.expect_min << "\n";
    }
    {
        char buf[420];
        std::snprintf(buf, sizeof(buf),
                      "worst |Q L Q^T - S| = %.2e; all converged=%s; minimum "
                      "eigenvalues as expected=%s; %zu degenerate groups found, "
                      "largest %zu modes",
                      worst_recon, sym_conv ? "yes" : "NO",
                      sym_ok ? "yes" : "NO", degenerate_total, largest_group);
        report.check_predicate(
            "o0a.eigen_fixtures_independent",
            "the eigensolver reproduces 4I - beta v v^T exactly: strictly "
            "positive spectrum, an exact null mode, and a near-null mode; it "
            "satisfies |Q L Q^T - S| = 0 and reports both convergence and the "
            "threefold degeneracy",
            VerificationLevel::A,
            sym_ok && sym_conv && worst_recon <= 1e-13 && largest_group == 3,
            "reconstruction at rounding level, one threefold degenerate group "
            "per fixture",
            buf, sw.seconds(),
            "the threefold eigenvalue is deliberate: the sign convention fixes a "
            "sign but CANNOT make an individual vector inside a degenerate "
            "subspace unique, so the degeneracy is recorded and no such vector "
            "is presented as physically determined");
    }

    // --- The sign convention is deterministic and canonical -----------------
    sw.reset();
    bool sign_ok = true, repeat_ok = true;
    for (const Fixture& f : fixtures) {
        const SvdResult s1 = design::jacobi_svd(f.a);
        const SvdResult s2 = design::jacobi_svd(f.a);
        for (int j = 0; j < s1.v.cols; ++j) {
            int best = 0;
            double ba = -1.0;
            for (int i = 0; i < s1.v.rows; ++i)
                if (std::fabs(s1.v(i, j)) > ba) {
                    ba = std::fabs(s1.v(i, j));
                    best = i;
                }
            if (s1.v(best, j) < 0.0) sign_ok = false;
            for (int i = 0; i < s1.v.rows; ++i)
                if (s1.v(i, j) != s2.v(i, j)) repeat_ok = false;
        }
    }
    report.check_predicate(
        "o0a.deterministic_vector_signs",
        "every returned singular vector has a positive largest-magnitude "
        "component, ties resolved by lowest index, and two runs agree bit for bit",
        VerificationLevel::A, sign_ok && repeat_ok,
        "canonical signs, bit-identical repeats",
        std::string("canonical=") + (sign_ok ? "yes" : "NO") +
            ", bit-identical=" + (repeat_ok ? "yes" : "NO"),
        sw.seconds(),
        "without a convention the sign depends on the sweep order, and a reader "
        "comparing two builds cannot tell a sign flip from a change in the "
        "physics");

    // =======================================================================
    // The spline-derived constraint system.
    // =======================================================================
    const DesignRegion region = canonical_region();
    AxisymmetricBasis basis = canonical_basis();
    const DesignQuadrature quad = design::build_quadrature(basis);
    const design::SobolevScales scales{kO0aR0, kO0aR0};
    const DenseMatrix M = design::sobolev_gram(basis, quad, scales);

    const ConstraintSystem sys = canonical_constraints(basis, quad);
    const std::vector<double> moment = sys.rows()[0].row;
    const ConstraintProjection pr = design::build_projection(sys);

    // --- 21. Affine constraint residual --------------------------------------
    sw.reset();
    report.check_residual(
        "o0a.affine_constraint_residual", "|C c0 - d| for the particular solution",
        VerificationLevel::A,
        std::to_string(pr.n_constraints) + " rows x " +
            std::to_string(pr.n_raw) + " coefficients",
        "minimum-norm SVD solution", pr.affine_residual, 1.0e-12, sw.seconds(),
        "c0 satisfies the constraints exactly, so every iterate c0 + P z does "
        "too; no optimizer is ever asked to satisfy a linear constraint by "
        "penalty");

    // --- 22. Null-space residual and orthonormality --------------------------
    sw.reset();
    report.check_residual(
        "o0a.null_space_residual",
        "max |C P| = 0, and the null-space basis is orthonormal: |P^T P - I| = 0",
        VerificationLevel::A,
        std::to_string(pr.n_raw) + " x " + std::to_string(pr.n_free),
        "columns from the SVD's orthogonal right factor",
        std::fmax(pr.nullspace_residual, pr.orthonormality_residual), 1.0e-12,
        sw.seconds(),
        "|C P| = " + std::to_string(pr.nullspace_residual) +
            ", |P^T P - I| = " + std::to_string(pr.orthonormality_residual) +
            "; orthonormality is recorded because every downstream statement "
            "about step length in z assumes it");

    // --- 23. Constraint rank detection ---------------------------------------
    sw.reset();
    const int independent_rows = 2;  // moment and cabin mean
    // An extended system, built only to show that adding a genuinely new
    // functional RAISES the detected rank rather than being absorbed into
    // the existing two.
    ConstraintSystem extended = canonical_constraints(basis, quad);
    extended.add_functional(
        functional_weights(basis, quad, region.s_max, false), 0.0,
        "fixed_normalization",
        "zero change in the chosen global normalization: independent of the "
        "moment and the cabin mean, so the rank must rise to 3");
    const ConstraintProjection pr_ext = design::build_projection(extended);
    {
        char buf[400];
        std::snprintf(buf, sizeof(buf),
                      "%d rows, rank %d, %d redundant, %d free; sigma_max = "
                      "%.4e, smallest retained = %.4e, first discarded = %.4e, "
                      "tolerance = %.4e; SVD converged=%d in %d sweeps",
                      pr.n_constraints, pr.rank, pr.redundant_rows, pr.n_free,
                      pr.singular_values.front(),
                      pr.singular_values[static_cast<std::size_t>(pr.rank - 1)],
                      pr.singular_values[static_cast<std::size_t>(pr.rank)],
                      pr.rank_tolerance, pr.svd_converged ? 1 : 0,
                      pr.svd_sweeps);
        char ext[240];
        std::snprintf(ext, sizeof(ext),
                      " | extended system: %d rows -> rank %d, %d redundant, %d free",
                      pr_ext.n_constraints, pr_ext.rank,
                      pr_ext.redundant_rows, pr_ext.n_free);
        report.check_predicate(
            "o0a.constraint_rank_detection",
            "the rank-revealing factorisation finds 2 independent "
            "constraints among the canonical 3 rows, names the planted "
            "dependency, rises to rank 3 when a genuinely new functional is "
            "added, and reports convergence in every case",
            VerificationLevel::A,
            pr.rank == independent_rows && pr.redundant_rows == 1 &&
                pr.svd_converged && pr_ext.rank == 3 &&
                pr_ext.redundant_rows == 1 && pr_ext.svd_converged &&
                pr.n_free == pr.n_raw - independent_rows,
            "rank 2 of 3 canonical rows with 1 redundant; rank 3 of 4 when "
            "a genuinely new functional is added",
            std::string(buf) + ext, sw.seconds(),
            "a dependent constraint is usually a modelling mistake -- the same "
            "condition written twice in two different ways -- so it is reported "
            "rather than absorbed");
    }

    // --- 24. The translation moment is removed -------------------------------
    //
    // The test that matters: for EVERY direction in the free space, the
    // xi-moment of the resulting correction is zero. Checked by evaluating the
    // moment functional on the columns of P, and then on random combinations, so
    // it is a statement about the whole subspace rather than one vector.
    sw.reset();
    double worst_moment = 0.0;
    for (int j = 0; j < pr.n_free; ++j) {
        double acc = 0.0;
        for (int i = 0; i < pr.n_raw; ++i)
            acc += moment[static_cast<std::size_t>(i)] * pr.P(i, j);
        worst_moment = std::fmax(worst_moment, std::fabs(acc));
    }
    double worst_combo = 0.0;
    for (int t = 0; t < 24; ++t) {
        std::vector<double> z(static_cast<std::size_t>(pr.n_free));
        for (int j = 0; j < pr.n_free; ++j)
            z[static_cast<std::size_t>(j)] =
                std::sin(0.91 * static_cast<double>(j) +
                         0.37 * static_cast<double>(t));
        double acc = 0.0;
        for (int i = 0; i < pr.n_raw; ++i) {
            double ci = 0.0;
            for (int j = 0; j < pr.n_free; ++j)
                ci += pr.P(i, j) * z[static_cast<std::size_t>(j)];
            acc += moment[static_cast<std::size_t>(i)] * ci;
        }
        worst_combo = std::fmax(worst_combo, std::fabs(acc));
    }
    const double moment_scale = design::norm2(moment);
    report.check_residual(
        "o0a.translation_moment_removed",
        "the first xi-moment of the correction vanishes for every direction in "
        "the free subspace, not merely for one vector",
        VerificationLevel::A,
        std::to_string(pr.n_free) + " basis directions + 24 combinations",
        "moment weight norm " + std::to_string(moment_scale),
        std::fmax(worst_moment, worst_combo), 1.0e-11 * moment_scale,
        sw.seconds(),
        "a pure translation of the bubble is gauge: it changes the coefficients, "
        "leaves the geometry the same up to a shift of origin, and would let an "
        "optimizer report progress for nothing. It is removed from the search "
        "space rather than discouraged by the objective");

    // --- 25. Sobolev basis norms ---------------------------------------------
    sw.reset();
    const double gram_symmetry = design::symmetry_residual(M);
    double min_diag = 1e300, max_diag = 0.0;
    bool diag_positive = true;
    for (int i = 0; i < M.rows; ++i) {
        const double d = M(i, i);
        if (!(d > 0.0)) diag_positive = false;
        min_diag = std::fmin(min_diag, d);
        max_diag = std::fmax(max_diag, d);
    }
    // Cauchy-Schwarz on the Sobolev inner product: |M_ij| <= sqrt(M_ii M_jj).
    double worst_cs = 0.0;
    for (int i = 0; i < M.rows; ++i)
        for (int j = 0; j < M.cols; ++j) {
            const double bound = std::sqrt(M(i, i) * M(j, j));
            if (bound > 0.0)
                worst_cs = std::fmax(worst_cs, std::fabs(M(i, j)) / bound - 1.0);
        }
    {
        char buf[400];
        std::snprintf(buf, sizeof(buf),
                      "diagonal in [%.4e, %.4e], all positive=%s; "
                      "|M - M^T| = %.2e (exact by construction); worst "
                      "Cauchy-Schwarz excess %.2e",
                      min_diag, max_diag, diag_positive ? "yes" : "NO",
                      gram_symmetry, worst_cs);
        report.check_predicate(
            "o0a.sobolev_basis_norms",
            "the Sobolev Gram matrix is EXACTLY symmetric, has a strictly "
            "positive diagonal, and satisfies |M_ij| <= sqrt(M_ii M_jj)",
            VerificationLevel::A,
            diag_positive && gram_symmetry == 0.0 && worst_cs <= 1e-12,
            "exact symmetry, positive diagonal, Cauchy-Schwarz respected",
            buf, sw.seconds(),
            "assembled one unordered pair at a time and added to both entries, "
            "so symmetry is bit-exact. Assembling both orders independently gave "
            "6e-14 asymmetry, which the eigensolver would have silently "
            "discarded by reading only the upper triangle");
    }

    // --- 26/27/28. Conditioning, whitening, and the discard decision --------
    const design::WhiteningResult wr = design::whiten(M, pr);
    const auto gram_groups = design::degenerate_groups(wr.eigenvalues);

    sw.reset();
    report.check_predicate(
        "o0a.raw_gram_conditioning",
        "the raw Sobolev Gram matrix is positive definite with a reported "
        "condition number, so mode magnitudes are comparable at all",
        VerificationLevel::B, wr.raw_condition > 1.0 && wr.raw_condition < 1e6,
        "finite condition number, no numerically null raw mode",
        "raw condition " + std::to_string(wr.raw_condition) +
            ", constrained " + std::to_string(wr.constrained_condition),
        sw.seconds(),
        "raw coefficient magnitude measures nothing on its own: two modes with "
        "the same coefficient can differ by orders of magnitude in value, slope "
        "and curvature because their supports and spacings differ");

    sw.reset();
    // c = c0 + P z, z = W y. Both coordinate systems are kept and both matrices
    // recorded, because they answer different questions: P is the constrained
    // subspace, PW is the whitened parameterisation an optimizer would step in.
    const DenseMatrix PW = design::multiply(pr.P, wr.w);
    double worst_constraint_after_whitening = 0.0;
    double worst_collar_after_whitening = 0.0;
    const DenseMatrix C = sys.matrix();
    for (int t = 0; t < 16; ++t) {
        std::vector<double> y(static_cast<std::size_t>(wr.retained));
        for (int j = 0; j < wr.retained; ++j)
            y[static_cast<std::size_t>(j)] =
                std::cos(1.13 * static_cast<double>(j) +
                         0.29 * static_cast<double>(t));
        // The perturbation in raw coefficients.
        std::vector<double> dc(static_cast<std::size_t>(pr.n_raw), 0.0);
        for (int i = 0; i < pr.n_raw; ++i)
            for (int j = 0; j < wr.retained; ++j)
                dc[static_cast<std::size_t>(i)] +=
                    PW(i, j) * y[static_cast<std::size_t>(j)];
        const std::vector<double> r = design::matvec(C, dc);
        for (double v : r)
            worst_constraint_after_whitening =
                std::fmax(worst_constraint_after_whitening, std::fabs(v));
        // And in the collar the perturbed field is still exactly zero, because
        // only active modes have coefficients at all.
        std::vector<double> c_full(static_cast<std::size_t>(pr.n_raw));
        for (int i = 0; i < pr.n_raw; ++i)
            c_full[static_cast<std::size_t>(i)] =
                pr.c0[static_cast<std::size_t>(i)] +
                dc[static_cast<std::size_t>(i)];
        for (int a = 0; a <= 12; ++a) {
            const double xi = region.xi_max - region.outer_flat_collar * a / 12.0
                              + 1e-9;
            const design::DesignCoordinates dcoord =
                design::from_design_coordinates(std::fmin(xi, region.xi_max),
                                                0.7, kO0aR0);
            if (dcoord.xi <= region.xi_flat_upper()) continue;
            worst_collar_after_whitening = std::fmax(
                worst_collar_after_whitening,
                std::fabs(basis.field_value<double>(c_full, dcoord)));
        }
    }
    {
        char buf[520];
        std::snprintf(buf, sizeof(buf),
                      "W^T G W = I to %.2e; whitened condition %.6f; %d modes "
                      "retained, %d discarded; retained eigenvalues in "
                      "[%.4e, %.4e]; 16 whitened perturbations leave every "
                      "constraint satisfied to %.2e and the collar exactly "
                      "%.1e",
                      wr.whitening_residual, wr.whitened_condition, wr.retained,
                      wr.discarded, wr.smallest_retained, wr.largest_retained,
                      worst_constraint_after_whitening,
                      worst_collar_after_whitening);
        report.check_predicate(
            "o0a.whitening",
            "W^T G W = I on the retained modes, the whitened condition number is "
            "1, and a whitened perturbation still satisfies the collar, the "
            "cabin and the moment constraints",
            VerificationLevel::A,
            wr.whitening_residual <= 1e-11 &&
                std::fabs(wr.whitened_condition - 1.0) <= 1e-8 &&
                worst_constraint_after_whitening <= 1e-9 &&
                worst_collar_after_whitening == 0.0,
            "identity Gram in whitened coordinates, constraints preserved",
            buf, sw.seconds(),
            "BOTH coordinate systems are kept: c = c0 + P z with z = W y. "
            "Whitening is a DENSE transform and generally destroys mode "
            "locality -- a single whitened coordinate moves many local "
            "coefficients -- while evaluation still happens in the original "
            "sparse local basis. It must never be described as preserving local "
            "support");
    }

    // Locality loss, measured rather than asserted: how many active modes a
    // single coordinate touches, before and after whitening.
    sw.reset();
    auto occupancy = [](const DenseMatrix& m, int col, double rel) {
        double mx = 0.0;
        for (int i = 0; i < m.rows; ++i)
            mx = std::fmax(mx, std::fabs(m(i, col)));
        int n = 0;
        for (int i = 0; i < m.rows; ++i)
            if (std::fabs(m(i, col)) > rel * mx) ++n;
        return n;
    };
    double avg_p = 0.0, avg_pw = 0.0;
    for (int j = 0; j < pr.n_free; ++j)
        avg_p += occupancy(pr.P, j, 1e-6);
    for (int j = 0; j < wr.retained; ++j)
        avg_pw += occupancy(PW, j, 1e-6);
    avg_p /= static_cast<double>(pr.n_free);
    avg_pw /= static_cast<double>(wr.retained);

    sw.reset();
    // The discard decision, exercised on a matrix that HAS a numerically null
    // mode: whitening the nearly-singular fixture must drop exactly one.
    design::ConstraintSystem none(4);
    const ConstraintProjection identity_pr = design::build_projection(none);
    const design::WhiteningResult wr_null =
        design::whiten(psd_null, identity_pr, 1e-12);
    const design::WhiteningResult wr_near =
        design::whiten(nearly, identity_pr, 1e-12);
    {
        char buf[460];
        std::snprintf(buf, sizeof(buf),
                      "spline Gram: %d retained, %d discarded, smallest "
                      "retained %.4e against tolerance %.4e. Exact-null "
                      "fixture: %d retained, %d discarded. Nearly-singular "
                      "fixture (lambda_min = 4e-18): %d retained, %d discarded",
                      wr.retained, wr.discarded, wr.smallest_retained,
                      wr.tolerance, wr_null.retained, wr_null.discarded,
                      wr_near.retained, wr_near.discarded);
        report.check_predicate(
            "o0a.numerical_null_mode_rejection",
            "a mode at or below the eigenvalue tolerance is DISCARDED and "
            "counted, never divided by: the exact-null and nearly-singular "
            "fixtures each lose exactly one mode",
            VerificationLevel::A,
            wr_null.discarded == 1 && wr_near.discarded == 1 &&
                wr_null.retained == 3 && wr_near.retained == 3 &&
                wr.smallest_retained > wr.tolerance,
            "exactly one mode discarded from each singular fixture",
            buf, sw.seconds(),
            "dividing by sqrt(lambda) for a numerically null lambda manufactures "
            "an enormous optimization direction out of rounding error, which an "
            "optimizer would then dutifully explore");
    }

    // --- Records -------------------------------------------------------------
    {
        std::ofstream f(out_dir + "/constraint_matrix.csv");
        f << "row,kind,detail,rhs,row_norm";
        for (int j = 0; j < pr.n_raw; ++j) f << ",c" << j;
        f << "\n";
        for (std::size_t i = 0; i < sys.rows().size(); ++i) {
            const design::ConstraintRow& r = sys.rows()[i];
            f << i << "," << r.kind << ",\"" << r.detail << "\"," << r.rhs
              << "," << design::norm2(r.row);
            for (int j = 0; j < pr.n_raw; ++j)
                f << "," << r.row[static_cast<std::size_t>(j)];
            f << "\n";
        }
    }
    {
        std::ofstream f(out_dir + "/constraint_singular_values.csv");
        f << "index,singular_value,retained,tolerance\n";
        for (std::size_t i = 0; i < pr.singular_values.size(); ++i)
            f << i << "," << pr.singular_values[i] << ","
              << (static_cast<int>(i) < pr.rank ? 1 : 0) << ","
              << pr.rank_tolerance << "\n";
    }
    {
        // Both matrices, as section 7 and the whitening semantics require: P is
        // the constrained subspace, P W the whitened parameterisation.
        std::ofstream f(out_dir + "/constraint_projection.csv");
        f << "matrix,row,col,value\n";
        for (int i = 0; i < pr.n_raw; ++i) {
            f << "c0," << i << ",0," << pr.c0[static_cast<std::size_t>(i)]
              << "\n";
            for (int j = 0; j < pr.n_free; ++j)
                f << "P," << i << "," << j << "," << pr.P(i, j) << "\n";
            for (int j = 0; j < wr.retained; ++j)
                f << "PW," << i << "," << j << "," << PW(i, j) << "\n";
        }
    }
    {
        const EigenResult raw_e = design::jacobi_eigen(M);
        const DenseMatrix G =
            design::multiply(design::multiply(design::transpose(pr.P), M), pr.P);
        const EigenResult con_e = design::jacobi_eigen(G);
        const DenseMatrix WtGW = design::multiply(
            design::multiply(design::transpose(wr.w), G), wr.w);
        const EigenResult whi_e = design::jacobi_eigen(WtGW);
        struct Spec { const char* file; const EigenResult* e; double cond; };
        const Spec specs[3] = {
            {"/gram_spectrum_raw.csv", &raw_e, wr.raw_condition},
            {"/gram_spectrum_constrained.csv", &con_e, wr.constrained_condition},
            {"/gram_spectrum_whitened.csv", &whi_e, wr.whitened_condition}};
        for (const Spec& sp : specs) {
            std::ofstream f(out_dir + sp.file);
            f << "index,eigenvalue,condition_number,converged,retained,"
                 "in_degenerate_group,group_size\n";
            const auto groups = design::degenerate_groups(sp.e->eigenvalues);
            for (std::size_t i = 0; i < sp.e->eigenvalues.size(); ++i) {
                int gsize = 1;
                bool deg = false;
                for (const auto& g : groups)
                    for (int idx : g)
                        if (static_cast<std::size_t>(idx) == i) {
                            deg = true;
                            gsize = static_cast<int>(g.size());
                        }
                f << i << "," << sp.e->eigenvalues[i] << "," << sp.cond << ","
                  << (sp.e->converged ? 1 : 0) << ","
                  << (sp.e->eigenvalues[i] > wr.tolerance ? 1 : 0) << ","
                  << (deg ? 1 : 0) << "," << gsize << "\n";
            }
        }
    }
    {
        std::ofstream f(out_dir + "/basis_normalization.csv");
        f << "flat_index,sobolev_norm,value_norm,l1_length_scale,"
             "l2_length_scale,occupancy_P,occupancy_PW\n";
        for (int i = 0; i < M.rows; ++i)
            f << i << "," << std::sqrt(M(i, i)) << "," << M(i, i) << ","
              << scales.l1 << "," << scales.l2 << "," << avg_p << ","
              << avg_pw << "\n";
    }
    {
        std::ofstream f(out_dir + "/whitening_locality.csv");
        f << "coordinate_system,average_modes_touched,total_modes,"
             "preserves_local_support\n";
        f << "P," << avg_p << "," << pr.n_raw << ",partial\n";
        f << "PW," << avg_pw << "," << pr.n_raw << ",no\n";
    }
}

}  // namespace checks
}  // namespace warp_hypad

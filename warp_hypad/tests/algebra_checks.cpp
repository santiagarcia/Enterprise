#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

#include "checks.hpp"
#include "reference/ref_linalg.hpp"
#include "warp_hypad/core/mat4.hpp"

namespace warp_hypad {
namespace checks {
namespace {

using reference::RefMat4;

Mat4<double> from_ref(const RefMat4& r) {
    Mat4<double> m;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) m(i, j) = r(i, j);
    return m;
}

Mat4<double> from_rows(const double (&v)[4][4]) {
    Mat4<double> m;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) m(i, j) = v[i][j];
    return m;
}

double max_diff_ref(const Mat4<double>& a, const RefMat4& b) {
    double worst = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            const double d = std::fabs(a(i, j) - b(i, j));
            if (d > worst) worst = d;
        }
    return worst;
}

// Seed is fixed and recorded in the report notes so any failure is reproducible.
constexpr std::uint64_t kSeed = 0x5741525046414354ULL;  // "WARPFACT"
constexpr int kRandomCases = 500;

struct ResidualSample {
    std::string family;
    double cond1;
    double residual_left;
    double residual_right;
    double det_rel_error;
};

}  // namespace

void run_algebra_checks(VerificationReport& report, const std::string& out_dir) {
    std::vector<ResidualSample> samples;

    // -----------------------------------------------------------------------
    // Matrix multiplication against a hand-calculated result.
    //
    // D = diag(1,2,3,4), B = the 1..16 matrix in row-major order.
    // D*B scales row i by D_ii; B*D scales column j by D_jj. The two differ,
    // so this single case also catches a row/column transposition -- which a
    // symmetric or identity test case could not.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const double dvals[4][4] = {{1, 0, 0, 0}, {0, 2, 0, 0}, {0, 0, 3, 0}, {0, 0, 0, 4}};
        const double bvals[4][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        const Mat4<double> d = from_rows(dvals);
        const Mat4<double> b = from_rows(bvals);

        const double db_expected[4][4] = {
            {1, 2, 3, 4}, {10, 12, 14, 16}, {27, 30, 33, 36}, {52, 56, 60, 64}};
        const double bd_expected[4][4] = {
            {1, 4, 9, 16}, {5, 12, 21, 32}, {9, 20, 33, 48}, {13, 28, 45, 64}};

        const Mat4<double> db = mat_mul(d, b);
        const Mat4<double> bd = mat_mul(b, d);

        double worst = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                worst = std::fmax(worst, std::fabs(db(i, j) - db_expected[i][j]));
                worst = std::fmax(worst, std::fabs(bd(i, j) - bd_expected[i][j]));
            }

        report.check_residual(
            "algebra.matmul_hand_calculated",
            "diag(1,2,3,4)*B scales rows; B*diag(1,2,3,4) scales columns",
            VerificationLevel::A, "n/a", "n/a", worst, 0.0, sw.seconds(),
            "asymmetric by construction, so it also detects a transposed "
            "row/column convention");
    }

    // -----------------------------------------------------------------------
    // Matrix multiplication against the independent reference.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        reference::MatrixGenerator gen(kSeed);
        double worst = 0.0;
        for (int n = 0; n < kRandomCases; ++n) {
            const RefMat4 ra = gen.general();
            const RefMat4 rb = gen.general();
            const Mat4<double> a = from_ref(ra);
            const Mat4<double> b = from_ref(rb);
            worst = std::fmax(worst, max_diff_ref(mat_mul(a, b), reference::ref_mul(ra, rb)));
        }
        std::ostringstream note;
        note << kRandomCases << " random matrices, seed 0x" << std::hex << kSeed;
        report.check_residual(
            "algebra.matmul_vs_independent",
            "mat_mul agrees with an independently written naive triple loop",
            VerificationLevel::D, "n/a", "entries in [-2, 2]", worst, 1e-12,
            sw.seconds(), note.str());
    }

    // -----------------------------------------------------------------------
    // Determinant against known analytic values.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;

        // Minkowski: det diag(-1,1,1,1) = -1. The single most important
        // determinant in this project.
        worst = std::fmax(worst, std::fabs(det(Mat4<double>::minkowski()) - (-1.0)));

        // Identity.
        worst = std::fmax(worst, std::fabs(det(Mat4<double>::identity()) - 1.0));

        // Diagonal: 2*3*5*7 = 210.
        const double diag[4][4] = {{2, 0, 0, 0}, {0, 3, 0, 0}, {0, 0, 5, 0}, {0, 0, 0, 7}};
        worst = std::fmax(worst, std::fabs(det(from_rows(diag)) - 210.0));

        // Upper triangular: product of the diagonal = 1*2*3*4 = 24.
        const double tri[4][4] = {{1, 9, -3, 5}, {0, 2, 7, 1}, {0, 0, 3, -2}, {0, 0, 0, 4}};
        worst = std::fmax(worst, std::fabs(det(from_rows(tri)) - 24.0));

        // Two identical rows: exactly singular, det = 0.
        const double dup[4][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {1, 2, 3, 4}, {9, 1, 4, 2}};
        worst = std::fmax(worst, std::fabs(det(from_rows(dup)) - 0.0));

        report.check_residual(
            "algebra.det_analytic",
            "det of Minkowski (-1), identity (1), diag(2,3,5,7) (210), an "
            "upper-triangular matrix (24) and a rank-deficient matrix (0)",
            VerificationLevel::B, "n/a", "n/a", worst, 1e-12, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Determinant against the independent recursive Laplace expansion.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        reference::MatrixGenerator gen(kSeed + 1);
        double worst_rel = 0.0;
        for (int n = 0; n < kRandomCases; ++n) {
            const RefMat4 r = gen.general();
            const double mine = det(from_ref(r));
            const double theirs = reference::ref_det(r);
            const double scale = std::fmax(1.0, std::fabs(theirs));
            worst_rel = std::fmax(worst_rel, std::fabs(mine - theirs) / scale);
        }
        report.check_residual(
            "algebra.det_vs_independent",
            "the paired-2x2-minor determinant agrees with recursive Laplace "
            "expansion along row 0",
            VerificationLevel::D, "n/a", "entries in [-2, 2]", worst_rel, 1e-12,
            sw.seconds(),
            "different summation order and grouping, so a shared bug is "
            "unlikely");
    }

    // -----------------------------------------------------------------------
    // Inversion: A*A^-1 = I and A^-1*A = I across four matrix families.
    //
    // Both orderings are tested because they are not the same statement for a
    // buggy implementation -- a transposed adjugate satisfies neither, but a
    // subtly wrong one can satisfy one and not the other.
    // -----------------------------------------------------------------------
    {
        struct Family {
            const char* name;
            const char* description;
        };
        const Family families[4] = {
            {"diagonal", "diagonal matrices with mixed signs"},
            {"symmetric", "general symmetric matrices"},
            {"lorentzian", "perturbed Lorentzian metrics, signature -+++"},
            {"general", "general nonsingular matrices"},
        };

        for (int f = 0; f < 4; ++f) {
            Stopwatch sw;
            reference::MatrixGenerator gen(kSeed + 10 + static_cast<std::uint64_t>(f));
            double worst_left = 0.0;
            double worst_right = 0.0;
            int rejected = 0;

            for (int n = 0; n < kRandomCases; ++n) {
                RefMat4 r;
                switch (f) {
                    case 0: r = gen.diagonal(); break;
                    case 1: r = gen.symmetric(); break;
                    case 2: r = gen.lorentzian_like(); break;
                    default: r = gen.general(); break;
                }
                const Mat4<double> a = from_ref(r);
                const InverseResult<double> inv = inverse(a);
                if (!inv.ok) { ++rejected; continue; }

                const double left = identity_residual(a, inv.inverse);
                const double right = identity_residual(inv.inverse, a);
                worst_left = std::fmax(worst_left, left);
                worst_right = std::fmax(worst_right, right);

                RefMat4 ref_inv;
                double det_rel = 0.0;
                if (reference::ref_inverse(r, ref_inv)) {
                    const double theirs = reference::ref_det(r);
                    const double scale = std::fmax(1.0, std::fabs(theirs));
                    det_rel = std::fabs(inv.det_real - theirs) / scale;
                }
                samples.push_back({families[f].name, inv.cond1, left, right, det_rel});
            }

            std::ostringstream note;
            note << kRandomCases << " matrices, " << rejected
                 << " rejected as singular before the residual test";

            report.check_residual(
                std::string("algebra.inverse_identity_left [") + families[f].name + "]",
                "A * A^-1 == I",
                VerificationLevel::A, "n/a", families[f].description, worst_left,
                1e-9, sw.seconds(), note.str());

            report.check_residual(
                std::string("algebra.inverse_identity_right [") + families[f].name + "]",
                "A^-1 * A == I",
                VerificationLevel::A, "n/a", families[f].description, worst_right,
                1e-9, sw.seconds(),
                "tested separately from the left identity on purpose");
        }
    }

    // -----------------------------------------------------------------------
    // Inversion agrees with Gauss-Jordan with partial pivoting.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        reference::MatrixGenerator gen(kSeed + 20);
        double worst = 0.0;
        int compared = 0;
        for (int n = 0; n < kRandomCases; ++n) {
            const RefMat4 r = gen.lorentzian_like();
            const InverseResult<double> mine = inverse(from_ref(r));
            RefMat4 theirs;
            if (!mine.ok || !reference::ref_inverse(r, theirs)) continue;
            worst = std::fmax(worst, max_diff_ref(mine.inverse, theirs));
            ++compared;
        }
        std::ostringstream note;
        note << compared << " matrices compared; adjugate vs Gauss-Jordan with "
                            "partial pivoting";
        report.check_residual(
            "algebra.inverse_vs_independent",
            "the adjugate inverse agrees with Gauss-Jordan elimination",
            VerificationLevel::D, "n/a", "perturbed Lorentzian metrics", worst,
            1e-9, sw.seconds(), note.str());
    }

    // -----------------------------------------------------------------------
    // Minkowski is its own inverse in the (-,+,+,+) signature. Exact.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const Mat4<double> eta = Mat4<double>::minkowski();
        const InverseResult<double> inv = inverse(eta);
        const double worst = inv.ok ? max_abs_diff(inv.inverse, eta) : 1.0;
        report.check_residual(
            "algebra.minkowski_self_inverse",
            "eta^{mu nu} == eta_{mu nu} == diag(-1,1,1,1)",
            VerificationLevel::B, "n/a", "signature -+++", worst, 0.0,
            sw.seconds(), "exact: every entry is +-1 or 0, so no rounding occurs");
    }

    // -----------------------------------------------------------------------
    // Singular matrices are REJECTED, not silently inverted.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        reference::MatrixGenerator gen(kSeed + 30);
        int accepted = 0;
        const int trials = 200;
        for (int n = 0; n < trials; ++n) {
            if (inverse(from_ref(gen.singular())).ok) ++accepted;
        }
        std::ostringstream computed;
        computed << accepted << " of " << trials << " exactly singular matrices "
                 << "were accepted";
        report.check_predicate(
            "algebra.singular_detected",
            "inverse() reports ok == false for exactly singular matrices "
            "(row 3 = 2*row 0 - 3*row 1)",
            VerificationLevel::A, accepted == 0, "0 accepted", computed.str(),
            sw.seconds(),
            "silently returning garbage here would corrupt every downstream "
            "curvature value");
    }

    // -----------------------------------------------------------------------
    // Nearly singular matrices: either rejected, or accepted WITH an
    // ill-conditioning flag. What must never happen is quiet acceptance of a
    // result that has lost all its significant digits.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        reference::MatrixGenerator gen(kSeed + 40);
        int quietly_accepted = 0;
        int rejected = 0;
        int flagged = 0;
        const int trials = 200;
        for (int n = 0; n < trials; ++n) {
            const RefMat4 r = gen.near_singular(1e-13);
            const InverseResult<double> inv = inverse(from_ref(r));
            if (!inv.ok) { ++rejected; continue; }
            if (inv.ill_conditioned) { ++flagged; continue; }
            // Accepted without a flag: only acceptable if the residual is
            // genuinely good.
            if (identity_residual(from_ref(r), inv.inverse) > 1e-6) ++quietly_accepted;
        }
        std::ostringstream computed;
        computed << rejected << " rejected, " << flagged << " flagged "
                 << "ill-conditioned, " << quietly_accepted
                 << " accepted quietly with a bad residual";
        report.check_predicate(
            "algebra.near_singular_reported",
            "a near-singular matrix is either rejected or flagged "
            "ill-conditioned; it is never quietly accepted with a bad residual",
            VerificationLevel::A, quietly_accepted == 0,
            "0 quietly accepted with a bad residual", computed.str(),
            sw.seconds(),
            "perturbation 1e-13 on a rank-deficient matrix");
    }

    // -----------------------------------------------------------------------
    // The ill-conditioning flag must actually fire somewhere.
    //
    // At a perturbation of 1e-13 the determinant test rejects every matrix, so
    // the check above passes without ever exercising the ill_conditioned path.
    // A flag that is never set is not a working flag. This sweep walks the
    // perturbation across the range where the behaviour transitions from
    // "accepted and well conditioned" to "rejected as singular", and requires
    // that the intermediate regime -- accepted but flagged -- is non-empty.
    //
    // It simultaneously re-checks, at every perturbation scale, that nothing is
    // ever quietly accepted with a bad residual.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const double epsilons[6] = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
        int total_flagged = 0;
        int total_quietly_bad = 0;
        double max_cond_accepted = 0.0;
        std::ostringstream breakdown;

        for (int e = 0; e < 6; ++e) {
            reference::MatrixGenerator gen(kSeed + 60 + static_cast<std::uint64_t>(e));
            int accepted = 0, flagged = 0, rejected = 0, quietly_bad = 0;
            double worst_unflagged_residual = 0.0;
            double cond_at_worst = 0.0;
            for (int n = 0; n < 100; ++n) {
                const RefMat4 r = gen.near_singular(epsilons[e]);
                const InverseResult<double> inv = inverse(from_ref(r));
                if (!inv.ok) { ++rejected; continue; }
                ++accepted;
                max_cond_accepted = std::fmax(max_cond_accepted, inv.cond1);
                if (inv.ill_conditioned) {
                    ++flagged;
                    continue;
                }
                const double res = identity_residual(from_ref(r), inv.inverse);
                if (res > worst_unflagged_residual) {
                    worst_unflagged_residual = res;
                    cond_at_worst = inv.cond1;
                }
                if (res > 1e-6) ++quietly_bad;
            }
            total_flagged += flagged;
            total_quietly_bad += quietly_bad;
            if (e) breakdown << "; ";
            breakdown << "eps=" << epsilons[e] << ": " << accepted << " ok ("
                      << flagged << " flagged, " << quietly_bad << " quietly bad, "
                      << "worst unflagged residual " << worst_unflagged_residual
                      << " at cond1 " << cond_at_worst << "), " << rejected
                      << " rejected";
        }

        report.check_predicate(
            "algebra.ill_conditioning_flag_fires",
            "across perturbations from 1e-4 to 1e-14 the ill_conditioned flag "
            "is set for at least one accepted matrix, and no matrix at any "
            "scale is quietly accepted with a bad residual",
            VerificationLevel::A,
            total_flagged > 0 && total_quietly_bad == 0,
            "at least one flagged, zero quietly accepted with a bad residual",
            breakdown.str(), sw.seconds(),
            "without this, the near-singular check above passes vacuously "
            "because the determinant test rejects everything first");
    }

    // -----------------------------------------------------------------------
    // Conditioning diagnostics are self-consistent: cond1 == ||A||*||A^-1||,
    // and is never below 1 (a true lower bound for any induced norm).
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        reference::MatrixGenerator gen(kSeed + 50);
        double worst = 0.0;
        int below_one = 0;
        for (int n = 0; n < kRandomCases; ++n) {
            const InverseResult<double> inv = inverse(from_ref(gen.general()));
            if (!inv.ok) continue;
            const double expect = inv.norm1_a * inv.norm1_inv;
            const double scale = std::fmax(1.0, expect);
            worst = std::fmax(worst, std::fabs(inv.cond1 - expect) / scale);
            if (inv.cond1 < 1.0 - 1e-12) ++below_one;
        }
        report.check_residual(
            "algebra.condition_number_consistent",
            "cond1 == ||A||_1 * ||A^-1||_1, and cond1 >= 1",
            VerificationLevel::A, "n/a", "general matrices",
            worst + static_cast<double>(below_one), 1e-14, sw.seconds(),
            "adds the count of cond1 < 1 violations to the residual");
    }

    // -----------------------------------------------------------------------
    // Export the residual/conditioning samples for the Python figures.
    // -----------------------------------------------------------------------
    {
        std::error_code ec;
        std::filesystem::create_directories(out_dir, ec);
        std::ofstream os(out_dir + "/inverse_residuals.csv", std::ios::binary);
        if (os) {
            os << "family,cond1,residual_left,residual_right,det_rel_error\n";
            for (const ResidualSample& s : samples) {
                os << s.family << ',' << s.cond1 << ',' << s.residual_left << ','
                   << s.residual_right << ',' << s.det_rel_error << '\n';
            }
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad

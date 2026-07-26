#include <cmath>

#include "checks.hpp"
#include "warp_hypad/core/tensor.hpp"

namespace warp_hypad {
namespace checks {
namespace {

// Diagonal matrices make every contraction hand-computable in closed form
// while still being sensitive to index placement: a diagonal factor picks out
// exactly one term of the sum, so a transposed or mis-slotted index changes
// WHICH term survives, not merely its value.
Mat4<double> diag(double a, double b, double c, double d) {
    Mat4<double> m;
    m(0, 0) = a; m(1, 1) = b; m(2, 2) = c; m(3, 3) = d;
    return m;
}

}  // namespace

void run_contraction_checks(VerificationReport& report) {
    // --- trace_mixed: M^mu_mu ---
    {
        Stopwatch sw;
        const double got = trace_mixed(diag(1.0, 2.0, 3.0, 4.0));
        report.check_value(
            "contract.trace_mixed",
            "trace of diag(1,2,3,4) == 10",
            VerificationLevel::A, "n/a", "n/a", 10.0, got, 0.0, sw.seconds());
    }
    {
        Stopwatch sw;
        // Minkowski trace in the (-,+,+,+) signature: -1 + 1 + 1 + 1 = 2.
        const double got = trace_mixed(Mat4<double>::minkowski());
        report.check_value(
            "contract.trace_minkowski",
            "trace of diag(-1,1,1,1) == 2",
            VerificationLevel::B, "n/a", "signature -+++", 2.0, got, 0.0,
            sw.seconds(),
            "a sign error in the signature would give 4 or -2 here");
    }

    // --- full_contract: A^{mu nu} B_{mu nu} ---
    {
        Stopwatch sw;
        // A(i,j) = i + j, B(i,j) = 1. Sum over all 16 entries of (i+j)
        // = 4*sum_i i + 4*sum_j j = 4*6 + 4*6 = 48.
        Mat4<double> a, b;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                a(i, j) = static_cast<double>(i + j);
                b(i, j) = 1.0;
            }
        report.check_value(
            "contract.full_contract",
            "sum_{mu,nu} (mu + nu) * 1 == 48",
            VerificationLevel::A, "n/a", "n/a", 48.0, full_contract(a, b), 0.0,
            sw.seconds());
    }
    {
        Stopwatch sw;
        // g^{mu nu} g_{mu nu} = delta^mu_mu = 4 in four dimensions. The single
        // most useful sanity contraction in the whole pipeline.
        const Mat4<double> eta = Mat4<double>::minkowski();
        report.check_value(
            "contract.metric_full_trace",
            "eta^{mu nu} eta_{mu nu} == 4 (dimension of spacetime)",
            VerificationLevel::A, "n/a", "signature -+++", 4.0,
            full_contract(eta, eta), 0.0, sw.seconds());
    }

    // --- raise_first: T^mu_nu = g^{mu alpha} T_{alpha nu} ---
    {
        Stopwatch sw;
        // Raising with eta then lowering with eta must return the original.
        const Mat4<double> eta = Mat4<double>::minkowski();
        Mat4<double> t;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) t(i, j) = 1.0 + i * 4.0 + j;

        const Mat4<double> roundtrip = raise_first(eta, raise_first(eta, t));
        report.check_residual(
            "contract.raise_lower_roundtrip",
            "eta^{mu a} eta_{a b} T^{b}{}_{nu} == T^{mu}{}_{nu}",
            VerificationLevel::A, "n/a", "n/a", max_abs_diff(roundtrip, t), 0.0,
            sw.seconds(),
            "exact: eta is its own inverse and every entry is +-1");
    }

    // --- contract_first: the Christoffel assembly step ---
    {
        Stopwatch sw;
        // a = diag(2,3,5,7), d(s,m,n) = 100s + 10m + n.
        // Because a is diagonal, out(r,m,n) = a(r,r) * d(r,m,n).
        // So out(1,2,3) = 3 * (100 + 20 + 3) = 3 * 123 = 369.
        const Mat4<double> a = diag(2.0, 3.0, 5.0, 7.0);
        Rank3<double> d;
        for (int s = 0; s < 4; ++s)
            for (int m = 0; m < 4; ++m)
                for (int n = 0; n < 4; ++n)
                    d(s, m, n) = 100.0 * s + 10.0 * m + n;

        const Rank3<double> out = contract_first(a, d);

        double worst = 0.0;
        for (int r = 0; r < 4; ++r)
            for (int m = 0; m < 4; ++m)
                for (int n = 0; n < 4; ++n) {
                    const double expect = a(r, r) * d(r, m, n);
                    worst = std::fmax(worst, std::fabs(out(r, m, n) - expect));
                }

        report.check_residual(
            "contract.contract_first_diagonal",
            "for diagonal a, (a . d)(r,m,n) == a(r,r) * d(r,m,n) at all 64 "
            "components; spot value out(1,2,3) == 369",
            VerificationLevel::A, "n/a", "a = diag(2,3,5,7)", worst, 0.0,
            sw.seconds(),
            "this is the g^{rho sigma} contraction that produces "
            "Gamma^rho_{mu nu}");
    }
    {
        Stopwatch sw;
        const Mat4<double> a = diag(2.0, 3.0, 5.0, 7.0);
        Rank3<double> d;
        for (int s = 0; s < 4; ++s)
            for (int m = 0; m < 4; ++m)
                for (int n = 0; n < 4; ++n)
                    d(s, m, n) = 100.0 * s + 10.0 * m + n;
        report.check_value(
            "contract.contract_first_spot_value",
            "out(1,2,3) == 3 * 123 == 369",
            VerificationLevel::A, "n/a", "a = diag(2,3,5,7)", 369.0,
            contract_first(a, d)(1, 2, 3), 0.0, sw.seconds(),
            "hand-calculated single component");
    }

    // --- contract_13: the Riemann -> Ricci step ---
    {
        Stopwatch sw;
        // t(a,b,c,d) = (a == c) ? (10b + d) : 0.
        // out(s,n) = sum_r t(r,s,r,n) = sum_r (10s + n) = 4 * (10s + n).
        // So out(2,3) = 4 * 23 = 92.
        Rank4<double> t;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d)
                        t(a, b, c, d) = (a == c) ? (10.0 * b + d) : 0.0;

        const Mat4<double> out = contract_13(t);

        double worst = 0.0;
        for (int s = 0; s < 4; ++s)
            for (int n = 0; n < 4; ++n) {
                const double expect = 4.0 * (10.0 * s + n);
                worst = std::fmax(worst, std::fabs(out(s, n) - expect));
            }

        report.check_residual(
            "contract.contract_13",
            "out(s,n) == sum_r t(r,s,r,n); with the constructed t this is "
            "4*(10s+n), so out(2,3) == 92",
            VerificationLevel::A, "n/a", "n/a", worst, 0.0, sw.seconds(),
            "this is R_{sigma nu} = R^rho{}_{sigma rho nu}");
    }
    {
        Stopwatch sw;
        // A rank-4 that is nonzero only where the contracted slots differ must
        // contract to exactly zero. This distinguishes contract_13 from the
        // other five possible index pairings.
        Rank4<double> t;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d)
                        t(a, b, c, d) = (a == c) ? 0.0 : (1.0 + a + b + c + d);

        report.check_residual(
            "contract.contract_13_slot_selectivity",
            "a rank-4 vanishing wherever slot 0 == slot 2 contracts to zero "
            "under contract_13",
            VerificationLevel::A, "n/a", "n/a", max_abs_real(contract_13(t)),
            0.0, sw.seconds(),
            "fails loudly if any other index pairing were contracted");
    }

    // --- lower_first ---
    {
        Stopwatch sw;
        const Mat4<double> g = diag(2.0, 3.0, 5.0, 7.0);
        Rank4<double> t;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d)
                        t(a, b, c, d) = 1.0 + a + b + c + d;

        const Rank4<double> out = lower_first(g, t);

        double worst = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    for (int d = 0; d < 4; ++d) {
                        const double expect = g(a, a) * t(a, b, c, d);
                        worst = std::fmax(worst, std::fabs(out(a, b, c, d) - expect));
                    }

        report.check_residual(
            "contract.lower_first_diagonal",
            "for diagonal g, (g . t)(a,b,c,d) == g(a,a) * t(a,b,c,d) at all "
            "256 components",
            VerificationLevel::A, "n/a", "g = diag(2,3,5,7)", worst, 0.0,
            sw.seconds(),
            "needed in Milestone 5 to test Riemann pair symmetry with all "
            "indices lowered");
    }

    // --- Symmetry residual helpers report what they claim to ---
    {
        Stopwatch sw;
        Rank3<double> sym;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    sym(a, b, c) = static_cast<double>(a) + b * c;  // symmetric in b,c

        Rank3<double> asym;
        asym(0, 1, 2) = 1.0;  // one antisymmetric entry, magnitude 1

        const double sym_res = lower_pair_symmetry_residual(sym);
        const double asym_res = lower_pair_symmetry_residual(asym);

        report.check_predicate(
            "contract.symmetry_detector_works",
            "lower_pair_symmetry_residual returns 0 for a b<->c symmetric "
            "tensor and 1 for a deliberately antisymmetric one",
            VerificationLevel::A,
            sym_res == 0.0 && std::fabs(asym_res - 1.0) < 1e-15,
            "0 for symmetric, 1 for the planted asymmetry",
            "symmetric: " + std::to_string(sym_res) +
                ", planted: " + std::to_string(asym_res),
            sw.seconds(),
            "a residual function that always returns 0 would make Milestone 4 "
            "pass vacuously");
    }

    // --- symmetry_residual on Mat4, same argument ---
    {
        Stopwatch sw;
        Mat4<double> sym;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) sym(i, j) = static_cast<double>(i * j + i + j);
        Mat4<double> asym = sym;
        asym(0, 2) += 3.0;

        const double s = symmetry_residual(sym);
        const double a = symmetry_residual(asym);

        report.check_predicate(
            "contract.mat4_symmetry_detector_works",
            "symmetry_residual returns 0 for a symmetric matrix and 3 for a "
            "planted asymmetry of 3",
            VerificationLevel::A, s == 0.0 && std::fabs(a - 3.0) < 1e-15,
            "0 for symmetric, 3 for the planted asymmetry",
            "symmetric: " + std::to_string(s) + ", planted: " + std::to_string(a),
            sw.seconds());
    }
}

}  // namespace checks
}  // namespace warp_hypad

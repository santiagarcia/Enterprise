#include <cmath>
#include <random>
#include <sstream>

#include "b1c_checks.hpp"
#include "program_b/effective_metric/hodge_dual.hpp"

// ---------------------------------------------------------------------------
// B1C step 1: four-dimensional Hodge duals.
//
// Everything downstream -- the Bel extraction (A7a-d), the magnetoelectric
// benchmark (B17), and the Z reconstruction (A11) -- is built on these, so they
// are verified in isolation first and nothing proceeds until they pass.
// ---------------------------------------------------------------------------

namespace program_b {
namespace checks {

using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace {

// A tensor carrying Z's antisymmetries: antisymmetric in (ab) and in (cd).
// Built from a fixed seed so a failure is reproducible.
Rank4<double> random_Z_symmetry(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> d(-1.5, 1.5);
    double raw[4][4][4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int e = 0; e < 4; ++e) raw[a][b][c][e] = d(rng);

    Rank4<double> Z;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int e = 0; e < 4; ++e)
                    Z(a, b, c, e) = 0.25 * (raw[a][b][c][e] - raw[b][a][c][e] -
                                            raw[a][b][e][c] + raw[b][a][e][c]);
    return Z;
}

// A non-trivial laboratory background: not Minkowski, not diagonal, so the
// sqrt(-g) factor and the index raising actually do something. Any
// hard-coded Minkowski Levi-Civita symbol would fail here.
Mat4<double> curved_lab() {
    Mat4<double> g;
    g(0, 0) = -1.7;
    g(1, 1) = 1.3;  g(2, 2) = 0.8;  g(3, 3) = 2.1;
    g(0, 1) = 0.25; g(1, 0) = 0.25;
    g(2, 3) = -0.4; g(3, 2) = -0.4;
    return g;
}

}  // namespace

void run_dual_checks(VerificationReport& report) {
    const Mat4<double> eta = Mat4<double>::minkowski();
    const Mat4<double> lab = curved_lab();
    const Rank4<double> Z = random_Z_symmetry(0x5A17C0DEULL);

    // -----------------------------------------------------------------------
    // Levi-Civita tensor normalisation, both index positions.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Mat4<double>& g : {eta, lab}) {
            const double root = std::sqrt(-warp_hypad::det(g));
            const Rank4<double> el = levi_civita_lower(g);
            const Rank4<double> eu = levi_civita_upper(g);
            worst = std::fmax(worst, std::fabs(el(0, 1, 2, 3) - root));
            worst = std::fmax(worst, std::fabs(eu(0, 1, 2, 3) + 1.0 / root));
            // Total antisymmetry: one odd swap must flip the sign.
            worst = std::fmax(worst, std::fabs(el(1, 0, 2, 3) + el(0, 1, 2, 3)));
            worst = std::fmax(worst, std::fabs(eu(0, 1, 3, 2) + eu(0, 1, 2, 3)));
            // Repeated index must vanish.
            worst = std::fmax(worst, std::fabs(el(0, 0, 2, 3)));
        }
        report.check_residual(
            "b1c.dual.levi_civita_normalisation",
            "eps_{0123} = +sqrt(-det g) and eps^{0123} = -1/sqrt(-det g), both "
            "totally antisymmetric, for Minkowski AND a curved non-diagonal lab",
            VerificationLevel::A, "n/a", "project convention, section 11", worst,
            1e-14, sw.seconds(),
            "the contravariant sign is forced by the (-,+,+,+) signature, not "
            "chosen; the curved lab would catch a hard-coded Minkowski symbol");
    }

    // -----------------------------------------------------------------------
    // eps_{abcd} eps^{cdef} = -2 (d_a^e d_b^f - d_a^f d_b^e)
    //
    // This identity IS the -1 in the double dual, so it is checked directly.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Mat4<double>& g : {eta, lab}) {
            const Rank4<double> el = levi_civita_lower(g);
            const Rank4<double> eu = levi_civita_upper(g);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int e = 0; e < 4; ++e)
                        for (int f = 0; f < 4; ++f) {
                            double acc = 0.0;
                            for (int c = 0; c < 4; ++c)
                                for (int d = 0; d < 4; ++d)
                                    acc += el(a, b, c, d) * eu(c, d, e, f);
                            const double want =
                                -2.0 * ((a == e ? 1.0 : 0.0) * (b == f ? 1.0 : 0.0) -
                                        (a == f ? 1.0 : 0.0) * (b == e ? 1.0 : 0.0));
                            worst = std::fmax(worst, std::fabs(acc - want));
                        }
        }
        report.check_residual(
            "b1c.dual.contraction_identity",
            "eps_{abcd} eps^{cdef} = -2(delta_a^e delta_b^f - delta_a^f delta_b^e)",
            VerificationLevel::A, "n/a", "Minkowski and curved lab", worst,
            1e-13, sw.seconds(),
            "the minus sign here is the source of **F = -F; checking it "
            "directly means the double-dual result is explained, not just "
            "observed");
    }

    // -----------------------------------------------------------------------
    // **F = -F for a two-form. The sign the correction task asks to confirm.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        std::mt19937_64 rng(12345);
        std::uniform_real_distribution<double> d(-1.0, 1.0);

        for (const Mat4<double>& g : {eta, lab}) {
            Mat4<double> F;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b) {
                    const double v = d(rng);
                    F(a, b) = v;
                    F(b, a) = -v;
                }
            const Mat4<double> FF = dual_2form(dual_2form(F, g), g);
            for (int i = 0; i < 16; ++i)
                worst = std::fmax(worst, std::fabs(FF.m[i] + F.m[i]));
        }
        report.check_residual(
            "b1c.dual.two_form_double_dual_is_minus_identity",
            "**F = -F for an arbitrary two-form, in four-dimensional Lorentzian "
            "signature",
            VerificationLevel::A, "n/a", "Minkowski and curved lab", worst,
            1e-13, sw.seconds(),
            "CONFIRMED SIGN: minus. In Riemannian signature it would be +F, so "
            "this pins the signature-dependent branch of the convention");
    }

    // -----------------------------------------------------------------------
    // Left dual twice, right dual twice, both = -Z.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst_l = 0.0, worst_r = 0.0;
        for (const Mat4<double>& g : {eta, lab}) {
            worst_l = std::fmax(worst_l,
                                rank4_max_abs_diff(dual_left(dual_left(Z, g), g),
                                                   rank4_negate(Z)));
            worst_r = std::fmax(worst_r,
                                rank4_max_abs_diff(dual_right(dual_right(Z, g), g),
                                                   rank4_negate(Z)));
        }
        report.check_residual(
            "b1c.dual.left_dual_squared_is_minus_identity",
            "*(*Z) = -Z", VerificationLevel::A, "n/a",
            "random Z with Z's antisymmetries; Minkowski and curved lab",
            worst_l, 1e-12, sw.seconds());
        report.check_residual(
            "b1c.dual.right_dual_squared_is_minus_identity",
            "(Z*)* = -Z", VerificationLevel::A, "n/a",
            "random Z with Z's antisymmetries; Minkowski and curved lab",
            worst_r, 1e-12, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Left and right duals act on disjoint index pairs, so they commute.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Mat4<double>& g : {eta, lab})
            worst = std::fmax(worst,
                              rank4_max_abs_diff(dual_right(dual_left(Z, g), g),
                                                 dual_left(dual_right(Z, g), g)));
        report.check_residual(
            "b1c.dual.left_right_commute",
            "*(Z*) = (*Z)*, since the two duals contract disjoint index pairs",
            VerificationLevel::A, "n/a", "Minkowski and curved lab", worst,
            1e-12, sw.seconds());
    }

    // Double dual squared: (*Z*) applied twice returns +Z, since each of the
    // two duals contributes -1.
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Mat4<double>& g : {eta, lab})
            worst = std::fmax(worst,
                              rank4_max_abs_diff(dual_both(dual_both(Z, g), g), Z));
        report.check_residual(
            "b1c.dual.double_dual_squared_is_identity",
            "(*Z*) applied twice returns +Z: (-1) from each index pair",
            VerificationLevel::A, "n/a", "Minkowski and curved lab", worst,
            1e-12, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Vacuum analytic case: the dual of the vacuum Z, computed in closed form.
    //
    // For Z^{abcd} = (1/2)(g^{ac}g^{bd} - g^{ad}g^{bc}), the left dual is
    // (1/2) eps^{ab}_{~~ef} Z^{efcd} = (1/2) eps^{abcd}, which is checkable
    // against the Levi-Civita tensor directly.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Mat4<double>& g : {eta, lab}) {
            const warp_hypad::InverseResult<double> gi = warp_hypad::inverse(g);
            Rank4<double> Zv;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int c = 0; c < 4; ++c)
                        for (int d = 0; d < 4; ++d)
                            Zv(a, b, c, d) = 0.5 * (gi.inverse(a, c) * gi.inverse(b, d) -
                                                    gi.inverse(a, d) * gi.inverse(b, c));
            const Rank4<double> got = dual_left(Zv, g);
            const Rank4<double> eu = levi_civita_upper(g);
            for (int i = 0; i < 256; ++i)
                worst = std::fmax(worst, std::fabs(got.t[i] - 0.5 * eu.t[i]));
        }
        report.check_residual(
            "b1c.dual.vacuum_analytic",
            "the left dual of the vacuum Z equals (1/2) eps^{abcd}",
            VerificationLevel::B, "n/a", "Minkowski and curved lab", worst,
            1e-13, sw.seconds(),
            "an independent closed form, not a self-consistency identity");
    }

    // -----------------------------------------------------------------------
    // Coordinate-transformation behaviour.
    //
    // eps is a genuine tensor once the sqrt(-g) is included, so the dual must
    // transform tensorially. A dual built from a bare permutation symbol would
    // fail this the moment det L != 1 -- which is why L below has det != 1.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        Mat4<double> L = Mat4<double>::identity();
        const double c = std::cos(0.5), s = std::sin(0.5);
        L(1, 1) = c;  L(1, 2) = -s;
        L(2, 1) = s;  L(2, 2) = c;
        L(3, 3) = 1.4;   // det L = 1.4, so the density factor is exercised
        L(0, 0) = 1.15;

        const warp_hypad::InverseResult<double> Li = warp_hypad::inverse(L);

        // g'_{ab} = (L^-1)^p_a (L^-1)^q_b g_{pq}
        Mat4<double> gp;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                double acc = 0.0;
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q)
                        acc += Li.inverse(p, a) * Li.inverse(q, b) * lab(p, q);
                gp(a, b) = acc;
            }

        // Z'^{abcd} = L^a_p L^b_q L^c_r L^d_s Z^{pqrs}
        auto push4 = [&](const Rank4<double>& T) {
            Rank4<double> o;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int cc = 0; cc < 4; ++cc)
                        for (int d = 0; d < 4; ++d) {
                            double acc = 0.0;
                            for (int p = 0; p < 4; ++p)
                                for (int q = 0; q < 4; ++q)
                                    for (int r = 0; r < 4; ++r)
                                        for (int t = 0; t < 4; ++t)
                                            acc += L(a, p) * L(b, q) * L(cc, r) *
                                                   L(d, t) * T(p, q, r, t);
                            o(a, b, cc, d) = acc;
                        }
            return o;
        };

        const double worst =
            rank4_max_abs_diff(dual_left(push4(Z), gp), push4(dual_left(Z, lab)));

        report.check_residual(
            "b1c.dual.transforms_tensorially",
            "*(L.Z) computed in the new coordinates equals L.(*Z) computed in "
            "the old ones",
            VerificationLevel::A, "n/a", "linear map with det L = 1.61", worst,
            1e-12, sw.seconds(),
            "det L != 1 deliberately: a dual built from a bare permutation "
            "symbol rather than the Levi-Civita TENSOR would fail exactly here");
    }
}

}  // namespace checks
}  // namespace program_b

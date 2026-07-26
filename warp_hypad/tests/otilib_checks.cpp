#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "otilib_checks.hpp"
#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/tensor.hpp"
#include "warp_hypad/sensitivity/otilib_traits.hpp"

namespace warp_hypad {
namespace checks {
namespace {

using sensitivity::extract_full_gradient;
using sensitivity::extract_full_hessian;
using sensitivity::extract_gradient;
using sensitivity::extract_hessian;
using sensitivity::OTIScalar;
using sensitivity::seed_parameter;

// Evaluation point for the composite tests. Chosen away from any special value
// (no zeros, no integer arguments) so a coincidental agreement is unlikely.
constexpr double kA = 1.30;
constexpr double kB = 0.70;
constexpr double kC = 0.45;

struct DerivRecord {
    std::string function;
    std::string quantity;
    double analytic;
    double oti;
    double abs_error;
};

std::vector<DerivRecord>& records() {
    static std::vector<DerivRecord> r;
    return r;
}

void record(const std::string& fn, const std::string& q, double analytic,
            double oti) {
    records().push_back({fn, q, analytic, oti, std::fabs(oti - analytic)});
}

// --- Elementary function checks -------------------------------------------
//
// Each seeds a single direction, applies one function, and compares the first
// and second derivative against the closed form. Second derivatives matter
// here: order-2 OTI propagation is the part that a first-order AD library would
// get wrong, so testing only gradients would not distinguish the two.

double check_unary(VerificationReport& report, const char* name,
                   OTIScalar (*f)(const OTIScalar&), double x, double d1,
                   double d2) {
    Stopwatch sw;
    const OTIScalar sx = seed_parameter(x, 0);
    const OTIScalar r = f(sx);

    const double g = extract_gradient(r, 0);
    const double h = extract_hessian(r, 0, 0);

    record(name, "d/dx", d1, g);
    record(name, "d2/dx2", d2, h);

    const double err = std::fmax(std::fabs(g - d1), std::fabs(h - d2));
    const double scale = std::fmax(1.0, std::fmax(std::fabs(d1), std::fabs(d2)));

    std::ostringstream note;
    note << "x=" << x << ", d/dx=" << d1 << ", d2/dx2=" << d2;
    report.check_residual(std::string("oti.elementary.") + name,
                          "first and second derivatives against the closed form",
                          VerificationLevel::B, "n/a", note.str(), err / scale,
                          1e-12, sw.seconds());
    return err;
}

// f(a,b,c) = a^2 sin(b) + exp(a c) -- the function named in the project spec.
template <typename S>
S composite(const S& a, const S& b, const S& c) {
    using std::exp;
    using std::sin;
    // Written with no reuse of a temporary as both input and output; see the
    // aliasing rule in otilib_scalar.hpp.
    const S a2 = a * a;
    const S sb = sin(b);
    const S term1 = a2 * sb;
    const S ac = a * c;
    const S term2 = exp(ac);
    return term1 + term2;
}

}  // namespace

void run_otilib_checks(VerificationReport& report, const std::string& out_dir) {
    // -----------------------------------------------------------------------
    // Elementary functions.
    // -----------------------------------------------------------------------
    {
        const double x = 1.3;
        check_unary(report, "sqrt", [](const OTIScalar& v) { using sensitivity::sqrt; return sqrt(v); },
                    x, 0.5 / std::sqrt(x), -0.25 / (x * std::sqrt(x)));
        check_unary(report, "exp", [](const OTIScalar& v) { using sensitivity::exp; return exp(v); },
                    x, std::exp(x), std::exp(x));
        check_unary(report, "log", [](const OTIScalar& v) { using sensitivity::log; return log(v); },
                    x, 1.0 / x, -1.0 / (x * x));
        check_unary(report, "sin", [](const OTIScalar& v) { using sensitivity::sin; return sin(v); },
                    x, std::cos(x), -std::sin(x));
        check_unary(report, "cos", [](const OTIScalar& v) { using sensitivity::cos; return cos(v); },
                    x, -std::sin(x), -std::cos(x));
        check_unary(report, "tan", [](const OTIScalar& v) { using sensitivity::tan; return tan(v); },
                    x, 1.0 / (std::cos(x) * std::cos(x)),
                    2.0 * std::tan(x) / (std::cos(x) * std::cos(x)));
        check_unary(report, "sinh", [](const OTIScalar& v) { using sensitivity::sinh; return sinh(v); },
                    x, std::cosh(x), std::sinh(x));
        check_unary(report, "cosh", [](const OTIScalar& v) { using sensitivity::cosh; return cosh(v); },
                    x, std::sinh(x), std::cosh(x));
        check_unary(report, "tanh", [](const OTIScalar& v) { using sensitivity::tanh; return tanh(v); },
                    x, 1.0 - std::tanh(x) * std::tanh(x),
                    -2.0 * std::tanh(x) * (1.0 - std::tanh(x) * std::tanh(x)));
    }

    // Arithmetic: division is the one worth stating explicitly, since OTI
    // division is implemented via pow(-1) rather than as a primitive.
    {
        Stopwatch sw;
        const OTIScalar a = seed_parameter(kA, 0);
        const OTIScalar b = sensitivity::make_constant(kB);
        const OTIScalar q = a / b;
        const double g = extract_gradient(q, 0);
        const double h = extract_hessian(q, 0, 0);
        record("divide", "d/da (a/b)", 1.0 / kB, g);
        record("divide", "d2/da2 (a/b)", 0.0, h);
        report.check_residual(
            "oti.elementary.divide",
            "d(a/b)/da == 1/b and d2(a/b)/da2 == 0",
            VerificationLevel::B, "n/a", "a=1.3, b=0.7 constant",
            std::fmax(std::fabs(g - 1.0 / kB), std::fabs(h)), 1e-12, sw.seconds(),
            "OTI division is pow(-1) composed with multiply, not a primitive");
    }

    // -----------------------------------------------------------------------
    // Composite function: full gradient and Hessian from ONE evaluation.
    // -----------------------------------------------------------------------
    const double f_val = kA * kA * std::sin(kB) + std::exp(kA * kC);
    const double gx[3] = {2 * kA * std::sin(kB) + kC * std::exp(kA * kC),
                          kA * kA * std::cos(kB),
                          kA * std::exp(kA * kC)};
    const double hx[3][3] = {
        {2 * std::sin(kB) + kC * kC * std::exp(kA * kC), 2 * kA * std::cos(kB),
         std::exp(kA * kC) * (1.0 + kA * kC)},
        {2 * kA * std::cos(kB), -kA * kA * std::sin(kB), 0.0},
        {std::exp(kA * kC) * (1.0 + kA * kC), 0.0, kA * kA * std::exp(kA * kC)}};

    OTIScalar f_oti;
    {
        Stopwatch sw;
        const sensitivity::SeededParameters p = sensitivity::seed_all(kA, kB, kC);
        f_oti = composite<OTIScalar>(p.velocity, p.radius, p.sigma);

        report.check_value("oti.composite.value",
                           "f(a,b,c) = a^2 sin(b) + exp(ac), real part",
                           VerificationLevel::B, "n/a",
                           "a=1.3, b=0.7, c=0.45", f_val, f_oti.real(), 1e-14,
                           sw.seconds());
    }
    {
        Stopwatch sw;
        const sensitivity::Gradient g = extract_full_gradient(f_oti);
        double worst = 0.0;
        for (int i = 0; i < 3; ++i) {
            worst = std::fmax(worst, std::fabs(g[i] - gx[i]));
            record("composite", std::string("df/dp") + std::to_string(i), gx[i], g[i]);
        }
        report.check_residual(
            "oti.composite.gradient",
            "all three first derivatives of a^2 sin(b) + exp(ac)",
            VerificationLevel::B, "n/a", "a=1.3, b=0.7, c=0.45", worst, 1e-13,
            sw.seconds(),
            "three directions seeded simultaneously, one evaluation");
    }
    {
        Stopwatch sw;
        const sensitivity::Hessian h = extract_full_hessian(f_oti);
        double worst = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                worst = std::fmax(worst, std::fabs(h(i, j) - hx[i][j]));
                if (i <= j)
                    record("composite",
                           "d2f/dp" + std::to_string(i) + "dp" + std::to_string(j),
                           hx[i][j], h(i, j));
            }
        report.check_residual(
            "oti.composite.hessian",
            "all six independent second derivatives of a^2 sin(b) + exp(ac)",
            VerificationLevel::B, "n/a", "a=1.3, b=0.7, c=0.45", worst, 1e-12,
            sw.seconds(),
            "reading the raw coefficients instead of get_deriv would fail this "
            "on the three diagonal entries only");
    }

    // The Hessian is symmetric by construction here (both (i,j) and (j,i) read
    // one packed coefficient), so symmetry alone proves nothing. Recorded
    // explicitly so nobody later mistakes it for evidence.
    {
        Stopwatch sw;
        const sensitivity::Hessian h = extract_full_hessian(f_oti);
        double worst = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                worst = std::fmax(worst, std::fabs(h(i, j) - h(j, i)));
        report.check_residual(
            "oti.hessian_symmetry", "H(i,j) == H(j,i)", VerificationLevel::A,
            "n/a", "n/a", worst, 0.0, sw.seconds(),
            "WEAK: symmetric by construction, since both entries read the same "
            "packed coefficient. Kept only to document that fact");
    }

    // -----------------------------------------------------------------------
    // Seed isolation: seeding one parameter must leave the others' channels 0.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (int seeded = 0; seeded < 3; ++seeded) {
            const OTIScalar a = (seeded == 0) ? seed_parameter(kA, 0) : OTIScalar(kA);
            const OTIScalar b = (seeded == 1) ? seed_parameter(kB, 1) : OTIScalar(kB);
            const OTIScalar c = (seeded == 2) ? seed_parameter(kC, 2) : OTIScalar(kC);
            const OTIScalar r = composite<OTIScalar>(a, b, c);
            for (int d = 0; d < 3; ++d) {
                if (d == seeded) continue;
                worst = std::fmax(worst, std::fabs(extract_gradient(r, d)));
            }
        }
        report.check_residual(
            "oti.seed_isolation",
            "seeding direction k leaves every other derivative channel exactly 0",
            VerificationLevel::A, "n/a", "3 runs, one seed each", worst, 0.0,
            sw.seconds(),
            "catches an off-by-one in the seed index, which would otherwise "
            "produce plausible but wrong sensitivities");
    }

    // Simultaneous seeding must equal separate runs.
    {
        Stopwatch sw;
        double worst = 0.0;
        for (int d = 0; d < 3; ++d) {
            const OTIScalar a = (d == 0) ? seed_parameter(kA, 0) : OTIScalar(kA);
            const OTIScalar b = (d == 1) ? seed_parameter(kB, 1) : OTIScalar(kB);
            const OTIScalar c = (d == 2) ? seed_parameter(kC, 2) : OTIScalar(kC);
            const OTIScalar single = composite<OTIScalar>(a, b, c);
            worst = std::fmax(worst, std::fabs(extract_gradient(single, d) -
                                               extract_gradient(f_oti, d)));
        }
        report.check_residual(
            "oti.simultaneous_vs_separate",
            "3 directions seeded at once give the same derivatives as 3 "
            "single-seed runs",
            VerificationLevel::A, "n/a", "n/a", worst, 0.0, sw.seconds(),
            "this is the property that makes gradient cost independent of "
            "parameter count");
    }

    // Zero-seed equivalence: with no seeds, the real part must equal double.
    {
        Stopwatch sw;
        const OTIScalar r = composite<OTIScalar>(OTIScalar(kA), OTIScalar(kB), OTIScalar(kC));
        const double d = composite<double>(kA, kB, kC);
        report.check_residual(
            "oti.zero_seed_matches_double",
            "with all seeds zero the OTI real part equals the double result",
            VerificationLevel::A, "n/a", "n/a", std::fabs(r.real() - d), 0.0,
            sw.seconds(), "exact equality required");
    }

    // -----------------------------------------------------------------------
    // Aliasing hazard: documented as a regression guard.
    //
    // This asserts the LIBRARY MISBEHAVES under aliasing, which is why the
    // wrapper never aliases. If a future OTIlib made aliasing safe this case
    // would fail, and that failure would be the signal to relax the wrapper.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        onumm3n2_t a, b, aliased;
        onumm3n2_set_r(kA, &a);
        a.e1 = 1.0;
        onumm3n2_set_r(kB, &b);
        b.e2 = 1.0;

        onumm3n2_t clean;
        onumm3n2_mul_oo_to(&a, &b, &clean);

        aliased = a;
        onumm3n2_mul_oo_to(&aliased, &b, &aliased);  // deliberate aliasing

        const bool differs = (aliased.e1 != clean.e1) || (aliased.e2 != clean.e2) ||
                             (aliased.e12 != clean.e12);
        report.check_predicate(
            "oti.aliasing_is_unsafe",
            "OTIlib *_to functions produce different results when the output "
            "aliases an input, confirming the wrapper's no-alias rule is needed",
            VerificationLevel::A, differs,
            "aliased result differs from the clean result",
            differs ? "differs, as expected" : "identical -- aliasing may now be safe",
            sw.seconds(),
            "if this ever FAILS, upstream may have made aliasing safe; revisit "
            "the wrapper rather than assuming a regression");
    }

    // -----------------------------------------------------------------------
    // Generic-kernel compatibility: the Milestone 1 algebra with OTIScalar.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        Mat4<OTIScalar> m;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                m(i, j) = scalar_cast<OTIScalar>((i == j) ? 2.0 + i : 0.25 * (i - j));
        const InverseResult<OTIScalar> inv = inverse(m);
        const double resid = inv.ok ? identity_residual(m, inv.inverse) : 1.0;
        report.check_residual(
            "oti.generic_matrix_inverse",
            "Mat4/inverse/identity_residual instantiate and run with OTIScalar",
            VerificationLevel::A, "n/a", "n/a", resid, 1e-12, sw.seconds(),
            "the geometry kernels are unmodified; only the scalar type changed");
    }

    // -----------------------------------------------------------------------
    // MANDATORY: derivative of a matrix inverse.
    //
    //     d(A^-1)/dp = -A^-1 (dA/dp) A^-1
    //
    // This is the identity that governs the inverse metric g^{mu nu}, which
    // every curvature quantity in the project flows through. It is a genuinely
    // independent check: the left side comes from OTI propagating through the
    // adjugate inverse, the right side from an analytic identity evaluated with
    // ordinary double arithmetic on separately computed pieces.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;

        // A(p) depends on all three parameters, with entries chosen so the
        // matrix stays well conditioned and every parameter reaches every row.
        auto build = [](const OTIScalar& p0, const OTIScalar& p1,
                        const OTIScalar& p2) {
            Mat4<OTIScalar> a;
            const OTIScalar one = scalar_cast<OTIScalar>(1.0);
            a(0, 0) = scalar_cast<OTIScalar>(-1.0) - p0 * p0;
            a(0, 1) = p0 * p1 * scalar_cast<OTIScalar>(0.1);
            a(1, 0) = a(0, 1);
            a(1, 1) = one + p1 * scalar_cast<OTIScalar>(0.2);
            a(2, 2) = one + p2 * scalar_cast<OTIScalar>(0.15);
            a(3, 3) = one + p0 * p2 * scalar_cast<OTIScalar>(0.05);
            a(0, 3) = p2 * scalar_cast<OTIScalar>(0.07);
            a(3, 0) = a(0, 3);
            a(1, 2) = p0 * scalar_cast<OTIScalar>(0.03);
            a(2, 1) = a(1, 2);
            return a;
        };

        const sensitivity::SeededParameters p = sensitivity::seed_all(kA, kB, kC);
        const Mat4<OTIScalar> A = build(p.velocity, p.radius, p.sigma);
        const InverseResult<OTIScalar> Ainv = inverse(A);

        double worst = 0.0;
        if (!Ainv.ok) {
            worst = 1.0;
        } else {
            // Real-valued A and A^-1 for the right-hand side.
            Mat4<double> Ar, Ainvr;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    Ar(i, j) = real_part(A(i, j));
                    Ainvr(i, j) = real_part(Ainv.inverse(i, j));
                }

            for (int d = 0; d < 3; ++d) {
                // dA/dp_d, extracted from the OTI entries of A itself.
                Mat4<double> dA;
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        dA(i, j) = extract_gradient(A(i, j), d);

                // rhs = -A^-1 dA A^-1
                const Mat4<double> rhs =
                    mat_scale(mat_mul(mat_mul(Ainvr, dA), Ainvr), -1.0);

                // lhs = d(A^-1)/dp_d, straight from OTI propagation.
                Mat4<double> lhs;
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        lhs(i, j) = extract_gradient(Ainv.inverse(i, j), d);

                worst = std::fmax(worst, max_abs_diff(lhs, rhs));

                record("matrix_inverse",
                       std::string("dAinv/dp") + std::to_string(d) + " max|lhs-rhs|",
                       0.0, max_abs_diff(lhs, rhs));
            }
        }

        report.check_residual(
            "oti.matrix_inverse_derivative_identity",
            "d(A^-1)/dp == -A^-1 (dA/dp) A^-1 for all three parameters",
            VerificationLevel::D, "n/a", "a=1.3, b=0.7, c=0.45", worst, 1e-11,
            sw.seconds(),
            "mandatory before OTI may be used for inverse metrics; the two "
            "sides come from independent routes (AD through the adjugate vs an "
            "analytic identity in double)");
    }

    // -----------------------------------------------------------------------
    // Write analytic_derivatives.csv
    // -----------------------------------------------------------------------
    {
        std::ofstream os(out_dir + "/analytic_derivatives.csv", std::ios::binary);
        if (os) {
            os << "function,quantity,analytic,oti,absolute_error\n";
            os.setf(std::ios::scientific);
            os.precision(12);
            for (const DerivRecord& r : records()) {
                os << r.function << ',' << r.quantity << ',' << r.analytic << ','
                   << r.oti << ',' << r.abs_error << '\n';
            }
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad

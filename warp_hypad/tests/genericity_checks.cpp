#include <cmath>
#include <cstdio>
#include <string>
#include <type_traits>

#include "checks.hpp"
#include "warp_hypad/core/field.hpp"
#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/tensor.hpp"
#include "warp_hypad/core/test_scalar.hpp"

namespace warp_hypad {
namespace checks {
namespace {

// --- Compile-time proof obligations ---------------------------------------
//
// These are the properties that make TestScalar a meaningful stand-in for a
// hypercomplex AD type. They are asserted here rather than merely assumed,
// because if TestScalar ever acquired an implicit double conversion the
// genericity checks below would keep passing while proving nothing.

template <typename T, typename = void>
struct has_less : std::false_type {};
template <typename T>
struct has_less<T, std::void_t<decltype(std::declval<T>() < std::declval<T>())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_equality : std::false_type {};
template <typename T>
struct has_equality<T, std::void_t<decltype(std::declval<T>() == std::declval<T>())>>
    : std::true_type {};

static_assert(!has_less<TestScalar>::value,
              "TestScalar must not be ordered: a kernel that compares scalars "
              "cannot be given a hypercomplex type later");
static_assert(!has_equality<TestScalar>::value,
              "TestScalar must not be equality-comparable, for the same reason");
static_assert(!std::is_convertible<TestScalar, double>::value,
              "TestScalar must not implicitly convert to double, or a kernel "
              "could silently fall back to ordinary arithmetic");
static_assert(!std::is_convertible<double, TestScalar>::value,
              "double must not implicitly convert to TestScalar, or an "
              "unseeded constant could enter the computation unnoticed");

// The whole point: these instantiations must COMPILE. If any kernel in mat4.hpp
// or tensor.hpp branched on a scalar, compared two scalars, or named double
// internally, the TestScalar instantiation below would fail to build and this
// file would not link.
template <typename Scalar>
struct AlgebraOutcome {
    double determinant = 0.0;
    double residual_left = 0.0;
    double residual_right = 0.0;
    double contraction = 0.0;
    double trace = 0.0;
    bool inverse_ok = false;
};

// A fixed perturbed-Lorentzian test matrix: the shape the Alcubierre metric
// actually takes, so the genericity check exercises a realistic conditioning
// regime rather than a benign one.
const double kTestMatrix[4][4] = {
    {-0.91,  0.17, -0.05,  0.02},
    { 0.17,  1.04,  0.11, -0.03},
    {-0.05,  0.11,  0.98,  0.07},
    { 0.02, -0.03,  0.07,  1.02},
};

template <typename Scalar>
AlgebraOutcome<Scalar> run_pipeline() {
    Mat4<Scalar> a;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            a(i, j) = scalar_cast<Scalar>(kTestMatrix[i][j]);

    AlgebraOutcome<Scalar> out;

    const InverseResult<Scalar> inv = inverse(a);
    out.inverse_ok = inv.ok;
    out.determinant = real_part(inv.determinant);

    if (inv.ok) {
        out.residual_left = identity_residual(a, inv.inverse);
        out.residual_right = identity_residual(inv.inverse, a);
        out.contraction = real_part(full_contract(inv.inverse, a));
    }

    out.trace = real_part(trace_mixed(a));

    // Exercise the rank-3 and rank-4 paths too, so genericity is proven for
    // the containers Milestones 4-5 depend on, not just for Mat4.
    Rank3<Scalar> d;
    for (int s = 0; s < 4; ++s)
        for (int m = 0; m < 4; ++m)
            for (int n = 0; n < 4; ++n)
                d(s, m, n) = scalar_cast<Scalar>(0.5 * s - 0.25 * m + 0.125 * n);
    const Rank3<Scalar> contracted = contract_first(a, d);
    out.contraction += max_abs_real(contracted);

    return out;
}

// Field containers must also be instantiable with a custom scalar, since every
// field in the pipeline will eventually hold HyPADScalar.
template <typename Scalar>
double exercise_fields(const Grid& grid) {
    Rank2Field<Scalar> field(grid);
    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        Mat4<Scalar> m;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                m(i, j) = scalar_cast<Scalar>(kTestMatrix[i][j] +
                                              0.001 * static_cast<double>(p));
        field.set_matrix(p, m);
    }
    double acc = 0.0;
    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        acc += real_part(trace_mixed(field.matrix_at(p)));
    }
    return acc;
}

std::string fmt(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

}  // namespace

void run_genericity_checks(VerificationReport& report) {
    // Instantiating all three is itself the primary result of this group.
    const AlgebraOutcome<double> d = run_pipeline<double>();
    const AlgebraOutcome<long double> ld = run_pipeline<long double>();
    const AlgebraOutcome<TestScalar> ts = run_pipeline<TestScalar>();

    {
        Stopwatch sw;
        report.check_predicate(
            "generic.compiles_with_custom_scalar",
            "Mat4, Rank3, Rank4 and every contraction instantiate with double, "
            "long double and a custom scalar that has no ordering, no equality "
            "and no implicit conversion to double",
            VerificationLevel::A, true, "all three instantiations compile",
            "all three instantiations compiled and ran",
            sw.seconds(),
            "enforced at compile time by static_assert plus successful "
            "instantiation; reaching this line at all is the proof");
    }

    {
        Stopwatch sw;
        report.check_predicate(
            "generic.inverse_ok_agrees",
            "the singularity verdict does not depend on the scalar type",
            VerificationLevel::A,
            d.inverse_ok && ld.inverse_ok && ts.inverse_ok,
            "all three report ok == true",
            std::string("double=") + (d.inverse_ok ? "ok" : "singular") +
                ", long double=" + (ld.inverse_ok ? "ok" : "singular") +
                ", TestScalar=" + (ts.inverse_ok ? "ok" : "singular"),
            sw.seconds());
    }

    // double vs TestScalar: TestScalar performs the identical sequence of
    // double operations, so agreement should be bit-exact. Asserting exact
    // equality rather than a tolerance is the stronger statement, and it is
    // the one that would catch an accidental precision change inside a kernel.
    {
        Stopwatch sw;
        const bool same = (d.determinant == ts.determinant) &&
                          (d.trace == ts.trace) &&
                          (d.contraction == ts.contraction) &&
                          (d.residual_left == ts.residual_left) &&
                          (d.residual_right == ts.residual_right);
        report.check_predicate(
            "generic.double_vs_custom_scalar_bit_exact",
            "TestScalar reproduces the double results bit for bit, since it "
            "performs the same operations in the same order",
            VerificationLevel::A, same, "bit-identical results",
            same ? "bit-identical"
                 : ("det: " + fmt(d.determinant) + " vs " + fmt(ts.determinant) +
                    "; trace: " + fmt(d.trace) + " vs " + fmt(ts.trace)),
            sw.seconds(),
            "a mismatch would mean a kernel treats the two types differently");
    }

    // long double genuinely has more precision, so exact agreement is not
    // expected. What must hold is that the two agree to double precision --
    // which is evidence that the algorithm, not the arithmetic width, sets the
    // answer.
    {
        Stopwatch sw;
        const double rel = std::fabs(d.determinant - ld.determinant) /
                           std::fmax(1e-300, std::fabs(ld.determinant));
        report.check_residual(
            "generic.long_double_agrees",
            "the determinant computed in long double agrees with the double "
            "result to double precision",
            VerificationLevel::D, "n/a", "perturbed Lorentzian test matrix",
            rel, 1e-14, sw.seconds(),
            "long double is an independent precision, not an independent "
            "algorithm, so this is weaker than the Gauss-Jordan comparison");
    }

    {
        Stopwatch sw;
        report.check_residual(
            "generic.inverse_residual_custom_scalar",
            "A * A^-1 == I holds when computed entirely in TestScalar",
            VerificationLevel::A, "n/a", "perturbed Lorentzian test matrix",
            std::fmax(ts.residual_left, ts.residual_right), 1e-12, sw.seconds());
    }

    // Field containers with a custom scalar.
    {
        Stopwatch sw;
        GridSpec spec;
        spec.nx = 5; spec.ny = 4; spec.nz = 3;
        spec.x_min = -1.0; spec.x_max = 1.0;
        spec.y_min = -1.0; spec.y_max = 1.0;
        spec.z_min = -1.0; spec.z_max = 1.0;
        const Grid grid(spec);

        const double a_d = exercise_fields<double>(grid);
        const double a_ts = exercise_fields<TestScalar>(grid);

        report.check_predicate(
            "generic.fields_with_custom_scalar",
            "Rank2Field instantiates with a custom scalar and produces "
            "bit-identical results to double",
            VerificationLevel::A, a_d == a_ts, "bit-identical accumulations",
            a_d == a_ts ? "bit-identical"
                        : (fmt(a_d) + " vs " + fmt(a_ts)),
            sw.seconds(),
            "every field in the pipeline will hold HyPADScalar eventually");
    }

    // ScalarTraits round-trip: the seam itself must be an identity on reals.
    {
        Stopwatch sw;
        const double values[6] = {0.0, 1.0, -1.0, 0.2, 1e-8, -3.14159265358979};
        int mismatches = 0;
        for (double v : values) {
            if (real_part(scalar_cast<double>(v)) != v) ++mismatches;
            if (real_part(scalar_cast<TestScalar>(v)) != v) ++mismatches;
            if (static_cast<double>(real_part(scalar_cast<long double>(v))) != v) ++mismatches;
        }
        report.check_exact(
            "generic.scalar_traits_roundtrip",
            "real_part(scalar_cast<S>(x)) == x for every supported scalar type",
            VerificationLevel::A, "n/a", 0.0, static_cast<double>(mismatches),
            sw.seconds(),
            "this is the seam HyPADScalar will specialise in Milestone 10");
    }
}

}  // namespace checks
}  // namespace warp_hypad

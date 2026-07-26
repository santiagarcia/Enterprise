#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "derivative_checks.hpp"
#include "warp_hypad/core/test_scalar.hpp"
#include "warp_hypad/geometry/finite_difference.hpp"

namespace warp_hypad {
namespace checks {
namespace {

using geometry::Axis;

// --- Error accounting ------------------------------------------------------
//
// Interior and boundary errors are accumulated separately, always. A single
// global norm hides the boundary, which is exactly where a warp-bubble figure
// shows artefacts, and averaging the two would let a first-order edge pass
// behind a second-order interior.
struct ErrorNorms {
    double interior_max = 0.0;
    double boundary_max = 0.0;
    double l1 = 0.0;
    double l2 = 0.0;
    double linf = 0.0;
    std::size_t interior_count = 0;
    std::size_t boundary_count = 0;

    void add(double err, bool boundary) {
        const double a = std::fabs(err);
        if (boundary) {
            if (a > boundary_max) boundary_max = a;
            ++boundary_count;
        } else {
            if (a > interior_max) interior_max = a;
            ++interior_count;
        }
        l1 += a;
        l2 += a * a;
        if (a > linf) linf = a;
    }

    void finalise() {
        const double n = static_cast<double>(interior_count + boundary_count);
        if (n > 0.0) {
            l1 /= n;
            l2 = std::sqrt(l2 / n);
        }
    }
};

// Analytic test field: value plus the derivatives we compare against.
struct TestField {
    std::string name;
    std::function<double(double, double, double)> f;
    std::function<double(double, double, double, int)> d1;      // axis
    std::function<double(double, double, double, int)> d2;      // axis
    std::function<double(double, double, double, int, int)> dmixed;
};

TestField constant_field() {
    return {"constant",
            [](double, double, double) { return 2.75; },
            [](double, double, double, int) { return 0.0; },
            [](double, double, double, int) { return 0.0; },
            [](double, double, double, int, int) { return 0.0; }};
}

TestField linear_field() {
    return {"linear",
            [](double x, double y, double z) { return 1.5 + 2.0 * x - 3.0 * y + 0.5 * z; },
            [](double, double, double, int a) {
                const double d[3] = {2.0, -3.0, 0.5};
                return d[a];
            },
            [](double, double, double, int) { return 0.0; },
            [](double, double, double, int, int) { return 0.0; }};
}

// x^4 + 2y^3 - 3z^2 + xyz -- the polynomial named in the project spec.
TestField polynomial_field() {
    return {
        "polynomial",
        [](double x, double y, double z) {
            return x * x * x * x + 2.0 * y * y * y - 3.0 * z * z + x * y * z;
        },
        [](double x, double y, double z, int a) {
            if (a == 0) return 4.0 * x * x * x + y * z;
            if (a == 1) return 6.0 * y * y + x * z;
            return -6.0 * z + x * y;
        },
        [](double x, double y, double z, int a) {
            (void)z;
            if (a == 0) return 12.0 * x * x;
            if (a == 1) return 12.0 * y;
            return -6.0;
        },
        [](double x, double y, double z, int a, int b) {
            if (a == b) return 0.0;  // not used; mixed pairs only
            const int s = a + b;
            if (s == 1) return z;  // xy
            if (s == 2) return y;  // xz
            return x;              // yz
        }};
}

// sin(Ax + px) cos(By + py) exp(Cz)
//
// The PHASES matter and are not decorative. With pure sin(Ax) on the symmetric
// domain [-1,1], the leading boundary truncation coefficient for d/dx is
// proportional to cos(Ax) evaluated at x = +-1 -- and cos is even, so BOTH
// faces carry the same coefficient. Choosing A = 1.7 makes cos(1.7) = -0.129 at
// both of them, the leading term nearly cancels, the subleading term dominates,
// and the observed boundary order comes out as noise (measured: 8.2, -3.3, 1.4
// across successive refinements).
//
// A phase breaks the symmetry so the two faces have different coefficients and
// the max over the boundary is governed by the larger one. With A = 1.1,
// px = 0.35: |cos| is 0.12 at x=+1 but 0.73 at x=-1, and |sin| is 0.99 at x=+1.
// So both the first- and second-derivative boundary tests have a healthy leading
// term somewhere on the boundary.
TestField trig_field() {
    const double A = 1.1, B = 1.7, C = 0.8;
    const double PX = 0.35, PY = 0.6;
    return {
        "trigonometric",
        [=](double x, double y, double z) {
            return std::sin(A * x + PX) * std::cos(B * y + PY) * std::exp(C * z);
        },
        [=](double x, double y, double z, int a) {
            const double sx = std::sin(A * x + PX), cx = std::cos(A * x + PX);
            const double cy = std::cos(B * y + PY), sy = std::sin(B * y + PY);
            const double ez = std::exp(C * z);
            if (a == 0) return A * cx * cy * ez;
            if (a == 1) return -B * sx * sy * ez;
            return C * sx * cy * ez;
        },
        [=](double x, double y, double z, int a) {
            const double sx = std::sin(A * x + PX);
            const double cy = std::cos(B * y + PY);
            const double ez = std::exp(C * z);
            if (a == 0) return -A * A * sx * cy * ez;
            if (a == 1) return -B * B * sx * cy * ez;
            return C * C * sx * cy * ez;
        },
        [=](double x, double y, double z, int a, int b) {
            const double sx = std::sin(A * x + PX), cx = std::cos(A * x + PX);
            const double cy = std::cos(B * y + PY), sy = std::sin(B * y + PY);
            const double ez = std::exp(C * z);
            const int s = a + b;
            if (s == 1) return -A * B * cx * sy * ez;  // xy
            if (s == 2) return A * C * cx * cy * ez;   // xz
            return -B * C * sx * sy * ez;              // yz
        }};
}

Grid make_grid(int n) {
    GridSpec s;
    s.nx = s.ny = s.nz = n;
    s.x_min = -1.0; s.x_max = 1.0;
    s.y_min = -1.0; s.y_max = 1.0;
    s.z_min = -1.0; s.z_max = 1.0;
    return Grid(s);
}

template <typename Scalar>
ScalarField<Scalar> sample(const Grid& g, const TestField& tf) {
    ScalarField<Scalar> f(g);
    for (int i = 0; i < g.nx(); ++i)
        for (int j = 0; j < g.ny(); ++j)
            for (int k = 0; k < g.nz(); ++k)
                f.at(i, j, k) =
                    scalar_cast<Scalar>(tf.f(g.x(i), g.y(j), g.z(k)));
    return f;
}

enum class Op { First, Second, Mixed };

ErrorNorms measure(const Grid& g, const TestField& tf, Op op, int a, int b) {
    const ScalarField<double> f = sample<double>(g, tf);
    ScalarField<double> num(g);
    switch (op) {
        case Op::First:
            num = geometry::first_derivative(f, g, static_cast<Axis>(a));
            break;
        case Op::Second:
            num = geometry::second_derivative(f, g, static_cast<Axis>(a));
            break;
        case Op::Mixed:
            num = geometry::mixed_derivative(f, g, static_cast<Axis>(a),
                                             static_cast<Axis>(b));
            break;
    }

    ErrorNorms e;
    for (int i = 0; i < g.nx(); ++i)
        for (int j = 0; j < g.ny(); ++j)
            for (int k = 0; k < g.nz(); ++k) {
                const std::size_t p = g.flat(i, j, k);
                const double x = g.x(i), y = g.y(j), z = g.z(k);
                double exact = 0.0;
                switch (op) {
                    case Op::First: exact = tf.d1(x, y, z, a); break;
                    case Op::Second: exact = tf.d2(x, y, z, a); break;
                    case Op::Mixed: exact = tf.dmixed(x, y, z, a, b); break;
                }
                // A point is "boundary" for this operator if it uses a one-sided
                // stencil on ANY axis the operator touches.
                bool bnd = geometry::is_axis_boundary(g, p, static_cast<Axis>(a));
                if (op == Op::Mixed)
                    bnd = bnd || geometry::is_axis_boundary(g, p, static_cast<Axis>(b));
                e.add(num.at(p) - exact, bnd);
            }
    e.finalise();
    return e;
}

struct ConvergenceRow {
    std::string field;
    std::string op;
    int n;
    double h;
    double interior_max;
    double boundary_max;
    double l1, l2, linf;
    double order_interior;  // NaN for the first resolution
    double order_boundary;
};

std::vector<ConvergenceRow>& rows() {
    static std::vector<ConvergenceRow> r;
    return r;
}

const int kResolutions[4] = {11, 21, 41, 81};

// Observed order between successive resolutions with h halving.
double observed_order(double e_coarse, double e_fine) {
    if (!(e_coarse > 0.0) || !(e_fine > 0.0)) return std::nan("");
    return std::log(e_coarse / e_fine) / std::log(2.0);
}

}  // namespace

void run_derivative_checks(VerificationReport& report, const std::string& out_dir) {
    const TestField fields[4] = {constant_field(), linear_field(),
                                 polynomial_field(), trig_field()};
    const char* axis_name[3] = {"x", "y", "z"};

    // -----------------------------------------------------------------------
    // Exactness on constant and linear fields.
    //
    // These are Level A, not convergence: a second-order stencil reproduces a
    // linear field's first derivative exactly, and annihilates a constant
    // exactly. Any nonzero result is a coefficient error, not truncation.
    // -----------------------------------------------------------------------
    {
        const Grid g = make_grid(21);
        for (int fi = 0; fi < 2; ++fi) {
            const TestField& tf = fields[fi];
            Stopwatch sw;
            double worst_d1 = 0.0, worst_d2 = 0.0;
            for (int a = 0; a < 3; ++a) {
                const ErrorNorms e1 = measure(g, tf, Op::First, a, a);
                const ErrorNorms e2 = measure(g, tf, Op::Second, a, a);
                worst_d1 = std::fmax(worst_d1, e1.linf);
                worst_d2 = std::fmax(worst_d2, e2.linf);
            }
            report.check_residual(
                "fd.exact." + tf.name + ".first",
                tf.name == "constant"
                    ? "d/dx of a constant is exactly 0 on every axis"
                    : "d/dx of a linear field is exact on every axis",
                VerificationLevel::A, "21x21x21", "n/a", worst_d1, 1e-13,
                sw.seconds(), "interior and boundary stencils together");
            report.check_residual(
                "fd.exact." + tf.name + ".second",
                "second derivative of a " + tf.name + " field is exactly 0",
                VerificationLevel::A, "21x21x21", "n/a", worst_d2, 1e-11,
                sw.seconds());
        }
    }

    // -----------------------------------------------------------------------
    // Grid convergence: four resolutions, interior and boundary separately.
    // -----------------------------------------------------------------------
    struct OpSpec {
        Op op;
        int a, b;
        std::string label;
    };
    std::vector<OpSpec> ops;
    for (int a = 0; a < 3; ++a) {
        ops.push_back({Op::First, a, a, std::string("d/d") + axis_name[a]});
        ops.push_back({Op::Second, a, a, std::string("d2/d") + axis_name[a] + "2"});
    }
    ops.push_back({Op::Mixed, 0, 1, "d2/dxdy"});
    ops.push_back({Op::Mixed, 0, 2, "d2/dxdz"});
    ops.push_back({Op::Mixed, 1, 2, "d2/dydz"});

    for (int fi = 2; fi < 4; ++fi) {  // polynomial and trigonometric
        const TestField& tf = fields[fi];
        for (const OpSpec& os : ops) {
            std::vector<ErrorNorms> per_res;
            std::vector<double> hs;
            for (int r = 0; r < 4; ++r) {
                const Grid g = make_grid(kResolutions[r]);
                per_res.push_back(measure(g, tf, os.op, os.a, os.b));
                hs.push_back(g.dx());
            }

            // A stencil of order 2 is EXACT whenever the derivative it truncates
            // against vanishes identically. For x^4 + 2y^3 - 3z^2 + xyz that is
            // most of the operators: d/dz truncates against d3/dz3, which is 0.
            // The measured error is then pure roundoff, which GROWS as h shrinks
            // (~eps/h for first derivatives, ~eps/h^2 for second), so the
            // "observed order" comes out negative and means nothing.
            //
            // Detect that regime on the COARSEST grid, where roundoff has not
            // yet grown and truncation would still be large if it existed at
            // all, and assert exactness instead of an order.
            const bool exact_regime = per_res[0].interior_max < 1e-10 &&
                                      per_res[0].boundary_max < 1e-10;

            double worst_interior_order = 1e9;
            double worst_boundary_order = 1e9;
            for (int r = 0; r < 4; ++r) {
                ConvergenceRow row;
                row.field = tf.name;
                row.op = os.label;
                row.n = kResolutions[r];
                row.h = hs[static_cast<std::size_t>(r)];
                row.interior_max = per_res[static_cast<std::size_t>(r)].interior_max;
                row.boundary_max = per_res[static_cast<std::size_t>(r)].boundary_max;
                row.l1 = per_res[static_cast<std::size_t>(r)].l1;
                row.l2 = per_res[static_cast<std::size_t>(r)].l2;
                row.linf = per_res[static_cast<std::size_t>(r)].linf;
                row.order_interior = std::nan("");
                row.order_boundary = std::nan("");
                if (r > 0) {
                    row.order_interior = observed_order(
                        per_res[static_cast<std::size_t>(r - 1)].interior_max,
                        per_res[static_cast<std::size_t>(r)].interior_max);
                    row.order_boundary = observed_order(
                        per_res[static_cast<std::size_t>(r - 1)].boundary_max,
                        per_res[static_cast<std::size_t>(r)].boundary_max);
                    // Judge convergence on the FINEST pair. That is where the
                    // asymptotic regime actually holds; coarser pairs still
                    // carry subleading terms and would understate the order for
                    // a scheme that is genuinely second order.
                    if (r == 3) {
                        worst_interior_order = row.order_interior;
                        worst_boundary_order = row.order_boundary;
                    }
                }
                rows().push_back(row);
            }

            Stopwatch sw;

            if (exact_regime) {
                // The operator annihilates this field exactly. Assert that,
                // rather than a convergence order that does not exist.
                const double worst_fine =
                    std::fmax(per_res[3].interior_max, per_res[3].boundary_max);
                report.check_residual(
                    "fd.exact_operator." + tf.name + "." + os.label,
                    "the stencil is exact for this field (the derivative it "
                    "truncates against vanishes identically), so the error stays "
                    "at roundoff on every grid",
                    VerificationLevel::A, "n = 11..81", tf.name, worst_fine,
                    1e-8, sw.seconds(),
                    "no convergence order is asserted here: the measured error "
                    "is pure roundoff and rises as ~eps/h^k under refinement");
                continue;
            }

            std::ostringstream note;
            note << "n = 11,21,41,81; finest-pair interior order "
                 << worst_interior_order << ", boundary " << worst_boundary_order;

            // Expect 2. Accept >= 1.85 to allow for the finite asymptotic
            // regime; anything below means the scheme is not second order.
            report.check_predicate(
                "fd.order.interior." + tf.name + "." + os.label,
                "observed interior convergence order approaches 2",
                VerificationLevel::E, worst_interior_order >= 1.85,
                ">= 1.85 (theoretical 2)",
                std::to_string(worst_interior_order), sw.seconds(), note.str());

            report.check_predicate(
                "fd.order.boundary." + tf.name + "." + os.label,
                "observed boundary convergence order approaches 2, i.e. the "
                "one-sided stencils are the same order as the interior ones",
                VerificationLevel::E, worst_boundary_order >= 1.85,
                ">= 1.85 (theoretical 2)",
                std::to_string(worst_boundary_order), sw.seconds(),
                "a three-point one-sided second derivative would score ~1 here");
        }
    }

    // -----------------------------------------------------------------------
    // Mixed-derivative commutativity.
    //
    // D_a and D_b touch disjoint index axes and each stencil depends only on
    // its own axis index, so they commute EXACTLY, not just to truncation
    // order. This checks arithmetic, and would catch a stencil whose choice
    // accidentally depended on another axis.
    // -----------------------------------------------------------------------
    {
        const Grid g = make_grid(21);
        const TestField& tf = fields[3];
        const ScalarField<double> f = sample<double>(g, tf);
        Stopwatch sw;
        double worst = 0.0;
        const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
        for (const auto& pr : pairs) {
            const ScalarField<double> ab = geometry::mixed_derivative(
                f, g, static_cast<Axis>(pr[0]), static_cast<Axis>(pr[1]));
            const ScalarField<double> ba = geometry::mixed_derivative(
                f, g, static_cast<Axis>(pr[1]), static_cast<Axis>(pr[0]));
            for (std::size_t p = 0; p < g.num_points(); ++p)
                worst = std::fmax(worst, std::fabs(ab.at(p) - ba.at(p)));
        }
        report.check_residual(
            "fd.mixed_commutes",
            "D_a(D_b f) == D_b(D_a f) at every point, including corners",
            VerificationLevel::A, "21x21x21", "trigonometric field", worst,
            1e-12, sw.seconds(),
            "exact, not to truncation order: the operators touch disjoint axes");
    }

    // -----------------------------------------------------------------------
    // Tensor-field differentiation, each component a different function.
    // -----------------------------------------------------------------------
    {
        const Grid g = make_grid(41);
        Stopwatch sw;

        // g_{mu nu}(x) = (mu+1) * sin(x) + (nu+1) * y^2 + (mu*nu+1) * z
        Rank2Field<double> field(g);
        for (int i = 0; i < g.nx(); ++i)
            for (int j = 0; j < g.ny(); ++j)
                for (int k = 0; k < g.nz(); ++k) {
                    const std::size_t p = g.flat(i, j, k);
                    const double x = g.x(i), y = g.y(j), z = g.z(k);
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu)
                            field.at(p, mu, nu) = (mu + 1) * std::sin(x) +
                                                  (nu + 1) * y * y +
                                                  (mu * nu + 1) * z;
                }

        const Rank3Field<double> grad =
            geometry::spacetime_gradient_frozen_time(field, g);

        double worst_interior = 0.0, worst_boundary = 0.0;
        double worst_time = 0.0;
        for (int i = 0; i < g.nx(); ++i)
            for (int j = 0; j < g.ny(); ++j)
                for (int k = 0; k < g.nz(); ++k) {
                    const std::size_t p = g.flat(i, j, k);
                    const double x = g.x(i), y = g.y(j);
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu) {
                            const double ex[3] = {(mu + 1) * std::cos(x),
                                                  (nu + 1) * 2.0 * y,
                                                  static_cast<double>(mu * nu + 1)};
                            for (int ax = 0; ax < 3; ++ax) {
                                const double err =
                                    grad.at(p, ax + 1, mu, nu) - ex[ax];
                                const bool bnd = geometry::is_axis_boundary(
                                    g, p, static_cast<Axis>(ax));
                                if (bnd)
                                    worst_boundary = std::fmax(worst_boundary, std::fabs(err));
                                else
                                    worst_interior = std::fmax(worst_interior, std::fabs(err));
                            }
                            worst_time = std::fmax(worst_time,
                                                   std::fabs(grad.at(p, 0, mu, nu)));
                        }
                }

        // Tolerances are the theoretical truncation bounds, not round numbers.
        // At n = 41 on [-1,1], h = 0.05. The largest third derivative present is
        // d3/dx3 of (mu+1)sin(x) with mu+1 = 4, so |f'''| <= 4.
        //   interior centred:   h^2/6 * 4 = 1.67e-3
        //   boundary one-sided: h^2/3 * 4 = 3.33e-3
        // Allowing 1.5x headroom over each keeps the check meaningful rather
        // than merely satisfiable.
        const double h = g.dx();
        const double interior_bound = 1.5 * (h * h / 6.0) * 4.0;
        const double boundary_bound = 1.5 * (h * h / 3.0) * 4.0;

        report.check_residual(
            "fd.tensor_field.interior",
            "every rank-2 component is differentiated against its own analytic "
            "derivative (interior points)",
            VerificationLevel::B, "41x41x41",
            "g_{mu nu} = (mu+1)sin x + (nu+1)y^2 + (mu nu+1)z", worst_interior,
            interior_bound, sw.seconds(),
            "16 components x 3 axes, each with a distinct analytic form; "
            "tolerance is 1.5x the leading truncation term h^2|f'''|/6");

        report.check_residual(
            "fd.tensor_field.boundary",
            "same, at boundary points using one-sided stencils",
            VerificationLevel::B, "41x41x41", "same field", worst_boundary,
            boundary_bound, sw.seconds(),
            "reported separately from the interior, always; the one-sided "
            "constant is 2x the centred one, hence the larger bound");

        report.check_residual(
            "fd.tensor_field.time_slot_zero",
            "with no time-derivative field supplied, d_0 g_{mu nu} is exactly 0",
            VerificationLevel::A, "41x41x41", "n/a", worst_time, 0.0,
            sw.seconds(),
            "the frozen-time assumption is explicit at the call site, never "
            "implicit in the operator");
    }

    // -----------------------------------------------------------------------
    // Scalar genericity: identical results for double, long double, TestScalar.
    // -----------------------------------------------------------------------
    {
        const Grid g = make_grid(21);
        const TestField& tf = fields[3];
        Stopwatch sw;

        const ScalarField<double> fd = sample<double>(g, tf);
        const ScalarField<TestScalar> ft = sample<TestScalar>(g, tf);
        const ScalarField<long double> fl = sample<long double>(g, tf);

        const ScalarField<double> dd = geometry::first_derivative(fd, g, Axis::X);
        const ScalarField<TestScalar> dt = geometry::first_derivative(ft, g, Axis::X);
        const ScalarField<long double> dl = geometry::first_derivative(fl, g, Axis::X);

        double worst_ts = 0.0, worst_ld = 0.0;
        for (std::size_t p = 0; p < g.num_points(); ++p) {
            worst_ts = std::fmax(worst_ts,
                                 std::fabs(dd.at(p) - real_part(dt.at(p))));
            worst_ld = std::fmax(worst_ld,
                                 std::fabs(dd.at(p) - static_cast<double>(dl.at(p))));
        }

        report.check_residual(
            "fd.generic.test_scalar_bit_exact",
            "the derivative operators give bit-identical results for double and "
            "a custom scalar with no ordering and no implicit double conversion",
            VerificationLevel::A, "21x21x21", "trigonometric field", worst_ts,
            0.0, sw.seconds(),
            "compiling at all is the genericity proof; this checks the values "
            "too");

        report.check_residual(
            "fd.generic.long_double_agrees",
            "long double agrees with double to double precision",
            VerificationLevel::D, "21x21x21", "trigonometric field", worst_ld,
            1e-12, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Grid validation refuses to differentiate an impossible grid.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        GridSpec s;
        s.nx = 3; s.ny = 21; s.nz = 21;
        s.x_min = -1; s.x_max = 1; s.y_min = -1; s.y_max = 1; s.z_min = -1; s.z_max = 1;
        const Grid small(s);
        const bool first_ok = geometry::validate_for_derivatives(small, false).empty();
        const bool second_rejected =
            !geometry::validate_for_derivatives(small, true).empty();
        report.check_predicate(
            "fd.grid_validation",
            "a 3-point axis is accepted for first derivatives and rejected for "
            "second derivatives, whose one-sided stencil needs 4 points",
            VerificationLevel::A, first_ok && second_rejected,
            "first accepted, second rejected",
            std::string("first ") + (first_ok ? "accepted" : "REJECTED") +
                ", second " + (second_rejected ? "rejected" : "ACCEPTED"),
            sw.seconds());
    }

    // --- Export convergence table ------------------------------------------
    {
        std::ofstream os(out_dir + "/convergence.csv", std::ios::binary);
        if (os) {
            os << "field,operator,n,h,interior_max,boundary_max,l1,l2,linf,"
                  "observed_order_interior,observed_order_boundary\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (const ConvergenceRow& r : rows()) {
                os << r.field << ',' << r.op << ',' << r.n << ',' << r.h << ','
                   << r.interior_max << ',' << r.boundary_max << ',' << r.l1 << ','
                   << r.l2 << ',' << r.linf << ',';
                if (std::isfinite(r.order_interior)) os << r.order_interior;
                os << ',';
                if (std::isfinite(r.order_boundary)) os << r.order_boundary;
                os << '\n';
            }
        }
    }

    // --- Export a 1D profile for the analytic-versus-numerical figure -------
    {
        const Grid g = make_grid(41);
        const TestField& tf = fields[3];
        const ScalarField<double> f = sample<double>(g, tf);
        const ScalarField<double> d1 = geometry::first_derivative(f, g, Axis::X);
        const ScalarField<double> d2 = geometry::second_derivative(f, g, Axis::X);

        const int jc = g.center_index_y();
        const int kc = g.center_index_z();

        std::ofstream os(out_dir + "/profile_x.csv", std::ios::binary);
        if (os && jc >= 0 && kc >= 0) {
            os << "x,f,d1_numeric,d1_analytic,d2_numeric,d2_analytic,"
                  "d1_error,d2_error,is_boundary\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int i = 0; i < g.nx(); ++i) {
                const std::size_t p = g.flat(i, jc, kc);
                const double x = g.x(i), y = g.y(jc), z = g.z(kc);
                const double a1 = tf.d1(x, y, z, 0);
                const double a2 = tf.d2(x, y, z, 0);
                os << x << ',' << f.at(p) << ',' << d1.at(p) << ',' << a1 << ','
                   << d2.at(p) << ',' << a2 << ',' << (d1.at(p) - a1) << ','
                   << (d2.at(p) - a2) << ',' << (i == 0 || i == g.nx() - 1 ? 1 : 0)
                   << '\n';
            }
        }
    }

    // --- Export an error map for the slice figure --------------------------
    {
        const Grid g = make_grid(41);
        const TestField& tf = fields[3];
        const ScalarField<double> f = sample<double>(g, tf);
        const ScalarField<double> d1 = geometry::first_derivative(f, g, Axis::X);
        const int kc = g.center_index_z();
        std::ofstream os(out_dir + "/error_map_xy.csv", std::ios::binary);
        if (os && kc >= 0) {
            os << "x,y,error\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int i = 0; i < g.nx(); ++i)
                for (int j = 0; j < g.ny(); ++j) {
                    const std::size_t p = g.flat(i, j, kc);
                    const double e = d1.at(p) - tf.d1(g.x(i), g.y(j), g.z(kc), 0);
                    os << g.x(i) << ',' << g.y(j) << ',' << e << '\n';
                }
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad

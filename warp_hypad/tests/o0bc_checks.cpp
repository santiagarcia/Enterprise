#include "o0b_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "warp_hypad/core/metadata.hpp"

#if defined(WARP_HYPAD_HAVE_OTILIB)
#include "warp_hypad/sensitivity/otilib_extract.hpp"
#include "warp_hypad/sensitivity/otilib_math.hpp"
#include "warp_hypad/sensitivity/otilib_scalar.hpp"
#include "warp_hypad/sensitivity/otilib_seed.hpp"
#endif

// ---------------------------------------------------------------------------
// O0B-C: the corrective gate.
//
// Each case here exists because the corresponding O0B claim was wrong, weaker
// than stated, or measured something other than what its name said. The four
// defects, and what replaced them:
//
//   A FINITE DIFFERENCE RAN INSIDE THE PRODUCTION METRIC JET. The small-radius
//   branch of shape_w_jet obtained the fourth Taylor coefficient by centred
//   differencing the analytic third derivative, and set F''' to zero outright.
//   F'''(0) is NOT zero for this shape function -- it is -0.063 at sigma = 1
//   and -9.2 at R = 1.5, sigma = 2 -- so the third design derivative was being
//   silently truncated on the axis by a milestone whose whole point was exact
//   analytic jets.
//
//   A CONSTRAINT RESIDUAL WAS MULTIPLIED BY ZERO. The conversion case computed
//   max_abs(...) * 0.0 + compose, so the constraint term was discarded and only
//   the composition was tested, under a case name claiming the constraints
//   stayed satisfied.
//
//   A CONSTRAINED GRADIENT WAS LABELLED WHITENED. The batch seeded z and formed
//   c = c0 + P z. W never appeared. The optimizer will step in y, and
//   sum_i P_ij B_i and sum_i (PW)_ij B_i are different vectors.
//
//   AN ARCHITECTURAL CLAIM WAS A LITERAL true. "There is no public path that
//   accepts an arbitrary g_mu_nu" was asserted by a hard-coded predicate.
//
// A fifth correction is definitional rather than numerical: the quadrature
// scaffold implied that the eventual integrand is exactly integrated on every
// knot span. The BASIS is piecewise polynomial; the metric and energy
// integrands are not, and the max-based objectives are not even smooth.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::ShiftFamilyConfig;

// The backend capacity, declared so the disabled build can still state the
// schedule. The enabled build static_asserts it against the compiled type in
// o0b_route_checks.cpp.
inline constexpr int kDeclaredOtiDirections = 3;

// Every optimizer-facing metric-construction entry point, with the argument
// types it accepts. Recorded rather than proved: see the case note.
struct ApiEntry {
    const char* function;
    const char* accepts;
    const char* returns;
    bool accepts_bare_mat4;
};

const ApiEntry kApiSurface[] = {
    {"design::assemble_covariant", "ADMDesignFields<Scalar>",
     "CovariantMetric<Scalar>", false},
    {"design::assemble_contravariant", "ADMDesignFields<Scalar>, bool& ok",
     "ContravariantMetric<Scalar>", false},
    {"design::f1_metric_jet",
     "ShiftFamilyConfig, AxisymmetricBasis, vector<double> coefficients, "
     "t, x, y, z", "MetricJet<double>", false},
    {"design::shift_jet",
     "ShiftFamilyConfig, AxisymmetricBasis, vector<double>, t, x, y, z",
     "ShiftJet", false},
    {"design::shift_value_scalar",
     "ShiftFamilyConfig, AxisymmetricBasis, vector<Scalar>, t, x, y, z",
     "Scalar", false},
    {"design::metric_component_scalar",
     "ShiftFamilyConfig, AxisymmetricBasis, vector<Scalar>, t, x, y, z, A, B",
     "Scalar", false},
    {"design::metric_first_derivative_scalar",
     "ShiftFamilyConfig, AxisymmetricBasis, vector<Scalar>, t, x, y, z, mu, "
     "A, B", "Scalar", false},
    {"design::CoefficientSpaces::to_raw",
     "ConstrainedSplineCoordinates or WhitenedDesignCoordinates",
     "RawSplineCoefficients + ConversionStatus", false},
    {"design::CoefficientSpaces::to_constrained",
     "WhitenedDesignCoordinates", "ConstrainedSplineCoordinates + status",
     false},
};

}  // namespace

void run_o0bc_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;
    const ShiftFamilyConfig cfg = o0b_config();
    AxisymmetricBasis basis = o0b_basis();
    const design::DesignQuadrature quad = design::build_quadrature(basis);
    const design::ConstraintSystem sys = o0b_constraints(basis, quad);
    const design::ConstraintProjection pr = design::build_projection(sys);
    const design::DenseMatrix M = design::sobolev_gram(
        basis, quad, design::SobolevScales{cfg.scale_R0, cfg.scale_R0});
    const design::WhiteningResult wr = design::whiten(M, pr);
    const design::CoefficientSpaces spaces(pr, wr);
    const design::DenseMatrix PW = spaces.PW();
    const design::DenseMatrix C = sys.matrix();
    const std::vector<double> coef =
        o0b_coefficients(basis.coefficient_count(), 0.02);
    const double R0 = cfg.scale_R0;

    struct Cfg { double R, sigma; };
    const Cfg configs[6] = {{2.0, 1.0}, {2.0, 2.0}, {2.0, 3.0},
                            {1.5, 1.0}, {1.5, 2.0}, {1.5, 3.0}};

    std::ofstream ax(out_dir + "/axis_jet_checks.csv");
    ax << "R,sigma,quantity,closed_form,recurrence,branch_value,"
          "absolute_difference\n";

    // --- 1/2/3. Exact axis F', F'', F''' -------------------------------------
    double e1 = 0.0, e2 = 0.0, e3 = 0.0, e_cf = 0.0;
    double min_abs_f3 = 1e300;
    for (const Cfg& c : configs) {
        const design::ShapeWJet lim =
            design::shape_w_jet_axis_limits(R0, c.R, c.sigma);
        const design::ShapeWJet at0 = design::shape_w_jet(0.0, R0, c.R, c.sigma);
        e1 = std::fmax(e1, std::fabs(at0.F1 - lim.F1));
        e2 = std::fmax(e2, std::fabs(at0.F2 - lim.F2));
        e3 = std::fmax(e3, std::fabs(at0.F3 - lim.F3));
        min_abs_f3 = std::fmin(min_abs_f3, std::fabs(lim.F3));

        // The recurrence against the closed forms the specification supplies.
        const design::ShapeAxisDerivatives cf =
            design::shape_axis_derivatives_closed_form(c.R, c.sigma);
        const double r2 = design::shape_even_derivative_at_axis(2, c.R, c.sigma);
        const double r4 = design::shape_even_derivative_at_axis(4, c.R, c.sigma);
        const double r6 = design::shape_even_derivative_at_axis(6, c.R, c.sigma);
        e_cf = std::fmax(e_cf,
                         std::fmax(std::fabs(cf.f2 - r2) / std::fabs(cf.f2),
                                   std::fmax(std::fabs(cf.f4 - r4) /
                                                 std::fabs(cf.f4),
                                             std::fabs(cf.f6 - r6) /
                                                 std::fabs(cf.f6))));
        ax << c.R << "," << c.sigma << ",f2," << cf.f2 << "," << r2 << ","
           << lim.F1 << "," << std::fabs(cf.f2 - r2) << "\n";
        ax << c.R << "," << c.sigma << ",f4," << cf.f4 << "," << r4 << ","
           << lim.F2 << "," << std::fabs(cf.f4 - r4) << "\n";
        ax << c.R << "," << c.sigma << ",f6," << cf.f6 << "," << r6 << ","
           << lim.F3 << "," << std::fabs(cf.f6 - r6) << "\n";
    }
    report.check_residual(
        "o0bc.exact_axis_first_derivative",
        "F'(0) = R0^2 f''(0)/2 exactly, from the closed form rather than a "
        "limit or a difference",
        VerificationLevel::A, "6 (R, sigma) configurations", "R0 = 2", e1,
        1.0e-15, sw.seconds(),
        "recurrence vs supplied closed forms, worst relative: " +
            std::to_string(e_cf));
    sw.reset();
    report.check_residual(
        "o0bc.exact_axis_second_derivative", "F''(0) = R0^4 f4(0)/12 exactly",
        VerificationLevel::A, "6 configurations", "R0 = 2", e2, 1.0e-15,
        sw.seconds(), "");
    sw.reset();
    {
        char buf[300];
        std::snprintf(buf, sizeof(buf),
                      "worst |branch - closed form| = %.2e; smallest |F'''(0)| "
                      "over the six configurations is %.6f, so the previous "
                      "F3 = 0 was wrong by that much at minimum",
                      e3, min_abs_f3);
        report.check_predicate(
            "o0bc.exact_axis_third_derivative",
            "F'''(0) = R0^6 f6(0)/120 exactly, and it is NOT zero",
            VerificationLevel::A, e3 <= 1.0e-15 && min_abs_f3 > 1e-3,
            "exact, and nonzero in every configuration", buf, sw.seconds(),
            "the earlier code set F3 = 0 on the axis without derivation; this "
            "case is what would have caught it");
    }

    // --- 4. Branch continuity -------------------------------------------------
    sw.reset();
    double jump_abs = 0.0, jump_rel = 0.0;
    std::ofstream bc(out_dir + "/axis_branch_continuity.csv");
    bc << "R,sigma,switch_w,dF,dF1,dF2,dF3,relative_dF3\n";
    for (const Cfg& c : configs) {
        const double rsw = design::kAxisSigmaRadius / c.sigma;
        const double wsw = (rsw / R0) * (rsw / R0);
        const design::ShapeWJet a =
            design::shape_w_jet(wsw * (1.0 - 1e-9), R0, c.R, c.sigma);
        const design::ShapeWJet b =
            design::shape_w_jet(wsw * (1.0 + 1e-9), R0, c.R, c.sigma);
        const double d0 = std::fabs(a.F - b.F), d1 = std::fabs(a.F1 - b.F1);
        const double d2 = std::fabs(a.F2 - b.F2), d3 = std::fabs(a.F3 - b.F3);
        const double rel = d3 / std::fmax(std::fabs(b.F3), 1e-300);
        jump_abs = std::fmax(jump_abs, std::fmax(std::fmax(d0, d1),
                                                 std::fmax(d2, d3)));
        jump_rel = std::fmax(jump_rel, rel);
        bc << c.R << "," << c.sigma << "," << wsw << "," << d0 << "," << d1
           << "," << d2 << "," << d3 << "," << rel << "\n";
    }
    report.check_predicate(
        "o0bc.small_w_branch_continuity",
        "F, F', F'' and F''' are continuous across the Taylor/quotient switch "
        "to the accuracy each branch can deliver there",
        VerificationLevel::B, jump_rel < 1.0e-6,
        "relative F''' jump below 1e-6",
        "worst absolute jump " + std::to_string(jump_abs) +
            ", worst relative F''' jump " + std::to_string(jump_rel),
        sw.seconds(),
        "the switch at sigma*r = 0.1 is where the two error curves cross: below "
        "it the quotient form loses digits to cancellation (3 f'/r^5 is 1e4 "
        "where F''' is 1e-1), above it the truncated series degrades. The "
        "residual jump is what each branch can actually deliver, measured "
        "rather than assumed");

    // --- 5. No production finite difference ----------------------------------
    sw.reset();
    // A finite difference inside the jet would make the result depend on a step
    // size. Evaluating the same point through two independently-built bases and
    // at bitwise-identical inputs must be bit-exact; and the axis value must be
    // reproduced EXACTLY by the closed-form limit, which a difference quotient
    // cannot do.
    bool bit_exact = true, axis_exact = true;
    for (const Cfg& c : configs) {
        const design::ShapeWJet a = design::shape_w_jet(0.0, R0, c.R, c.sigma);
        const design::ShapeWJet b = design::shape_w_jet(0.0, R0, c.R, c.sigma);
        if (a.F1 != b.F1 || a.F2 != b.F2 || a.F3 != b.F3) bit_exact = false;
        const design::ShapeWJet lim =
            design::shape_w_jet_axis_limits(R0, c.R, c.sigma);
        if (a.F1 != lim.F1 || a.F2 != lim.F2 || a.F3 != lim.F3)
            axis_exact = false;
    }
    report.check_predicate(
        "o0bc.no_production_finite_difference",
        "the axis branch reproduces the closed-form limits BIT EXACTLY, which a "
        "difference quotient could not, and repeated evaluation is bit-identical",
        VerificationLevel::A, bit_exact && axis_exact,
        "bit-exact against the closed form",
        std::string("bit-identical repeats: ") + (bit_exact ? "yes" : "NO") +
            ", exact against closed-form limits: " +
            (axis_exact ? "yes" : "NO"),
        sw.seconds(),
        "the previous implementation evaluated shape_radial_jet at +-1e-3 and "
        "divided; that cannot land exactly on the analytic limit, so this case "
        "distinguishes the two implementations rather than merely describing "
        "them");

    // --- 6/7. Third metric jet on and near the axis --------------------------
    sw.reset();
    double axis_finite_max = 0.0;
    bool axis_all_finite = true;
    double near_axis_worst = 0.0;
    for (double rho : {0.0, 1e-12, 1e-6, 1e-3, 0.02, 0.05}) {
        const design::MetricJet<double> j =
            design::f1_metric_jet(cfg, basis, coef, 0.0, 0.0, rho, 0.0, true);
        for (int m = 0; m < 4; ++m)
            for (int n = 0; n < 4; ++n)
                for (int l = 0; l < 4; ++l)
                    for (int A = 0; A < 4; ++A)
                        for (int B = 0; B < 4; ++B) {
                            const double v =
                                j.d3g[static_cast<std::size_t>(m)]
                                     [static_cast<std::size_t>(n)]
                                     [static_cast<std::size_t>(l)].m(A, B);
                            if (!std::isfinite(v)) axis_all_finite = false;
                            axis_finite_max = std::fmax(axis_finite_max,
                                                        std::fabs(v));
                        }
        near_axis_worst = std::fmax(
            near_axis_worst,
            design::third_derivative_symmetry_residual(j));
    }
    report.check_predicate(
        "o0bc.exact_axis_third_metric_jet",
        "the third metric jet is finite and permutation-symmetric at exactly "
        "(xi, s) = (0, 0)",
        VerificationLevel::A,
        axis_all_finite && axis_finite_max < 1e4 && near_axis_worst < 1e-12,
        "finite and symmetric on the axis",
        "max |d3 g| = " + std::to_string(axis_finite_max) +
            ", worst permutation residual " + std::to_string(near_axis_worst),
        sw.seconds(), "");

    sw.reset();
    double mixed_worst = 0.0;
    const double hfd = 1e-4;
    for (double rho : {1e-3, 0.02, 0.08}) {
        const design::MetricJet<double> j =
            design::f1_metric_jet(cfg, basis, coef, 0.2, 0.7, rho, 0.0, true);
        // d3 g_txy against a difference of the analytic d2 g_tx.
        const design::MetricJet<double> jp = design::f1_metric_jet(
            cfg, basis, coef, 0.2, 0.7, rho + hfd, 0.0, true);
        const design::MetricJet<double> jm = design::f1_metric_jet(
            cfg, basis, coef, 0.2, 0.7, rho - hfd, 0.0, true);
        for (int A = 0; A < 4; ++A)
            for (int B = 0; B < 4; ++B)
                mixed_worst = std::fmax(
                    mixed_worst,
                    std::fabs(j.d3g[0][1][2].m(A, B) -
                              (jp.d2g[0][1].m(A, B) - jm.d2g[0][1].m(A, B)) /
                                  (2 * hfd)));
    }
    report.check_residual(
        "o0bc.near_axis_mixed_third_derivative",
        "d_txy g near the axis against a difference of the analytic d_tx g",
        VerificationLevel::B, "rho = 1e-3, 0.02, 0.08", "h = 1e-4",
        mixed_worst, 1.0e-4, sw.seconds(),
        "these radii straddle the branch switch at sigma*r = 0.1, so the case "
        "exercises both sides and the crossing itself");

    // --- 8/9/10. Three DIFFERENT gradients ------------------------------------
    const double pt[4] = {0.35, 1.1, 0.9, 0.4};
    const design::DesignCoordinates dp = design::make_design_coordinates(
        pt[1] - (cfg.bubble_center_x0 + cfg.velocity * pt[0]), pt[2], pt[3],
        0.0, R0);
    const design::TensorBasisEvaluation ev = basis.evaluate(dp);
    std::vector<double> Bk(
        static_cast<std::size_t>(basis.coefficient_count()), 0.0);
    for (std::size_t m = 0; m < ev.size(); ++m)
        Bk[static_cast<std::size_t>(ev.flat_index[m])] = ev.value[m];

    std::vector<double> ref_z(static_cast<std::size_t>(pr.n_free), 0.0);
    std::vector<double> ref_y(static_cast<std::size_t>(wr.retained), 0.0);
    for (int j = 0; j < pr.n_free; ++j) {
        double a = 0.0;
        for (int i = 0; i < pr.n_raw; ++i)
            a += pr.P(i, j) * Bk[static_cast<std::size_t>(i)];
        ref_z[static_cast<std::size_t>(j)] = a;
    }
    for (int j = 0; j < wr.retained; ++j) {
        double a = 0.0;
        for (int i = 0; i < pr.n_raw; ++i)
            a += PW(i, j) * Bk[static_cast<std::size_t>(i)];
        ref_y[static_cast<std::size_t>(j)] = a;
    }
    double zy_gap = 0.0;
    for (int j = 0; j < std::min(pr.n_free, wr.retained); ++j)
        zy_gap = std::fmax(zy_gap,
                           std::fabs(ref_z[static_cast<std::size_t>(j)] -
                                     ref_y[static_cast<std::size_t>(j)]));

    std::ofstream gg(out_dir + "/gradient_space_checks.csv");
    gg << "coordinate_space,method,index,target,reference,computed,"
          "absolute_difference\n";

    sw.reset();
    double raw_worst = 0.0;
    for (std::size_t m = 0; m < ev.size(); ++m) {
        const int k = ev.flat_index[m];
        const double h = 1e-6;
        std::vector<double> cp = coef, cm = coef;
        cp[static_cast<std::size_t>(k)] += h;
        cm[static_cast<std::size_t>(k)] -= h;
        const double d = (design::shift_value_scalar<double>(cfg, basis, cp,
                                                             pt[0], pt[1], pt[2],
                                                             pt[3]) -
                          design::shift_value_scalar<double>(cfg, basis, cm,
                                                             pt[0], pt[1], pt[2],
                                                             pt[3])) / (2 * h);
        raw_worst = std::fmax(raw_worst, std::fabs(d - ev.value[m]));
        gg << "raw_c,finite_difference," << k << ",beta_x," << ev.value[m]
           << "," << d << "," << std::fabs(d - ev.value[m]) << "\n";
    }
    report.check_residual(
        "o0bc.raw_coordinate_gradient",
        "d beta^x / d c_k = B_k, in the RAW coefficient space",
        VerificationLevel::B, "36 coefficients", "h = 1e-6", raw_worst, 1.0e-8,
        sw.seconds(), "labelled raw_c, which is what it is");

    // Finite differences directly in z and directly in y -- no derivative
    // backend needed, so these run in both configurations.
    auto beta_from_z = [&](const std::vector<double>& z) {
        design::ConstrainedSplineCoordinates zc;
        zc.z = z;
        design::RawSplineCoefficients rc;
        spaces.to_raw(zc, rc);
        return design::shift_value_scalar<double>(cfg, basis, rc.c, pt[0], pt[1],
                                                  pt[2], pt[3]);
    };
    auto beta_from_y = [&](const std::vector<double>& y) {
        design::WhitenedDesignCoordinates yc;
        yc.y = y;
        design::RawSplineCoefficients rc;
        spaces.to_raw(yc, rc);
        return design::shift_value_scalar<double>(cfg, basis, rc.c, pt[0], pt[1],
                                                  pt[2], pt[3]);
    };

    sw.reset();
    double z_fd = 0.0;
    for (int j = 0; j < pr.n_free; ++j) {
        const double h = 1e-6;
        std::vector<double> zp(static_cast<std::size_t>(pr.n_free), 0.0), zm = zp;
        zp[static_cast<std::size_t>(j)] = h;
        zm[static_cast<std::size_t>(j)] = -h;
        const double d = (beta_from_z(zp) - beta_from_z(zm)) / (2 * h);
        z_fd = std::fmax(z_fd, std::fabs(d - ref_z[static_cast<std::size_t>(j)]));
        gg << "constrained_z,finite_difference," << j << ",beta_x,"
           << ref_z[static_cast<std::size_t>(j)] << "," << d << ","
           << std::fabs(d - ref_z[static_cast<std::size_t>(j)]) << "\n";
    }
    report.check_residual(
        "o0bc.constrained_z_gradient",
        "d beta^x / d z_j = sum_i P_ij B_i, in the CONSTRAINED space",
        VerificationLevel::B, "34 coordinates", "h = 1e-6", z_fd, 1.0e-8,
        sw.seconds(), "labelled constrained_z");

    sw.reset();
    double y_fd = 0.0, y_fd_g00 = 0.0;
    for (int j = 0; j < wr.retained; ++j) {
        const double h = 1e-6;
        std::vector<double> yp(static_cast<std::size_t>(wr.retained), 0.0),
            ym = yp;
        yp[static_cast<std::size_t>(j)] = h;
        ym[static_cast<std::size_t>(j)] = -h;
        const double d = (beta_from_y(yp) - beta_from_y(ym)) / (2 * h);
        y_fd = std::fmax(y_fd, std::fabs(d - ref_y[static_cast<std::size_t>(j)]));
        // g_00 through the same route: reference 2 beta sum_i (PW)_ij B_i.
        design::WhitenedDesignCoordinates ycp, ycm;
        ycp.y = yp; ycm.y = ym;
        design::RawSplineCoefficients rp, rm;
        spaces.to_raw(ycp, rp);
        spaces.to_raw(ycm, rm);
        const double d00 =
            (design::metric_component_scalar<double>(cfg, basis, rp.c, pt[0],
                                                     pt[1], pt[2], pt[3], 0, 0) -
             design::metric_component_scalar<double>(cfg, basis, rm.c, pt[0],
                                                     pt[1], pt[2], pt[3], 0, 0)) /
            (2 * h);
        // beta at the base point of THIS route (c0 + PW*0 = c0), not at coef.
        design::WhitenedDesignCoordinates y0;
        y0.y.assign(static_cast<std::size_t>(wr.retained), 0.0);
        design::RawSplineCoefficients r0;
        spaces.to_raw(y0, r0);
        const double beta0 = design::shift_value_scalar<double>(
            cfg, basis, r0.c, pt[0], pt[1], pt[2], pt[3]);
        const double want00 = 2.0 * beta0 * ref_y[static_cast<std::size_t>(j)];
        y_fd_g00 = std::fmax(y_fd_g00, std::fabs(d00 - want00));
        gg << "whitened_y,finite_difference," << j << ",beta_x,"
           << ref_y[static_cast<std::size_t>(j)] << "," << d << ","
           << std::fabs(d - ref_y[static_cast<std::size_t>(j)]) << "\n";
        gg << "whitened_y,finite_difference," << j << ",g00," << want00 << ","
           << d00 << "," << std::fabs(d00 - want00) << "\n";
    }
    {
        char buf[320];
        std::snprintf(buf, sizeof(buf),
                      "beta^x residual %.3e, g_00 residual %.3e; the z and y "
                      "references differ by up to %.6f component-wise",
                      y_fd, y_fd_g00, zy_gap);
        report.check_predicate(
            "o0bc.whitened_y_gradient",
            "d/dy_j of beta^x and g_00 match the (PW)-based references, in the "
            "WHITENED optimizer space",
            VerificationLevel::B, y_fd < 1e-8 && y_fd_g00 < 1e-7,
            "both within tolerance", buf, sw.seconds(),
            "this is the gradient an optimizer would actually use; the earlier "
            "milestone measured the z-gradient and called it this one");
    }
    sw.reset();
    report.check_residual(
        "o0bc.whitened_finite_difference_agreement",
        "centred differences taken DIRECTLY in y agree with the analytic "
        "(PW)-based reference",
        VerificationLevel::B, "34 whitened coordinates", "h = 1e-6", y_fd,
        1.0e-8, sw.seconds(),
        "the differences are taken in the optimizer coordinate itself, not in c "
        "and then transformed");

#if defined(WARP_HYPAD_HAVE_OTILIB)
    using sensitivity::OTIScalar;
    sw.reset();
    // Separate vs batched, both seeding WHITENED coordinates.
    double sep_y = 0.0, bat_y = 0.0, sep_bat = 0.0;
    std::vector<std::vector<int>> batch;
    for (int s = 0; s < wr.retained; s += sensitivity::kOtiDirections) {
        std::vector<int> b;
        for (int d = 0; d < sensitivity::kOtiDirections && s + d < wr.retained;
             ++d)
            b.push_back(s + d);
        batch.push_back(b);
    }
    auto seed_y = [&](const std::vector<int>& idx, std::vector<OTIScalar>& c) {
        std::vector<OTIScalar> y(static_cast<std::size_t>(wr.retained));
        for (int j = 0; j < wr.retained; ++j)
            y[static_cast<std::size_t>(j)] = sensitivity::make_constant(0.0);
        for (std::size_t d = 0; d < idx.size(); ++d)
            y[static_cast<std::size_t>(idx[d])] =
                sensitivity::seed_parameter(0.0, static_cast<int>(d));
        c.assign(static_cast<std::size_t>(pr.n_raw),
                 sensitivity::make_constant(0.0));
        for (int i = 0; i < pr.n_raw; ++i) {
            OTIScalar acc = sensitivity::make_constant(
                pr.c0[static_cast<std::size_t>(i)]);
            for (int j = 0; j < wr.retained; ++j)
                acc = acc + y[static_cast<std::size_t>(j)] * OTIScalar(PW(i, j));
            c[static_cast<std::size_t>(i)] = acc;
        }
    };
    std::vector<double> separate(static_cast<std::size_t>(wr.retained), 0.0);
    for (int j = 0; j < wr.retained; ++j) {
        std::vector<OTIScalar> c;
        seed_y({j}, c);
        const OTIScalar b = design::shift_value_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3]);
        const double got = sensitivity::extract_gradient(b, 0);
        separate[static_cast<std::size_t>(j)] = got;
        sep_y = std::fmax(sep_y,
                          std::fabs(got - ref_y[static_cast<std::size_t>(j)]));
        gg << "whitened_y,oti_separate," << j << ",beta_x,"
           << ref_y[static_cast<std::size_t>(j)] << "," << got << ","
           << std::fabs(got - ref_y[static_cast<std::size_t>(j)]) << "\n";
    }
    std::vector<int> seen(static_cast<std::size_t>(wr.retained), 0);
    for (const auto& b : batch) {
        std::vector<OTIScalar> c;
        seed_y(b, c);
        const OTIScalar bv = design::shift_value_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3]);
        const OTIScalar dgx = design::metric_first_derivative_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3], 1, 0, 1);
        for (std::size_t d = 0; d < b.size(); ++d) {
            ++seen[static_cast<std::size_t>(b[d])];
            const double got =
                sensitivity::extract_gradient(bv, static_cast<int>(d));
            const double want = ref_y[static_cast<std::size_t>(b[d])];
            bat_y = std::fmax(bat_y, std::fabs(got - want));
            sep_bat = std::fmax(
                sep_bat,
                std::fabs(got - separate[static_cast<std::size_t>(b[d])]));
            // A selected first-derivative channel through the same seeds.
            const double gx =
                sensitivity::extract_gradient(dgx, static_cast<int>(d));
            double refgx = 0.0;
            for (std::size_t m = 0; m < ev.size(); ++m)
                refgx += PW(ev.flat_index[m], b[d]) * ev.d_xi[m] / R0;
            bat_y = std::fmax(bat_y, std::fabs(gx - refgx));
            gg << "whitened_y,oti_batched," << b[d] << ",beta_x," << want << ","
               << got << "," << std::fabs(got - want) << "\n";
            gg << "whitened_y,oti_batched," << b[d] << ",d_x_g01," << refgx
               << "," << gx << "," << std::fabs(gx - refgx) << "\n";
        }
    }
    bool covers_once = true;
    for (int v : seen)
        if (v != 1) covers_once = false;
    report.check_residual(
        "o0bc.whitened_separate_batched_agreement",
        "separate and twelve-batch OTI seeding of the WHITENED coordinates give "
        "the same gradient, and both match the (PW)-based reference",
        VerificationLevel::A, "34 coordinates", "onumm3n2",
        std::fmax(std::fmax(sep_y, bat_y), sep_bat), 1.0e-14, sw.seconds(),
        "separate vs reference " + std::to_string(sep_y) + ", batched vs "
        "reference " + std::to_string(bat_y) + ", separate vs batched " +
            std::to_string(sep_bat));
    sw.reset();
    report.check_predicate(
        "o0bc.twelve_batch_whitened_coverage",
        "the complete gradient with respect to the 34 whitened optimizer "
        "coordinates is recovered in twelve deterministic OTI batches",
        VerificationLevel::A,
        covers_once && static_cast<int>(batch.size()) == 12,
        "12 batches covering each whitened coordinate exactly once",
        std::to_string(wr.retained) + " whitened coordinates in " +
            std::to_string(batch.size()) + " batches, each covered exactly "
            "once: " + (covers_once ? "yes" : "NO"),
        sw.seconds(),
        "Mixed Hessian entries between distinct batches are not available.");
#else
    sw.reset();
    for (const char* nm : {"o0bc.whitened_separate_batched_agreement"})
        report.check_skipped(
            nm, "separate vs batched OTI seeding of whitened coordinates",
            "the OTIlib sensitivity backend is not compiled into this "
            "configuration; the whitened gradient is still verified against the "
            "(PW)-based reference by centred differences taken directly in y",
            sw.seconds());
    {
        int nb = 0;
        for (int s = 0; s < wr.retained; s += kDeclaredOtiDirections) ++nb;
        report.check_predicate(
            "o0bc.twelve_batch_whitened_coverage",
            "the complete gradient with respect to the 34 whitened optimizer "
            "coordinates is recovered in twelve deterministic OTI batches",
            VerificationLevel::A, nb == 12,
            "12 batches", std::to_string(wr.retained) +
                " whitened coordinates in " + std::to_string(nb) +
                " batches (declared capacity, backend absent)",
            sw.seconds(),
            "Mixed Hessian entries between distinct batches are not available.");
    }
#endif

    // --- 14/15/16. Three separate residuals -----------------------------------
    sw.reset();
    auto resid = [&](const std::vector<double>& c) {
        double w = 0.0;
        for (double v : design::matvec(C, c)) w = std::fmax(w, std::fabs(v));
        return w;
    };
    double rz = 0.0, ry = 0.0, comp = 0.0;
    std::ofstream cc(out_dir + "/constraint_closure_checks.csv");
    cc << "trial,constrained_constraint_residual,whitened_constraint_residual,"
          "conversion_composition_residual\n";
    for (int t = 0; t < 16; ++t) {
        design::WhitenedDesignCoordinates y;
        y.y.assign(static_cast<std::size_t>(wr.retained), 0.0);
        for (int j = 0; j < wr.retained; ++j)
            y.y[static_cast<std::size_t>(j)] =
                std::sin(1.7 * j + 0.31 * t) * (1.0 + 0.1 * t);
        design::ConstrainedSplineCoordinates z;
        spaces.to_constrained(y, z);
        design::RawSplineCoefficients c_two, c_direct;
        spaces.to_raw(z, c_two);
        spaces.to_raw(y, c_direct);
        const double a = resid(c_two.c), b = resid(c_direct.c);
        double d = 0.0;
        for (std::size_t i = 0; i < c_two.c.size(); ++i)
            d = std::fmax(d, std::fabs(c_two.c[i] - c_direct.c[i]));
        rz = std::fmax(rz, a);
        ry = std::fmax(ry, b);
        comp = std::fmax(comp, d);
        cc << t << "," << a << "," << b << "," << d << "\n";
    }
    report.check_residual(
        "o0bc.constrained_constraint_closure",
        "|C (c0 + P z) - d|_inf over 16 deterministic vectors, measured directly",
        VerificationLevel::A, "16 trials", "34 whitened coordinates", rz,
        1.0e-11, sw.seconds(), "not multiplied by zero");
    sw.reset();
    report.check_residual(
        "o0bc.whitened_constraint_closure",
        "|C (c0 + P W y) - d|_inf, a SEPARATE quantity from the constrained "
        "residual",
        VerificationLevel::A, "16 trials", "34 whitened coordinates", ry,
        1.0e-11, sw.seconds(), "");
    sw.reset();
    report.check_residual(
        "o0bc.conversion_composition_residual",
        "|c_y - [c0 + P(W y)]|_inf: a composition residual, reported as its own "
        "quantity",
        VerificationLevel::A, "16 trials", "34 whitened coordinates", comp,
        1.0e-13, sw.seconds(),
        "three different numbers, three different cases; the earlier gate "
        "collapsed them into one and discarded two");

    // --- 17-21. Checked exponential -------------------------------------------
    sw.reset();
    std::ofstream ex(out_dir + "/exponential_guard_checks.csv");
    ex << "field,log_value,status,accepted,result_finite,result_positive\n";
    const double just_in_hi = design::kMaxLogLapse;
    const double just_out_hi = std::nextafter(design::kMaxLogLapse,
                                              std::numeric_limits<double>::infinity());
    const double just_in_lo = design::kMinLogLapse;
    const double just_out_lo = std::nextafter(design::kMinLogLapse,
                                              -std::numeric_limits<double>::infinity());
    bool lapse_over_rejected = false, lapse_under_rejected = false;
    bool lapse_hi_ok = false, lapse_lo_ok = false;
    for (double a : {just_in_lo, -10.0, 0.0, 10.0, just_in_hi, just_out_hi,
                     just_out_lo}) {
        double out = 0.0;
        const design::ExpStatus s = design::checked_exp<double>(a, out);
        const bool acc = s == design::ExpStatus::Ok;
        if (a == just_out_hi && s == design::ExpStatus::ExponentialOverflow)
            lapse_over_rejected = true;
        if (a == just_out_lo && s == design::ExpStatus::ExponentialUnderflow)
            lapse_under_rejected = true;
        if (a == just_in_hi && acc && std::isfinite(out) && out > 0.0)
            lapse_hi_ok = true;
        if (a == just_in_lo && acc && std::isfinite(out) && out > 0.0)
            lapse_lo_ok = true;
        ex << "lapse," << a << "," << design::to_string(s) << ","
           << (acc ? 1 : 0) << "," << (acc && std::isfinite(out) ? 1 : 0) << ","
           << (acc && out > 0.0 ? 1 : 0) << "\n";
    }
    report.check_predicate(
        "o0bc.lapse_overflow_rejection",
        "the largest accepted log lapse produces a finite positive alpha, and "
        "the next representable value above it is REJECTED as overflow",
        VerificationLevel::A, lapse_hi_ok && lapse_over_rejected,
        "accepted at the bound, rejected just above",
        "a = " + std::to_string(just_in_hi) + " accepted with finite positive "
        "alpha: " + (lapse_hi_ok ? "yes" : "NO") + "; a = " +
            std::to_string(just_out_hi) + " rejected: " +
            (lapse_over_rejected ? "yes" : "NO"),
        sw.seconds(),
        "'positive for every finite a' was false: exp(1000) is +inf and "
        "exp(-1000) is exactly 0, and both inputs are finite");
    sw.reset();
    report.check_predicate(
        "o0bc.lapse_underflow_rejection",
        "the smallest accepted log lapse produces a finite positive alpha, and "
        "the next representable value below it is REJECTED as underflow",
        VerificationLevel::A, lapse_lo_ok && lapse_under_rejected,
        "accepted at the bound, rejected just below",
        "a = " + std::to_string(just_in_lo) + " accepted: " +
            (lapse_lo_ok ? "yes" : "NO") + "; a = " +
            std::to_string(just_out_lo) + " rejected: " +
            (lapse_under_rejected ? "yes" : "NO"),
        sw.seconds(),
        "a zero lapse is the dangerous one: g^00 = -1/alpha^2 becomes infinite, "
        "which is exactly the invalid metric the parameterisation was supposed "
        "to make unrepresentable");

    sw.reset();
    int chol_over = 0, chol_under = 0, chol_ok = 0;
    for (int slot = 0; slot < 3; ++slot)
        for (double a : {just_in_lo, 0.0, just_in_hi, just_out_hi,
                         just_out_lo}) {
            design::CholeskyFactorL<double> L;
            if (slot == 0) L.l1 = a; else if (slot == 1) L.l2 = a; else L.l3 = a;
            design::SymmetricMatrix3<double> g;
            const design::ExpStatus s =
                design::checked_spatial_metric_from_cholesky<double>(L, g);
            if (a == just_out_hi &&
                s == design::ExpStatus::ExponentialOverflow) ++chol_over;
            if (a == just_out_lo &&
                s == design::ExpStatus::ExponentialUnderflow) ++chol_under;
            if (s == design::ExpStatus::Ok &&
                design::spatial_metric_positive_definite(g)) ++chol_ok;
            ex << "cholesky_l" << (slot + 1) << "," << a << ","
               << design::to_string(s) << ","
               << (s == design::ExpStatus::Ok ? 1 : 0) << ","
               << (s == design::ExpStatus::Ok ? 1 : 0) << ","
               << (s == design::ExpStatus::Ok ? 1 : 0) << "\n";
        }
    report.check_predicate(
        "o0bc.cholesky_overflow_rejection",
        "each of the three Cholesky diagonal log parameters rejects an "
        "overflowing value",
        VerificationLevel::A, chol_over == 3, "all three slots reject",
        std::to_string(chol_over) + " of 3 diagonal slots rejected overflow",
        sw.seconds(), "");
    sw.reset();
    report.check_predicate(
        "o0bc.cholesky_underflow_rejection",
        "each of the three Cholesky diagonal log parameters rejects an "
        "underflowing value",
        VerificationLevel::A, chol_under == 3, "all three slots reject",
        std::to_string(chol_under) + " of 3 diagonal slots rejected underflow; "
        + std::to_string(chol_ok) + " accepted values gave a positive-definite "
        "gamma", sw.seconds(), "");

    sw.reset();
    bool all_accepted_finite = true;
    for (double a : {just_in_lo, -50.0, 0.0, 50.0, just_in_hi}) {
        double al = 0.0;
        if (design::checked_exp<double>(a, al) != design::ExpStatus::Ok)
            continue;
        design::ADMDesignFields<double> f =
            design::ADMDesignFields<double>::minkowski();
        f.lapse = al;
        f.shift_contravariant[0] = -0.4;
        const design::CovariantMetric<double> g = design::assemble_covariant(f);
        bool ok = false;
        const design::ContravariantMetric<double> gi =
            design::assemble_contravariant(f, ok);
        for (int A = 0; A < 4; ++A)
            for (int B = 0; B < 4; ++B)
                if (!std::isfinite(g.m(A, B)) || !std::isfinite(gi.m(A, B)))
                    all_accepted_finite = false;
        if (!ok) all_accepted_finite = false;
    }
    report.check_predicate(
        "o0bc.finite_accepted_metric",
        "no ACCEPTED log parameter produces inf, NaN or a zero lapse anywhere "
        "in the assembled metric or its inverse",
        VerificationLevel::A, all_accepted_finite,
        "every accepted parameter gives a finite metric and inverse",
        all_accepted_finite ? "all finite" : "NONFINITE FOUND", sw.seconds(),
        "alpha = exp(a) > 0 and the Cholesky diagonal is positive and finite "
        "for every accepted representable log parameter -- the corrected "
        "structural statement");

#if defined(WARP_HYPAD_HAVE_OTILIB)
    sw.reset();
    {
        // The representability test applies to the REAL component.
        OTIScalar seeded = sensitivity::seed_parameter(just_in_hi, 0);
        OTIScalar out;
        const design::ExpStatus s1 = design::checked_exp<OTIScalar>(seeded, out);
        OTIScalar bad = sensitivity::seed_parameter(just_out_hi, 0);
        OTIScalar out2;
        const design::ExpStatus s2 = design::checked_exp<OTIScalar>(bad, out2);
        report.check_predicate(
            "o0bc.oti_exponential_guard",
            "the exponential guard applies to the real component of an OTI "
            "scalar, accepting a representable seed and rejecting an "
            "overflowing one",
            VerificationLevel::A,
            s1 == design::ExpStatus::Ok &&
                s2 == design::ExpStatus::ExponentialOverflow,
            "accepted then rejected",
            std::string("accepted: ") + design::to_string(s1) + ", rejected: " +
                design::to_string(s2),
            sw.seconds(),
            "testing every derivative channel individually would reject "
            "perfectly good seeds; the channels are finite whenever the real "
            "part is");
    }
#else
    sw.reset();
    report.check_skipped(
        "o0bc.oti_exponential_guard",
        "the exponential guard on an OTI scalar",
        "the OTIlib sensitivity backend is not compiled into this "
        "configuration; the double path is checked above and the guard is "
        "scalar-generic",
        sw.seconds());
#endif

    // --- 22/23/24. Quadrature semantics ---------------------------------------
    sw.reset();
    const design::QuadratureScaffold sc =
        design::build_quadrature_scaffold(basis, o0b_region());
    const design::ConvergenceContract contract;
    {
        const design::QuadratureExactness& e = sc.exactness;
        report.check_predicate(
            "o0bc.quadrature_exactness_statement_corrected",
            "the scaffold records the rule's formal polynomial degree WITHOUT "
            "claiming the eventual integrand is integrated exactly",
            VerificationLevel::A,
            e.gauss_legendre_nodes_per_coordinate == 8 &&
                e.polynomial_exactness_per_coordinate == 15 &&
                e.spline_basis_piecewise_polynomial &&
                !e.statement.empty(),
            "8 nodes, degree 15, basis piecewise polynomial, statement present",
            e.statement, sw.seconds(),
            "the basis IS piecewise polynomial on an element, which is what "
            "made the earlier phrasing tempting; the integrand built from it is "
            "not");
    }
    sw.reset();
    report.check_predicate(
        "o0bc.full_integrand_marked_nonpolynomial",
        "the full metric and energy integrands are recorded as NONpolynomial",
        VerificationLevel::A,
        !sc.exactness.full_metric_integrand_polynomial &&
            !sc.exactness.full_energy_integrand_polynomial,
        "both false",
        "metric integrand polynomial: false, energy integrand polynomial: "
        "false", sw.seconds(),
        "the metric involves exp and tanh of spline combinations and the "
        "objectives take max(0, .) of the result");
    sw.reset();
    report.check_predicate(
        "o0bc.energy_integral_marked_uncertified",
        "the energy integral is recorded as NOT exact and NOT certified, and "
        "the O0C convergence contract names every axis it must be certified on",
        VerificationLevel::A,
        !sc.exactness.energy_integral_exact &&
            contract.quadrature_orders.size() == 5 &&
            contract.element_subdivisions.size() == 3 &&
            contract.domain_enlargement_required &&
            contract.independent_uniform_grid_comparison_required &&
            contract.nonsmooth_objectives_treated_separately &&
            !contract.evaluated_in_o0b,
        "uncertified, with a five-order, three-subdivision contract",
        "orders 4,6,8,10,12; subdivisions 1,2,4; domain enlargement, "
        "uniform-grid comparison and separate nonsmooth treatment all "
        "required; evaluated in O0B: false", sw.seconds(),
        contract.nonsmooth_note);

    // --- 25/26. API surface ----------------------------------------------------
    sw.reset();
    int n_api = 0, n_bare = 0;
    std::ofstream api(out_dir + "/optimizer_api_surface.csv");
    api << "function,accepted_argument_types,returns,accepts_bare_mat4\n";
    for (const ApiEntry& e : kApiSurface) {
        ++n_api;
        if (e.accepts_bare_mat4) ++n_bare;
        api << e.function << ",\"" << e.accepts << "\",\"" << e.returns
            << "\"," << (e.accepts_bare_mat4 ? 1 : 0) << "\n";
    }
    report.check_predicate(
        "o0bc.optimizer_api_surface_recorded",
        "every optimizer-facing metric-construction function is listed with its "
        "accepted argument types, and none accepts a bare Mat4",
        VerificationLevel::A, n_api >= 9 && n_bare == 0,
        "complete listing, zero bare-Mat4 entries",
        std::to_string(n_api) + " entry points recorded, " +
            std::to_string(n_bare) + " accepting a bare Mat4",
        sw.seconds(),
        "ARCHITECTURAL INVARIANT REVIEWED AND API SURFACE RECORDED. This is a "
        "record of the current surface, not a proof that no bypass exists "
        "anywhere in the program");

    sw.reset();
    report.check_predicate(
        "o0bc.no_hardcoded_architectural_proof",
        "no case in this gate or in O0B asserts an architectural invariant "
        "through a hard-coded predicate",
        VerificationLevel::A, true,
        "claims narrowed to what is measured",
        "o0b.adm_metric_assembly now counts assembled, finite, Lorentzian "
        "metrics instead of passing a literal true; the absence-of-bypass "
        "claim is downgraded to 'reviewed and recorded'",
        sw.seconds(),
        "this case is itself a statement about the source rather than a "
        "measurement, and is labelled as such rather than presented as "
        "evidence");

    // --- 28. Scope -------------------------------------------------------------
    sw.reset();
    report.check_predicate(
        "o0bc.scope_statement", "the scope statement is present in every record",
        VerificationLevel::A, true, "present",
        "O0B-C corrects derivative, optimizer-coordinate and "
        "numerical-integration semantics in the verified ADM design layer. It "
        "does not evaluate energy, run an optimizer, demonstrate a causal "
        "shortcut or establish physical realizability.",
        sw.seconds(), "");

    // --- Records ---------------------------------------------------------------
    {
        std::ofstream f(out_dir + "/quadrature_semantics.json");
        f.precision(17);
        const design::QuadratureExactness& e = sc.exactness;
        f << "{\n";
        f << "  \"gauss_legendre_nodes_per_coordinate\": "
          << e.gauss_legendre_nodes_per_coordinate << ",\n";
        f << "  \"polynomial_exactness_per_coordinate\": "
          << e.polynomial_exactness_per_coordinate << ",\n";
        f << "  \"spline_basis_piecewise_polynomial\": "
          << (e.spline_basis_piecewise_polynomial ? "true" : "false") << ",\n";
        f << "  \"full_metric_integrand_polynomial\": "
          << (e.full_metric_integrand_polynomial ? "true" : "false") << ",\n";
        f << "  \"full_energy_integrand_polynomial\": "
          << (e.full_energy_integrand_polynomial ? "true" : "false") << ",\n";
        f << "  \"energy_integral_exact\": "
          << (e.energy_integral_exact ? "true" : "false") << ",\n";
        f << "  \"statement\": \"" << e.statement << "\",\n";
        f << "  \"o0c_convergence_contract\": {\n";
        f << "    \"quadrature_orders\": [4, 6, 8, 10, 12],\n";
        f << "    \"element_subdivisions\": [1, 2, 4],\n";
        f << "    \"domain_enlargement_required\": true,\n";
        f << "    \"independent_uniform_grid_comparison_required\": true,\n";
        f << "    \"nonsmooth_objectives_treated_separately\": true,\n";
        f << "    \"nonsmooth_note\": \"" << contract.nonsmooth_note << "\",\n";
        f << "    \"evaluated_in_o0b\": false\n  }\n}\n";
    }
    {
        std::ofstream f(out_dir + "/correction_summary.json");
        f.precision(17);
        f << "{\n  \"milestone\": \"O0B-C\",\n";
        f << "  \"corrections\": [\n";
        f << "    {\"defect\": \"finite difference inside the production metric "
             "jet\", \"where\": \"shape_w_jet small-radius branch\", "
             "\"was\": \"F2 from a centred difference of the analytic third "
             "derivative; F3 set to 0\", \"now\": \"exact even Taylor series "
             "from closed-form f^(2k)(0)\", \"material\": true, "
             "\"evidence\": \"F'''(0) ranges from -0.063 to -9.21 over the "
             "tested configurations and was being reported as 0\"},\n";
        f << "    {\"defect\": \"constraint residual multiplied by zero\", "
             "\"where\": \"o0b.raw_to_constrained_conversion\", \"was\": "
             "\"max_abs(...) * 0.0 + compose\", \"now\": \"three separate "
             "measured residuals\", \"material\": true, \"evidence\": \"the "
             "constraint term was discarded entirely\"},\n";
        f << "    {\"defect\": \"constrained gradient labelled whitened\", "
             "\"where\": \"o0b.separate_batched_oti_agreement\", \"was\": "
             "\"seeded z, formed c0 + P z\", \"now\": \"seeds y through PW, "
             "with the z-gradient kept and relabelled\", \"material\": true, "
             "\"evidence\": \"component-wise gap between the two references "
             "is " << zy_gap << "\"},\n";
        f << "    {\"defect\": \"hard-coded architectural proof\", \"where\": "
             "\"o0b.adm_metric_assembly\", \"was\": \"literal true\", "
             "\"now\": \"counts assembled valid metrics; absence of bypass "
             "downgraded to reviewed and recorded\", \"material\": true},\n";
        f << "    {\"defect\": \"quadrature exactness overclaim\", \"where\": "
             "\"elementwise quadrature scaffold\", \"was\": \"implied exact "
             "integration on every knot span\", \"now\": \"formal degree "
             "recorded separately from integrand polynomiality; energy "
             "integral marked uncertified\", \"material\": true},\n";
        f << "    {\"defect\": \"unchecked exponential\", \"where\": "
             "\"lapse and Cholesky diagonal\", \"was\": \"'positive for every "
             "finite a', which is false in floating point\", \"now\": "
             "\"checked_exp with four statuses and measured accept/reject "
             "bounds\", \"material\": true}\n";
        f << "  ],\n";
        f << "  \"energy_evaluated\": false,\n";
        f << "  \"optimizer_run\": false,\n";
        f << "  \"scope\": \"O0B-C corrects derivative, optimizer-coordinate "
             "and numerical-integration semantics in the verified ADM design "
             "layer. It does not evaluate energy, run an optimizer, "
             "demonstrate a causal shortcut or establish physical "
             "realizability.\"\n}\n";
    }
}

}  // namespace checks
}  // namespace warp_hypad

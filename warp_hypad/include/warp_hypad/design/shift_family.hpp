#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "warp_hypad/design/adm_design.hpp"
#include "warp_hypad/design/axisymmetric_basis.hpp"
#include "warp_hypad/metrics/metric.hpp"

// ---------------------------------------------------------------------------
// O0B-F1: the shift-only design family, and its analytic derivative jet.
//
//     alpha = 1,   gamma_ij = delta_ij,   beta^y = beta^z = 0,
//     beta^x = beta^x_baseline + delta beta^x
//
// THE BASELINE IS THE PRODUCTION SHAPE FUNCTION, NOT THE F0 PROJECTION.
// F0 is a representation test that O0A concluded is not resolution-certified;
// building on it would fold a 1.4-magnitude second-derivative error into every
// metric this milestone produces. The baseline here is
// metrics::shape_function, the same analytic tanh form the A3-A6 gates verify,
// and the F1 correction is added SEPARATELY on top of it.
//
// EVERYTHING IS A FUNCTION OF (xi, s), INCLUDING THE BASELINE.
//
// This is the observation that makes the whole jet regular. The baseline
// depends on the radius through r_s^2 = (x-x_c)^2 + y^2 + z^2, and in design
// coordinates
//
//     w = r_s^2 / R0^2 = xi^2 + s
//
// which is a POLYNOMIAL in xi and LINEAR in s. So every derivative of the
// baseline with respect to the design coordinates is a polynomial combination
// of F', F'', F''' where F(w) = f(R0 sqrt(w)) -- there is no 1/rho anywhere in
// the chain, exactly as in O0A. The only place the square root survives is
// inside F itself, and f is EVEN in r, so F and all its w-derivatives have
// finite limits on the axis. Those limits are supplied in closed form rather
// than approached numerically; see radial_jet().
//
// PRODUCTION SUPPLIES f AND f' ONLY. metrics::shape_function and
// metrics::shape_function_radial_derivative are verified from A3 onward, but
// nothing above first order exists there, and a third-order metric jet needs
// f''' . They are differentiated here in closed form from the same tanh
// expression, and radial_jet_matches_production() checks the two overlapping
// orders against production so the extension cannot drift from the function it
// claims to extend.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// f and its first three radial derivatives, from the tanh form.
//
//   f    = (t+ - t-) / d,                    t± = tanh(sigma (r ± R))
//   f'   = sigma [ (1-t+^2) - (1-t-^2) ] / d
//   f''  = sigma^2 [ -2 t+(1-t+^2) + 2 t-(1-t-^2) ] / d
//   f''' = sigma^3 [ -2(1-t+^2)(1-3t+^2) + 2(1-t-^2)(1-3t-^2) ] / d
//
// using (tanh)' = 1-t^2, and d/dr(1-t^2) = -2 t sigma (1-t^2).
struct ShapeRadialJet {
    double f = 1.0;
    double d1 = 0.0;   // df/dr
    double d2 = 0.0;   // d2f/dr2
    double d3 = 0.0;   // d3f/dr3
};

inline ShapeRadialJet shape_radial_jet(double r, double R, double sigma) {
    ShapeRadialJet j;
    const double tp = std::tanh(sigma * (r + R));
    const double tm = std::tanh(sigma * (r - R));
    const double d = 2.0 * std::tanh(sigma * R);
    const double sp = 1.0 - tp * tp;   // sech^2 at r+R
    const double sm = 1.0 - tm * tm;
    j.f = (tp - tm) / d;
    j.d1 = sigma * (sp - sm) / d;
    j.d2 = sigma * sigma * (-2.0 * tp * sp + 2.0 * tm * sm) / d;
    j.d3 = sigma * sigma * sigma *
           (-2.0 * sp * (1.0 - 3.0 * tp * tp) +
             2.0 * sm * (1.0 - 3.0 * tm * tm)) / d;
    return j;
}

// F(w) = f(R0 sqrt(w)) and its first three w-derivatives, where w = r^2/R0^2.
//
// With r = R0 sqrt(w) and dr/dw = R0/(2 sqrt(w)) = R0^2/(2r):
//
//     F'   = (R0^2/2)  f'/r
//     F''  = (R0^4/4)  (f''/r^2 - f'/r^3)
//     F''' = (R0^6/8)  (f'''/r^3 - 3 f''/r^4 + 3 f'/r^5)
//
// Each bracket is finite on the axis because f is EVEN in r: f' = O(r),
// f'' - f'/r = O(r^2), and so on. Near the axis the quotients cancel to
// rounding, so an EXACT EVEN TAYLOR SERIES is used instead. No finite
// difference executes anywhere in this path.
//
// THE AXIS COEFFICIENTS ARE CLOSED FORM, NOT DIFFERENCED. f is even, so
//
//     f(r)  = sum_k c_2k r^2k,        c_2k = f^(2k)(0) / (2k)!
//     F(w)  = sum_k c_2k R0^2k w^k
//
// and therefore F^(k)(0) = k! c_2k R0^2k, i.e.
//
//     F'(0)   = R0^2 f''(0)  / 2
//     F''(0)  = R0^4 f''''(0)/ 12
//     F'''(0) = R0^6 f6(0)   / 120
//
// The even derivatives of f at r = 0 follow from
// f^(n)(0) = sigma^n g^(n)(sigma R) / T with g = tanh and T = tanh(sigma R),
// because g^(n) is ODD for even n. Expanding g^(n) as a polynomial in tanh via
// the recurrence p -> p'(t)(1-t^2) gives the coefficients below; they were
// checked against that recurrence symbolically and against a high-precision
// numerical derivative.
//
// AN EARLIER VERSION SET F3 = 0 HERE AND OBTAINED F2 FROM A CENTRED DIFFERENCE
// OF THE ANALYTIC THIRD DERIVATIVE. Both were wrong. F'''(0) is not zero for
// this shape function -- it is -0.0634 at sigma = 1, -0.306 at sigma = 3 and
// -9.21 at R = 1.5, sigma = 2 -- so the third design derivative was being
// silently truncated on the axis, and a finite difference was running inside a
// production metric path that claims to be analytic throughout.
//
// WHERE THE SWITCH GOES, AND WHY IT IS NOT ARBITRARY.
//
// The quotient branch above suffers catastrophic cancellation as r -> 0: at
// sigma r = 0.02 the term 3 f'/r^5 is of order 1e4 while F''' is of order 1e-1,
// so five or six digits are lost. The Taylor branch has the opposite behaviour,
// degrading as w grows. Measured against a high-order reference series, the two
// error curves cross near sigma r = 0.1 for every (R, sigma) in use, where each
// branch is accurate to about 1e-7 relative or better. That is the switch, and
// the gate measures the residual jump rather than assuming it is negligible.
inline constexpr double kAxisSigmaRadius = 0.1;
inline constexpr int kAxisTaylorTerms = 6;   // through w^6, i.e. f^(12)(0)

// The n-th derivative of tanh at u, as a polynomial in T = tanh(u).
//
// p_1(t) = 1 - t^2 and p_{n+1}(t) = p_n'(t) (1 - t^2). Generated rather than
// transcribed: six hand-copied coefficient lists is six chances to lose a digit,
// and this repository has already lost one that way. The explicit closed forms
// for f'', f'''' and f^(6) are retained in the gate as an INDEPENDENT check
// that this recurrence reproduces them.
inline double tanh_derivative_poly(int n, double T) {
    constexpr int kMax = 20;
    double p[kMax + 3] = {0.0};
    int deg = 2;
    p[0] = 1.0; p[2] = -1.0;                 // 1 - t^2
    for (int k = 1; k < n; ++k) {
        double d[kMax + 3] = {0.0};
        for (int i = 1; i <= deg; ++i) d[i - 1] = static_cast<double>(i) * p[i];
        const int ddeg = deg - 1;
        double q[kMax + 3] = {0.0};
        for (int i = 0; i <= ddeg; ++i) { q[i] += d[i]; q[i + 2] -= d[i]; }
        deg = ddeg + 2;
        for (int i = 0; i <= deg; ++i) p[i] = q[i];
    }
    double acc = 0.0;
    for (int i = deg; i >= 0; --i) acc = acc * T + p[i];
    return acc;
}

// f^(n)(0) for EVEN n.
//
// f(r) = [g(sigma(r+R)) - g(sigma(r-R))] / (2 tanh(sigma R)) with g = tanh, so
// f^(n)(0) = sigma^n [g^(n)(sigma R) - g^(n)(-sigma R)] / (2T). For even n,
// g^(n) is ODD, so the bracket is 2 g^(n)(sigma R) and
//
//     f^(n)(0) = sigma^n g^(n)(sigma R) / T.
inline double shape_even_derivative_at_axis(int n, double R, double sigma) {
    const double T = std::tanh(sigma * R);
    return std::pow(sigma, n) * tanh_derivative_poly(n, T) / T;
}

// The explicit closed forms the specification supplies, kept so the gate can
// check the recurrence against them rather than trusting either alone.
struct ShapeAxisDerivatives {
    double f2 = 0.0, f4 = 0.0, f6 = 0.0;
};

inline ShapeAxisDerivatives shape_axis_derivatives_closed_form(double R,
                                                               double sigma) {
    const double T = std::tanh(sigma * R);
    const double T2 = T * T, T4 = T2 * T2, T6 = T4 * T2;
    const double s2 = sigma * sigma, s4 = s2 * s2, s6 = s4 * s2;
    ShapeAxisDerivatives d;
    d.f2 = -2.0 * s2 * (1.0 - T2);
    d.f4 = s4 * (16.0 - 40.0 * T2 + 24.0 * T4);
    d.f6 = s6 * (-272.0 + 1232.0 * T2 - 1680.0 * T4 + 720.0 * T6);
    return d;
}

// a_k = R0^(2k) f^(2k)(0) / (2k)!, the even Taylor coefficients of F in w.
inline double shape_w_taylor_coefficient(int k, double R0, double R,
                                         double sigma) {
    if (k == 0) return 1.0;                  // f(0) = 1 exactly
    double fact = 1.0;
    for (int i = 2; i <= 2 * k; ++i) fact *= static_cast<double>(i);
    return std::pow(R0, 2 * k) *
           shape_even_derivative_at_axis(2 * k, R, sigma) / fact;
}

struct ShapeWJet {
    double F = 1.0;
    double F1 = 0.0;
    double F2 = 0.0;
    double F3 = 0.0;
    bool axis_branch = false;
};

inline ShapeWJet shape_w_jet(double w, double R0, double R, double sigma) {
    ShapeWJet j;
    const double r = R0 * std::sqrt(std::fmax(w, 0.0));
    if (sigma * r < kAxisSigmaRadius) {
        // The exact even Taylor series in w, carried to w^kAxisTaylorTerms --
        // three orders beyond the third derivative -- so that F, F', F'' and
        // F''' are all continuous across the branch to rounding.
        //
        //   F    = sum_k       a_k w^k
        //   F'   = sum_k     k a_k w^(k-1)
        //   F''  = sum_k k(k-1) a_k w^(k-2)
        //   F''' = sum_k k(k-1)(k-2) a_k w^(k-3)
        //
        // with a_k = R0^2k f^(2k)(0) / (2k)!, evaluated by Horner below.
        double a[kAxisTaylorTerms + 1];
        for (int k = 0; k <= kAxisTaylorTerms; ++k)
            a[k] = shape_w_taylor_coefficient(k, R0, R, sigma);
        double F = 0.0, F1 = 0.0, F2 = 0.0, F3 = 0.0;
        for (int k = kAxisTaylorTerms; k >= 0; --k) {
            F = F * w + a[k];
            if (k >= 1) F1 = F1 * w + static_cast<double>(k) * a[k];
            if (k >= 2) F2 = F2 * w + static_cast<double>(k * (k - 1)) * a[k];
            if (k >= 3)
                F3 = F3 * w + static_cast<double>(k * (k - 1) * (k - 2)) * a[k];
        }
        j.F = F;
        j.F1 = F1;
        j.F2 = F2;
        j.F3 = F3;
        j.axis_branch = true;
        return j;
    }
    const ShapeRadialJet s = shape_radial_jet(r, R, sigma);
    const double r2 = r * r, r3 = r2 * r, r4 = r3 * r, r5 = r4 * r;
    const double R2 = R0 * R0;
    j.F = s.f;
    j.F1 = 0.5 * R2 * s.d1 / r;
    j.F2 = 0.25 * R2 * R2 * (s.d2 / r2 - s.d1 / r3);
    j.F3 = 0.125 * R2 * R2 * R2 *
           (s.d3 / r3 - 3.0 * s.d2 / r4 + 3.0 * s.d1 / r5);
    return j;
}

// The exact axis limits, for the gate to check directly rather than through the
// branch.
inline ShapeWJet shape_w_jet_axis_limits(double R0, double R, double sigma) {
    // F^(k)(0) = k! a_k. Equivalently F'(0) = R0^2 f''(0)/2,
    // F''(0) = R0^4 f4(0)/12 and F'''(0) = R0^6 f6(0)/120, which is the form the
    // specification states and the gate checks against the closed forms.
    ShapeWJet j;
    j.F = 1.0;
    j.F1 = shape_w_taylor_coefficient(1, R0, R, sigma);
    j.F2 = 2.0 * shape_w_taylor_coefficient(2, R0, R, sigma);
    j.F3 = 6.0 * shape_w_taylor_coefficient(3, R0, R, sigma);
    j.axis_branch = true;
    return j;
}

// The overlapping orders, checked against production so the closed-form
// extension cannot drift from the function it extends.
inline double radial_jet_matches_production(double r, double R, double sigma) {
    metrics::AlcubierreParameters<double> p(1.0, R, sigma);
    const ShapeRadialJet j = shape_radial_jet(r, R, sigma);
    const double f_prod = metrics::shape_function(r, p);
    const double d1_prod = metrics::shape_function_radial_derivative(r, p);
    return std::fmax(std::fabs(j.f - f_prod), std::fabs(j.d1 - d1_prod));
}

// --- Design-coordinate derivative table --------------------------------------

// D[p][q] = d^p/dxi^p d^q/ds^q of a scalar design field, for p + q <= 3.
// Storage is a 4x4 table with the p+q > 3 entries left zero and never read.
struct DesignDerivativeTable {
    double d[4][4] = {{0}};
    double& at(int p, int q) { return d[p][q]; }
    double at(int p, int q) const { return d[p][q]; }
};

// The baseline beta^x = -v_s f, expressed through w = xi^2 + s.
//
//   dw/dxi = 2 xi,  d2w/dxi2 = 2,  dw/ds = 1, everything else zero
//
// so, writing F for the w-jet,
//
//   q_xi    = 2 xi F'                 q_s     = F'
//   q_xixi  = 2 F' + 4 xi^2 F''       q_xis   = 2 xi F''      q_ss  = F''
//   q_xixixi= 12 xi F'' + 8 xi^3 F''' q_xixis = 2 F'' + 4 xi^2 F'''
//   q_xiss  = 2 xi F'''               q_sss   = F'''
inline DesignDerivativeTable baseline_shift_table(double xi, double s,
                                                  double v_s, double R0,
                                                  double R, double sigma) {
    const double w = xi * xi + s;
    const ShapeWJet F = shape_w_jet(w, R0, R, sigma);
    const double k = -v_s;   // beta^x = -v_s f: the sign read from production
    DesignDerivativeTable t;
    t.at(0, 0) = k * F.F;
    t.at(1, 0) = k * (2.0 * xi * F.F1);
    t.at(0, 1) = k * F.F1;
    t.at(2, 0) = k * (2.0 * F.F1 + 4.0 * xi * xi * F.F2);
    t.at(1, 1) = k * (2.0 * xi * F.F2);
    t.at(0, 2) = k * F.F2;
    t.at(3, 0) = k * (12.0 * xi * F.F2 + 8.0 * xi * xi * xi * F.F3);
    t.at(2, 1) = k * (2.0 * F.F2 + 4.0 * xi * xi * F.F3);
    t.at(1, 2) = k * (2.0 * xi * F.F3);
    t.at(0, 3) = k * F.F3;
    return t;
}

// The F1 correction, from the O0A basis. Only the design-coordinate channels
// are taken; the Cartesian chain below is shared with the baseline, so both
// travel the same verified route.
inline DesignDerivativeTable correction_shift_table(
    const AxisymmetricBasis& basis, const std::vector<double>& coefficients,
    const DesignCoordinates& p) {
    DesignDerivativeTable t;
    const TensorBasisEvaluation e = basis.evaluate(p, /*want_third=*/true);
    for (std::size_t m = 0; m < e.size(); ++m) {
        const double c = coefficients[static_cast<std::size_t>(e.flat_index[m])];
        t.at(0, 0) += c * e.value[m];
        t.at(1, 0) += c * e.d_xi[m];
        t.at(0, 1) += c * e.d_s[m];
        t.at(2, 0) += c * e.d_xi_xi[m];
        t.at(1, 1) += c * e.d_xi_s[m];
        t.at(0, 2) += c * e.d_s_s[m];
        t.at(3, 0) += c * e.d3_xi[m];
        t.at(0, 3) += c * e.d3_s[m];
    }
    // The two MIXED third derivatives are not among O0A's exported channels, so
    // they are formed here from O0A's own univariate kernel -- the same
    // functions, at the same point, one order higher in the other variable.
    // Nothing about O0A is modified or re-derived.
    for (std::size_t m = 0; m < e.size(); ++m) {
        const double c = coefficients[static_cast<std::size_t>(e.flat_index[m])];
        const double Nxi2 =
            basis.axial().derivative_of(e.index_xi[m], p.xi, 2);
        const double Nxi1 =
            basis.axial().derivative_of(e.index_xi[m], p.xi, 1);
        const double Ms1 =
            basis.radial().derivative_of(e.index_s[m], p.radial_squared, 1);
        const double Ms2 =
            basis.radial().derivative_of(e.index_s[m], p.radial_squared, 2);
        t.at(2, 1) += c * Nxi2 * Ms1;
        t.at(1, 2) += c * Nxi1 * Ms2;
    }
    return t;
}

// --- The spacetime chain rule -------------------------------------------------
//
// xi = (x - x_c(t))/R0 with x_c = x_c0 + v_s t, and s = (y^2+z^2)/R0^2. So
//
//     d/dt -> (-v_s/R0) d/dxi        d/dx -> (1/R0) d/dxi
//     d/dy -> a d/ds                 d/dz -> b d/ds,     a = 2y/R0^2,
//                                                        b = 2z/R0^2,
//     d2s/dy2 = d2s/dz2 = c = 2/R0^2,  d2s/dydz = 0, higher s-derivatives zero.
//
// An index from {t, x} therefore RAISES THE xi ORDER and multiplies by its own
// constant; an index from {y, z} acts through s and picks up the a, b, c terms.
// That factorisation is what lets one routine serve every channel up to third
// order without a separate hand-derived formula per component.

struct ChainFactors {
    double kt = 0.0;   // -v_s/R0
    double kx = 0.0;   // 1/R0
    double a = 0.0;    // ds/dy
    double b = 0.0;    // ds/dz
    double c = 0.0;    // d2s/dy2 = d2s/dz2
};

inline ChainFactors chain_factors(double y, double z, double v_s, double R0) {
    ChainFactors f;
    f.kt = -v_s / R0;
    f.kx = 1.0 / R0;
    f.a = 2.0 * y / (R0 * R0);
    f.b = 2.0 * z / (R0 * R0);
    f.c = 2.0 / (R0 * R0);
    return f;
}

// One spacetime derivative of a design field, for a multi-index of length 0..3.
// `idx` holds spacetime indices (0=t, 1=x, 2=y, 3=z); order within the
// multi-index is irrelevant, which is exactly the symmetry the gate checks.
inline double chained_derivative(const DesignDerivativeTable& T,
                                 const ChainFactors& F, const int* idx,
                                 int n) {
    // Split the multi-index into the xi-acting part {t,x} and the s-acting part
    // {y,z}, accumulating the constant from each xi-acting index.
    int n_xi = 0;
    double k = 1.0;
    int n_y = 0, n_z = 0;
    for (int i = 0; i < n; ++i) {
        if (idx[i] == kTimeIndex) { ++n_xi; k *= F.kt; }
        else if (idx[i] == kIndexX) { ++n_xi; k *= F.kx; }
        else if (idx[i] == kIndexY) ++n_y;
        else ++n_z;
    }
    const int n_s = n_y + n_z;
    if (n_xi + n_s > 3) return 0.0;

    // The s-chain applied to (d^n_xi/dxi^n_xi) q. Derived once here rather than
    // per component:
    //   0 s-indices : q
    //   1           : q_s * g          g = a or b
    //   2           : q_ss * g1 g2 + q_s * (c if the two indices match, else 0)
    //   3           : q_sss * g1 g2 g3
    //                 + q_ss * ( c * the g of the index NOT in the matching pair,
    //                            summed over the three pairings )
    const double q0 = T.at(n_xi, 0);
    const double q1 = (n_xi + 1 <= 3) ? T.at(n_xi, 1) : 0.0;
    const double q2 = (n_xi + 2 <= 3) ? T.at(n_xi, 2) : 0.0;
    const double q3 = (n_xi + 3 <= 3) ? T.at(n_xi, 3) : 0.0;

    if (n_s == 0) return k * q0;
    if (n_s == 1) return k * q1 * (n_y == 1 ? F.a : F.b);
    if (n_s == 2) {
        const double g1 = (n_y >= 1) ? F.a : F.b;
        const double g2 = (n_y == 2) ? F.a : F.b;
        // d2s/dy dz = 0, so the q_s term survives only for yy or zz.
        const double same = (n_y == 2 || n_z == 2) ? F.c : 0.0;
        return k * (q2 * g1 * g2 + q1 * same);
    }
    // n_s == 3
    double g1, g2, g3;
    if (n_y == 3) { g1 = g2 = g3 = F.a; }
    else if (n_y == 2) { g1 = g2 = F.a; g3 = F.b; }
    else if (n_y == 1) { g1 = F.a; g2 = g3 = F.b; }
    else { g1 = g2 = g3 = F.b; }
    // The q_ss term: sum over the three ways of pairing two of the indices,
    // each contributing c only when that pair is yy or zz, times the g of the
    // remaining index.
    double pair_term = 0.0;
    if (n_y == 3) pair_term = 3.0 * F.c * F.a;
    else if (n_y == 2) pair_term = F.c * F.b;              // the yy pair only
    else if (n_y == 1) pair_term = F.c * F.a;              // the zz pair only
    else pair_term = 3.0 * F.c * F.b;
    return k * (q3 * g1 * g2 * g3 + q2 * pair_term);
}

// --- The F1 family -----------------------------------------------------------

struct ShiftFamilyConfig {
    double velocity = 0.6;
    double bubble_radius = 2.0;
    double wall_sigma = 1.0;
    double scale_R0 = 2.0;
    double bubble_center_x0 = 0.0;
};

// beta^x and its full spacetime derivative jet, baseline plus correction.
struct ShiftJet {
    double value = 0.0;
    double d1[4] = {0, 0, 0, 0};
    double d2[4][4] = {{0}};
    double d3[4][4][4] = {{{0}}};
};

inline ShiftJet shift_jet(const ShiftFamilyConfig& cfg,
                          const AxisymmetricBasis& basis,
                          const std::vector<double>& coefficients,
                          double t, double x, double y, double z,
                          bool with_correction = true) {
    const double R0 = cfg.scale_R0;
    const double xc = cfg.bubble_center_x0 + cfg.velocity * t;
    const DesignCoordinates p =
        make_design_coordinates(x - xc, y, z, 0.0, R0);

    DesignDerivativeTable T = baseline_shift_table(
        p.xi, p.radial_squared, cfg.velocity, R0, cfg.bubble_radius,
        cfg.wall_sigma);
    if (with_correction && !coefficients.empty()) {
        const DesignDerivativeTable C =
            correction_shift_table(basis, coefficients, p);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) T.at(i, j) += C.at(i, j);
    }

    const ChainFactors F = chain_factors(y, z, cfg.velocity, R0);
    ShiftJet j;
    j.value = T.at(0, 0);
    for (int m = 0; m < 4; ++m) {
        const int one[1] = {m};
        j.d1[m] = chained_derivative(T, F, one, 1);
        for (int n = 0; n < 4; ++n) {
            const int two[2] = {m, n};
            j.d2[m][n] = chained_derivative(T, F, two, 2);
            for (int l = 0; l < 4; ++l) {
                const int three[3] = {m, n, l};
                j.d3[m][n][l] = chained_derivative(T, F, three, 3);
            }
        }
    }
    return j;
}

// --- Scalar-generic coefficient path -----------------------------------------
//
// For coefficient derivatives the geometry stays double and only the
// COEFFICIENTS carry channels, exactly as in O0A. The baseline has no
// coefficient dependence at all, so it enters as a constant; the correction is
// the linear form sum_k c_k B_k, and the metric components are built from it.
//
// beta^x is LINEAR in c, so d beta/d c_k = B_k exactly. g_01 = beta^x is
// likewise linear. g_00 = -1 + (beta^x)^2 is QUADRATIC, so d g_00/d c_k =
// 2 beta^x B_k -- still an exact reference, and a more searching one than a
// linear identity because it mixes the baseline into the derivative.
template <typename Scalar>
Scalar shift_value_scalar(const ShiftFamilyConfig& cfg,
                          const AxisymmetricBasis& basis,
                          const std::vector<Scalar>& coefficients,
                          double t, double x, double y, double z) {
    const double R0 = cfg.scale_R0;
    const double xc = cfg.bubble_center_x0 + cfg.velocity * t;
    const DesignCoordinates p = make_design_coordinates(x - xc, y, z, 0.0, R0);
    const double w = p.xi * p.xi + p.radial_squared;
    const ShapeWJet F = shape_w_jet(w, R0, cfg.bubble_radius, cfg.wall_sigma);
    Scalar acc = scalar_cast<Scalar>(-cfg.velocity * F.F);
    if (coefficients.empty()) return acc;
    const TensorBasisEvaluation e = basis.evaluate(p);
    for (std::size_t m = 0; m < e.size(); ++m)
        acc = acc + coefficients[static_cast<std::size_t>(e.flat_index[m])] *
                        scalar_cast<Scalar>(e.value[m]);
    return acc;
}

// One covariant metric component, scalar-generic in the coefficients.
template <typename Scalar>
Scalar metric_component_scalar(const ShiftFamilyConfig& cfg,
                               const AxisymmetricBasis& basis,
                               const std::vector<Scalar>& coefficients,
                               double t, double x, double y, double z,
                               int A, int B) {
    const Scalar b =
        shift_value_scalar<Scalar>(cfg, basis, coefficients, t, x, y, z);
    if (A == 0 && B == 0)
        return scalar_cast<Scalar>(-1.0) + b * b;
    if ((A == 0 && B == 1) || (A == 1 && B == 0)) return b;
    if (A == B) return scalar_cast<Scalar>(1.0);
    return scalar_cast<Scalar>(0.0);
}

// d_mu of one covariant component, scalar-generic. Only the channels the
// gradient cases need; the full jet stays double, because a jet over an OTI
// scalar would multiply the working set by the channel count for no gain.
template <typename Scalar>
Scalar metric_first_derivative_scalar(const ShiftFamilyConfig& cfg,
                                      const AxisymmetricBasis& basis,
                                      const std::vector<Scalar>& coefficients,
                                      double t, double x, double y, double z,
                                      int mu, int A, int B) {
    const double R0 = cfg.scale_R0;
    const double xc = cfg.bubble_center_x0 + cfg.velocity * t;
    const DesignCoordinates p = make_design_coordinates(x - xc, y, z, 0.0, R0);
    const double w = p.xi * p.xi + p.radial_squared;
    const ShapeWJet F = shape_w_jet(w, R0, cfg.bubble_radius, cfg.wall_sigma);
    const ChainFactors CF = chain_factors(y, z, cfg.velocity, R0);

    // The baseline's design-derivative table, as constants.
    DesignDerivativeTable T = baseline_shift_table(
        p.xi, p.radial_squared, cfg.velocity, R0, cfg.bubble_radius,
        cfg.wall_sigma);
    (void)F;
    const int one[1] = {mu};
    Scalar dbeta = scalar_cast<Scalar>(chained_derivative(T, CF, one, 1));
    Scalar beta = scalar_cast<Scalar>(T.at(0, 0));

    if (!coefficients.empty()) {
        const TensorBasisEvaluation e = basis.evaluate(p, true);
        for (std::size_t m = 0; m < e.size(); ++m) {
            const Scalar& c =
                coefficients[static_cast<std::size_t>(e.flat_index[m])];
            DesignDerivativeTable Bt;
            Bt.at(0, 0) = e.value[m];
            Bt.at(1, 0) = e.d_xi[m];
            Bt.at(0, 1) = e.d_s[m];
            beta = beta + c * scalar_cast<Scalar>(e.value[m]);
            dbeta = dbeta +
                    c * scalar_cast<Scalar>(chained_derivative(Bt, CF, one, 1));
        }
    }
    if (A == 0 && B == 0) return scalar_cast<Scalar>(2.0) * beta * dbeta;
    if ((A == 0 && B == 1) || (A == 1 && B == 0)) return dbeta;
    return scalar_cast<Scalar>(0.0);
}

// The metric jet of the F1 family. With alpha = 1 and gamma = I the assembly
// collapses to
//
//     g_00 = -1 + (beta^x)^2,   g_01 = beta^x,   g_ij = delta_ij
//
// so every derivative follows from the shift jet by the product rule on
// (beta^x)^2. Note g_01 = beta^x here because gamma = I makes beta_x = beta^x;
// the two are kept distinct in adm_design.hpp precisely because that
// coincidence does not survive a nontrivial spatial metric.
inline MetricJet<double> f1_metric_jet(const ShiftFamilyConfig& cfg,
                                       const AxisymmetricBasis& basis,
                                       const std::vector<double>& coefficients,
                                       double t, double x, double y, double z,
                                       bool with_correction = true) {
    const ShiftJet b =
        shift_jet(cfg, basis, coefficients, t, x, y, z, with_correction);

    ADMDesignFields<double> adm = ADMDesignFields<double>::minkowski();
    adm.shift_contravariant[0] = b.value;

    MetricJet<double> jet;
    jet.g_cov = assemble_covariant(adm);
    bool ok = false;
    jet.g_contra = assemble_contravariant(adm, ok);
    jet.valid = ok;
    jet.rejection = ok ? MetricRejection::None
                       : MetricRejection::SingularSpatialMetric;

    auto zero = []() {
        CovariantMetric<double> m;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) m.m(i, j) = 0.0;
        return m;
    };
    for (int m = 0; m < 4; ++m) {
        jet.dg[static_cast<std::size_t>(m)] = zero();
        // d g_00 = 2 beta d beta ; d g_01 = d beta ; spatial block constant.
        jet.dg[static_cast<std::size_t>(m)].m(0, 0) = 2.0 * b.value * b.d1[m];
        jet.dg[static_cast<std::size_t>(m)].m(0, 1) = b.d1[m];
        jet.dg[static_cast<std::size_t>(m)].m(1, 0) = b.d1[m];
        for (int n = 0; n < 4; ++n) {
            jet.d2g[static_cast<std::size_t>(m)][static_cast<std::size_t>(n)] =
                zero();
            auto& d2 = jet.d2g[static_cast<std::size_t>(m)]
                              [static_cast<std::size_t>(n)];
            d2.m(0, 0) = 2.0 * (b.d1[m] * b.d1[n] + b.value * b.d2[m][n]);
            d2.m(0, 1) = b.d2[m][n];
            d2.m(1, 0) = b.d2[m][n];
            for (int l = 0; l < 4; ++l) {
                jet.d3g[static_cast<std::size_t>(m)][static_cast<std::size_t>(n)]
                       [static_cast<std::size_t>(l)] = zero();
                auto& d3 = jet.d3g[static_cast<std::size_t>(m)]
                                  [static_cast<std::size_t>(n)]
                                  [static_cast<std::size_t>(l)];
                d3.m(0, 0) = 2.0 * (b.d2[m][n] * b.d1[l] +
                                    b.d1[m] * b.d2[n][l] +
                                    b.d1[n] * b.d2[m][l] +
                                    b.value * b.d3[m][n][l]);
                d3.m(0, 1) = b.d3[m][n][l];
                d3.m(1, 0) = b.d3[m][n][l];
            }
        }
    }
    return jet;
}

}  // namespace design
}  // namespace warp_hypad

#pragma once

#include <algorithm>
#include <vector>

#include "program_b/inverse_design/complex_scalar.hpp"

// ---------------------------------------------------------------------------
// Robust real-and-complex root finding for a quartic.
//
// WHY THE PREVIOUS SCAN WAS WRONG. B4 located roots by sampling the polynomial
// and bracketing SIGN CHANGES. A root of EVEN multiplicity does not change the
// sign, so a doubled root is invisible to that method -- and the doubled root is
// precisely the metric-induced case, the one the physics cares about most. The
// B4 factorization test passed only because it compared lambda against a
// separately computed cone rather than asking the root finder for the branches.
// The scan also could not see complex roots at all, so evanescent branches were
// silently absent rather than reported.
//
// THE METHOD. Durand-Kerner (Weierstrass) simultaneous iteration:
//
//     z_i  <-  z_i - p(z_i) / ( a_n prod_{j != i} (z_i - z_j) )
//
// It converges to ALL roots at once, real and complex alike, from deterministic
// distinct complex seeds. Repeated roots converge more slowly -- the attainable
// accuracy for a root of multiplicity m is about eps^(1/m), which for a double
// root is about 1e-8 -- so every root is returned with its residual and a
// conditioning warning rather than a bare number.
//
// Multiplicity is estimated by CLUSTERING: roots closer together than the
// accuracy their own residuals imply are one root of higher multiplicity. That is
// the honest reading of a numerically split double root, and case 4 of the gate
// exercises a deliberately near-split pair.
//
// Real parts only, and deliberately: root finding involves comparisons,
// clustering and sorting throughout. It is a DIAGNOSTIC and branch-selection
// step, and it is kept out of any differentiated path.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

struct PolynomialRoot {
    double real = 0.0;
    double imaginary = 0.0;
    double residual = 0.0;        // |p(z)|, normalised by the coefficient scale
    int multiplicity = 1;
    bool is_real = false;
    bool ill_conditioned = false; // multiplicity > 1: accuracy ~ eps^(1/m)
};

struct RootSet {
    std::vector<PolynomialRoot> roots;     // all four, complex included
    std::vector<PolynomialRoot> real_roots;
    int real_count = 0;
    int complex_pair_count = 0;
    int distinct_count = 0;
    double worst_residual = 0.0;
    bool converged = false;
    int iterations = 0;
    std::string warning = "none";
};

namespace detail {

using RC = Complex<double>;

inline RC poly_eval(const double c[5], const RC& z) {
    // Horner, degree 4.
    RC acc = complex_from_real<double>(c[4]);
    for (int i = 3; i >= 0; --i)
        acc = acc * z + complex_from_real<double>(c[i]);
    return acc;
}

inline RC poly_derivative(const double c[5], const RC& z) {
    const double d[4] = {c[1], 2.0 * c[2], 3.0 * c[3], 4.0 * c[4]};
    RC acc = complex_from_real<double>(d[3]);
    for (int i = 2; i >= 0; --i)
        acc = acc * z + complex_from_real<double>(d[i]);
    return acc;
}

}  // namespace detail

// Roots of c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4.
inline RootSet quartic_roots(const double c_in[5]) {
    using detail::RC;
    RootSet out;

    double c[5];
    double scale = 0.0;
    for (int i = 0; i < 5; ++i) scale = std::fmax(scale, std::fabs(c_in[i]));
    if (scale <= 0.0) {
        out.warning = "the polynomial is identically zero";
        return out;
    }
    for (int i = 0; i < 5; ++i) c[i] = c_in[i] / scale;

    // Degree reduction if the leading coefficients vanish. Reported, not hidden.
    int degree = 4;
    while (degree > 0 && std::fabs(c[degree]) < 1e-14) --degree;
    if (degree < 1) {
        out.warning = "the polynomial has no roots (degree 0 after reduction)";
        return out;
    }
    if (degree < 4)
        out.warning = "leading coefficient(s) negligible: degree reduced to " +
                      std::to_string(degree);

    // Deterministic distinct seeds on a spiral, so a rerun is bit-identical.
    std::vector<RC> z(static_cast<std::size_t>(degree));
    const RC seed(0.4, 0.9);
    RC p = complex_from_real<double>(1.0);
    for (int i = 0; i < degree; ++i) {
        p = p * seed;
        z[static_cast<std::size_t>(i)] = p;
    }

    // Durand-Kerner.
    const int max_iter = 2000;
    for (out.iterations = 0; out.iterations < max_iter; ++out.iterations) {
        double move = 0.0;
        for (int i = 0; i < degree; ++i) {
            RC denom = complex_from_real<double>(c[degree]);
            for (int j = 0; j < degree; ++j)
                if (j != i) denom = denom * (z[i] - z[j]);
            const double dn = norm_squared(denom);
            if (dn < 1e-300) continue;   // coincident seeds; the next sweep fixes it
            const RC step = detail::poly_eval(c, z[i]) / denom;
            z[i] = z[i] - step;
            move = std::fmax(move, std::sqrt(norm_squared(step)));
        }
        if (move < 1e-15) { out.converged = true; break; }
    }

    // Newton polish, which sharpens simple roots and does no harm to multiple
    // ones.
    for (int i = 0; i < degree; ++i)
        for (int it = 0; it < 60; ++it) {
            const RC f = detail::poly_eval(c, z[i]);
            const RC fp = detail::poly_derivative(c, z[i]);
            if (norm_squared(fp) < 1e-300) break;
            const RC step = f / fp;
            z[i] = z[i] - step;
            if (norm_squared(step) < 1e-32) break;
        }

    // Sort by real part, then imaginary, so the output is deterministic.
    std::sort(z.begin(), z.end(), [](const RC& a, const RC& b) {
        if (a.re != b.re) return a.re < b.re;
        return a.im < b.im;
    });

    for (int i = 0; i < degree; ++i) {
        PolynomialRoot r;
        r.real = z[i].re;
        r.imaginary = z[i].im;
        r.residual = std::sqrt(norm_squared(detail::poly_eval(c, z[i])));
        out.roots.push_back(r);
        out.worst_residual = std::fmax(out.worst_residual, r.residual);
    }

    // --- Multiplicity by clustering ------------------------------------------
    //
    // The tolerance is set by the residuals themselves: a root pair separated by
    // less than the accuracy their own residuals imply cannot be resolved, and
    // calling them one double root is the honest reading.
    const double cluster_tol =
        std::fmax(1e-7, 4.0 * std::sqrt(std::fmax(out.worst_residual, 1e-32)));

    std::vector<int> label(out.roots.size(), -1);
    int clusters = 0;
    for (std::size_t i = 0; i < out.roots.size(); ++i) {
        if (label[i] >= 0) continue;
        label[i] = clusters;
        for (std::size_t j = i + 1; j < out.roots.size(); ++j) {
            if (label[j] >= 0) continue;
            const double dr = out.roots[i].real - out.roots[j].real;
            const double di = out.roots[i].imaginary - out.roots[j].imaginary;
            if (std::sqrt(dr * dr + di * di) < cluster_tol) label[j] = clusters;
        }
        ++clusters;
    }
    out.distinct_count = clusters;
    for (std::size_t i = 0; i < out.roots.size(); ++i) {
        int m = 0;
        for (std::size_t j = 0; j < out.roots.size(); ++j)
            if (label[j] == label[i]) ++m;
        out.roots[i].multiplicity = m;
        out.roots[i].ill_conditioned = m > 1;
    }

    // --- Real versus complex -------------------------------------------------
    //
    // A root counts as real when its imaginary part is below the accuracy its
    // multiplicity permits, NOT below a fixed epsilon: a double real root
    // computed to 1e-8 legitimately carries an imaginary part of that size.
    for (PolynomialRoot& r : out.roots) {
        const double imag_tol =
            cluster_tol * static_cast<double>(r.multiplicity);
        r.is_real = std::fabs(r.imaginary) <= imag_tol;
        if (r.is_real) {
            ++out.real_count;
            out.real_roots.push_back(r);
        }
    }
    out.complex_pair_count = (degree - out.real_count) / 2;

    if (out.worst_residual > 1e-6 && out.warning == "none")
        out.warning = "large residual: the polynomial may be badly scaled or "
                      "have a high-multiplicity root";
    else if (out.distinct_count < degree && out.warning == "none")
        out.warning = "repeated root detected: attainable accuracy is about "
                      "eps^(1/m), so about 1e-8 for a double root";
    return out;
}

// Rebuild the coefficients from the roots, as an independent check that the set
// returned really is the root set of the polynomial handed in.
inline void coefficients_from_roots(const RootSet& rs, double leading,
                                    double c_out[5]) {
    using detail::RC;
    std::vector<RC> poly;
    poly.push_back(complex_from_real<double>(1.0));
    for (const PolynomialRoot& r : rs.roots) {
        const RC root(r.real, r.imaginary);
        std::vector<RC> next(poly.size() + 1);
        for (std::size_t i = 0; i < poly.size(); ++i) {
            next[i + 1] = next[i + 1] + poly[i];
            next[i] = next[i] - poly[i] * root;
        }
        poly.swap(next);
    }
    for (int i = 0; i < 5; ++i) c_out[i] = 0.0;
    for (std::size_t i = 0; i < poly.size() && i < 5; ++i)
        c_out[i] = leading * poly[i].re;
}

}  // namespace inverse_design
}  // namespace program_b

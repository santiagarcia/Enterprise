#pragma once

#include <cmath>
#include <string>

#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// A paired-real complex scalar, templated on the underlying real scalar.
//
// WHY NOT std::complex. std::complex<T> has undefined behaviour for T outside
// {float, double, long double}, and in practice its operator/ and abs() use
// comparisons and scaling tricks that branch on the value -- exactly what a
// hypercomplex AD type cannot tolerate. This type is written from the
// definitions instead, with no comparison on a Scalar anywhere.
//
// WHY NOT ASK OTIlib FOR A COMPLEX AD SCALAR. OTIlib does ship a complex module,
// but the adapter in this project deliberately builds only the scalar onumm3n2
// path, and nothing outside warp_hypad/sensitivity/ may name an oti/ symbol.
// Pairing two REAL AD scalars gives complex arithmetic whose derivative channels
// flow through the real and imaginary parts independently, which is precisely
// what a real-valued objective differentiated with respect to real physical
// parameters needs. No assumption is made about upstream complex support.
//
//     Complex<double>      plain complex arithmetic
//     Complex<OTIScalar>   complex arithmetic carrying derivative channels
//
// CONJUGATION AND ADJOINT, kept distinct because for a lossy medium they are
// different statements:
//
//     conj(z)              complex conjugate of a scalar
//     hermitian_adjoint(M) conjugate transpose  -- energy / passivity
//     transpose(M)         plain transpose      -- Lorentz reciprocity
//
// Reciprocity is the UNCONJUGATED transpose. Using the adjoint there would be a
// physics error that only shows up once the medium is lossy, which is exactly
// when B3 starts.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

using warp_hypad::real_part;
using warp_hypad::scalar_cast;

template <typename Scalar>
struct Complex {
    Scalar re;
    Scalar im;

    Complex() : re(scalar_cast<Scalar>(0.0)), im(scalar_cast<Scalar>(0.0)) {}
    explicit Complex(const Scalar& r)
        : re(r), im(scalar_cast<Scalar>(0.0)) {}
    Complex(const Scalar& r, const Scalar& i) : re(r), im(i) {}
};

template <typename Scalar>
Complex<Scalar> complex_from_real(double r) {
    return Complex<Scalar>(scalar_cast<Scalar>(r), scalar_cast<Scalar>(0.0));
}
template <typename Scalar>
Complex<Scalar> complex_i() {
    return Complex<Scalar>(scalar_cast<Scalar>(0.0), scalar_cast<Scalar>(1.0));
}

template <typename Scalar>
Complex<Scalar> operator+(const Complex<Scalar>& a, const Complex<Scalar>& b) {
    return Complex<Scalar>(a.re + b.re, a.im + b.im);
}
template <typename Scalar>
Complex<Scalar> operator-(const Complex<Scalar>& a, const Complex<Scalar>& b) {
    return Complex<Scalar>(a.re - b.re, a.im - b.im);
}
template <typename Scalar>
Complex<Scalar> operator-(const Complex<Scalar>& a) {
    return Complex<Scalar>(-a.re, -a.im);
}
template <typename Scalar>
Complex<Scalar> operator*(const Complex<Scalar>& a, const Complex<Scalar>& b) {
    return Complex<Scalar>(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
template <typename Scalar>
Complex<Scalar> operator*(const Complex<Scalar>& a, const Scalar& s) {
    return Complex<Scalar>(a.re * s, a.im * s);
}

// Division by the textbook formula z*w-bar / |w|^2. No Smith scaling, because
// that branches on a magnitude comparison; the models here keep denominators
// well away from zero and the gate measures how far.
template <typename Scalar>
Complex<Scalar> operator/(const Complex<Scalar>& a, const Complex<Scalar>& b) {
    const Scalar d = b.re * b.re + b.im * b.im;
    return Complex<Scalar>((a.re * b.re + a.im * b.im) / d,
                           (a.im * b.re - a.re * b.im) / d);
}

template <typename Scalar>
Complex<Scalar> conj(const Complex<Scalar>& a) {
    return Complex<Scalar>(a.re, -a.im);
}

// |z|^2, real-valued and differentiable. Preferred over abs() everywhere an
// objective is built, because sqrt at the origin is not differentiable.
template <typename Scalar>
Scalar norm_squared(const Complex<Scalar>& a) {
    return a.re * a.re + a.im * a.im;
}

// Principal square root, written so no branch is taken on a sign.
//
//     sqrt(z)_re = sqrt( (|z| + re)/2 )
//     sqrt(z)_im = im / (2 sqrt(z)_re)
//
// The second line carries the sign of im automatically, which the textbook
// sign(im)*sqrt((|z|-re)/2) form does by an explicit branch. Found by ADL, so
// the generic `using std::sqrt; sqrt(x)` in the parameter transforms picks this
// up when Scalar is Complex -- which is what lets the whole dispersion model be
// instantiated at Complex<double> for complex-step differentiation.
template <typename Scalar>
Complex<Scalar> sqrt(const Complex<Scalar>& a) {
    using std::sqrt;
    const Scalar r = sqrt(a.re * a.re + a.im * a.im);
    const Scalar sr = sqrt(scalar_cast<Scalar>(0.5) * (r + a.re));
    return Complex<Scalar>(sr, a.im / (scalar_cast<Scalar>(2.0) * sr));
}

template <typename Scalar>
Complex<Scalar> reciprocal(const Complex<Scalar>& a) {
    const Scalar d = a.re * a.re + a.im * a.im;
    return Complex<Scalar>(a.re / d, -a.im / d);
}

// --- Complex 3x3, the working size for the constitutive blocks ---------------

template <typename Scalar>
struct CMat3 {
    Complex<Scalar> m[9];

    Complex<Scalar>& operator()(int r, int c) { return m[3 * r + c]; }
    const Complex<Scalar>& operator()(int r, int c) const { return m[3 * r + c]; }

    static CMat3 zero() { return CMat3(); }
    static CMat3 identity() {
        CMat3 out;
        for (int i = 0; i < 3; ++i) out(i, i) = complex_from_real<Scalar>(1.0);
        return out;
    }
    static CMat3 scaled_identity(const Complex<Scalar>& z) {
        CMat3 out;
        for (int i = 0; i < 3; ++i) out(i, i) = z;
        return out;
    }
};

template <typename Scalar>
CMat3<Scalar> operator+(const CMat3<Scalar>& a, const CMat3<Scalar>& b) {
    CMat3<Scalar> o;
    for (int i = 0; i < 9; ++i) o.m[i] = a.m[i] + b.m[i];
    return o;
}
template <typename Scalar>
CMat3<Scalar> operator-(const CMat3<Scalar>& a, const CMat3<Scalar>& b) {
    CMat3<Scalar> o;
    for (int i = 0; i < 9; ++i) o.m[i] = a.m[i] - b.m[i];
    return o;
}
template <typename Scalar>
CMat3<Scalar> operator*(const CMat3<Scalar>& a, const CMat3<Scalar>& b) {
    CMat3<Scalar> o;
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) {
            Complex<Scalar> s;
            for (int j = 0; j < 3; ++j) s = s + a(i, j) * b(j, k);
            o(i, k) = s;
        }
    return o;
}
template <typename Scalar>
CMat3<Scalar> scale(const CMat3<Scalar>& a, const Complex<Scalar>& z) {
    CMat3<Scalar> o;
    for (int i = 0; i < 9; ++i) o.m[i] = a.m[i] * z;
    return o;
}

// PLAIN transpose. This is the one Lorentz reciprocity uses.
template <typename Scalar>
CMat3<Scalar> transpose(const CMat3<Scalar>& a) {
    CMat3<Scalar> o;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) o(i, j) = a(j, i);
    return o;
}

// CONJUGATE transpose. This is the one energy and passivity use.
template <typename Scalar>
CMat3<Scalar> hermitian_adjoint(const CMat3<Scalar>& a) {
    CMat3<Scalar> o;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) o(i, j) = conj(a(j, i));
    return o;
}

// Closed-form 3x3 inverse by the adjugate. Branch-free straight-line
// arithmetic, for the same reason Mat4::inverse uses the cofactor form: an
// AD scalar cannot tolerate a pivot selection.
template <typename Scalar>
CMat3<Scalar> invert(const CMat3<Scalar>& a, Complex<Scalar>& determinant) {
    const Complex<Scalar> c00 = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
    const Complex<Scalar> c01 = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
    const Complex<Scalar> c02 = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);
    const Complex<Scalar> c10 = a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2);
    const Complex<Scalar> c11 = a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0);
    const Complex<Scalar> c12 = a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1);
    const Complex<Scalar> c20 = a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1);
    const Complex<Scalar> c21 = a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2);
    const Complex<Scalar> c22 = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);

    determinant = a(0, 0) * c00 + a(0, 1) * c01 + a(0, 2) * c02;
    const Complex<Scalar> inv_det = reciprocal(determinant);

    CMat3<Scalar> o;
    o(0, 0) = c00 * inv_det;  o(0, 1) = c10 * inv_det;  o(0, 2) = c20 * inv_det;
    o(1, 0) = c01 * inv_det;  o(1, 1) = c11 * inv_det;  o(1, 2) = c21 * inv_det;
    o(2, 0) = c02 * inv_det;  o(2, 1) = c12 * inv_det;  o(2, 2) = c22 * inv_det;
    return o;
}

// --- Residual helpers, on real parts -----------------------------------------

template <typename Scalar>
double cmat3_max_abs(const CMat3<Scalar>& a) {
    double w = 0.0;
    for (int i = 0; i < 9; ++i) {
        w = std::fmax(w, std::fabs(real_part(a.m[i].re)));
        w = std::fmax(w, std::fabs(real_part(a.m[i].im)));
    }
    return w;
}

template <typename Scalar>
double cmat3_max_diff(const CMat3<Scalar>& a, const CMat3<Scalar>& b) {
    double w = 0.0;
    for (int i = 0; i < 9; ++i) {
        w = std::fmax(w, std::fabs(real_part(a.m[i].re) - real_part(b.m[i].re)));
        w = std::fmax(w, std::fabs(real_part(a.m[i].im) - real_part(b.m[i].im)));
    }
    return w;
}

// max |A - A^T|, the reciprocity residual for a complex block.
template <typename Scalar>
double transpose_residual(const CMat3<Scalar>& a) {
    return cmat3_max_diff(a, transpose(a));
}

// max |A - A^dagger|, the hermiticity residual.
template <typename Scalar>
double hermiticity_residual(const CMat3<Scalar>& a) {
    return cmat3_max_diff(a, hermitian_adjoint(a));
}

template <typename Scalar>
double imaginary_magnitude(const CMat3<Scalar>& a) {
    double w = 0.0;
    for (int i = 0; i < 9; ++i)
        w = std::fmax(w, std::fabs(real_part(a.m[i].im)));
    return w;
}

}  // namespace inverse_design
}  // namespace program_b

// ---------------------------------------------------------------------------
// Complex<double> as a Scalar in its own right.
//
// This is what makes COMPLEX-STEP differentiation possible. Instantiating the
// dispersion model at Scalar = Complex<double> gives quantities of type
// Complex<Complex<double>>: the INNER complex unit is the physical one that
// carries loss, the OUTER one is an infinitesimal perturbation. Setting a
// physical parameter to (p, h) and reading the outer imaginary part divided by h
// gives the derivative with no subtractive cancellation at all -- exact to
// machine precision for arbitrarily small h, which a centred difference can
// never be.
//
// real() returns the real part, which is the right choice for the trait's stated
// purpose: comparisons, tolerances and reporting. It is never used inside a
// physics expression.
// ---------------------------------------------------------------------------

namespace warp_hypad {

template <>
struct ScalarTraits<program_b::inverse_design::Complex<double>> {
    static double real(const program_b::inverse_design::Complex<double>& x) {
        return x.re;
    }
    static program_b::inverse_design::Complex<double> from_double(double x) {
        return program_b::inverse_design::Complex<double>(x, 0.0);
    }
    static constexpr bool is_differential = true;
    static constexpr int derivative_channels = 1;
};

}  // namespace warp_hypad

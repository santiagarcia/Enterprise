#pragma once

#include <cmath>

#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// A deliberately minimal custom scalar type whose only job is to fail to
// compile if a kernel cheats.
//
// TestScalar provides exactly the operations a HyPAD-style hypercomplex scalar
// can be expected to provide, and NOTHING else. In particular it deliberately
// omits:
//
//   * implicit conversion to or from double, so a kernel cannot silently drop
//     into ordinary arithmetic;
//   * operator<, >, <=, >=, ==, !=, so a kernel cannot branch on a value;
//   * fabs/abs, which is a branch in disguise and is non-differentiable at 0.
//
// If Milestone 1's algebra compiles and runs correctly with TestScalar, then it
// contains no hidden dependency on double and no hidden ordering decision --
// which is precisely the property Milestone 10 needs when HyPADScalar is
// dropped in. This is a compile-time proof obligation, not a runtime test, and
// it is why the type carries no derivative channels: adding them would test
// something else.
// ---------------------------------------------------------------------------

namespace warp_hypad {

class TestScalar {
public:
    TestScalar() : v_(0.0) {}
    explicit TestScalar(double v) : v_(v) {}

    // Explicit accessor only -- there is no operator double().
    double value() const { return v_; }

    friend TestScalar operator+(const TestScalar& a, const TestScalar& b) {
        return TestScalar(a.v_ + b.v_);
    }
    friend TestScalar operator-(const TestScalar& a, const TestScalar& b) {
        return TestScalar(a.v_ - b.v_);
    }
    friend TestScalar operator*(const TestScalar& a, const TestScalar& b) {
        return TestScalar(a.v_ * b.v_);
    }
    friend TestScalar operator/(const TestScalar& a, const TestScalar& b) {
        return TestScalar(a.v_ / b.v_);
    }
    friend TestScalar operator-(const TestScalar& a) { return TestScalar(-a.v_); }

    TestScalar& operator+=(const TestScalar& o) { v_ += o.v_; return *this; }
    TestScalar& operator-=(const TestScalar& o) { v_ -= o.v_; return *this; }
    TestScalar& operator*=(const TestScalar& o) { v_ *= o.v_; return *this; }
    TestScalar& operator/=(const TestScalar& o) { v_ /= o.v_; return *this; }

private:
    double v_;
};

// Smooth elementary functions, found by ADL from generic kernels. These are the
// ones the metrics in Milestones 3-4 actually need. Each has a well-defined
// derivative everywhere on its domain, which is the criterion for inclusion.
inline TestScalar sqrt(const TestScalar& x) { return TestScalar(std::sqrt(x.value())); }
inline TestScalar exp(const TestScalar& x) { return TestScalar(std::exp(x.value())); }
inline TestScalar log(const TestScalar& x) { return TestScalar(std::log(x.value())); }
inline TestScalar tanh(const TestScalar& x) { return TestScalar(std::tanh(x.value())); }
inline TestScalar sinh(const TestScalar& x) { return TestScalar(std::sinh(x.value())); }
inline TestScalar cosh(const TestScalar& x) { return TestScalar(std::cosh(x.value())); }
inline TestScalar sin(const TestScalar& x) { return TestScalar(std::sin(x.value())); }
inline TestScalar cos(const TestScalar& x) { return TestScalar(std::cos(x.value())); }
inline TestScalar pow(const TestScalar& x, double p) {
    return TestScalar(std::pow(x.value(), p));
}

template <>
struct ScalarTraits<TestScalar> {
    static double real(const TestScalar& x) { return x.value(); }
    static TestScalar from_double(double x) { return TestScalar(x); }
    static constexpr bool is_differential = false;
    static constexpr int derivative_channels = 0;
};

}  // namespace warp_hypad

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Univariate quintic B-spline basis with analytic derivatives to third order.
//
// DEGREE FIVE IS NOT A STYLE CHOICE. The metric this basis will eventually
// parameterise is differentiated twice to get curvature and three times for the
// contracted-Bianchi diagnostic. A degree-p spline with simple interior knots is
// C^(p-1); at p = 5 that is C4, so the second derivative is C2 and the third is
// C1. A cubic basis would give a merely continuous second derivative, and the
// curvature built from it would show knot-line artefacts that look like physics.
//
// THE MINIMUM BASIS COUNT. A clamped degree-5 knot vector needs p+1 = 6 repeats
// at each end, so a basis of size n has n + 6 knots and n - 6 simple interior
// knots. A family with n < 6 has NO interior knots at all -- it is a single
// Bernstein patch, with no local support and nothing for a spline to do. n = 6
// is the smallest admissible size and is checked, not assumed.
//
// DERIVATIVES ARE RECURRENCES, NOT DIFFERENCES. The standard triangular
// recurrence (the NURBS Book, algorithm A2.3) produces the value and all
// derivatives up to a requested order for the p+1 locally active functions in
// one pass, to full precision, with no step size to choose. Finite differences
// appear in this project only as an INDEPENDENT verification path, never as the
// production route -- a basis whose derivatives were themselves differenced
// would have its own truncation error entangled with the metric's.
//
// LOCAL SUPPORT. N_i is nonzero only on [knot_i, knot_{i+p+1}), so at any point
// at most p+1 = 6 functions are active. `evaluate` returns exactly those, and
// returns an EMPTY evaluation outside the parameter interval -- not small
// numbers, not clamped edge values. Section 6 of the O0A specification turns
// that property into a hard guarantee about the flat collar, and it only holds
// if "outside the support" means exactly zero.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

inline constexpr int kQuinticDegree = 5;
inline constexpr int kMinQuinticBasisCount = kQuinticDegree + 1;  // 6

// Only the locally active functions, with their global indices.
struct BasisEvaluation1D {
    std::vector<int> active_indices;
    std::vector<double> value;
    std::vector<double> first;
    std::vector<double> second;
    std::vector<double> third;

    std::size_t size() const { return active_indices.size(); }
    bool empty() const { return active_indices.empty(); }
};

class QuinticBSplineBasis1D {
public:
    // A clamped uniform basis of `n_basis` functions on [a, b].
    //
    // `valid()` is false when n_basis is below the degree-5 minimum; the
    // constructor does not throw, because the gate must be able to CHECK the
    // rejection rather than crash on it.
    QuinticBSplineBasis1D(double a, double b, int n_basis) {
        a_ = a;
        b_ = b;
        n_basis_ = n_basis;
        if (n_basis < kMinQuinticBasisCount || !(b > a)) {
            reason_ = (n_basis < kMinQuinticBasisCount)
                          ? "degree 5 needs at least 6 basis functions per "
                            "active direction"
                          : "the parameter interval must have positive length";
            return;
        }
        const int p = kQuinticDegree;
        const int n_interior = n_basis - p - 1;  // simple interior knots
        knots_.reserve(static_cast<std::size_t>(n_basis + p + 1));
        for (int k = 0; k <= p; ++k) knots_.push_back(a);
        for (int k = 1; k <= n_interior; ++k)
            knots_.push_back(a + (b - a) * static_cast<double>(k) /
                                     static_cast<double>(n_interior + 1));
        for (int k = 0; k <= p; ++k) knots_.push_back(b);
        valid_ = true;
    }

    // An explicit knot vector, for the repeated-interior-knot case. The size
    // relation n_basis = knots.size() - p - 1 is what defines the basis count.
    explicit QuinticBSplineBasis1D(std::vector<double> knots) {
        const int p = kQuinticDegree;
        if (static_cast<int>(knots.size()) < 2 * (p + 1)) {
            reason_ = "a clamped degree-5 knot vector needs at least 12 knots";
            return;
        }
        for (std::size_t i = 1; i < knots.size(); ++i) {
            if (knots[i] < knots[i - 1]) {
                reason_ = "knots must be non-decreasing";
                return;
            }
        }
        knots_ = std::move(knots);
        n_basis_ = static_cast<int>(knots_.size()) - p - 1;
        a_ = knots_[static_cast<std::size_t>(p)];
        b_ = knots_[static_cast<std::size_t>(n_basis_)];
        if (n_basis_ < kMinQuinticBasisCount) {
            reason_ = "degree 5 needs at least 6 basis functions";
            return;
        }
        if (!(b_ > a_)) {
            reason_ = "the parameter interval must have positive length";
            return;
        }
        valid_ = true;
    }

    int degree() const { return kQuinticDegree; }
    int size() const { return n_basis_; }
    const std::vector<double>& knots() const { return knots_; }
    double lower() const { return a_; }
    double upper() const { return b_; }
    bool valid() const { return valid_; }
    const std::string& rejection_reason() const { return reason_; }

    // The support of basis function i: [knots_i, knots_{i+p+1}].
    double support_lower(int i) const {
        return knots_[static_cast<std::size_t>(i)];
    }
    double support_upper(int i) const {
        return knots_[static_cast<std::size_t>(i + kQuinticDegree + 1)];
    }

    // The multiplicity of the first and last knot values. A CLAMPED end has
    // exactly degree+1 = 6; a uniform/open end has 1. "Exactly" matters in both
    // directions: 7 repeats would mean an interior knot had collided with the
    // endpoint, shrinking the space without saying so.
    int end_multiplicity_lower(double tol = 1e-12) const {
        int m = 0;
        for (double t : knots_)
            if (std::fabs(t - knots_.front()) <= tol) ++m;
        return m;
    }
    int end_multiplicity_upper(double tol = 1e-12) const {
        int m = 0;
        for (double t : knots_)
            if (std::fabs(t - knots_.back()) <= tol) ++m;
        return m;
    }

    // The declared semantics of the knot vector, read back from the vector
    // itself rather than remembered from the constructor -- so a vector built by
    // one path and described by another cannot disagree with its description.
    enum class Semantics { Clamped, LowerClampedOpenUpper, Uniform, Irregular };

    Semantics semantics(double tol = 1e-12) const {
        const int lo = end_multiplicity_lower(tol);
        const int hi = end_multiplicity_upper(tol);
        const int want = kQuinticDegree + 1;
        if (lo == want && hi == want) return Semantics::Clamped;
        if (lo == want && hi == 1) return Semantics::LowerClampedOpenUpper;
        if (lo == 1 && hi == 1) return Semantics::Uniform;
        return Semantics::Irregular;
    }

    static const char* to_string(Semantics s) {
        switch (s) {
            case Semantics::Clamped: return "clamped";
            case Semantics::LowerClampedOpenUpper:
                return "lower_clamped_open_upper";
            case Semantics::Uniform: return "uniform";
            case Semantics::Irregular: return "irregular";
        }
        return "unknown";
    }

    // Uniform interior spacing, or 0 when the spacing is not uniform. Used to
    // confirm a lattice is the lattice it claims to be.
    double uniform_spacing(double tol = 1e-12) const {
        double h = 0.0;
        for (std::size_t i = 1; i < knots_.size(); ++i) {
            const double d = knots_[i] - knots_[i - 1];
            if (d <= tol) continue;
            if (h == 0.0)
                h = d;
            else if (std::fabs(d - h) > 1e-9 * h)
                return 0.0;
        }
        return h;
    }

    // The number of distinct interior knots, and the multiplicity of each. A
    // multiplicity m interior knot leaves C^(p-m) continuity there.
    int interior_knot_multiplicity(double x, double tol = 1e-12) const {
        int m = 0;
        for (std::size_t i = static_cast<std::size_t>(kQuinticDegree + 1);
             i + static_cast<std::size_t>(kQuinticDegree + 1) < knots_.size();
             ++i)
            if (std::fabs(knots_[i] - x) <= tol) ++m;
        return m;
    }

    std::vector<double> distinct_interior_knots(double tol = 1e-12) const {
        std::vector<double> out;
        for (std::size_t i = static_cast<std::size_t>(kQuinticDegree + 1);
             i + static_cast<std::size_t>(kQuinticDegree + 1) < knots_.size();
             ++i) {
            if (out.empty() || std::fabs(knots_[i] - out.back()) > tol)
                out.push_back(knots_[i]);
        }
        return out;
    }

    // The knot span containing x, or -1 when x is outside [a, b].
    //
    // The upper end is INCLUDED and mapped to the last nonempty span, which is
    // the standard clamped convention: without it, x = b would report no active
    // functions even though N_{n-1}(b) = 1.
    int find_span(double x) const {
        if (!valid_) return -1;
        const int p = kQuinticDegree;
        const int n = n_basis_ - 1;
        if (x < a_ || x > b_) return -1;
        if (x >= knots_[static_cast<std::size_t>(n + 1)]) {
            // Walk back over any zero-length spans at the clamped end.
            int i = n;
            while (i > p && knots_[static_cast<std::size_t>(i)] ==
                                knots_[static_cast<std::size_t>(i + 1)])
                --i;
            return i;
        }
        int low = p, high = n + 1, mid = (low + high) / 2;
        while (x < knots_[static_cast<std::size_t>(mid)] ||
               x >= knots_[static_cast<std::size_t>(mid + 1)]) {
            if (x < knots_[static_cast<std::size_t>(mid)])
                high = mid;
            else
                low = mid;
            mid = (low + high) / 2;
        }
        return mid;
    }

    // Value and first three derivatives of the locally active functions.
    //
    // Returns an EMPTY evaluation outside [a, b]. That is the exactness the
    // flat-collar guarantee rests on: a function whose support does not contain
    // the point is not in the returned list at all, so no rounding can leak a
    // contribution.
    BasisEvaluation1D evaluate(double coordinate) const {
        BasisEvaluation1D out;
        const int span = find_span(coordinate);
        if (span < 0) return out;

        constexpr int p = kQuinticDegree;
        constexpr int kOrders = 3;  // first, second, third

        // ndu holds the value triangle; its lower triangle carries the knot
        // differences the derivative pass reuses.
        double ndu[p + 1][p + 1];
        double left[p + 1], right[p + 1];
        ndu[0][0] = 1.0;
        for (int j = 1; j <= p; ++j) {
            left[j] = coordinate -
                      knots_[static_cast<std::size_t>(span + 1 - j)];
            right[j] = knots_[static_cast<std::size_t>(span + j)] - coordinate;
            double saved = 0.0;
            for (int r = 0; r < j; ++r) {
                ndu[j][r] = right[r + 1] + left[j - r];
                // A zero denominator means a repeated knot has collapsed this
                // sub-interval. The B-spline convention is 0/0 = 0, and using it
                // here is what makes a multiple interior knot evaluate correctly
                // rather than produce a NaN.
                const double temp = (ndu[j][r] != 0.0)
                                        ? ndu[r][j - 1] / ndu[j][r]
                                        : 0.0;
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        double ders[kOrders + 1][p + 1];
        for (int j = 0; j <= p; ++j) ders[0][j] = ndu[j][p];

        double aa[2][p + 1];
        for (int r = 0; r <= p; ++r) {
            int s1 = 0, s2 = 1;
            aa[0][0] = 1.0;
            for (int k = 1; k <= kOrders; ++k) {
                double d = 0.0;
                const int rk = r - k, pk = p - k;
                if (r >= k) {
                    const double den = ndu[pk + 1][rk];
                    aa[s2][0] = (den != 0.0) ? aa[s1][0] / den : 0.0;
                    d = aa[s2][0] * ndu[rk][pk];
                }
                const int j1 = (rk >= -1) ? 1 : -rk;
                const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
                for (int j = j1; j <= j2; ++j) {
                    const double den = ndu[pk + 1][rk + j];
                    aa[s2][j] = (den != 0.0)
                                    ? (aa[s1][j] - aa[s1][j - 1]) / den
                                    : 0.0;
                    d += aa[s2][j] * ndu[rk + j][pk];
                }
                if (r <= pk) {
                    const double den = ndu[pk + 1][r];
                    aa[s2][k] = (den != 0.0) ? -aa[s1][k - 1] / den : 0.0;
                    d += aa[s2][k] * ndu[r][pk];
                }
                ders[k][r] = d;
                const int tmp = s1;
                s1 = s2;
                s2 = tmp;
            }
        }

        // The recurrence returns derivatives scaled by the falling factorial.
        int factor = p;
        for (int k = 1; k <= kOrders; ++k) {
            for (int j = 0; j <= p; ++j)
                ders[k][j] *= static_cast<double>(factor);
            factor *= (p - k);
        }

        out.active_indices.reserve(p + 1);
        out.value.reserve(p + 1);
        out.first.reserve(p + 1);
        out.second.reserve(p + 1);
        out.third.reserve(p + 1);
        for (int j = 0; j <= p; ++j) {
            out.active_indices.push_back(span - p + j);
            out.value.push_back(ders[0][j]);
            out.first.push_back(ders[1][j]);
            out.second.push_back(ders[2][j]);
            out.third.push_back(ders[3][j]);
        }
        return out;
    }

    // The value of one basis function, for diagnostics and for the exact-zero
    // check. Returns 0.0 when the point is outside that function's support --
    // exactly zero, by construction, not by cancellation.
    double value_of(int index, double coordinate) const {
        const BasisEvaluation1D e = evaluate(coordinate);
        for (std::size_t k = 0; k < e.size(); ++k)
            if (e.active_indices[k] == index) return e.value[k];
        return 0.0;
    }

    double derivative_of(int index, double coordinate, int order) const {
        const BasisEvaluation1D e = evaluate(coordinate);
        for (std::size_t k = 0; k < e.size(); ++k) {
            if (e.active_indices[k] != index) continue;
            switch (order) {
                case 0: return e.value[k];
                case 1: return e.first[k];
                case 2: return e.second[k];
                case 3: return e.third[k];
                default: return 0.0;
            }
        }
        return 0.0;
    }

private:
    double a_ = 0.0, b_ = 1.0;
    int n_basis_ = 0;
    std::vector<double> knots_;
    bool valid_ = false;
    std::string reason_;
};

}  // namespace design
}  // namespace warp_hypad

#include "warp_hypad/core/grid.hpp"

#include <cmath>
#include <sstream>

namespace warp_hypad {
namespace {

// x_i = x_min*(1-t) + x_max*t with t = i/(n-1). Exact at both endpoints and,
// for odd n, at the midpoint. See the header for why this form is preferred.
std::vector<double> make_axis(int n, double lo, double hi) {
    std::vector<double> out(static_cast<std::size_t>(n));
    if (n == 1) {
        out[0] = lo;
        return out;
    }
    const double inv = 1.0 / static_cast<double>(n - 1);
    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) * inv;
        out[static_cast<std::size_t>(i)] = lo * (1.0 - t) + hi * t;
    }
    // Pin the endpoints, and the midpoint for odd n. The interpolation form
    // already lands on these for many bounds, but i*inv is not exactly 0.5 or
    // 1.0 for every n (1/(n-1) is generally inexact in binary), so pinning
    // turns "usually exact" into an unconditional guarantee that the Milestone 1
    // gate can assert with ==.
    out[0] = lo;
    out[static_cast<std::size_t>(n - 1)] = hi;
    if (n % 2 == 1) out[static_cast<std::size_t>((n - 1) / 2)] = 0.5 * (lo + hi);
    return out;
}

bool check_axis(int n, double lo, double hi, const char* name,
                std::string& error) {
    std::ostringstream os;
    if (n < 1) {
        os << "grid." << name << " must be >= 1, got " << n;
        error = os.str();
        return false;
    }
    if (!(std::isfinite(lo) && std::isfinite(hi))) {
        os << "grid bounds for " << name << " must be finite";
        error = os.str();
        return false;
    }
    if (n > 1 && !(hi > lo)) {
        os << "grid bounds for " << name << " require max > min, got min=" << lo
           << " max=" << hi;
        error = os.str();
        return false;
    }
    return true;
}

}  // namespace

bool validate_grid_spec(const GridSpec& spec, std::string& error) {
    if (!check_axis(spec.nx, spec.x_min, spec.x_max, "nx/x", error)) return false;
    if (!check_axis(spec.ny, spec.y_min, spec.y_max, "ny/y", error)) return false;
    if (!check_axis(spec.nz, spec.z_min, spec.z_max, "nz/z", error)) return false;

    // Guard against an allocation that would fail later in a less obvious place.
    const double points = static_cast<double>(spec.nx) *
                          static_cast<double>(spec.ny) *
                          static_cast<double>(spec.nz);
    if (points > 2.0e9) {
        error = "grid has more than 2e9 points; refusing to allocate";
        return false;
    }
    return true;
}

Grid::Grid(const GridSpec& spec) : spec_(spec) {
    xs_ = make_axis(spec_.nx, spec_.x_min, spec_.x_max);
    ys_ = make_axis(spec_.ny, spec_.y_min, spec_.y_max);
    zs_ = make_axis(spec_.nz, spec_.z_min, spec_.z_max);

    dx_ = (spec_.nx > 1)
              ? (spec_.x_max - spec_.x_min) / static_cast<double>(spec_.nx - 1)
              : 0.0;
    dy_ = (spec_.ny > 1)
              ? (spec_.y_max - spec_.y_min) / static_cast<double>(spec_.ny - 1)
              : 0.0;
    dz_ = (spec_.nz > 1)
              ? (spec_.z_max - spec_.z_min) / static_cast<double>(spec_.nz - 1)
              : 0.0;

    // A degenerate axis contributes a factor of 1, not 0, so that a single-slice
    // grid still integrates sensibly over the non-degenerate directions.
    cell_volume_ = (dx_ > 0.0 ? dx_ : 1.0) * (dy_ > 0.0 ? dy_ : 1.0) *
                   (dz_ > 0.0 ? dz_ : 1.0);
}

}  // namespace warp_hypad

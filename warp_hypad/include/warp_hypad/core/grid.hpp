#pragma once

#include <cstddef>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Regular Cartesian grid over a rectangular box in the spatial coordinates
// (x, y, z). Time is not gridded: the metrics implemented through Milestone 4
// are evaluated on a single time slice, and any temporal derivative that is
// assumed to vanish is stated explicitly at the point of use.
//
// Coordinates are stored explicitly rather than recomputed, so that every
// consumer sees bit-identical values.
//
// Coordinate generation uses the two-endpoint interpolation form
//
//     x_i = x_min * (1 - t) + x_max * t,     t = i / (n - 1)
//
// rather than the accumulate form x_min + i*dx, which drifts by up to n ulp and
// generally misses x_max outright. The two endpoints and, for odd n, the
// midpoint are then pinned to their exact values.
//
// The consequence is a deliberate trade: the three pinned coordinates are
// exactly right, and the *spacing* is therefore uniform only to within a few
// ulp rather than bit-exactly. The Milestone 1 gate reflects this honestly -- it
// asserts endpoints and midpoint with ==, and uniformity at 1e-12 relative.
// ---------------------------------------------------------------------------

namespace warp_hypad {

struct GridSpec {
    int nx = 41;
    int ny = 41;
    int nz = 41;
    double x_min = -10.0;
    double x_max = 10.0;
    double y_min = -10.0;
    double y_max = 10.0;
    double z_min = -10.0;
    double z_max = 10.0;
};

// Returns false and fills `error` if the spec cannot produce a usable grid.
bool validate_grid_spec(const GridSpec& spec, std::string& error);

class Grid {
public:
    explicit Grid(const GridSpec& spec);

    const GridSpec& spec() const { return spec_; }

    int nx() const { return spec_.nx; }
    int ny() const { return spec_.ny; }
    int nz() const { return spec_.nz; }

    std::size_t num_points() const {
        return static_cast<std::size_t>(spec_.nx) *
               static_cast<std::size_t>(spec_.ny) *
               static_cast<std::size_t>(spec_.nz);
    }

    // Uniform spacing along each axis. Zero when the corresponding dimension
    // has a single point (a degenerate but legal configuration, useful for
    // 2D-like slices).
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double dz() const { return dz_; }

    // Volume element for the midpoint quadrature rule used by the objective.
    double cell_volume() const { return cell_volume_; }

    double x(int i) const { return xs_[static_cast<std::size_t>(i)]; }
    double y(int j) const { return ys_[static_cast<std::size_t>(j)]; }
    double z(int k) const { return zs_[static_cast<std::size_t>(k)]; }

    const std::vector<double>& xs() const { return xs_; }
    const std::vector<double>& ys() const { return ys_; }
    const std::vector<double>& zs() const { return zs_; }

    // Flat index. z varies fastest, then y, then x, so a walk along +z is
    // contiguous. See docs/memory_layout.md for why this pairs with the
    // point-major component layout in field.hpp.
    std::size_t flat(int i, int j, int k) const {
        return (static_cast<std::size_t>(i) * static_cast<std::size_t>(spec_.ny) +
                static_cast<std::size_t>(j)) *
                   static_cast<std::size_t>(spec_.nz) +
               static_cast<std::size_t>(k);
    }

    void unflatten(std::size_t f, int& i, int& j, int& k) const {
        const std::size_t nz = static_cast<std::size_t>(spec_.nz);
        const std::size_t ny = static_cast<std::size_t>(spec_.ny);
        k = static_cast<int>(f % nz);
        const std::size_t rest = f / nz;
        j = static_cast<int>(rest % ny);
        i = static_cast<int>(rest / ny);
    }

    // Stride in flat-index units for a step of +1 along each axis. Finite
    // difference stencils index by these.
    std::size_t stride_x() const {
        return static_cast<std::size_t>(spec_.ny) *
               static_cast<std::size_t>(spec_.nz);
    }
    std::size_t stride_y() const { return static_cast<std::size_t>(spec_.nz); }
    std::size_t stride_z() const { return 1; }

    bool is_boundary(int i, int j, int k) const {
        return i == 0 || i == spec_.nx - 1 || j == 0 || j == spec_.ny - 1 ||
               k == 0 || k == spec_.nz - 1;
    }

    // Index of the central point along each axis, defined only for odd extents.
    // Returns -1 when the extent is even, so callers must handle it rather than
    // silently picking one side of the centre.
    int center_index_x() const { return (spec_.nx % 2 == 1) ? spec_.nx / 2 : -1; }
    int center_index_y() const { return (spec_.ny % 2 == 1) ? spec_.ny / 2 : -1; }
    int center_index_z() const { return (spec_.nz % 2 == 1) ? spec_.nz / 2 : -1; }

private:
    GridSpec spec_;
    std::vector<double> xs_, ys_, zs_;
    double dx_ = 0.0, dy_ = 0.0, dz_ = 0.0;
    double cell_volume_ = 0.0;
};

}  // namespace warp_hypad

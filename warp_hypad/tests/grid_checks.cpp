#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include "checks.hpp"
#include "warp_hypad/core/grid.hpp"

namespace warp_hypad {
namespace checks {
namespace {

std::string grid_label(const GridSpec& s) {
    std::ostringstream os;
    os << s.nx << "x" << s.ny << "x" << s.nz;
    return os.str();
}

// Two grids. The cube would hide an x/y/z mix-up because every axis is
// identical; the asymmetric one has distinct extents AND distinct bounds on
// every axis, so any transposition in flat(), unflatten() or the stride
// accessors shows up immediately.
const GridSpec& cube_spec() {
    static const GridSpec s = [] {
        GridSpec g;
        g.nx = 41; g.ny = 41; g.nz = 41;
        g.x_min = -10.0; g.x_max = 10.0;
        g.y_min = -10.0; g.y_max = 10.0;
        g.z_min = -10.0; g.z_max = 10.0;
        return g;
    }();
    return s;
}

const GridSpec& asymmetric_spec() {
    static const GridSpec s = [] {
        GridSpec g;
        g.nx = 17; g.ny = 9; g.nz = 5;
        g.x_min = -3.0; g.x_max = 7.0;
        g.y_min = 0.0;  g.y_max = 4.0;
        g.z_min = -1.0; g.z_max = 1.0;
        return g;
    }();
    return s;
}

// Largest relative deviation of consecutive spacings from the nominal spacing.
double spacing_deviation(const std::vector<double>& axis, double nominal) {
    if (axis.size() < 2 || nominal == 0.0) return 0.0;
    double worst = 0.0;
    for (std::size_t i = 0; i + 1 < axis.size(); ++i) {
        const double gap = axis[i + 1] - axis[i];
        const double rel = std::fabs(gap - nominal) / std::fabs(nominal);
        if (rel > worst) worst = rel;
    }
    return worst;
}

}  // namespace

void run_grid_checks(VerificationReport& report) {
    const GridSpec specs[2] = {cube_spec(), asymmetric_spec()};

    for (const GridSpec& spec : specs) {
        const Grid grid(spec);
        const std::string label = grid_label(spec);

        // --- Endpoints, exactly ---
        {
            Stopwatch sw;
            int mismatches = 0;
            if (grid.x(0) != spec.x_min) ++mismatches;
            if (grid.y(0) != spec.y_min) ++mismatches;
            if (grid.z(0) != spec.z_min) ++mismatches;
            report.check_exact(
                "grid.lower_bound_exact [" + label + "]",
                "x(0) == x_min for every axis",
                VerificationLevel::A, label, 0.0, mismatches, sw.seconds(),
                "coordinate generation pins the lower endpoint");
        }
        {
            Stopwatch sw;
            int mismatches = 0;
            if (grid.x(spec.nx - 1) != spec.x_max) ++mismatches;
            if (grid.y(spec.ny - 1) != spec.y_max) ++mismatches;
            if (grid.z(spec.nz - 1) != spec.z_max) ++mismatches;
            report.check_exact(
                "grid.upper_bound_exact [" + label + "]",
                "x(n-1) == x_max for every axis",
                VerificationLevel::A, label, 0.0, mismatches, sw.seconds(),
                "the accumulate form x_min + i*dx would fail this");
        }

        // --- Uniform spacing ---
        {
            Stopwatch sw;
            double worst = spacing_deviation(grid.xs(), grid.dx());
            worst = std::fmax(worst, spacing_deviation(grid.ys(), grid.dy()));
            worst = std::fmax(worst, spacing_deviation(grid.zs(), grid.dz()));
            report.check_residual(
                "grid.spacing_uniform [" + label + "]",
                "max_i |(x_{i+1} - x_i) - dx| / dx == 0",
                VerificationLevel::A, label, "n/a", worst, 1e-12, sw.seconds(),
                "tolerance, not exact: pinning three coordinates costs a few ulp "
                "of spacing uniformity by design");
        }

        // --- Centre coordinate, exactly, for odd extents ---
        {
            Stopwatch sw;
            int mismatches = 0;
            int checked = 0;
            if (grid.center_index_x() >= 0) {
                ++checked;
                if (grid.x(grid.center_index_x()) != 0.5 * (spec.x_min + spec.x_max)) ++mismatches;
            }
            if (grid.center_index_y() >= 0) {
                ++checked;
                if (grid.y(grid.center_index_y()) != 0.5 * (spec.y_min + spec.y_max)) ++mismatches;
            }
            if (grid.center_index_z() >= 0) {
                ++checked;
                if (grid.z(grid.center_index_z()) != 0.5 * (spec.z_min + spec.z_max)) ++mismatches;
            }
            std::ostringstream note;
            note << checked << " of 3 axes have odd extent and were checked";
            report.check_exact(
                "grid.center_coordinate_exact [" + label + "]",
                "x((n-1)/2) == (x_min + x_max)/2 for odd n",
                VerificationLevel::A, label, 0.0, mismatches, sw.seconds(),
                note.str());
        }

        // --- Index round-trip, every point, both directions ---
        {
            Stopwatch sw;
            long long mismatches = 0;
            for (int i = 0; i < spec.nx; ++i) {
                for (int j = 0; j < spec.ny; ++j) {
                    for (int k = 0; k < spec.nz; ++k) {
                        const std::size_t f = grid.flat(i, j, k);
                        int bi = -1, bj = -1, bk = -1;
                        grid.unflatten(f, bi, bj, bk);
                        if (bi != i || bj != j || bk != k) ++mismatches;
                    }
                }
            }
            report.check_exact(
                "grid.index_roundtrip_ijk [" + label + "]",
                "(i,j,k) -> flat -> (i,j,k) is the identity at every point",
                VerificationLevel::A, label, 0.0,
                static_cast<double>(mismatches), sw.seconds(),
                "exhaustive over all grid points");
        }
        {
            Stopwatch sw;
            long long mismatches = 0;
            const std::size_t n = grid.num_points();
            for (std::size_t f = 0; f < n; ++f) {
                int i = -1, j = -1, k = -1;
                grid.unflatten(f, i, j, k);
                if (grid.flat(i, j, k) != f) ++mismatches;
            }
            report.check_exact(
                "grid.index_roundtrip_flat [" + label + "]",
                "flat -> (i,j,k) -> flat is the identity for every flat index",
                VerificationLevel::A, label, 0.0,
                static_cast<double>(mismatches), sw.seconds(),
                "exhaustive over all flat indices");
        }

        // --- The flat map is a bijection onto [0, N) ---
        {
            Stopwatch sw;
            const std::size_t n = grid.num_points();
            std::vector<unsigned char> seen(n, 0);
            long long out_of_range = 0;
            long long collisions = 0;
            for (int i = 0; i < spec.nx; ++i) {
                for (int j = 0; j < spec.ny; ++j) {
                    for (int k = 0; k < spec.nz; ++k) {
                        const std::size_t f = grid.flat(i, j, k);
                        if (f >= n) { ++out_of_range; continue; }
                        if (seen[f]) ++collisions;
                        seen[f] = 1;
                    }
                }
            }
            long long unvisited = 0;
            for (std::size_t f = 0; f < n; ++f) if (!seen[f]) ++unvisited;

            report.check_exact(
                "grid.flat_index_bijective [" + label + "]",
                "flat() is a bijection from the index box onto [0, N)",
                VerificationLevel::A, label, 0.0,
                static_cast<double>(out_of_range + collisions + unvisited),
                sw.seconds(),
                "counts out-of-range + collisions + unvisited slots");
        }

        // --- Strides agree with flat() ---
        {
            Stopwatch sw;
            long long mismatches = 0;
            for (int i = 0; i + 1 < spec.nx; ++i)
                if (grid.flat(i + 1, 0, 0) - grid.flat(i, 0, 0) != grid.stride_x()) ++mismatches;
            for (int j = 0; j + 1 < spec.ny; ++j)
                if (grid.flat(0, j + 1, 0) - grid.flat(0, j, 0) != grid.stride_y()) ++mismatches;
            for (int k = 0; k + 1 < spec.nz; ++k)
                if (grid.flat(0, 0, k + 1) - grid.flat(0, 0, k) != grid.stride_z()) ++mismatches;
            report.check_exact(
                "grid.strides_match_flat [" + label + "]",
                "stride_a() == flat(...+1...) - flat(...) for each axis",
                VerificationLevel::A, label, 0.0,
                static_cast<double>(mismatches), sw.seconds(),
                "finite-difference stencils in Milestone 2 index by these");
        }

        // --- Derived quantities ---
        {
            Stopwatch sw;
            const double expected_dx =
                (spec.x_max - spec.x_min) / static_cast<double>(spec.nx - 1);
            report.check_value(
                "grid.dx_value [" + label + "]",
                "dx == (x_max - x_min) / (nx - 1)",
                VerificationLevel::A, label, "n/a", expected_dx, grid.dx(),
                1e-15, sw.seconds());
        }
        {
            Stopwatch sw;
            const double expected_vol = grid.dx() * grid.dy() * grid.dz();
            report.check_value(
                "grid.cell_volume [" + label + "]",
                "cell_volume == dx * dy * dz",
                VerificationLevel::A, label, "n/a", expected_vol,
                grid.cell_volume(), 1e-15, sw.seconds());
        }
        {
            Stopwatch sw;
            const double expected_n = static_cast<double>(spec.nx) *
                                      static_cast<double>(spec.ny) *
                                      static_cast<double>(spec.nz);
            report.check_exact(
                "grid.num_points [" + label + "]",
                "num_points == nx * ny * nz",
                VerificationLevel::A, label, expected_n,
                static_cast<double>(grid.num_points()), sw.seconds());
        }
    }

    // --- Spec validation rejects malformed input ---
    {
        Stopwatch sw;
        std::string err;
        bool all_rejected = true;

        GridSpec bad = cube_spec();
        bad.nx = 0;
        if (validate_grid_spec(bad, err)) all_rejected = false;

        bad = cube_spec();
        bad.x_max = bad.x_min - 1.0;
        if (validate_grid_spec(bad, err)) all_rejected = false;

        bad = cube_spec();
        bad.y_max = std::numeric_limits<double>::quiet_NaN();
        if (validate_grid_spec(bad, err)) all_rejected = false;

        bad = cube_spec();
        bad.nx = 2000; bad.ny = 2000; bad.nz = 2000;  // 8e9 points
        if (validate_grid_spec(bad, err)) all_rejected = false;

        report.check_predicate(
            "grid.validation_rejects_malformed",
            "validate_grid_spec rejects n<1, inverted bounds, non-finite "
            "bounds and unallocatable sizes",
            VerificationLevel::A, all_rejected,
            "all four malformed specs rejected",
            all_rejected ? "all four rejected" : "at least one was accepted",
            sw.seconds());
    }

    // A well-formed spec must still be accepted -- otherwise the check above
    // would pass trivially with a validator that rejects everything.
    {
        Stopwatch sw;
        std::string err;
        const bool ok = validate_grid_spec(cube_spec(), err) &&
                        validate_grid_spec(asymmetric_spec(), err);
        report.check_predicate(
            "grid.validation_accepts_valid",
            "validate_grid_spec accepts both reference specs",
            VerificationLevel::A, ok, "both accepted",
            ok ? "both accepted" : ("rejected: " + err), sw.seconds(),
            "guards against a validator that rejects everything");
    }
}

}  // namespace checks
}  // namespace warp_hypad

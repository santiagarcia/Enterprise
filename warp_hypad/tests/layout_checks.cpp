#include <sstream>

#include "checks.hpp"
#include "warp_hypad/core/field.hpp"

namespace warp_hypad {
namespace checks {
namespace {

// Every component carries a fingerprint of its own full index. If any two
// slots alias, or any index is transposed anywhere in the addressing chain,
// the value read back is a *different* index tuple rather than merely a wrong
// number -- so the failure identifies the bug instead of just reporting it.
//
//     tag(i,j,k,mu,nu) = 1e6*i + 1e4*j + 1e2*k + 10*mu + nu
//
// With extents below 100 and mu,nu <= 3 every tag is a distinct integer below
// 1e8, exactly representable in double, so equality here is bit-exact and not
// a tolerance question.
double tag2(int i, int j, int k, int mu, int nu) {
    return 1e6 * i + 1e4 * j + 1e2 * k + 10.0 * mu + nu;
}

double tag3(int i, int j, int k, int a, int b, int c) {
    return 1e6 * i + 1e4 * j + 1e2 * k + 16.0 * a + 4.0 * b + c;
}

// The rank-4 component index reaches 64*3+16*3+4*3+3 = 255, which overflows a
// 1e2 stride for k and would make tags collide across different k. Wider
// strides keep every tag distinct; the largest value is about 1.01e10, still
// far below 2^53, so equality stays bit-exact.
double tag4(int i, int j, int k, int a, int b, int c, int d) {
    return 1e9 * i + 1e7 * j + 1e5 * k + 64.0 * a + 16.0 * b + 4.0 * c + d;
}

const GridSpec& layout_spec() {
    // Small and asymmetric: exhaustive rank-4 coverage is 256 components per
    // point, so 11*7*5 = 385 points is 98,560 components -- enough to be
    // exhaustive, small enough to stay instant.
    static const GridSpec s = [] {
        GridSpec g;
        g.nx = 11; g.ny = 7; g.nz = 5;
        g.x_min = -1.0; g.x_max = 1.0;
        g.y_min = -2.0; g.y_max = 2.0;
        g.z_min = -3.0; g.z_max = 3.0;
        return g;
    }();
    return s;
}

std::string label() {
    std::ostringstream os;
    const GridSpec& s = layout_spec();
    os << s.nx << "x" << s.ny << "x" << s.nz;
    return os.str();
}

}  // namespace

void run_layout_checks(VerificationReport& report) {
    const Grid grid(layout_spec());
    const GridSpec& s = layout_spec();
    const std::string lbl = label();

    // --- Rank 2: write every component, read every component ---
    {
        Stopwatch sw;
        Rank2Field<double> field(grid);
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu)
                            field.at(p, mu, nu) = tag2(i, j, k, mu, nu);
                }

        long long mismatches = 0;
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu)
                            if (field.at(p, mu, nu) != tag2(i, j, k, mu, nu)) ++mismatches;
                }

        report.check_exact(
            "layout.rank2_roundtrip",
            "every rank-2 component stores and retrieves its own index tag "
            "exactly",
            VerificationLevel::A, lbl, 0.0, static_cast<double>(mismatches),
            sw.seconds(), "exhaustive: 385 points x 16 components");
    }

    // --- Rank 2: the raw buffer really is point-major, 16*p + 4*mu + nu ---
    //
    // The check above would pass for ANY self-consistent addressing scheme.
    // This one pins the documented layout, which is what an external reader of
    // an exported field dump actually depends on.
    {
        Stopwatch sw;
        Rank2Field<double> field(grid);
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu)
                            field.at(p, mu, nu) = tag2(i, j, k, mu, nu);
                }

        long long mismatches = 0;
        const std::vector<double>& raw = field.raw();
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu) {
                            const std::size_t off =
                                16 * p + static_cast<std::size_t>(4 * mu + nu);
                            if (raw[off] != tag2(i, j, k, mu, nu)) ++mismatches;
                        }
                }

        report.check_exact(
            "layout.rank2_documented_offset",
            "raw()[16*p + 4*mu + nu] holds component (mu,nu) at point p",
            VerificationLevel::A, lbl, 0.0, static_cast<double>(mismatches),
            sw.seconds(),
            "pins the on-disk layout documented in metadata.json");
    }

    // --- Rank 3 ---
    {
        Stopwatch sw;
        Rank3Field<double> field(grid);
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b)
                            for (int c = 0; c < 4; ++c)
                                field.at(p, a, b, c) = tag3(i, j, k, a, b, c);
                }

        long long mismatches = 0;
        const std::vector<double>& raw = field.raw();
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b)
                            for (int c = 0; c < 4; ++c) {
                                const double want = tag3(i, j, k, a, b, c);
                                if (field.at(p, a, b, c) != want) ++mismatches;
                                const std::size_t off =
                                    64 * p + static_cast<std::size_t>(16 * a + 4 * b + c);
                                if (raw[off] != want) ++mismatches;
                            }
                }

        report.check_exact(
            "layout.rank3_roundtrip_and_offset",
            "rank-3 accessor and raw()[64*p + 16*a + 4*b + c] agree with the "
            "index tag everywhere",
            VerificationLevel::A, lbl, 0.0, static_cast<double>(mismatches),
            sw.seconds(), "exhaustive: 385 points x 64 components");
    }

    // --- Rank 4 ---
    {
        Stopwatch sw;
        Rank4Field<double> field(grid);
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b)
                            for (int c = 0; c < 4; ++c)
                                for (int d = 0; d < 4; ++d)
                                    field.at(p, a, b, c, d) = tag4(i, j, k, a, b, c, d);
                }

        long long mismatches = 0;
        const std::vector<double>& raw = field.raw();
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b)
                            for (int c = 0; c < 4; ++c)
                                for (int d = 0; d < 4; ++d) {
                                    const double want = tag4(i, j, k, a, b, c, d);
                                    if (field.at(p, a, b, c, d) != want) ++mismatches;
                                    const std::size_t off =
                                        256 * p + static_cast<std::size_t>(
                                                      64 * a + 16 * b + 4 * c + d);
                                    if (raw[off] != want) ++mismatches;
                                }
                }

        report.check_exact(
            "layout.rank4_roundtrip_and_offset",
            "rank-4 accessor and raw()[256*p + 64*a + 16*b + 4*c + d] agree "
            "with the index tag everywhere",
            VerificationLevel::A, lbl, 0.0, static_cast<double>(mismatches),
            sw.seconds(), "exhaustive: 385 points x 256 components");
    }

    // --- Scalar field ---
    {
        Stopwatch sw;
        ScalarField<double> field(grid);
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k)
                    field.at(i, j, k) = 1e6 * i + 1e4 * j + 1e2 * k;

        long long mismatches = 0;
        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const double want = 1e6 * i + 1e4 * j + 1e2 * k;
                    if (field.at(i, j, k) != want) ++mismatches;
                    if (field.at(grid.flat(i, j, k)) != want) ++mismatches;
                }

        report.check_exact(
            "layout.scalar_field_roundtrip",
            "ScalarField (i,j,k) and flat accessors agree and store exactly",
            VerificationLevel::A, lbl, 0.0, static_cast<double>(mismatches),
            sw.seconds());
    }

    // --- Mat4 <-> Rank2Field conversion preserves index placement ---
    //
    // set_matrix/matrix_at copy 16 contiguous values. If Mat4's internal order
    // ever diverged from the field's, this would catch it -- and that would be
    // a silent transpose of every metric in the pipeline.
    {
        Stopwatch sw;
        Rank2Field<double> field(grid);
        long long mismatches = 0;

        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    Mat4<double> m;
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu)
                            m(mu, nu) = tag2(i, j, k, mu, nu);
                    field.set_matrix(p, m);
                }

        for (int i = 0; i < s.nx; ++i)
            for (int j = 0; j < s.ny; ++j)
                for (int k = 0; k < s.nz; ++k) {
                    const std::size_t p = grid.flat(i, j, k);
                    const Mat4<double> m = field.matrix_at(p);
                    for (int mu = 0; mu < 4; ++mu)
                        for (int nu = 0; nu < 4; ++nu) {
                            const double want = tag2(i, j, k, mu, nu);
                            if (m(mu, nu) != want) ++mismatches;
                            if (field.at(p, mu, nu) != want) ++mismatches;
                        }
                }

        report.check_exact(
            "layout.mat4_field_conversion",
            "set_matrix followed by matrix_at is the identity, with no "
            "transposition",
            VerificationLevel::A, lbl, 0.0, static_cast<double>(mismatches),
            sw.seconds(),
            "a mismatch here would silently transpose every metric");
    }

    // --- Rank3 tensor <-> field conversion ---
    {
        Stopwatch sw;
        Rank3Field<double> field(grid);
        long long mismatches = 0;

        for (std::size_t p = 0; p < grid.num_points(); ++p) {
            Rank3<double> t;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int c = 0; c < 4; ++c)
                        t(a, b, c) = static_cast<double>(p) * 100.0 + 16 * a + 4 * b + c;
            field.set_tensor(p, t);
        }
        for (std::size_t p = 0; p < grid.num_points(); ++p) {
            const Rank3<double> t = field.tensor_at(p);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int c = 0; c < 4; ++c) {
                        const double want =
                            static_cast<double>(p) * 100.0 + 16 * a + 4 * b + c;
                        if (t(a, b, c) != want) ++mismatches;
                    }
        }

        report.check_exact(
            "layout.rank3_tensor_conversion",
            "set_tensor followed by tensor_at is the identity for rank 3",
            VerificationLevel::A, lbl, 0.0, static_cast<double>(mismatches),
            sw.seconds());
    }
}

}  // namespace checks
}  // namespace warp_hypad

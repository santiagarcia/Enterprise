#pragma once

#include <cassert>
#include <cstddef>
#include <vector>

#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"
#include "warp_hypad/core/tensor.hpp"

// ---------------------------------------------------------------------------
// Contiguous storage for scalar and tensor fields over a Grid.
//
// Layout is POINT-MAJOR: all components at one grid point are adjacent, and the
// component index varies fastest.
//
//     rank-2 field:  data[16*p + 4*mu + nu]
//     rank-3 field:  data[64*p + 16*a + 4*b + c]
//     rank-4 field:  data[256*p + 64*a + 16*b + 4*c + d]
//
// where p = grid.flat(i, j, k).
//
// Why point-major rather than component-major (data[comp*num_points + p]):
//
//   The dominant cost in this pipeline is pointwise 4x4 algebra -- inverting
//   g_{mu nu}, contracting Christoffels, assembling Ricci. Those touch all 16
//   (or 64, or 256) components of a single point together and do O(4^3)-O(4^4)
//   work per point. Point-major makes each of those a dense sequential read.
//
//   Finite differencing wants the opposite: the same component at neighbouring
//   points. Under point-major that is a strided read, stride 16*grid_stride.
//   But differencing is only O(1) work per component per point, so it is the
//   cheaper side of the trade, and the access pattern is still regular and
//   prefetchable.
//
//   This choice becomes more lopsided in favour of point-major once HyPADScalar
//   replaces double, because the pointwise algebra grows by the number of
//   derivative channels while the stencil arithmetic does not change shape.
//
// This decision is revisited in Milestone 9 with the runtime-versus-grid-size
// measurements, not before.
// ---------------------------------------------------------------------------

namespace warp_hypad {

template <typename Scalar>
class ScalarField {
public:
    ScalarField() = default;

    explicit ScalarField(const Grid& grid)
        : grid_(&grid), data_(grid.num_points(), scalar_cast<Scalar>(0.0)) {}

    const Grid& grid() const { return *grid_; }
    std::size_t num_points() const { return data_.size(); }

    Scalar& at(std::size_t p) {
        assert(p < data_.size());
        return data_[p];
    }
    const Scalar& at(std::size_t p) const {
        assert(p < data_.size());
        return data_[p];
    }

    Scalar& at(int i, int j, int k) { return at(grid_->flat(i, j, k)); }
    const Scalar& at(int i, int j, int k) const {
        return at(grid_->flat(i, j, k));
    }

    std::vector<Scalar>& raw() { return data_; }
    const std::vector<Scalar>& raw() const { return data_; }

private:
    const Grid* grid_ = nullptr;
    std::vector<Scalar> data_;
};

template <typename Scalar>
class Rank2Field {
public:
    static constexpr int kComponents = 16;

    Rank2Field() = default;

    explicit Rank2Field(const Grid& grid)
        : grid_(&grid),
          data_(grid.num_points() * kComponents, scalar_cast<Scalar>(0.0)) {}

    const Grid& grid() const { return *grid_; }
    std::size_t num_points() const { return data_.size() / kComponents; }

    Scalar& at(std::size_t p, int mu, int nu) {
        const std::size_t idx = p * kComponents + static_cast<std::size_t>(4 * mu + nu);
        assert(idx < data_.size());
        return data_[idx];
    }
    const Scalar& at(std::size_t p, int mu, int nu) const {
        const std::size_t idx = p * kComponents + static_cast<std::size_t>(4 * mu + nu);
        assert(idx < data_.size());
        return data_[idx];
    }

    Mat4<Scalar> matrix_at(std::size_t p) const {
        Mat4<Scalar> m;
        const std::size_t base = p * kComponents;
        for (int c = 0; c < kComponents; ++c) m.m[c] = data_[base + static_cast<std::size_t>(c)];
        return m;
    }

    void set_matrix(std::size_t p, const Mat4<Scalar>& m) {
        const std::size_t base = p * kComponents;
        for (int c = 0; c < kComponents; ++c) data_[base + static_cast<std::size_t>(c)] = m.m[c];
    }

    std::vector<Scalar>& raw() { return data_; }
    const std::vector<Scalar>& raw() const { return data_; }

private:
    const Grid* grid_ = nullptr;
    std::vector<Scalar> data_;
};

template <typename Scalar>
class Rank3Field {
public:
    static constexpr int kComponents = 64;

    Rank3Field() = default;

    explicit Rank3Field(const Grid& grid)
        : grid_(&grid),
          data_(grid.num_points() * kComponents, scalar_cast<Scalar>(0.0)) {}

    const Grid& grid() const { return *grid_; }
    std::size_t num_points() const { return data_.size() / kComponents; }

    Scalar& at(std::size_t p, int a, int b, int c) {
        const std::size_t idx =
            p * kComponents + static_cast<std::size_t>(16 * a + 4 * b + c);
        assert(idx < data_.size());
        return data_[idx];
    }
    const Scalar& at(std::size_t p, int a, int b, int c) const {
        const std::size_t idx =
            p * kComponents + static_cast<std::size_t>(16 * a + 4 * b + c);
        assert(idx < data_.size());
        return data_[idx];
    }

    Rank3<Scalar> tensor_at(std::size_t p) const {
        Rank3<Scalar> t;
        const std::size_t base = p * kComponents;
        for (int c = 0; c < kComponents; ++c) t.t[c] = data_[base + static_cast<std::size_t>(c)];
        return t;
    }

    void set_tensor(std::size_t p, const Rank3<Scalar>& t) {
        const std::size_t base = p * kComponents;
        for (int c = 0; c < kComponents; ++c) data_[base + static_cast<std::size_t>(c)] = t.t[c];
    }

    std::vector<Scalar>& raw() { return data_; }
    const std::vector<Scalar>& raw() const { return data_; }

private:
    const Grid* grid_ = nullptr;
    std::vector<Scalar> data_;
};

template <typename Scalar>
class Rank4Field {
public:
    static constexpr int kComponents = 256;

    Rank4Field() = default;

    explicit Rank4Field(const Grid& grid)
        : grid_(&grid),
          data_(grid.num_points() * kComponents, scalar_cast<Scalar>(0.0)) {}

    const Grid& grid() const { return *grid_; }
    std::size_t num_points() const { return data_.size() / kComponents; }

    Scalar& at(std::size_t p, int a, int b, int c, int d) {
        const std::size_t idx =
            p * kComponents +
            static_cast<std::size_t>(64 * a + 16 * b + 4 * c + d);
        assert(idx < data_.size());
        return data_[idx];
    }
    const Scalar& at(std::size_t p, int a, int b, int c, int d) const {
        const std::size_t idx =
            p * kComponents +
            static_cast<std::size_t>(64 * a + 16 * b + 4 * c + d);
        assert(idx < data_.size());
        return data_[idx];
    }

    Rank4<Scalar> tensor_at(std::size_t p) const {
        Rank4<Scalar> t;
        const std::size_t base = p * kComponents;
        for (int c = 0; c < kComponents; ++c) t.t[c] = data_[base + static_cast<std::size_t>(c)];
        return t;
    }

    void set_tensor(std::size_t p, const Rank4<Scalar>& t) {
        const std::size_t base = p * kComponents;
        for (int c = 0; c < kComponents; ++c) data_[base + static_cast<std::size_t>(c)] = t.t[c];
    }

    std::vector<Scalar>& raw() { return data_; }
    const std::vector<Scalar>& raw() const { return data_; }

private:
    const Grid* grid_ = nullptr;
    std::vector<Scalar> data_;
};

}  // namespace warp_hypad

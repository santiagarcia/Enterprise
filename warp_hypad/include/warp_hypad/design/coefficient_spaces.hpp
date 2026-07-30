#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "warp_hypad/design/constraint_projection.hpp"

// ---------------------------------------------------------------------------
// O0B Phase 4: the three coefficient spaces, as DISTINCT TYPES.
//
//     c = c0 + P z        z = W y        c = c0 + P W y
//
// raw c            one number per active tensor product; local support
// constrained z    an unconstrained coordinate on the null space of C
// whitened y       a coordinate in which the Sobolev Gram is the identity
//
// All three are vectors of doubles, all three have plausible lengths, and two
// of them (z and y) have the SAME length whenever nothing is discarded. A bare
// std::vector<double> therefore cannot say which one it is, and passing a
// whitened vector where a raw one was wanted produces a perfectly finite,
// entirely wrong metric. The wrappers below make that a compile error instead.
//
// Conversions run in ONE DIRECTION ONLY -- y -> z -> c. The reverse requires a
// pseudo-inverse and is not needed by any O0B consumer, so it is not provided:
// an inverse that existed would eventually be used to "round-trip" a raw vector
// that has no exact preimage, and the discrepancy would surface as an
// unexplained residual much later.
//
// WHITENED INDICES HAVE NO PHYSICAL MODE IDENTITY. The whitened Gram spectrum
// is the identity, so EVERY orthonormal rotation of the whitened space is an
// equally valid choice of axes -- the spectrum is maximally degenerate. There
// is no sense in which whitened coordinate 7 is "the front-wall mode", and
// nothing here or downstream may label one that way unless a separate
// localization convention is constructed and verified. The degeneracy is
// recorded in the manifest so the claim cannot be made by accident.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// FNV-1a 64-bit over the bytes of a double vector, for the transform checksum.
// 0xCBF29CE484222325 is the offset basis; written in hex because the decimal
// form lost a digit once already in this repository.
inline std::uint64_t checksum_of(const std::vector<double>& v) {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (double d : v) {
        const unsigned char* b = reinterpret_cast<const unsigned char*>(&d);
        for (std::size_t i = 0; i < sizeof(double); ++i) {
            h ^= static_cast<std::uint64_t>(b[i]);
            h *= 0x100000001B3ULL;
        }
    }
    return h;
}

inline std::uint64_t checksum_of(const DenseMatrix& m) {
    return checksum_of(m.a);
}

inline std::string hex_checksum(std::uint64_t v) {
    static const char* digits = "0123456789abcdef";
    std::string s(16, '0');
    for (int i = 15; i >= 0; --i) { s[static_cast<std::size_t>(i)] =
        digits[v & 0xF]; v >>= 4; }
    return s;
}

// --- The three spaces --------------------------------------------------------

struct WhitenedDesignCoordinates {
    std::vector<double> y;
    static const char* space_id() { return "whitened_y"; }
    std::size_t size() const { return y.size(); }
};

struct ConstrainedSplineCoordinates {
    std::vector<double> z;
    static const char* space_id() { return "constrained_z"; }
    std::size_t size() const { return z.size(); }
};

struct RawSplineCoefficients {
    std::vector<double> c;
    static const char* space_id() { return "raw_c"; }
    std::size_t size() const { return c.size(); }
};

// Why a conversion can fail, reported rather than asserted so the gate can
// exercise each rejection.
enum class ConversionStatus {
    Ok,
    DimensionMismatch,
    NonFiniteInput,
};

inline const char* to_string(ConversionStatus s) {
    switch (s) {
        case ConversionStatus::Ok: return "ok";
        case ConversionStatus::DimensionMismatch: return "dimension_mismatch";
        case ConversionStatus::NonFiniteInput: return "nonfinite_input";
    }
    return "unknown";
}

inline bool all_finite(const std::vector<double>& v) {
    for (double d : v)
        if (!(d == d) || d == std::numeric_limits<double>::infinity() ||
            d == -std::numeric_limits<double>::infinity())
            return false;
    return true;
}

// The one object that owns the transforms and their recorded metadata.
class CoefficientSpaces {
public:
    CoefficientSpaces(const ConstraintProjection& projection,
                      const WhiteningResult& whitening)
        : pr_(projection), wr_(whitening) {
        pw_ = multiply(pr_.P, wr_.w);
    }

    int raw_dimension() const { return pr_.n_raw; }
    int constraint_rank() const { return pr_.rank; }
    int constrained_dimension() const { return pr_.n_free; }
    int whitened_dimension() const { return wr_.retained; }

    const DenseMatrix& P() const { return pr_.P; }
    const DenseMatrix& W() const { return wr_.w; }
    const DenseMatrix& PW() const { return pw_; }
    const std::vector<double>& c0() const { return pr_.c0; }

    // z = W y
    ConversionStatus to_constrained(const WhitenedDesignCoordinates& in,
                                    ConstrainedSplineCoordinates& out) const {
        if (static_cast<int>(in.y.size()) != wr_.retained)
            return ConversionStatus::DimensionMismatch;
        if (!all_finite(in.y)) return ConversionStatus::NonFiniteInput;
        out.z.assign(static_cast<std::size_t>(pr_.n_free), 0.0);
        for (int i = 0; i < pr_.n_free; ++i) {
            double acc = 0.0;
            for (int j = 0; j < wr_.retained; ++j)
                acc += wr_.w(i, j) * in.y[static_cast<std::size_t>(j)];
            out.z[static_cast<std::size_t>(i)] = acc;
        }
        return ConversionStatus::Ok;
    }

    // c = c0 + P z
    ConversionStatus to_raw(const ConstrainedSplineCoordinates& in,
                            RawSplineCoefficients& out) const {
        if (static_cast<int>(in.z.size()) != pr_.n_free)
            return ConversionStatus::DimensionMismatch;
        if (!all_finite(in.z)) return ConversionStatus::NonFiniteInput;
        out.c.assign(static_cast<std::size_t>(pr_.n_raw), 0.0);
        for (int i = 0; i < pr_.n_raw; ++i) {
            double acc = pr_.c0[static_cast<std::size_t>(i)];
            for (int j = 0; j < pr_.n_free; ++j)
                acc += pr_.P(i, j) * in.z[static_cast<std::size_t>(j)];
            out.c[static_cast<std::size_t>(i)] = acc;
        }
        return ConversionStatus::Ok;
    }

    // c = c0 + P W y, the composed route. Checked against the two-step route by
    // the gate: the composition is what an optimizer would actually call, and a
    // discrepancy between it and P(W y) would mean the recorded PW is not the
    // product it claims to be.
    ConversionStatus to_raw(const WhitenedDesignCoordinates& in,
                            RawSplineCoefficients& out) const {
        if (static_cast<int>(in.y.size()) != wr_.retained)
            return ConversionStatus::DimensionMismatch;
        if (!all_finite(in.y)) return ConversionStatus::NonFiniteInput;
        out.c.assign(static_cast<std::size_t>(pr_.n_raw), 0.0);
        for (int i = 0; i < pr_.n_raw; ++i) {
            double acc = pr_.c0[static_cast<std::size_t>(i)];
            for (int j = 0; j < wr_.retained; ++j)
                acc += pw_(i, j) * in.y[static_cast<std::size_t>(j)];
            out.c[static_cast<std::size_t>(i)] = acc;
        }
        return ConversionStatus::Ok;
    }

    // Everything the manifest has to record about the transforms.
    struct Record {
        int raw_dimension;
        int constraint_rank;
        int constrained_dimension;
        int whitened_dimension;
        std::string p_checksum;
        std::string w_checksum;
        std::string pw_checksum;
        std::string c0_checksum;
        std::string whitened_identity_note;
    };

    Record record() const {
        Record r;
        r.raw_dimension = pr_.n_raw;
        r.constraint_rank = pr_.rank;
        r.constrained_dimension = pr_.n_free;
        r.whitened_dimension = wr_.retained;
        r.p_checksum = hex_checksum(checksum_of(pr_.P));
        r.w_checksum = hex_checksum(checksum_of(wr_.w));
        r.pw_checksum = hex_checksum(checksum_of(pw_));
        r.c0_checksum = hex_checksum(checksum_of(pr_.c0));
        r.whitened_identity_note =
            "The whitened Gram is the identity, so its spectrum is maximally "
            "degenerate and every orthonormal rotation of the whitened axes is "
            "equally valid. A whitened coordinate index therefore has NO unique "
            "physical mode identity, and must not be labelled as a named "
            "geometric mode without a separate, verified localization "
            "convention.";
        return r;
    }

private:
    ConstraintProjection pr_;
    WhiteningResult wr_;
    DenseMatrix pw_;
};

}  // namespace design
}  // namespace warp_hypad

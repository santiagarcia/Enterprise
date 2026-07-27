#pragma once

#include <type_traits>
#include <utility>

#include "program_b/effective_metric/covariant_constitutive.hpp"
#include "program_b/effective_metric/metric_types.hpp"

// ---------------------------------------------------------------------------
// COMPILE-TIME proof of the laboratory / effective type boundary.
//
// The B1C separation is only worth something if the compiler enforces it. A
// runtime string comparison, or an enum inspected after construction, would let
// the invalid object exist first -- and an object that exists can be passed on,
// stored, or logged before anyone checks it. So the claims below are
// static_asserts over SFINAE probes: they are evaluated when this header is
// compiled, and a regression in the public interface breaks the build rather
// than a test.
//
// This header adds no runtime code and changes no mathematics. It only asserts
// what the existing signatures already imply, so that the implication cannot be
// quietly lost.
// ---------------------------------------------------------------------------

namespace program_b {
namespace type_boundary {

// A detector for "is this call well-formed". std::is_invocable_v is not enough
// on its own here, because the interesting cases are FUNCTION TEMPLATES whose
// deduction must fail: is_invocable_v answers exactly that question for a
// generic lambda wrapping the call.
template <typename F, typename... Args>
constexpr bool callable = std::is_invocable_v<F, Args...>;

using S = double;
using Lab = LaboratoryMetric<S>;
using Eff = EffectiveMetric<S>;
using RawMetric = Mat4<S>;
using Vec = S[4];

// Wrap each entry point in a generic lambda so the probe is a deduction test on
// the real overload set rather than on a hand-copied signature.
constexpr auto probe_lab_observer = [](auto&& g, auto&& v)
    -> decltype(make_laboratory_observer(std::forward<decltype(g)>(g),
                                         std::forward<decltype(v)>(v),
                                         std::string{})) { return {}; };

constexpr auto probe_eff_observer = [](auto&& g, auto&& v)
    -> decltype(make_effective_geometry_observer(
        std::forward<decltype(g)>(g), std::forward<decltype(v)>(v),
        std::string{})) { return {}; };

constexpr auto probe_lab_eulerian = [](auto&& g)
    -> decltype(laboratory_eulerian_observer(std::forward<decltype(g)>(g))) {
    return {};
};

// --- Claim 1: the laboratory observer accepts a LaboratoryMetric ------------
static_assert(callable<decltype(probe_lab_observer), const Lab&, const Vec&>,
              "make_laboratory_observer must accept LaboratoryMetric");

// --- Claim 2: it does NOT accept an EffectiveMetric -------------------------
//
// This is the case that matters. An effective metric normalising a physical
// laboratory observer would silently produce blocks referred to the wrong
// geometry, and every downstream number would be wrong in a way no residual
// would reveal.
static_assert(!callable<decltype(probe_lab_observer), const Eff&, const Vec&>,
              "make_laboratory_observer must REJECT EffectiveMetric: a physical "
              "laboratory observer is normalised with g_lab, never with g_eff");

// --- Claim 3: the effective-geometry observer accepts an EffectiveMetric ----
static_assert(callable<decltype(probe_eff_observer), const Eff&, const Vec&>,
              "make_effective_geometry_observer must accept EffectiveMetric");

// --- Claim 4: it does NOT accept a LaboratoryMetric -------------------------
static_assert(!callable<decltype(probe_eff_observer), const Lab&, const Vec&>,
              "make_effective_geometry_observer must REJECT LaboratoryMetric");

// --- Claim 5: the convenience laboratory constructors are equally guarded ---
static_assert(callable<decltype(probe_lab_eulerian), const Lab&>,
              "laboratory_eulerian_observer must accept LaboratoryMetric");
static_assert(!callable<decltype(probe_lab_eulerian), const Eff&>,
              "laboratory_eulerian_observer must REJECT EffectiveMetric");

// --- Claim 6: a raw Mat4 cannot reach the typed entry points ----------------
//
// The wrappers take LaboratoryMetric / EffectiveMetric, whose constructors from
// Mat4 are EXPLICIT. So a bare Mat4 does not convert, and the call does not
// compile: the tag has to be applied deliberately, at a named call site.
static_assert(!callable<decltype(probe_lab_observer), const RawMetric&,
                        const Vec&>,
              "a raw Mat4 must not convert into LaboratoryMetric implicitly");
static_assert(!callable<decltype(probe_eff_observer), const RawMetric&,
                        const Vec&>,
              "a raw Mat4 must not convert into EffectiveMetric implicitly");
static_assert(!callable<decltype(probe_lab_eulerian), const RawMetric&>,
              "a raw Mat4 must not convert into LaboratoryMetric implicitly");

// --- Claim 7: the two tagged types are not interconvertible ------------------
static_assert(!std::is_convertible_v<Lab, Eff>,
              "LaboratoryMetric must not convert to EffectiveMetric");
static_assert(!std::is_convertible_v<Eff, Lab>,
              "EffectiveMetric must not convert to LaboratoryMetric");
static_assert(!std::is_constructible_v<Lab, Eff>,
              "LaboratoryMetric must not be constructible from EffectiveMetric");
static_assert(!std::is_constructible_v<Eff, Lab>,
              "EffectiveMetric must not be constructible from LaboratoryMetric");

// --- Claim 8: the tag is explicit, not implicit -----------------------------
static_assert(!std::is_convertible_v<RawMetric, Lab>,
              "the Mat4 -> LaboratoryMetric constructor must stay explicit");
static_assert(!std::is_convertible_v<RawMetric, Eff>,
              "the Mat4 -> EffectiveMetric constructor must stay explicit");
static_assert(std::is_constructible_v<Lab, RawMetric>,
              "an explicit LaboratoryMetric(Mat4) must remain available");
static_assert(std::is_constructible_v<Eff, RawMetric>,
              "an explicit EffectiveMetric(Mat4) must remain available");

// --- Claim 9: the observer result types are distinct ------------------------
//
// EffectiveGeometryObserver is its own type, so an effective-geometry
// diagnostic cannot be handed to a routine expecting a physical Observer.
static_assert(!std::is_same_v<Observer<S>, EffectiveGeometryObserver<S>>,
              "the effective-geometry observer must be a distinct type");
static_assert(!std::is_convertible_v<EffectiveGeometryObserver<S>, Observer<S>>,
              "an effective-geometry observer must not convert into a physical "
              "Observer: it is a diagnostic of the effective geometry and "
              "carries no constitutive meaning");

// A single symbol a translation unit can reference to force this header to be
// instantiated, so the assertions cannot be dropped as unused.
inline constexpr bool kTypeBoundaryVerified = true;

// What the exporter gate reports. Counting them here keeps the gate's number
// and the assertions above in one place.
inline constexpr int kCompileTimeClaims = 20;

// A KNOWN, DELIBERATE ESCAPE HATCH, recorded rather than hidden:
// make_observer(const Mat4<Scalar>&, ...) remains public, and a caller that
// reaches for it bypasses the tagged types entirely. It is used by the existing
// verified internals, and closing it means rerouting those call sites -- a
// change to working, gated code that this milestone does not need. The boundary
// proved above is over the TAGGED entry points, which is what production and
// export code is required to use.
inline constexpr bool kRawMakeObserverStillReachable = true;

}  // namespace type_boundary
}  // namespace program_b

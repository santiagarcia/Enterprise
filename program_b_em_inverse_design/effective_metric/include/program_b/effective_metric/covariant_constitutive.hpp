#pragma once

#include <algorithm>
#include <string>

#include "program_b/constitutive/constitutive.hpp"
#include "program_b/effective_metric/metric_types.hpp"
#include "warp_hypad/core/tensor.hpp"

// ---------------------------------------------------------------------------
// COVARIANT constitutive tensor, and its observer-dependent decomposition.
//
// This file exists because B1's original conclusion was WRONG. B1 computed the
// Plebanski 3+1 map, found it singular where g_00 >= 0, and reported "no
// equivalent medium exists". That conflated two different statements:
//
//   (a) the 3+1 decomposition RELATIVE TO THE COORDINATE-TIME OBSERVER is
//       singular there -- true, and provable;
//   (b) no covariant constitutive representation exists there -- NOT proved,
//       and in fact false for the Alcubierre metric.
//
// The Alcubierre metric has det g = -1 and signature (-,+,+,+) EVERYWHERE, for
// every v_s including v_s = 5, because its (t,x) block has determinant -1 and
// therefore always exactly one negative eigenvalue. A Lorentzian nondegenerate
// metric always admits Z^{abcd}. What fails at v_s f = 1 is that d_t stops
// being timelike -- an ERGOSURFACE, structurally identical to Kerr's.
//
// This is a known result: susceptibility eigenvalues diverge on the
// ergosurface while the covariant tensor does not.
//
// PRIMARY SOURCE (verified):
//   S. Schuster and M. Visser, "Effective metrics and a fully covariant
//   description of constitutive tensors in electrodynamics",
//   Physical Review D 96, 124019 (2017). DOI 10.1103/PhysRevD.96.124019.
//   arXiv:1706.06280. Abstract verified against the arXiv record; that paper
//   introduces Z^{abcd} and identifies eps, mu^-1 and zeta as the three pieces
//   of its Bel decomposition with respect to an observer four-velocity -- the
//   structure this file implements.
//
// CONVENTIONS USED HERE (see docs/metric_constitutive_conventions.md section 10):
//   signature (-,+,+,+); Levi-Civita symbol with eps_{0123} = +1 in the
//   orthonormal frame; Heaviside-Lorentz with c = 1; Z carries four UPPER
//   indices; the constitutive relation is H^{ab} = Z^{abcd} F_{cd}.
//
//   The factor is pinned, not assumed: with the definition below, Minkowski
//   gives Z^{abcd} F_{cd} = F^{ab}, i.e. exact vacuum. Verified as
//   b1c.covariant_minkowski_is_vacuum.
// ---------------------------------------------------------------------------

namespace program_b {

using warp_hypad::Rank4;

// The five distinguishable outcomes. Deliberately NOT collapsed into one
// "no medium" verdict -- that collapse is the error this file corrects.
enum class MetricMediumStatus {
    // g_eff Lorentzian, Z^{abcd} exists, chosen observer admissible, and the
    // 3+1 tensors are finite and well conditioned.
    CovariantMetricMediumValid,

    // Z^{abcd} exists and is finite. The CHOSEN OBSERVER is not admissible
    // here (not timelike / not normalisable). Says nothing about the metric or
    // about other observers. This is the ergosurface case.
    SelectedObserverDecompositionInvalid,

    // Z^{abcd} exists and the observer is admissible, but the resulting 3+1
    // material tensors diverge or are too ill conditioned to be meaningful.
    ConstitutiveParametersSingular,

    // The response cannot come from any single metric (skewon/axion/
    // birefringent). A statement about the CONSTITUTIVE TENSOR, not the
    // observer.
    NonMetricConstitutiveTensor,

    // g_eff itself is degenerate or not Lorentzian. The only status that
    // genuinely means "no metric medium here".
    MetricDegenerateOrNonLorentzian,
};

const char* to_string(MetricMediumStatus s);
const char* describe(MetricMediumStatus s);

// --- The covariant constitutive tensor ---------------------------------------
//
//   Z^{abcd} = (1/2) sqrt(det g_eff / det g_lab)
//              [ g_eff^{ac} g_eff^{bd} - g_eff^{ad} g_eff^{bc} ]
//
// Both determinants are negative for Lorentzian metrics, so their ratio is
// positive and the square root is real. The density factor is what makes Z a
// tensor DENSITY of the right weight relative to the laboratory background.
//
// Computed with no reference whatsoever to any 3+1 split. That independence is
// the whole point.

template <typename Scalar>
struct CovariantConstitutive {
    Rank4<Scalar> Z;
    Scalar density_factor;      // sqrt(det g_eff / det g_lab)
    double det_eff_real = 0.0;
    double det_lab_real = 0.0;
    bool ok = false;
    MetricMediumStatus status = MetricMediumStatus::MetricDegenerateOrNonLorentzian;
    std::string message;
};

template <typename Scalar>
CovariantConstitutive<Scalar> covariant_constitutive(
    const Mat4<Scalar>& g_eff, const Mat4<Scalar>& g_lab) {
    using std::sqrt;
    CovariantConstitutive<Scalar> out;

    const warp_hypad::InverseResult<Scalar> ie = warp_hypad::inverse(g_eff);
    const warp_hypad::InverseResult<Scalar> il = warp_hypad::inverse(g_lab);
    out.det_eff_real = ie.det_real;
    out.det_lab_real = il.det_real;

    if (!ie.ok || !il.ok) {
        out.message = "metric singular: cannot invert g_eff or g_lab.";
        return out;
    }
    if (!(ie.det_real < 0.0) || !(il.det_real < 0.0)) {
        out.message =
            "det g >= 0: not Lorentzian in the (-,+,+,+) signature. This is the "
            "ONLY condition under which no covariant metric medium exists.";
        return out;
    }

    out.density_factor = sqrt(ie.determinant / il.determinant);
    const Scalar half = scalar_cast<Scalar>(0.5);

    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    out.Z(a, b, c, d) =
                        half * out.density_factor *
                        (ie.inverse(a, c) * ie.inverse(b, d) -
                         ie.inverse(a, d) * ie.inverse(b, c));

    out.ok = true;
    out.status = MetricMediumStatus::CovariantMetricMediumValid;
    out.message = "covariant constitutive tensor exists and is finite.";
    return out;
}

// Type-safe entry point. The raw-Mat4 overload above cannot tell which argument
// is which; this one cannot be called wrongly. New code uses this form.
template <typename Scalar>
CovariantConstitutive<Scalar> covariant_constitutive(
    const EffectiveMetric<Scalar>& g_eff, const LaboratoryMetric<Scalar>& g_lab) {
    return covariant_constitutive(g_eff.value, g_lab.value);
}

// --- Observers ----------------------------------------------------------------

template <typename Scalar>
struct Observer {
    Scalar V[4];              // contravariant four-velocity V^a
    std::string name;

    // g_lab_{ab} V^a V^b, which must equal -1.
    double norm_real = 0.0;
    bool admissible = false;
    std::string reason;
};

// Normalise a candidate four-velocity against the LABORATORY metric.
// Admissible iff g_lab(V,V) < 0, i.e. V is timelike. NOT assumed proportional
// to d_t.
template <typename Scalar>
Observer<Scalar> make_observer(const Mat4<Scalar>& g_lab, const Scalar V_in[4],
                               const std::string& name) {
    using std::sqrt;
    Observer<Scalar> o;
    o.name = name;

    Scalar n = scalar_cast<Scalar>(0.0);
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) n = n + g_lab(a, b) * V_in[a] * V_in[b];
    o.norm_real = real_part(n);

    if (!(o.norm_real < 0.0)) {
        o.reason =
            "g_lab(V,V) >= 0: the candidate four-velocity is null or spacelike, "
            "so no observer moves along it. This is a statement about the "
            "OBSERVER, not about the metric or the medium.";
        for (int a = 0; a < 4; ++a) o.V[a] = V_in[a];
        return o;
    }

    const Scalar s = scalar_cast<Scalar>(1.0) / sqrt(-n);
    for (int a = 0; a < 4; ++a) o.V[a] = V_in[a] * s;
    o.admissible = true;
    o.reason = "timelike and normalised to g_lab(V,V) = -1.";
    return o;
}

// The coordinate-time observer, V ~ d_t. Admissible only where g_lab_00 < 0 --
// this is the observer whose failure B1 originally mistook for a failure of the
// medium.
template <typename Scalar>
Observer<Scalar> coordinate_time_observer(const Mat4<Scalar>& g_lab) {
    Scalar V[4] = {scalar_cast<Scalar>(1.0), scalar_cast<Scalar>(0.0),
                   scalar_cast<Scalar>(0.0), scalar_cast<Scalar>(0.0)};
    return make_observer(g_lab, V, "coordinate_time");
}

// The Eulerian observer: normal to the constant-t slices, V_a ~ -N delta^0_a,
// equivalently V^a ~ g^{a0} / sqrt(-g^{00}). Timelike wherever the SLICING is
// spacelike (g^{00} < 0), which for a warp metric holds even inside the
// ergoregion -- that is precisely why it is the useful alternative.
template <typename Scalar>
Observer<Scalar> eulerian_observer(const Mat4<Scalar>& g_lab) {
    Observer<Scalar> bad;
    bad.name = "eulerian";
    const warp_hypad::InverseResult<Scalar> gi = warp_hypad::inverse(g_lab);
    if (!gi.ok) {
        bad.reason = "g_lab is singular.";
        for (int a = 0; a < 4; ++a) bad.V[a] = scalar_cast<Scalar>(0.0);
        return bad;
    }
    if (!(real_part(gi.inverse(0, 0)) < 0.0)) {
        bad.reason =
            "g^{00} >= 0: the constant-t slicing is not spacelike here, so no "
            "Eulerian observer exists either.";
        for (int a = 0; a < 4; ++a) bad.V[a] = scalar_cast<Scalar>(0.0);
        return bad;
    }
    Scalar V[4];
    for (int a = 0; a < 4; ++a) V[a] = gi.inverse(a, 0);
    return make_observer(g_lab, V, "eulerian");
}

// A boosted inertial laboratory observer: V ~ gamma(d_t + beta d_x), for a
// FLAT laboratory background.
template <typename Scalar>
Observer<Scalar> boosted_observer(const Mat4<Scalar>& g_lab, double beta) {
    Scalar V[4] = {scalar_cast<Scalar>(1.0), scalar_cast<Scalar>(beta),
                   scalar_cast<Scalar>(0.0), scalar_cast<Scalar>(0.0)};
    return make_observer(g_lab, V, "boosted_beta_" + std::to_string(beta));
}

// --- The two observer APIs, which must not be confused ------------------------
//
// A PHYSICAL laboratory observer is normalised with g_lab. It is the observer
// the constitutive blocks eps, mu^-1, zeta, zeta+ belong to; it is the V that
// enters (A7a)-(A7d) and (A11); and it is the only kind of observer for which
// "the medium responds this way" is a statement about an experiment.
//
// A DIAGNOSTIC observer of the effective geometry is normalised with g_eff. It
// answers questions about the analogue spacetime -- is this surface an
// ergosurface, is that worldline timelike in the effective geometry -- and it
// has NO constitutive interpretation whatsoever. Feeding one into the
// decomposition would silently produce blocks belonging to no observer at all.
//
// The two return DIFFERENT TYPES, so the substitution cannot compile.

// V normalised against the physical laboratory background.
template <typename Scalar>
Observer<Scalar> make_laboratory_observer(const LaboratoryMetric<Scalar>& g_lab,
                                          const Scalar V_in[4],
                                          const std::string& name) {
    return make_observer(g_lab.value, V_in, name);
}

// The three named laboratory observers, in typed form. The raw-Mat4 versions
// above are the primitives; new code should reach for these.
template <typename Scalar>
Observer<Scalar> laboratory_coordinate_time_observer(
    const LaboratoryMetric<Scalar>& g_lab) {
    return coordinate_time_observer(g_lab.value);
}

template <typename Scalar>
Observer<Scalar> laboratory_eulerian_observer(const LaboratoryMetric<Scalar>& g_lab) {
    return eulerian_observer(g_lab.value);
}

template <typename Scalar>
Observer<Scalar> laboratory_boosted_observer(const LaboratoryMetric<Scalar>& g_lab,
                                             double beta) {
    return boosted_observer(g_lab.value, beta);
}

// A four-velocity normalised against g_eff. Structurally identical to Observer
// and deliberately NOT convertible to it.
template <typename Scalar>
struct EffectiveGeometryObserver {
    Scalar V[4];              // contravariant four-velocity in the EFFECTIVE geometry
    std::string name;

    // g_eff_{ab} V^a V^b, which must equal -1.
    double norm_real = 0.0;
    bool admissible = false;
    std::string reason;
};

template <typename Scalar>
EffectiveGeometryObserver<Scalar> make_effective_geometry_observer(
    const EffectiveMetric<Scalar>& g_eff, const Scalar V_in[4],
    const std::string& name) {
    const Observer<Scalar> o = make_observer(g_eff.value, V_in, name);
    EffectiveGeometryObserver<Scalar> e;
    for (int a = 0; a < 4; ++a) e.V[a] = o.V[a];
    e.name = name;
    e.norm_real = o.norm_real;
    e.admissible = o.admissible;
    e.reason = o.admissible
                   ? "timelike in g_eff and normalised to g_eff(V,V) = -1. This "
                     "is a DIAGNOSTIC of the effective geometry and carries no "
                     "constitutive meaning."
                   : o.reason;
    return e;
}

// --- Conditioning diagnostics -------------------------------------------------
//
// Small numeric helpers shared by the covariant path and, for now, by the
// quarantined legacy decomposition in legacy_unverified.hpp.

// Symmetric 3x3 eigenvalues by cyclic Jacobi, ascending. Real parts only:
// eigenvalues are used for conditioning diagnostics, which are comparisons.
inline void symmetric_eigenvalues_3x3(const double a_in[3][3], double out[3]) {
    double a[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) a[i][j] = 0.5 * (a_in[i][j] + a_in[j][i]);

    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j) off = std::fmax(off, std::fabs(a[i][j]));
        double scale = 0.0;
        for (int i = 0; i < 3; ++i) scale = std::fmax(scale, std::fabs(a[i][i]));
        if (off <= 1e-18 * std::fmax(1.0, scale)) break;

        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::fabs(a[p][q]) <= 1e-300) continue;
                const double th = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double tt = (th >= 0.0 ? 1.0 : -1.0) /
                                  (std::fabs(th) + std::sqrt(th * th + 1.0));
                const double c = 1.0 / std::sqrt(tt * tt + 1.0);
                const double s = tt * c;
                for (int k = 0; k < 3; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
    }
    for (int i = 0; i < 3; ++i) out[i] = a[i][i];
    std::sort(out, out + 3);
}


}  // namespace program_b

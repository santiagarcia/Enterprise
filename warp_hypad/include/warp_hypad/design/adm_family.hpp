#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "warp_hypad/design/jet3.hpp"
#include "warp_hypad/design/energy_objectives.hpp"
#include "warp_hypad/design/shift_family.hpp"
#include "warp_hypad/design/wall_family.hpp"

// ---------------------------------------------------------------------------
// F2A: the first EXPANDED ADM family.
//
//     alpha    = exp( a(xi, s) )                 lapse correction
//     beta^x   = beta_base(xi, s) + Bx(xi, s)    axial shift correction
//     beta^rho = (rho/R0) Q(xi, s)               radial shift correction
//     gamma    = I                               unchanged in F2A
//
// WHY THESE THREE, IN THIS ORDER. O1 established that the shift-only family
// F1W buys 0.33 per cent of integrated negative energy and pays for it in
// curvature and in the NEC. That is a statement about the FAMILY, not about the
// objective: with alpha frozen at 1 and gamma frozen at I, the only thing a
// correction can do to rho is reshape one shift component. The lapse is the
// obvious next degree of freedom because it enters g_00 directly, and the
// radial shift because an axisymmetric wall has a radial direction the axial
// shift cannot reach.
//
// AXIS REGULARITY IS STRUCTURAL, NOT ENFORCED.
//
//   a and Bx are SCALARS in (xi, s). Since s = (y^2+z^2)/R0^2, any function of
//   s is automatically even in rho with vanishing first radial derivative --
//   O0A's whole reason for using s rather than rho.
//
//   beta^rho must VANISH on the axis, and it does so algebraically rather than
//   by a limit: the Cartesian components are
//
//       beta^y = (y/R0) Q(xi, s),   beta^z = (z/R0) Q(xi, s)
//
//   so the rho in beta^rho = (rho/R0) Q cancels against the y/rho of the
//   projection. There is no division by rho anywhere and no small-rho branch.
//   This is exactly O0A's RadialVector family, which has been present and
//   unused since that milestone.
//
// THE PROTECTED REGIONS STAY STRUCTURAL. All three fields are expanded on the
// SAME support-excluded F1W basis, so a correction to the lapse or the radial
// shift is as incapable of touching the cabin or the collar as the axial one
// was. Delta alpha, Delta beta^x and Delta beta^rho are all identically zero
// there, with all their derivatives, for every coordinate vector.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// Which fields a family switches on. The hierarchy is explicit so an escalation
// is a data change rather than a new code path.
struct AdmFamilySpec {
    bool lapse = false;          // a(xi, s)
    bool axial_shift = true;     // Bx(xi, s)
    bool radial_shift = false;   // Q(xi, s), through the eta prefactor
    bool spatial_diagonal = false;   // F2B
    bool spatial_cross = false;      // F2C
    const char* name = "F1W";

    int field_count() const {
        return (lapse ? 1 : 0) + (axial_shift ? 1 : 0) +
               (radial_shift ? 1 : 0) + (spatial_diagonal ? 2 : 0) +
               (spatial_cross ? 1 : 0);
    }
};

inline AdmFamilySpec f1w_spec() {
    AdmFamilySpec s;
    s.lapse = false; s.axial_shift = true; s.radial_shift = false;
    s.name = "F1W";
    return s;
}

inline AdmFamilySpec f2a_spec() {
    AdmFamilySpec s;
    s.lapse = true; s.axial_shift = true; s.radial_shift = true;
    s.name = "F2A";
    return s;
}

inline AdmFamilySpec f2b_spec() {
    AdmFamilySpec s = f2a_spec();
    s.spatial_diagonal = true;
    s.name = "F2B";
    return s;
}

// One scalar design field's spacetime jet, from the verified O0A/O0B route:
// the (xi, s) derivative table, then chained_derivative for every multi-index.
// Nothing new is derived here; the jet is assembled from calls.
template <typename Scalar>
Jet3<Scalar> scalar_field_jet(const AxisymmetricBasis& basis,
                              const std::vector<Scalar>& coefficients,
                              const DesignCoordinates& p,
                              const ChainFactors& F) {
    Jet3<Scalar> j;
    if (coefficients.empty()) return j;
    const TensorBasisEvaluation e = basis.evaluate(p, /*want_third=*/true);
    for (std::size_t m = 0; m < e.size(); ++m) {
        const Scalar& c =
            coefficients[static_cast<std::size_t>(e.flat_index[m])];
        const DesignDerivativeTable Bt = single_mode_table(basis, e, p, m);
        j.v = j.v + c * scalar_cast<Scalar>(Bt.at(0, 0));
        for (int i = 0; i < 4; ++i) {
            const int one[1] = {i};
            j.d1[i] = j.d1[i] +
                      c * scalar_cast<Scalar>(chained_derivative(Bt, F, one, 1));
            for (int k = 0; k < 4; ++k) {
                const int two[2] = {i, k};
                j.d2[i][k] =
                    j.d2[i][k] +
                    c * scalar_cast<Scalar>(chained_derivative(Bt, F, two, 2));
                for (int l = 0; l < 4; ++l) {
                    const int three[3] = {i, k, l};
                    j.d3[i][k][l] =
                        j.d3[i][k][l] +
                        c * scalar_cast<Scalar>(
                                chained_derivative(Bt, F, three, 3));
                }
            }
        }
    }
    return j;
}

// The Alcubierre baseline shift, as a jet. Same table, no coefficients.
template <typename Scalar>
Jet3<Scalar> baseline_shift_jet(const ShiftFamilyConfig& cfg,
                                const DesignCoordinates& p,
                                const ChainFactors& F) {
    const DesignDerivativeTable T = baseline_shift_table(
        p.xi, p.radial_squared, cfg.velocity, cfg.scale_R0, cfg.bubble_radius,
        cfg.wall_sigma);
    Jet3<Scalar> j;
    j.v = scalar_cast<Scalar>(T.at(0, 0));
    for (int i = 0; i < 4; ++i) {
        const int one[1] = {i};
        j.d1[i] = scalar_cast<Scalar>(chained_derivative(T, F, one, 1));
        for (int k = 0; k < 4; ++k) {
            const int two[2] = {i, k};
            j.d2[i][k] = scalar_cast<Scalar>(chained_derivative(T, F, two, 2));
            for (int l = 0; l < 4; ++l) {
                const int three[3] = {i, k, l};
                j.d3[i][k][l] =
                    scalar_cast<Scalar>(chained_derivative(T, F, three, 3));
            }
        }
    }
    return j;
}

// The ADM fields of a family, each as a jet, before assembly.
template <typename Scalar>
struct AdmFieldJets {
    Jet3<Scalar> log_lapse;    // a
    Jet3<Scalar> lapse;        // exp(a)
    Jet3<Scalar> beta[3];      // beta^x, beta^y, beta^z, CARTESIAN, index up
};

// THE FIELDS ARRIVE ALREADY PROLONGED, AND THAT IS NOT A CONVENIENCE.
//
// Each active field is handed in as a PARENT-LENGTH coefficient vector whose
// support-excluded entries are exact zeros -- the same prolongation F1W uses.
// The first version of this sliced a flat vector into parent-length blocks
// directly, which silently did two wrong things at once: it made every design
// vector 272 coefficients per field instead of 36, and it filled the excluded
// modes, quietly destroying the structural protected-region lock that the whole
// family depends on. Taking prolonged vectors makes that impossible to express.
template <typename Scalar>
struct AdmFieldCoefficients {
    std::vector<Scalar> log_lapse;      // parent-length, zeros outside support
    std::vector<Scalar> axial_shift;
    std::vector<Scalar> radial_shift;
};

template <typename Scalar>
AdmFieldJets<Scalar> adm_field_jets(const AdmFamilySpec& spec,
                                    const ShiftFamilyConfig& cfg,
                                    const AxisymmetricBasis& basis,
                                    const AdmFieldCoefficients<Scalar>& c,
                                    double t, double x, double y, double z) {
    const double R0 = cfg.scale_R0;
    const double xc = cfg.bubble_center_x0 + cfg.velocity * t;
    const DesignCoordinates p = make_design_coordinates(x - xc, y, z, 0.0, R0);
    const ChainFactors F = chain_factors(y, z, cfg.velocity, R0);

    AdmFieldJets<Scalar> out;

    // The baseline axial shift is always present; corrections are additive.
    out.beta[0] = baseline_shift_jet<Scalar>(cfg, p, F);

    if (spec.lapse && !c.log_lapse.empty())
        out.log_lapse = scalar_field_jet<Scalar>(basis, c.log_lapse, p, F);
    out.lapse = exp_jet(out.log_lapse);   // exp(0) = 1 when the lapse is off

    if (spec.axial_shift && !c.axial_shift.empty())
        out.beta[0] = out.beta[0] +
                      scalar_field_jet<Scalar>(basis, c.axial_shift, p, F);

    if (spec.radial_shift && !c.radial_shift.empty()) {
        // beta^y = (y/R0) Q,  beta^z = (z/R0) Q. The eta prefactor is applied
        // as a JET PRODUCT with the coordinate, so the vanishing on the axis
        // and every derivative of it come out of the arithmetic rather than
        // out of a hand-written limit.
        const Jet3<Scalar> Q =
            scalar_field_jet<Scalar>(basis, c.radial_shift, p, F);
        const Jet3<Scalar> jy = coordinate_jet<Scalar>(2, y);
        const Jet3<Scalar> jz = coordinate_jet<Scalar>(3, z);
        out.beta[1] = (jy * Q) * (1.0 / R0);
        out.beta[2] = (jz * Q) * (1.0 / R0);
    }
    return out;
}

// The covariant metric and its first three derivatives, for gamma = I.
//
//     g_00 = -alpha^2 + beta^i beta^i
//     g_0i = beta^i
//     g_ij = delta_ij
//
// Every entry is a Jet3 expression, so the derivative orders are produced by
// the arithmetic rather than transcribed.
template <typename Scalar>
MetricJet<Scalar> adm_metric_jet(const AdmFamilySpec& spec,
                                 const ShiftFamilyConfig& cfg,
                                 const AxisymmetricBasis& basis,
                                 const AdmFieldCoefficients<Scalar>& c,
                                 double t, double x, double y, double z) {
    const AdmFieldJets<Scalar> f =
        adm_field_jets<Scalar>(spec, cfg, basis, c, t, x, y, z);

    Jet3<Scalar> g00 = -(f.lapse * f.lapse);
    for (int i = 0; i < 3; ++i) g00 = g00 + f.beta[i] * f.beta[i];

    MetricJet<Scalar> jet;
    auto put = [&](CovariantMetric<Scalar>& m, int order, int A, int B, int D) {
        m.m(0, 0) = order == 0   ? g00.v
                    : order == 1 ? g00.d1[A]
                    : order == 2 ? g00.d2[A][B]
                                 : g00.d3[A][B][D];
        for (int i = 0; i < 3; ++i) {
            const Scalar b = order == 0   ? f.beta[i].v
                             : order == 1 ? f.beta[i].d1[A]
                             : order == 2 ? f.beta[i].d2[A][B]
                                          : f.beta[i].d3[A][B][D];
            m.m(0, i + 1) = b;
            m.m(i + 1, 0) = b;
        }
        // gamma = I is constant, so only the value carries the identity.
        if (order == 0)
            for (int i = 0; i < 3; ++i) m.m(i + 1, i + 1) = scalar_cast<Scalar>(1.0);
    };

    put(jet.g_cov, 0, 0, 0, 0);
    for (int a = 0; a < 4; ++a) {
        put(jet.dg[static_cast<std::size_t>(a)], 1, a, 0, 0);
        for (int b = 0; b < 4; ++b) {
            put(jet.d2g[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)],
                2, a, b, 0);
            for (int d = 0; d < 4; ++d)
                put(jet.d3g[static_cast<std::size_t>(a)]
                           [static_cast<std::size_t>(b)]
                           [static_cast<std::size_t>(d)], 3, a, b, d);
        }
    }

    ADMDesignFields<Scalar> adm;
    adm.lapse = f.lapse.v;
    for (int i = 0; i < 3; ++i)
        adm.shift_contravariant[static_cast<std::size_t>(i)] = f.beta[i].v;
    bool ok = false;
    jet.g_contra = assemble_contravariant(adm, ok);
    jet.valid = ok;
    jet.rejection =
        ok ? MetricRejection::None : MetricRejection::SingularSpatialMetric;
    return jet;
}

}  // namespace design
}  // namespace warp_hypad

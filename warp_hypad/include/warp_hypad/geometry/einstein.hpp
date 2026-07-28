#pragma once

#include <cmath>
#include <string>

#include "warp_hypad/core/field.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"
#include "warp_hypad/core/tensor.hpp"
#include "warp_hypad/geometry/curvature.hpp"

// ---------------------------------------------------------------------------
// Milestone A5: the Einstein tensor and the stress-energy INFERRED from it.
//
//     R_{mu nu}, R, g_{mu nu}  ->  G_{mu nu} = R_{mu nu} - (1/2) g_{mu nu} R
//                              ->  T_{mu nu}^inferred
//
// --- On the word "inferred" --------------------------------------------------
//
// The tensor this file produces is obtained by inserting a CHOSEN geometry into
// Einstein's equation and reading off what source would have to be present. It
// is not a material model, not a physical source, not a fabricated source and
// not a realizable source. No mechanism that generates it has been demonstrated,
// and A5 does not attempt one. The naming here is deliberate and the gate checks
// it: every field carrying that meaning is called `inferred_stress_energy`.
//
// --- On symmetrisation -------------------------------------------------------
//
// A4 measured that Ricci symmetry is an O(h^2) truncation statement rather than
// a roundoff identity, because its derivation passes through the product rule
// and a difference operator does not satisfy it. G_{mu nu} inherits exactly that
// asymmetry.
//
// So the raw tensor is kept. Symmetrising in place would destroy the only
// direct measure of how far the discretisation is from the continuum, and a
// large antisymmetric part would then be invisible rather than reported. The
// symmetric projection is available, separately, explicitly labelled, and the
// antisymmetry residual is exported alongside it so a caller can see what the
// projection threw away.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

// How a cosmological constant enters. Recorded on every result: the three
// readings below are numerically different tensors, and silently choosing one
// is how a de Sitter spacetime comes to be reported as containing matter.
enum class CosmologicalConstantTreatment {
    // G_{mu nu} alone. The pure geometric object. For de Sitter this is
    // NONZERO: G = -(3/L^2) g.
    GeometricEinsteinOnly,

    // G_{mu nu} + Lambda g_{mu nu}, kept on the geometry side of the field
    // equation. For de Sitter with Lambda = 3/L^2 this vanishes, which is the
    // statement that de Sitter is a VACUUM solution with a cosmological
    // constant -- it contains no matter.
    LambdaOnGeometricSide,

    // Lambda moved to the source side and reported as part of the inferred
    // stress-energy: T = (1/8pi)(G + Lambda g). A legitimate choice, but a
    // different physical statement, and never the default.
    LambdaAsEffectiveStressEnergy
};

inline const char* to_string(CosmologicalConstantTreatment t) {
    switch (t) {
        case CosmologicalConstantTreatment::GeometricEinsteinOnly:
            return "geometric_einstein_only";
        case CosmologicalConstantTreatment::LambdaOnGeometricSide:
            return "lambda_on_geometric_side";
        case CosmologicalConstantTreatment::LambdaAsEffectiveStressEnergy:
            return "lambda_as_effective_stress_energy";
    }
    return "unknown";
}

struct CosmologicalConstantPolicy {
    CosmologicalConstantTreatment treatment =
        CosmologicalConstantTreatment::GeometricEinsteinOnly;
    double lambda = 0.0;
    std::string justification;

    // The warp-drive default. Lambda = 0 unless a model explicitly says
    // otherwise; a nonzero Lambda smuggled in as a "small correction" would
    // change what the inferred source means.
    static CosmologicalConstantPolicy none() {
        CosmologicalConstantPolicy p;
        p.treatment = CosmologicalConstantTreatment::GeometricEinsteinOnly;
        p.lambda = 0.0;
        p.justification = "warp-drive metrics carry no cosmological constant";
        return p;
    }
};

template <typename Scalar>
struct EinsteinFields {
    // The directly computed tensor. NEVER overwritten by its own symmetric
    // projection.
    Rank2Field<Scalar> raw_einstein;

    // (1/2)(G_{mu nu} + G_{nu mu}). The physical symmetric projection used
    // downstream, exactly symmetric by construction.
    Rank2Field<Scalar> symmetric_einstein;

    ScalarField<Scalar> trace;   // G^mu_mu, which must equal -R in 4 dimensions
    ScalarField<Scalar> square;  // G_{mu nu} G^{mu nu}

    // nabla_mu G^{mu nu}. Zero in the continuum by the contracted Bianchi
    // identity; O(h^2) here, and NOT an exact finite-difference identity.
    Rank1Field<Scalar> divergence;

    // T^inferred = (1/8pi) * (G [+ Lambda g, per policy]), geometrized units.
    Rank2Field<Scalar> inferred_stress_energy;

    double antisymmetry_residual = 0.0;   // max |G_{mu nu} - G_{nu mu}|
    double divergence_residual = 0.0;     // max |nabla_mu G^{mu nu}|
    double divergence_relative = 0.0;     // the same, scaled by max |G|
    double einstein_scale = 0.0;          // max |G_{mu nu}| over the grid
    double trace_identity_residual = 0.0; // max |G^mu_mu + R|

    CosmologicalConstantPolicy lambda_policy;

    bool ok = false;
    std::string refusal_reason;
};

// G_{mu nu} = R_{mu nu} - (1/2) g_{mu nu} R, pointwise.
template <typename Scalar>
Mat4<Scalar> einstein_at(const Mat4<Scalar>& g, const Mat4<Scalar>& ricci,
                         const Scalar& ricci_scalar) {
    const Scalar half = scalar_cast<Scalar>(0.5);
    Mat4<Scalar> out;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            out(mu, nu) = ricci(mu, nu) - half * g(mu, nu) * ricci_scalar;
    return out;
}

namespace detail {

// G^{mu nu} = g^{mu a} G_{a b} g^{b nu}.
template <typename Scalar>
Rank2Field<Scalar> raise_both_field(const Grid& grid,
                                    const Rank2Field<Scalar>& g_inv,
                                    const Rank2Field<Scalar>& t_down) {
    Rank2Field<Scalar> out(grid);
    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        const Mat4<Scalar> gi = g_inv.matrix_at(p);
        out.set_matrix(p, mat_mul(mat_mul(gi, t_down.matrix_at(p)), gi));
    }
    return out;
}

}  // namespace detail

// Assemble G, its trace, square and inferred stress-energy from an already
// computed curvature. Divergence is filled by the driver below, which needs
// more than one time slice.
template <typename Scalar>
EinsteinFields<Scalar> assemble_einstein(const Grid& grid,
                                         const CurvatureFields<Scalar>& c,
                                         const CosmologicalConstantPolicy& lp) {
    EinsteinFields<Scalar> e;
    e.lambda_policy = lp;

    if (!c.ok) {
        e.ok = false;
        e.refusal_reason = "curvature was refused: " + c.refusal_reason;
        return e;
    }

    e.raw_einstein = Rank2Field<Scalar>(grid);
    e.symmetric_einstein = Rank2Field<Scalar>(grid);
    e.trace = ScalarField<Scalar>(grid);
    e.square = ScalarField<Scalar>(grid);
    e.inferred_stress_energy = Rank2Field<Scalar>(grid);

    const Scalar half = scalar_cast<Scalar>(0.5);
    const Scalar inv8pi = scalar_cast<Scalar>(1.0 / (8.0 * 3.14159265358979323846));
    const Scalar lambda = scalar_cast<Scalar>(lp.lambda);

    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        const Mat4<Scalar> gm = c.g.matrix_at(p);
        const Mat4<Scalar> gi = c.g_inv.matrix_at(p);
        const Mat4<Scalar> G = einstein_at(gm, c.ricci.matrix_at(p),
                                           c.ricci_scalar.at(p));
        e.raw_einstein.set_matrix(p, G);

        Mat4<Scalar> Gs;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                Gs(mu, nu) = half * (G(mu, nu) + G(nu, mu));
        e.symmetric_einstein.set_matrix(p, Gs);

        const double asym = symmetry_residual(G);
        if (asym > e.antisymmetry_residual) e.antisymmetry_residual = asym;

        const double scale = max_abs_real(G);
        if (scale > e.einstein_scale) e.einstein_scale = scale;

        // G^mu_mu = g^{mu nu} G_{mu nu}. The symmetric projection is used, so
        // the identity is tested against the tensor that actually goes
        // downstream rather than against one that carries a spurious
        // antisymmetric part.
        const Scalar tr = full_contract(gi, Gs);
        e.trace.at(p) = tr;
        e.trace_identity_residual = std::fmax(
            e.trace_identity_residual,
            std::fabs(real_part(tr) + real_part(c.ricci_scalar.at(p))));

        const Mat4<Scalar> G_up = mat_mul(mat_mul(gi, Gs), gi);
        e.square.at(p) = full_contract(G_up, Gs);

        // The inferred source. What enters depends entirely on the declared
        // policy, which is why the policy travels with the result.
        Mat4<Scalar> source = Gs;
        if (lp.treatment ==
            CosmologicalConstantTreatment::LambdaAsEffectiveStressEnergy) {
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    source(mu, nu) = source(mu, nu) + lambda * gm(mu, nu);
        }
        e.inferred_stress_energy.set_matrix(p, mat_scale(source, inv8pi));
    }

    e.ok = true;
    return e;
}

// G_{mu nu} + Lambda g_{mu nu}, the geometric-side reading. Separate function
// rather than a flag on the assembly, because it answers a different question:
// "is this spacetime a vacuum solution WITH a cosmological constant", not "what
// source does it contain".
template <typename Scalar>
Rank2Field<Scalar> einstein_plus_lambda_g(const Grid& grid,
                                          const EinsteinFields<Scalar>& e,
                                          const Rank2Field<Scalar>& g,
                                          double lambda) {
    Rank2Field<Scalar> out(grid);
    const Scalar lam = scalar_cast<Scalar>(lambda);
    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        Mat4<Scalar> m = e.symmetric_einstein.matrix_at(p);
        const Mat4<Scalar> gm = g.matrix_at(p);
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                m(mu, nu) = m(mu, nu) + lam * gm(mu, nu);
        out.set_matrix(p, m);
    }
    return out;
}

// --- The driver --------------------------------------------------------------
//
// Computes curvature, assembles G, and forms the contracted Bianchi residual
//
//     nabla_mu G^{mu nu} = d_mu G^{mu nu}
//                        + Gamma^mu_{mu lambda} G^{lambda nu}
//                        + Gamma^nu_{mu lambda} G^{mu lambda}
//
// The temporal slot of d_mu G^{mu nu} is the expensive part: it needs G on
// neighbouring time slices, so an AnalyticResample policy triggers TWO extra
// full curvature computations. That cost is why the Bianchi check runs on
// modest grids, and it is paid rather than avoided because the alternative --
// assuming d_0 G = 0 -- is precisely the frozen-time error A4 refuses.
template <typename Scalar, typename SampleMetric>
EinsteinFields<Scalar> compute_einstein(const Grid& grid, SampleMetric&& sample,
                                        double t, const TemporalPolicy& policy,
                                        const CosmologicalConstantPolicy& lp) {
    const CurvatureFields<Scalar> c =
        compute_curvature<Scalar>(grid, sample, t, policy);

    EinsteinFields<Scalar> e = assemble_einstein(grid, c, lp);
    if (!e.ok) return e;

    const Rank2Field<Scalar> G_up =
        detail::raise_both_field(grid, c.g_inv, e.symmetric_einstein);

    // d_0 G^{mu nu}, by the same policy that governed the metric.
    Rank3Field<Scalar> d0;
    bool have_d0 = false;
    if (policy.kind == TemporalPolicyKind::AnalyticResample) {
        const double h = policy.h;
        // Each neighbouring slice is computed and reduced to G^{mu nu} inside
        // its own scope, so only ONE extra CurvatureFields is alive at a time.
        // Holding both would put three rank-4 fields in memory simultaneously --
        // 580 MB at 33^3 before anything else.
        Rank2Field<Scalar> up_minus, up_plus;
        bool okm = false, okp = false;
        {
            const CurvatureFields<Scalar> cm =
                compute_curvature<Scalar>(grid, sample, t - h, policy);
            if (cm.ok) {
                const EinsteinFields<Scalar> em = assemble_einstein(grid, cm, lp);
                up_minus = detail::raise_both_field(grid, cm.g_inv,
                                                    em.symmetric_einstein);
                okm = true;
            }
        }
        {
            const CurvatureFields<Scalar> cp =
                compute_curvature<Scalar>(grid, sample, t + h, policy);
            if (cp.ok) {
                const EinsteinFields<Scalar> ep = assemble_einstein(grid, cp, lp);
                up_plus = detail::raise_both_field(grid, cp.g_inv,
                                                   ep.symmetric_einstein);
                okp = true;
            }
        }
        if (okm && okp) {
            d0 = detail::centred_time_slot(grid, up_minus, up_plus, h);
            have_d0 = true;
        }
    } else {
        // ExactlyStatic: the metric does not mention t, and neither does
        // anything built pointwise from it. Already verified by the staticity
        // probe inside compute_curvature.
        d0 = Rank3Field<Scalar>(grid);
        have_d0 = true;
    }

    const Rank3Field<Scalar> dG =
        spacetime_gradient<Scalar>(G_up, grid, have_d0 ? &d0 : nullptr);

    e.divergence = Rank1Field<Scalar>(grid);
    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        const Rank3<Scalar> gam = c.gamma.tensor_at(p);
        const Mat4<Scalar> Gu = G_up.matrix_at(p);
        for (int nu = 0; nu < 4; ++nu) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int mu = 0; mu < 4; ++mu) acc = acc + dG.at(p, mu, mu, nu);
            for (int mu = 0; mu < 4; ++mu)
                for (int lam = 0; lam < 4; ++lam)
                    acc = acc + gam(mu, mu, lam) * Gu(lam, nu) +
                          gam(nu, mu, lam) * Gu(mu, lam);
            e.divergence.at(p, nu) = acc;
        }
    }

    // Reported over the interior only. The divergence costs a THIRD nested
    // derivative pass, so the one-sided shell is contaminated three deep and a
    // whole-grid maximum would report the boundary treatment rather than the
    // identity.
    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        if (is_shell_point(grid, p, 3)) continue;
        for (int nu = 0; nu < 4; ++nu) {
            const double a = real_abs(e.divergence.at(p, nu));
            if (a > e.divergence_residual) e.divergence_residual = a;
        }
    }
    e.divergence_relative =
        e.einstein_scale > 0.0 ? e.divergence_residual / e.einstein_scale : 0.0;

    return e;
}

}  // namespace geometry
}  // namespace warp_hypad

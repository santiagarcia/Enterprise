#include "einstein_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "warp_hypad/geometry/einstein.hpp"
#include "warp_hypad/metrics/analytic_metrics.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// A5C: cosmological-constant source semantics.
//
// A5 shipped an enum member named `LambdaAsEffectiveStressEnergy` that computed
// (G + Lambda g)/(8 pi). That is the MATTER source. The effective
// cosmological-constant source is -Lambda g/(8 pi), which has the opposite sign
// and, for de Sitter, is the only nonzero one of the two.
//
// The bug was not a typo. A single policy-selected "the source" field cannot
// express the decomposition at all, because there are three distinct tensors:
//
//     T^matter = (G + Lambda g)/(8 pi)      what non-Lambda matter must supply
//     T^Lambda = -Lambda g/(8 pi)           what Lambda contributes as a source
//     T^total  = G/(8 pi)                   their sum, always
//
// This group verifies the four identities that define that decomposition, plus
// the two limits where the distinction is sharpest: de Sitter, where T^matter
// vanishes and T^Lambda does not, and Lambda = 0, where the reverse holds.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using geometry::CosmologicalConstantPolicy;

constexpr double kPi = 3.14159265358979323846;
constexpr double kThetaMin = 1.0;
constexpr double kThetaMax = kPi - 1.0;

// The same grid the A5 de Sitter group uses. The four algebraic identities
// below would hold on any grid -- they are exact -- but the two cases that
// assert T^matter VANISHES are truncation statements about G, so they must be
// evaluated where A5 already measured that truncation rather than on a coarser
// grid whose larger error would then need a looser tolerance.
Grid desitter_grid() {
    GridSpec s;
    s.nx = 33; s.x_min = 0.5;        s.x_max = 1.5;
    s.ny = 33; s.y_min = kThetaMin;  s.y_max = kThetaMax;
    s.nz = 5;  s.z_min = 0.0;        s.z_max = 1.0;
    return Grid(s);
}

Grid alcubierre_grid() {
    GridSpec s;
    s.nx = s.ny = s.nz = 13;
    s.x_min = s.y_min = s.z_min = -4.0;
    s.x_max = s.y_max = s.z_max = 4.0;
    return Grid(s);
}

template <typename MetricT>
auto sampler_for(const Grid& grid, const MetricT& m,
                 const typename MetricT::Parameters& p) {
    return [&grid, &m, &p](double t) {
        return metrics::build_metric_field<double>(grid, m, p, t);
    };
}

geometry::TemporalPolicy static_policy(const char* which) {
    return geometry::TemporalPolicy::exactly_static(
        std::string(which) + " has no t dependence in this chart", 1.0);
}

// Worst componentwise difference between two rank-2 fields over the grid.
double field_diff(const Grid& g, const Rank2Field<double>& a,
                  const Rank2Field<double>& b) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p)
        worst = std::fmax(worst, max_abs_diff(a.matrix_at(p), b.matrix_at(p)));
    return worst;
}

double field_max(const Grid& g, const Rank2Field<double>& a) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p)
        worst = std::fmax(worst, max_abs_real(a.matrix_at(p)));
    return worst;
}

// Interior only, shell depth 2, as everywhere else that a truncation-limited
// quantity is bounded. The one-sided shell is an order of magnitude worse and
// reporting it here would say more about the boundary stencil than about the
// decomposition.
double field_max_interior(const Grid& g, const Rank2Field<double>& a) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        if (geometry::is_shell_point(g, p, 2)) continue;
        worst = std::fmax(worst, max_abs_real(a.matrix_at(p)));
    }
    return worst;
}

double field_diff_interior(const Grid& g, const Rank2Field<double>& a,
                           const Rank2Field<double>& b) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        if (geometry::is_shell_point(g, p, 2)) continue;
        worst = std::fmax(worst, max_abs_diff(a.matrix_at(p), b.matrix_at(p)));
    }
    return worst;
}

}  // namespace

void run_lambda_source_checks(VerificationReport& report,
                              const std::string& out_dir) {
    Stopwatch sw;

    // --- The four defining identities, on a NONZERO Lambda ------------------
    //
    // Checked on de Sitter with the physically correct Lambda, because that is
    // where every term is nonzero and a sign error cannot hide.
    const double L = 2.0;
    const double lambda = 3.0 / (L * L);

    const Grid g = desitter_grid();
    metrics::DeSitterStaticMetric<double> m;
    metrics::DeSitterParameters<double> p(L);

    sw.reset();
    const auto c = geometry::compute_curvature<double>(
        g, sampler_for(g, m, p), 0.0, static_policy("static de Sitter"));
    auto e = geometry::assemble_einstein(
        g, c, CosmologicalConstantPolicy::with_lambda(
                  lambda, "de Sitter with L = 2, Lambda = 3/L^2"));

    const double inv8pi = 1.0 / (8.0 * kPi);
    const std::string params = "de Sitter L=2, Lambda=0.75";

    // 8 pi T^matter = G + Lambda g
    sw.reset();
    double worst = 0.0;
    for (std::size_t q = 0; q < g.num_points(); ++q) {
        const Mat4<double> Tm = e.source.matter_stress_energy.matrix_at(q);
        const Mat4<double> G = e.symmetric_einstein.matrix_at(q);
        const Mat4<double> gm = c.g.matrix_at(q);
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                worst = std::fmax(worst,
                                  std::fabs(8.0 * kPi * Tm(a, b) -
                                            (G(a, b) + lambda * gm(a, b))));
    }
    report.check_residual("lambda.matter_source_identity",
                          "8 pi T^matter = G + Lambda g",
                          VerificationLevel::A, "33x33x5", params, worst,
                          1.0e-14, sw.seconds(),
                          "the identity A5 named correctly and computed as if "
                          "it were the Lambda source");

    // 8 pi T^Lambda = -Lambda g
    sw.reset();
    worst = 0.0;
    for (std::size_t q = 0; q < g.num_points(); ++q) {
        const Mat4<double> Tl =
            e.source.lambda_effective_stress_energy.matrix_at(q);
        const Mat4<double> gm = c.g.matrix_at(q);
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                worst = std::fmax(
                    worst, std::fabs(8.0 * kPi * Tl(a, b) + lambda * gm(a, b)));
    }
    report.check_residual("lambda.effective_lambda_source_identity",
                          "8 pi T^Lambda = -Lambda g",
                          VerificationLevel::A, "33x33x5", params, worst,
                          1.0e-14, sw.seconds(),
                          "NOTE THE SIGN. This is the tensor A5's enum claimed "
                          "to produce; it differs from what A5 computed by more "
                          "than a sign, since that was G + Lambda g");

    // 8 pi T^total = G
    sw.reset();
    worst = 0.0;
    for (std::size_t q = 0; q < g.num_points(); ++q) {
        const Mat4<double> Tt =
            e.source.total_effective_stress_energy.matrix_at(q);
        const Mat4<double> G = e.symmetric_einstein.matrix_at(q);
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                worst = std::fmax(worst,
                                  std::fabs(8.0 * kPi * Tt(a, b) - G(a, b)));
    }
    report.check_residual("lambda.total_source_identity",
                          "8 pi T^total = G", VerificationLevel::A, "17x17x5",
                          params, worst, 1.0e-14, sw.seconds(),
                          "independent of Lambda, as it must be: moving Lambda "
                          "between the two sides of the field equation is "
                          "bookkeeping, not physics");

    // T^total = T^matter + T^Lambda
    sw.reset();
    worst = 0.0;
    for (std::size_t q = 0; q < g.num_points(); ++q) {
        const Mat4<double> Tt =
            e.source.total_effective_stress_energy.matrix_at(q);
        const Mat4<double> Tm = e.source.matter_stress_energy.matrix_at(q);
        const Mat4<double> Tl =
            e.source.lambda_effective_stress_energy.matrix_at(q);
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                worst = std::fmax(worst,
                                  std::fabs(Tt(a, b) - Tm(a, b) - Tl(a, b)));
    }
    report.check_residual("lambda.decomposition_closure",
                          "T^total = T^matter + T^Lambda",
                          VerificationLevel::A, "33x33x5", params, worst,
                          1.0e-15, sw.seconds(),
                          "the closure that makes the decomposition a "
                          "decomposition rather than three unrelated fields");

    // --- de Sitter: the case that separates matter from Lambda ---------------
    sw.reset();
    const double matter_max = field_max_interior(g, e.source.matter_stress_energy);
    report.check_residual(
        "lambda.desitter_matter_source_vanishes",
        "T^matter = 0 for de Sitter with Lambda = 3/L^2: it is a VACUUM",
        VerificationLevel::B, "33x33x5", params, matter_max, 1.5e-3,
        sw.seconds(),
        "the physical content of the correction. A5's mislabelled field "
        "computed exactly this quantity while calling it the Lambda source, so "
        "a reader would have concluded de Sitter's Lambda source was zero");

    sw.reset();
    const double lambda_max =
        field_max(g, e.source.lambda_effective_stress_energy);
    report.check_predicate(
        "lambda.desitter_lambda_source_nonzero",
        "T^Lambda != 0 for de Sitter",
        VerificationLevel::A, lambda_max > 1.0e-3, "max |T^Lambda| > 1e-3",
        "max |T^Lambda| = " + std::to_string(lambda_max), sw.seconds(),
        "expected Lambda/(8 pi) times the metric scale, i.e. about " +
            std::to_string(lambda * inv8pi));

    sw.reset();
    const double tot_vs_lam = field_diff_interior(
        g, e.source.total_effective_stress_energy,
        e.source.lambda_effective_stress_energy);
    report.check_residual(
        "lambda.desitter_total_equals_lambda_source",
        "T^total = T^Lambda for de Sitter, because T^matter vanishes",
        VerificationLevel::B, "33x33x5", params, tot_vs_lam, 1.5e-3,
        sw.seconds(), "");

    // --- Lambda = 0: the reverse limit --------------------------------------
    {
        const Grid ga = alcubierre_grid();
        metrics::AlcubierreMetric<double> am(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> ap(0.5, 2.0, 1.0);
        sw.reset();
        const auto ca = geometry::compute_curvature<double>(
            ga, sampler_for(ga, am, ap), 0.0,
            geometry::TemporalPolicy::analytic_resample(ga.dx() / 0.5,
                                                        "stencil aligned"));
        auto ea = geometry::assemble_einstein(ga, ca,
                                              CosmologicalConstantPolicy::none());

        const double lam_max =
            field_max(ga, ea.source.lambda_effective_stress_energy);
        report.check_exact(
            "lambda.zero_lambda_source_vanishes",
            "T^Lambda = 0 exactly when Lambda = 0",
            VerificationLevel::A, "13x13x13", 0.0, lam_max, sw.seconds(),
            "multiplication by a literal zero, so exact rather than small");

        sw.reset();
        const double m_vs_t = field_diff(ga, ea.source.matter_stress_energy,
                                         ea.source.total_effective_stress_energy);
        report.check_exact(
            "lambda.zero_lambda_matter_equals_total",
            "T^matter = T^total = G/(8 pi) when Lambda = 0",
            VerificationLevel::A, "13x13x13", 0.0, m_vs_t, sw.seconds(),
            "the warp-drive case: with no cosmological constant the "
            "decomposition collapses to the single tensor A5 reported, which is "
            "why the A5 numerical results are unaffected by this correction");
    }

    // --- The labels match the arithmetic ------------------------------------
    sw.reset();
    const std::string conv = e.source.field_equation_convention;
    report.check_predicate(
        "lambda.labels_match_arithmetic",
        "the recorded convention states the equation actually implemented",
        VerificationLevel::A,
        conv.find("G_{mu nu} + Lambda g_{mu nu} = 8 pi T^matter") !=
                std::string::npos &&
            conv.find("T^Lambda = -Lambda g/(8 pi)") != std::string::npos,
        "convention string names the implemented decomposition", conv,
        sw.seconds(),
        "the A5C requirement in one sentence: no name may assert one physical "
        "interpretation while the arithmetic implements another");

    sw.reset();
    report.check_predicate(
        "lambda.no_ambiguous_policy_remains",
        "the public API exposes no single policy-selected 'the source' field",
        VerificationLevel::A,
        e.source.lambda_on_geometric_side && e.source.lambda == lambda,
        "SourceDecomposition with four separately named tensors",
        "lambda = " + std::to_string(e.source.lambda) +
            ", lambda_on_geometric_side = true",
        sw.seconds(),
        "CosmologicalConstantTreatment and its LambdaAsEffectiveStressEnergy "
        "member are removed, not renamed: an enum cannot express a "
        "decomposition, and a caller that could select one tensor would still "
        "have to guess which one it received");

    // --- Determinism ---------------------------------------------------------
    sw.reset();
    auto e2 = geometry::assemble_einstein(
        g, c, CosmologicalConstantPolicy::with_lambda(lambda, "repeat"));
    const double drift =
        field_diff(g, e.source.matter_stress_energy,
                   e2.source.matter_stress_energy) +
        field_diff(g, e.source.lambda_effective_stress_energy,
                   e2.source.lambda_effective_stress_energy);
    report.check_exact("lambda.deterministic_records",
                       "two identical runs give bit-identical decompositions",
                       VerificationLevel::A, "33x33x5", 0.0, drift,
                       sw.seconds(), "");

    // --- Record --------------------------------------------------------------
    std::ofstream f(out_dir + "/lambda_source_decomposition.csv");
    f << "r,G00,T_matter_00,T_lambda_00,T_total_00,"
         "T_matter_plus_lambda_00,closure_residual\n";
    const int jm = g.center_index_y(), km = g.center_index_z();
    for (int i = 0; i < g.nx(); ++i) {
        const std::size_t q = g.flat(i, jm, km);
        const double G00 = e.symmetric_einstein.matrix_at(q)(0, 0);
        const double tm = e.source.matter_stress_energy.matrix_at(q)(0, 0);
        const double tl =
            e.source.lambda_effective_stress_energy.matrix_at(q)(0, 0);
        const double tt =
            e.source.total_effective_stress_energy.matrix_at(q)(0, 0);
        f << g.x(i) << "," << G00 << "," << tm << "," << tl << "," << tt << ","
          << (tm + tl) << "," << std::fabs(tt - tm - tl) << "\n";
    }
}

}  // namespace checks
}  // namespace warp_hypad

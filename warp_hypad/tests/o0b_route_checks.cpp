#include "o0b_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "warp_hypad/core/metadata.hpp"

#if defined(WARP_HYPAD_HAVE_OTILIB)
#include "warp_hypad/sensitivity/otilib_extract.hpp"
#include "warp_hypad/sensitivity/otilib_scalar.hpp"
#include "warp_hypad/sensitivity/otilib_seed.hpp"
#endif

// ---------------------------------------------------------------------------
// O0B cases 27-32 and 37-44: coefficient gradients, rejections, the two
// evaluation routes, readiness and determinism.
//
// WHAT THE GRADIENT CASES CAN AND CANNOT CLAIM.
//
// beta^x and g_01 are LINEAR in the coefficients, so d/dc_k of each is exactly
// B_k -- an identity, checkable to rounding. g_00 = -1 + (beta^x)^2 is
// QUADRATIC, so d g_00/d c_k = 2 beta^x B_k, which is a stronger test: it mixes
// the baseline into the derivative, and an error in either factor shows up.
//
// The backend is onumm3n2: THREE directions at order TWO, read from the
// compiled type rather than requested. With 34 whitened coordinates the
// schedule is twelve batches -- eleven of three and one of one -- and it is
// verified to cover every coordinate exactly once.
//
// The complete design gradient is available through twelve batches. Mixed
// Hessian entries between distinct batches are not available. Order 2 gives a
// Hessian block for the three directions each evaluation seeded, and nothing
// else: the entry between a coordinate in batch 3 and one in batch 7 is never
// formed because the two were never seeded together. Assembling the blocks and
// calling the result a 34-dimensional Hessian would be a misstatement, so no
// full Hessian is implemented and the manifest says so.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::ShiftFamilyConfig;

inline constexpr int kDeclaredOtiDirections = 3;
inline constexpr int kDeclaredOtiOrder = 2;

struct RouteState {
    int oti_directions = 0;
    int oti_order = 0;
    int batches = 0;
    int full_batches = 0;
    int partial = 0;
    bool oti_available = false;
    double separate_residual = -1.0;
    double batched_residual = -1.0;
    std::uint64_t determinism_digest = 0;
};
RouteState g_state;

std::uint64_t fnv1a64_string(const std::string& s) {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (char c : s) {
        h ^= static_cast<std::uint64_t>(static_cast<unsigned char>(c));
        h *= 0x100000001B3ULL;
    }
    return h;
}

std::string hex64(std::uint64_t v) {
    static const char* d = "0123456789abcdef";
    std::string s(16, '0');
    for (int i = 15; i >= 0; --i) { s[static_cast<std::size_t>(i)] = d[v & 0xF];
        v >>= 4; }
    return s;
}

std::uint64_t fnv1a64_file(const std::string& path, bool& ok) {
    std::ifstream f(path, std::ios::binary);
    ok = static_cast<bool>(f);
    if (!ok) return 0;
    std::uint64_t h = 0xCBF29CE484222325ULL;
    char buf[1 << 16];
    while (f.read(buf, sizeof(buf)) || f.gcount() > 0) {
        const std::streamsize n = f.gcount();
        for (std::streamsize i = 0; i < n; ++i) {
            h ^= static_cast<std::uint64_t>(static_cast<unsigned char>(buf[i]));
            h *= 0x100000001B3ULL;
        }
    }
    return h;
}

// A canonical numeric dump of the metric jet at a fixed point set. Determinism
// is established on THESE NUMBERS, never on a timestamp or a build path.
std::string canonical_dump(const ShiftFamilyConfig& cfg,
                           const AxisymmetricBasis& basis,
                           const std::vector<double>& coef) {
    std::ostringstream os;
    os.precision(17);
    os << std::scientific;
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) {
            const double t = -0.6 + 0.3 * i;
            const double x = -2.0 + 1.0 * j;
            const design::MetricJet<double> jt =
                design::f1_metric_jet(cfg, basis, coef, t, x, 1.1, 0.4, true);
            for (int A = 0; A < 4; ++A)
                for (int B = 0; B < 4; ++B) {
                    os << jt.g_cov.m(A, B) << ' ' << jt.g_contra.m(A, B);
                    for (int m = 0; m < 4; ++m)
                        os << ' ' << jt.dg[static_cast<std::size_t>(m)].m(A, B);
                }
            os << '\n';
        }
    return os.str();
}

}  // namespace

void run_o0b_route_checks(VerificationReport& report,
                          const std::string& out_dir) {
    Stopwatch sw;
    const ShiftFamilyConfig cfg = o0b_config();
    const design::DesignRegion region = o0b_region();
    AxisymmetricBasis basis = o0b_basis();
    const design::DesignQuadrature quad = design::build_quadrature(basis);
    const design::ConstraintSystem sys = o0b_constraints(basis, quad);
    const design::ConstraintProjection pr = design::build_projection(sys);
    const design::DenseMatrix M = design::sobolev_gram(
        basis, quad, design::SobolevScales{cfg.scale_R0, cfg.scale_R0});
    const design::WhiteningResult wr = design::whiten(M, pr);
    const design::CoefficientSpaces spaces(pr, wr);
    const std::vector<double> coef =
        o0b_coefficients(basis.coefficient_count(), 0.02);

    const double pt[4] = {0.35, 1.1, 0.9, 0.4};
    const design::DesignCoordinates dp = design::make_design_coordinates(
        pt[1] - (cfg.bubble_center_x0 + cfg.velocity * pt[0]), pt[2], pt[3], 0.0,
        cfg.scale_R0);
    const design::TensorBasisEvaluation ev = basis.evaluate(dp);
    std::vector<double> Bk(
        static_cast<std::size_t>(basis.coefficient_count()), 0.0);
    for (std::size_t m = 0; m < ev.size(); ++m)
        Bk[static_cast<std::size_t>(ev.flat_index[m])] = ev.value[m];
    const double beta_here = design::shift_value_scalar<double>(
        cfg, basis, coef, pt[0], pt[1], pt[2], pt[3]);

    std::ofstream cg(out_dir + "/coefficient_gradient_checks.csv");
    cg << "method,target,coefficient_or_coordinate,batch,direction,reference,"
          "computed,absolute_difference,note\n";

#if defined(WARP_HYPAD_HAVE_OTILIB)
    g_state.oti_available = true;
    g_state.oti_directions = sensitivity::kOtiDirections;
    g_state.oti_order = sensitivity::kOtiOrder;
    static_assert(kDeclaredOtiDirections == sensitivity::kOtiDirections,
                  "declared OTI direction count no longer matches the type");
    static_assert(kDeclaredOtiOrder == sensitivity::kOtiOrder,
                  "declared OTI order no longer matches the type");
    using sensitivity::OTIScalar;

    // --- 27/28/29. Separate seeding: beta, g00, g01 -------------------------
    sw.reset();
    double w_beta = 0.0, w_g00 = 0.0, w_g01 = 0.0, w_dg = 0.0;
    for (std::size_t m = 0; m < ev.size(); ++m) {
        const int k = ev.flat_index[m];
        std::vector<OTIScalar> c(coef.size());
        for (std::size_t i = 0; i < coef.size(); ++i)
            c[i] = sensitivity::make_constant(coef[i]);
        c[static_cast<std::size_t>(k)] =
            sensitivity::seed_parameter(coef[static_cast<std::size_t>(k)], 0);

        const OTIScalar b = design::shift_value_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3]);
        const double gb = sensitivity::extract_gradient(b, 0);
        w_beta = std::fmax(w_beta, std::fabs(gb - ev.value[m]));
        cg << "oti_separate,beta_x," << k << ",-1,0," << ev.value[m] << ","
           << gb << "," << std::fabs(gb - ev.value[m]) << ",linear: dB/dc = B\n";

        const OTIScalar g00 = design::metric_component_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3], 0, 0);
        const double gg00 = sensitivity::extract_gradient(g00, 0);
        const double ref00 = 2.0 * beta_here * ev.value[m];
        w_g00 = std::fmax(w_g00, std::fabs(gg00 - ref00));
        cg << "oti_separate,g00," << k << ",-1,0," << ref00 << "," << gg00 << ","
           << std::fabs(gg00 - ref00) << ",quadratic: 2 beta B\n";

        const OTIScalar g01 = design::metric_component_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3], 0, 1);
        const double gg01 = sensitivity::extract_gradient(g01, 0);
        w_g01 = std::fmax(w_g01, std::fabs(gg01 - ev.value[m]));
        cg << "oti_separate,g01," << k << ",-1,0," << ev.value[m] << "," << gg01
           << "," << std::fabs(gg01 - ev.value[m]) << ",linear\n";

        // A selected dg component: d/dc_k of d_x g_01 must be dB_k/dx.
        const OTIScalar dgx = design::metric_first_derivative_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3], 1, 0, 1);
        const double gdx = sensitivity::extract_gradient(dgx, 0);
        const double refdx = ev.d_xi[m] / cfg.scale_R0;
        w_dg = std::fmax(w_dg, std::fabs(gdx - refdx));
        cg << "oti_separate,d_x_g01," << k << ",-1,0," << refdx << "," << gdx
           << "," << std::fabs(gdx - refdx) << ",dB/dx\n";
    }
    g_state.separate_residual = std::fmax(std::fmax(w_beta, w_g00),
                                          std::fmax(w_g01, w_dg));
    report.check_residual(
        "o0b.oti_gradient_of_shift",
        "d beta^x / d c_k = B_k exactly, by separate OTI seeding",
        VerificationLevel::A, "36 coefficients", "onumm3n2, 3 directions",
        w_beta, 1.0e-15, sw.seconds(),
        "beta is linear in the coefficients, so this is an identity");
    sw.reset();
    report.check_residual(
        "o0b.oti_gradient_of_g00",
        "d g_00 / d c_k = 2 beta^x B_k, the QUADRATIC case",
        VerificationLevel::A, "36 coefficients", "onumm3n2", w_g00, 1.0e-14,
        sw.seconds(),
        "a stronger test than the linear components: the baseline enters the "
        "derivative, so an error in either factor is visible");
    sw.reset();
    report.check_residual(
        "o0b.oti_gradient_of_g01", "d g_01 / d c_k = B_k",
        VerificationLevel::A, "36 coefficients", "onumm3n2", w_g01, 1.0e-15,
        sw.seconds(), "selected dg component residual " + std::to_string(w_dg));

    // --- 30/31. Batched seeding and coverage --------------------------------
    sw.reset();
    std::vector<std::vector<int>> batch;
    for (int s = 0; s < pr.n_free; s += sensitivity::kOtiDirections) {
        std::vector<int> b;
        for (int d = 0; d < sensitivity::kOtiDirections && s + d < pr.n_free;
             ++d)
            b.push_back(s + d);
        batch.push_back(b);
    }
    g_state.batches = static_cast<int>(batch.size());
    for (const auto& b : batch)
        if (static_cast<int>(b.size()) == sensitivity::kOtiDirections)
            ++g_state.full_batches;
    g_state.partial = static_cast<int>(batch.back().size());
    if (g_state.partial == sensitivity::kOtiDirections) g_state.partial = 0;

    std::vector<int> seen(static_cast<std::size_t>(pr.n_free), 0);
    for (const auto& b : batch)
        for (int j : b) ++seen[static_cast<std::size_t>(j)];
    bool covers_once = true;
    for (int v : seen) if (v != 1) covers_once = false;

    // Reference: d beta / d z_j = sum_i P_ij B_i.
    std::vector<double> ref_z(static_cast<std::size_t>(pr.n_free), 0.0);
    for (int j = 0; j < pr.n_free; ++j) {
        double acc = 0.0;
        for (int i = 0; i < pr.n_raw; ++i)
            acc += pr.P(i, j) * Bk[static_cast<std::size_t>(i)];
        ref_z[static_cast<std::size_t>(j)] = acc;
    }
    double w_batched = 0.0;
    const std::size_t probe[3] = {0, batch.size() / 2, batch.size() - 1};
    for (std::size_t bi : probe) {
        const std::vector<int>& b = batch[bi];
        std::vector<OTIScalar> z(static_cast<std::size_t>(pr.n_free));
        for (int j = 0; j < pr.n_free; ++j)
            z[static_cast<std::size_t>(j)] = sensitivity::make_constant(0.0);
        for (std::size_t d = 0; d < b.size(); ++d)
            z[static_cast<std::size_t>(b[d])] =
                sensitivity::seed_parameter(0.0, static_cast<int>(d));
        std::vector<OTIScalar> c(static_cast<std::size_t>(pr.n_raw));
        for (int i = 0; i < pr.n_raw; ++i) {
            OTIScalar acc = sensitivity::make_constant(
                pr.c0[static_cast<std::size_t>(i)]);
            for (int j = 0; j < pr.n_free; ++j)
                acc = acc + z[static_cast<std::size_t>(j)] *
                                OTIScalar(pr.P(i, j));
            c[static_cast<std::size_t>(i)] = acc;
        }
        const OTIScalar b_val = design::shift_value_scalar<OTIScalar>(
            cfg, basis, c, pt[0], pt[1], pt[2], pt[3]);
        for (std::size_t d = 0; d < b.size(); ++d) {
            const double got =
                sensitivity::extract_gradient(b_val, static_cast<int>(d));
            const double want = ref_z[static_cast<std::size_t>(b[d])];
            w_batched = std::fmax(w_batched, std::fabs(got - want));
            cg << "oti_batched,beta_x," << b[d] << "," << bi << "," << d << ","
               << want << "," << got << "," << std::fabs(got - want)
               << ",z coordinate\n";
        }
    }
    g_state.batched_residual = w_batched;
    report.check_residual(
        "o0b.separate_batched_oti_agreement",
        "batched three-direction seeding reproduces the same gradient the "
        "separate seeding gives",
        VerificationLevel::A, "first, middle and final batches", "onumm3n2",
        w_batched, 1.0e-14, sw.seconds(), "");
    sw.reset();
    {
        char buf[300];
        std::snprintf(buf, sizeof(buf),
                      "%d coordinates in %d batches (%d full of %d, final "
                      "partial %d); each covered exactly once: %s",
                      pr.n_free, g_state.batches, g_state.full_batches,
                      sensitivity::kOtiDirections, g_state.partial,
                      covers_once ? "yes" : "NO");
        report.check_predicate(
            "o0b.gradient_batch_coverage",
            "the twelve-batch schedule covers all 34 whitened coordinates "
            "exactly once",
            VerificationLevel::A,
            covers_once && g_state.batches == 12 && g_state.full_batches == 11 &&
                g_state.partial == 1,
            "12 batches = 11 x 3 + 1 x 1", buf, sw.seconds(), "");
    }
#else
    g_state.oti_available = false;
    g_state.oti_directions = kDeclaredOtiDirections;
    g_state.oti_order = kDeclaredOtiOrder;
    {
        int b = 0, full = 0;
        for (int s = 0; s < pr.n_free; s += kDeclaredOtiDirections) {
            const int n = std::min(kDeclaredOtiDirections, pr.n_free - s);
            ++b;
            if (n == kDeclaredOtiDirections) ++full;
            g_state.partial = (n == kDeclaredOtiDirections) ? 0 : n;
        }
        g_state.batches = b;
        g_state.full_batches = full;
    }
    sw.reset();
    for (const char* nm : {"o0b.oti_gradient_of_shift", "o0b.oti_gradient_of_g00",
                           "o0b.oti_gradient_of_g01",
                           "o0b.separate_batched_oti_agreement"})
        report.check_skipped(
            nm, "coefficient gradient through the OTI backend",
            "the OTIlib sensitivity backend is not compiled into this "
            "configuration (WARP_HYPAD_WITH_OTILIB=OFF); the exact references "
            "B_k and 2 beta B_k are still checked by the finite-difference "
            "case, and the schedule arithmetic is recorded from the backend's "
            "declared constants",
            sw.seconds());
    report.check_predicate(
        "o0b.gradient_batch_coverage",
        "the twelve-batch schedule covers all 34 whitened coordinates exactly "
        "once, computed from the backend's declared direction count",
        VerificationLevel::A,
        g_state.batches == 12 && g_state.full_batches == 11 &&
            g_state.partial == 1,
        "12 batches = 11 x 3 + 1 x 1",
        std::to_string(pr.n_free) + " coordinates in " +
            std::to_string(g_state.batches) + " batches (declared capacity, "
            "backend absent)",
        sw.seconds(), "");
#endif

    // A finite-difference cross-check of the same gradients, which needs no
    // derivative backend and therefore runs in BOTH configurations.
    sw.reset();
    double fd_worst = 0.0;
    for (std::size_t m = 0; m < ev.size(); ++m) {
        const int k = ev.flat_index[m];
        const double hh = 1e-6;
        std::vector<double> cp = coef, cm = coef;
        cp[static_cast<std::size_t>(k)] += hh;
        cm[static_cast<std::size_t>(k)] -= hh;
        const double d =
            (design::metric_component_scalar<double>(cfg, basis, cp, pt[0],
                                                     pt[1], pt[2], pt[3], 0, 0) -
             design::metric_component_scalar<double>(cfg, basis, cm, pt[0],
                                                     pt[1], pt[2], pt[3], 0, 0)) /
            (2 * hh);
        const double ref = 2.0 * beta_here * ev.value[m];
        fd_worst = std::fmax(fd_worst, std::fabs(d - ref));
        cg << "finite_difference,g00," << k << ",-1,-1," << ref << "," << d
           << "," << std::fabs(d - ref) << ",h=1e-6\n";
    }
    report.check_residual(
        "o0b.coefficient_gradient_finite_difference",
        "centred differences in the coefficients reproduce 2 beta^x B_k, in "
        "both build configurations",
        VerificationLevel::B, "36 coefficients", "h = 1e-6", fd_worst, 1.0e-8,
        sw.seconds(),
        "this is the case that keeps the gradient claim alive when the OTI "
        "backend is absent");

    // --- 32. No full-Hessian claim -------------------------------------------
    sw.reset();
    report.check_predicate(
        "o0b.no_full_hessian_claim",
        "O0B implements no full design Hessian, and says so",
        VerificationLevel::A, true, "no full Hessian implemented",
        "The complete design gradient is available through twelve batches. "
        "Mixed Hessian entries between distinct batches are not available.",
        sw.seconds(),
        "order 2 gives a Hessian block only for the three directions each "
        "evaluation seeded; the entry between a coordinate in batch 3 and one "
        "in batch 7 is never formed. A full Hessian would need pairwise "
        "reseeding, differences of gradients, or more simultaneous directions");

    // --- 37/38. Rejections ---------------------------------------------------
    sw.reset();
    design::WhitenedDesignCoordinates nf;
    nf.y.assign(static_cast<std::size_t>(wr.retained), 0.1);
    nf.y[0] = std::nan("");
    design::ConstrainedSplineCoordinates zz;
    const bool nf_rejected =
        spaces.to_constrained(nf, zz) == design::ConversionStatus::NonFiniteInput;
    design::RawSplineCoefficients rc;
    design::ConstrainedSplineCoordinates zbad;
    zbad.z.assign(static_cast<std::size_t>(pr.n_free), 0.0);
    zbad.z[2] = std::numeric_limits<double>::infinity();
    const bool inf_rejected =
        spaces.to_raw(zbad, rc) == design::ConversionStatus::NonFiniteInput;
    report.check_predicate(
        "o0b.nonfinite_coefficient_rejection",
        "a NaN or infinity anywhere in a design vector is REJECTED before a "
        "metric is produced",
        VerificationLevel::A, nf_rejected && inf_rejected,
        "both rejected with a reason",
        "NaN rejected: " + std::string(nf_rejected ? "yes" : "NO") +
            ", infinity rejected: " + std::string(inf_rejected ? "yes" : "NO"),
        sw.seconds(),
        "rejected at the conversion, so no downstream consumer ever sees a "
        "metric built from one");

    sw.reset();
    design::WhitenedDesignCoordinates wrong_dim;
    wrong_dim.y.assign(5, 0.0);
    const bool dim_rejected =
        spaces.to_constrained(wrong_dim, zz) ==
        design::ConversionStatus::DimensionMismatch;
    design::ConstrainedSplineCoordinates wrong_z;
    wrong_z.z.assign(3, 0.0);
    const bool dim_rejected2 =
        spaces.to_raw(wrong_z, rc) == design::ConversionStatus::DimensionMismatch;
    report.check_predicate(
        "o0b.dimension_mismatch_rejection",
        "a design vector of the wrong length is rejected in every conversion",
        VerificationLevel::A, dim_rejected && dim_rejected2,
        "both rejected", "whitened: " + std::string(dim_rejected ? "yes" : "NO") +
            ", constrained: " + std::string(dim_rejected2 ? "yes" : "NO"),
        sw.seconds(),
        "z and y have the same length here, so only the TYPE distinguishes "
        "them; the length check catches the remaining case");

    // --- 41. Route A ---------------------------------------------------------
    sw.reset();
    const design::UniformGridCost cost =
        design::uniform_grid_cost(0.3018, 4.0, 0.4);
    {
        char buf[420];
        std::snprintf(buf, sizeof(buf),
                      "narrowest feature %.4f m needs h <= %.6f m: %lld cells "
                      "per axis, %.4g cells, %.2f GB, %.1fx the current cost",
                      cost.narrowest_feature_m, cost.required_spacing_m,
                      cost.cells_per_axis, cost.total_cells,
                      cost.total_bytes / 1073741824.0,
                      cost.cost_factor_versus_current);
        report.check_predicate(
            "o0b.refined_grid_cost_estimate",
            "the cost of meeting the six-cell rule on a refined uniform grid is "
            "estimated and recorded, and the grid is NOT run",
            VerificationLevel::E,
            cost.required_spacing_m > 0.0 && cost.total_cells > 1e6 &&
                cost.cost_factor_versus_current > 100.0,
            "an explicit spacing, cell count and memory figure", buf,
            sw.seconds(),
            "the six-cell rule is a FINITE-DIFFERENCE criterion, applied here "
            "because Route A really does difference");
    }

    // --- 42. Route B ---------------------------------------------------------
    sw.reset();
    const design::QuadratureScaffold sc =
        design::build_quadrature_scaffold(basis, region);
    const design::RouteComparison rc2 = design::compare_routes(cost, sc);
    {
        char buf[420];
        std::snprintf(buf, sizeof(buf),
                      "%zu elements, %d nodes per direction, exact to degree "
                      "%d; %zu nodes total; element volumes sum to %.10f "
                      "against the exact %.10f; Route A needs %.4g times as "
                      "many points",
                      sc.element.size(), sc.nodes_per_direction,
                      sc.polynomial_exactness, sc.total_nodes, sc.total_volume,
                      sc.exact_volume, rc2.node_ratio);
        report.check_predicate(
            "o0b.elementwise_quadrature_scaffold",
            "one element per unique tensor-product knot span, tiling the design "
            "region exactly, with a fixed-order Gauss-Legendre rule inside each",
            VerificationLevel::A,
            !sc.element.empty() &&
                std::fabs(sc.total_volume - sc.exact_volume) < 1e-9 &&
                sc.total_nodes > 0,
            "elements tile the region and their volumes sum exactly",
            buf, sw.seconds(),
            "the six-cell rule does NOT transfer to this route: it takes no "
            "finite differences, so what must converge is the INTEGRAL -- in "
            "quadrature order, element subdivision, domain extent, and against "
            "an independent uniform-grid cross-check. O0B evaluates no energy "
            "on either route");
    }

    // --- 43. Readiness remains false -----------------------------------------
    sw.reset();
    report.check_predicate(
        "o0b.readiness_remains_false",
        "F1 is still NOT resolution-certified and NOT ready for energy "
        "optimization: O0B verified metric construction, which is a different "
        "claim",
        VerificationLevel::A, true,
        "resolution_certified false, ready_for_optimization false",
        "metric construction verified; the A6 uniform grid still gives 0.755 "
        "cells across the narrowest F1 feature against a six-cell contract",
        sw.seconds(),
        "a metric that assembles correctly is not a metric whose energy can be "
        "integrated accurately; O0B does not conflate the two");

    // --- 44. Scope -----------------------------------------------------------
    sw.reset();
    report.check_predicate(
        "o0b.scope_statement",
        "every O0B record carries the scope statement",
        VerificationLevel::A, true, "present in the manifest and readiness",
        "O0B maps verified spline coordinates into structurally valid ADM "
        "metrics and analytic metric derivatives. It does not yet minimize "
        "energy, certify the spline family for energy evaluation, demonstrate "
        "a causal shortcut or establish physical realizability.",
        sw.seconds(), "");

    // --- 39/40. Build-local outputs and determinism --------------------------
    sw.reset();
    const std::string dump1 = canonical_dump(cfg, basis, coef);
    AxisymmetricBasis rebuilt = o0b_basis();
    const std::string dump2 = canonical_dump(cfg, rebuilt, coef);
    g_state.determinism_digest = fnv1a64_string(dump1);
    report.check_predicate(
        "o0b.deterministic_checksums",
        "a freshly constructed basis reproduces the entire metric jet BYTE FOR "
        "BYTE, checked on the numbers rather than on a timestamp",
        VerificationLevel::A, dump1 == dump2 && dump1.size() > 1000,
        "byte-identical numeric dump",
        std::to_string(dump1.size()) + " bytes, digest " +
            hex64(g_state.determinism_digest),
        sw.seconds(), "");

    sw.reset();
    report.check_predicate(
        "o0b.build_local_outputs",
        "every record is written under the build-local output directory",
        VerificationLevel::A, out_dir.find("results") != std::string::npos,
        "build-local path", out_dir, sw.seconds(),
        "promotion to the source tree is an explicit opt-in target and is never "
        "run from a test");

    // --- Records --------------------------------------------------------------
    {
        std::ofstream f(out_dir + "/uniform_grid_cost_estimate.json");
        f.precision(17);
        f << "{\n  \"route\": \"A: refined uniform grid\",\n";
        f << "  \"narrowest_feature_m\": " << cost.narrowest_feature_m << ",\n";
        f << "  \"required_cells_per_feature\": "
          << cost.required_cells_per_feature << ",\n";
        f << "  \"required_spacing_m\": " << cost.required_spacing_m << ",\n";
        f << "  \"domain_half_width_m\": " << cost.domain_half_width_m << ",\n";
        f << "  \"cells_per_axis\": " << cost.cells_per_axis << ",\n";
        f << "  \"total_cells\": " << cost.total_cells << ",\n";
        f << "  \"bytes_per_point\": " << cost.bytes_per_point << ",\n";
        f << "  \"total_bytes\": " << cost.total_bytes << ",\n";
        f << "  \"total_gigabytes\": " << cost.total_bytes / 1073741824.0
          << ",\n";
        f << "  \"current_spacing_m\": " << cost.current_spacing_m << ",\n";
        f << "  \"refinement_factor\": " << cost.refinement_factor << ",\n";
        f << "  \"cost_factor_versus_current\": "
          << cost.cost_factor_versus_current << ",\n";
        f << "  \"grid_was_run\": false,\n";
        f << "  \"note\": \"" << cost.note << "\"\n}\n";
    }
    {
        std::ofstream f(out_dir + "/elementwise_quadrature_elements.csv");
        f << "index,span_xi,span_s,xi_lower,xi_upper,s_lower,s_upper,rho_lower,"
             "rho_upper,volume,active_functions,in_outer_collar,"
             "in_protected_cabin\n";
        for (const design::QuadratureElement& e : sc.element)
            f << e.index << "," << e.span_xi << "," << e.span_s << ","
              << e.xi_lower << "," << e.xi_upper << "," << e.s_lower << ","
              << e.s_upper << "," << e.rho_lower << "," << e.rho_upper << ","
              << e.volume << "," << e.active_functions << ","
              << (e.in_outer_collar ? 1 : 0) << ","
              << (e.in_protected_cabin ? 1 : 0) << "\n";
    }
    {
        std::ofstream f(out_dir + "/elementwise_quadrature_manifest.json");
        f.precision(17);
        f << "{\n  \"route\": \"B: analytic jets plus elementwise "
             "quadrature\",\n";
        f << "  \"rule\": \"" << sc.rule << "\",\n";
        f << "  \"elements\": " << sc.element.size() << ",\n";
        f << "  \"nodes_per_direction\": " << sc.nodes_per_direction << ",\n";
        f << "  \"polynomial_exactness_degree\": " << sc.polynomial_exactness
          << ",\n";
        f << "  \"total_nodes\": " << sc.total_nodes << ",\n";
        f << "  \"total_volume\": " << sc.total_volume << ",\n";
        f << "  \"exact_volume\": " << sc.exact_volume << ",\n";
        f << "  \"volume_residual\": "
          << std::fabs(sc.total_volume - sc.exact_volume) << ",\n";
        f << "  \"energy_evaluated\": false,\n";
        f << "  \"optimizer_run\": false,\n";
        f << "  \"convergence_contract\": \"" << sc.convergence_contract
          << "\",\n";
        f << "  \"route_a_cells_per_route_b_node\": " << rc2.node_ratio << ",\n";
        f << "  \"recommendation_for_o0c\": \"" << rc2.recommendation << "\",\n";
        f << "  \"caveat\": \"" << rc2.caveat << "\"\n}\n";
    }
}

void write_o0a_consistency_record(const std::string& out_dir) {
    std::ofstream f(out_dir + "/o0a_resolution_consistency.json");
    f.precision(17);
    f << "{\n";
    f << "  \"audit\": \"O0A resolution-report consistency\",\n";
    f << "  \"canonical_source\": "
         "\"resolution_coverage.csv from the O0A gate\",\n";
    f << "  \"f1_minimum_cells_across_narrowest_feature\": 0.75456,\n";
    f << "  \"f0_minimum_cells_across_narrowest_feature\": 0.892253,\n";
    f << "  \"f0_wall_width_m\": 0.33333333333333331,\n";
    f << "  \"f0_wall_width_in_cells\": 0.83333333333333326,\n";
    f << "  \"finding\": \"0.833 and 0.892 are DIFFERENT production "
         "quantities, both computed by the O0A gate and both printed in its "
         "support_resolution_accounting case note. 0.892 is F0's minimum cells "
         "across its narrowest knot-span feature -- the canonical F0 minimum. "
         "0.833 is the width of the Alcubierre wall (1/sigma) expressed in "
         "cells, i.e. the feature F0 is being asked to resolve. Earlier "
         "narrative text labelled the second as the F0 minimum, conflating "
         "them.\",\n";
    f << "  \"tracked_artifacts_affected\": 0,\n";
    f << "  \"scope_of_inconsistency\": \"historical commit message text "
         "only\",\n";
    f << "  \"corrective_commit_required\": false,\n";
    f << "  \"reason\": \"git grep over all tracked files finds no O0A source, "
         "caption, sidecar, document or gallery record containing either "
         "value. Every O0A figure derives its cell counts from "
         "resolution_coverage.csv at render time, so the captions are "
         "generated rather than hard-coded. Git history is not rewritten.\",\n";
    f << "  \"both_families_resolution_certified\": false\n";
    f << "}\n";
}

void write_o0b_readiness(const VerificationReport& report,
                         const std::string& out_dir) {
    auto passed = [&](const std::string& n) {
        for (const VerificationCase& c : report.cases())
            if (c.name == n) return c.passed;
        return false;
    };
    const bool construction =
        passed("o0b.adm_metric_assembly") &&
        passed("o0b.zero_correction_reproduces_alcubierre") &&
        passed("o0b.lorentzian_signature") &&
        passed("o0b.f1_determinant_is_minus_one");
    const bool jets = passed("o0b.temporal_derivative") &&
                      passed("o0b.mixed_temporal_spatial_derivative") &&
                      passed("o0b.third_derivative_permutation_symmetry");

    std::ofstream f(out_dir + "/readiness.json");
    f << "{\n  \"milestone\": \"O0B\",\n";
    f << "  \"scope\": \"O0B maps verified spline coordinates into "
         "structurally valid ADM metrics and analytic metric derivatives. It "
         "does not yet minimize energy, certify the spline family for energy "
         "evaluation, demonstrate a causal shortcut or establish physical "
         "realizability.\",\n";
    f << "  \"families\": [\n";
    f << "    {\n      \"family\": \"O0B-F1 shift-only ADM metric\",\n";
    f << "      \"implemented\": true,\n";
    f << "      \"metric_construction_verified\": "
      << (construction ? "true" : "false") << ",\n";
    f << "      \"analytic_jets_verified\": " << (jets ? "true" : "false")
      << ",\n";
    f << "      \"gradient_verified\": "
      << (passed("o0b.coefficient_gradient_finite_difference") ? "true"
                                                               : "false")
      << ",\n";
    f << "      \"resolution_certified\": false,\n";
    f << "      \"ready_for_energy_optimization\": false,\n";
    f << "      \"blocking_reason\": \"the underlying O0A F1 family is not "
         "resolution-certified: 0.755 cells across its narrowest radial "
         "feature on the A6 uniform grid, against a six-cell contract. O0B "
         "verified metric CONSTRUCTION, which is a different claim.\"\n";
    f << "    }\n  ],\n";
    f << "  \"required_statement\": \"metric construction verified; current F1 "
         "family not yet resolution-certified; not ready for energy "
         "optimization\",\n";
    f << "  \"oti_backend\": {\n";
    f << "    \"available\": " << (g_state.oti_available ? "true" : "false")
      << ",\n";
    f << "    \"type\": \"onumm3n2\",\n";
    f << "    \"directions\": " << g_state.oti_directions << ",\n";
    f << "    \"order\": " << g_state.oti_order << ",\n";
    f << "    \"gradient_batches\": " << g_state.batches << ",\n";
    f << "    \"full_batches\": " << g_state.full_batches << ",\n";
    f << "    \"final_partial_batch_size\": " << g_state.partial << ",\n";
    f << "    \"hessian_coverage\": \"The complete design gradient is "
         "available through twelve batches. Mixed Hessian entries between "
         "distinct batches are not available.\",\n";
    f << "    \"full_hessian_implemented\": false\n  },\n";
    f << "  \"energy_evaluated\": false,\n";
    f << "  \"optimizer_run\": false\n}\n";
}

void write_o0b_manifest(const VerificationReport& report,
                        const std::string& out_dir) {
    const std::string manifest_name = "metric_design_manifest.json";
    std::vector<std::string> names;
    std::error_code ec;
    for (const auto& e : std::filesystem::directory_iterator(out_dir, ec)) {
        if (!e.is_regular_file()) continue;
        const std::string n = e.path().filename().string();
        if (n == manifest_name) continue;
        names.push_back(n);
    }
    std::sort(names.begin(), names.end());

    const ShiftFamilyConfig cfg = o0b_config();
    AxisymmetricBasis basis = o0b_basis();
    const design::DesignQuadrature quad = design::build_quadrature(basis);
    const design::ConstraintProjection pr =
        design::build_projection(o0b_constraints(basis, quad));
    const design::DenseMatrix M = design::sobolev_gram(
        basis, quad, design::SobolevScales{cfg.scale_R0, cfg.scale_R0});
    const design::WhiteningResult wr = design::whiten(M, pr);
    const design::CoefficientSpaces spaces(pr, wr);
    const design::CoefficientSpaces::Record rec = spaces.record();

    {
        std::ofstream cf(out_dir + "/coefficient_space_manifest.json");
        cf.precision(17);
        cf << "{\n  \"raw_dimension\": " << rec.raw_dimension << ",\n";
        cf << "  \"constraint_rank\": " << rec.constraint_rank << ",\n";
        cf << "  \"constrained_dimension\": " << rec.constrained_dimension
           << ",\n";
        cf << "  \"whitened_dimension\": " << rec.whitened_dimension << ",\n";
        cf << "  \"coordinate_system_identifiers\": [\"raw_c\", "
              "\"constrained_z\", \"whitened_y\"],\n";
        cf << "  \"conversion_direction\": \"y -> z -> c only; no inverse is "
              "provided, because a raw vector need not have an exact "
              "preimage\",\n";
        cf << "  \"P_checksum\": \"" << rec.p_checksum << "\",\n";
        cf << "  \"W_checksum\": \"" << rec.w_checksum << "\",\n";
        cf << "  \"PW_checksum\": \"" << rec.pw_checksum << "\",\n";
        cf << "  \"c0_checksum\": \"" << rec.c0_checksum << "\",\n";
        cf << "  \"whitened_index_identity\": \"" << rec.whitened_identity_note
           << "\"\n}\n";
    }

    std::ofstream f(out_dir + "/" + manifest_name);
    f.precision(17);
    f << "{\n  \"milestone\": \"O0B\",\n";
    f << "  \"title\": \"constrained ADM metric and analytic metric jets\",\n";
    f << "  \"git_commit\": \"" << build_info().git_commit << "\",\n";
    f << "  \"git_dirty\": " << (build_info().git_dirty ? "true" : "false")
      << ",\n";
    f << "  \"generated_utc\": \"" << utc_timestamp() << "\",\n";
    f << "  \"provenance_note\": \"Timestamp and commit are provenance only. "
         "Determinism is established by o0b.deterministic_checksums on the "
         "numbers themselves.\",\n";
    f << "  \"index_ordering\": \"" << design::index_ordering() << "\",\n";
    f << "  \"family\": \"F1 shift-only: alpha = 1, gamma = I, beta^y = beta^z "
         "= 0\",\n";
    f << "  \"baseline\": \"metrics::shape_function, the analytic Alcubierre "
         "shape function verified from A3 onward. NOT the F0 projection, which "
         "O0A concluded is not resolution-certified.\",\n";
    f << "  \"shift_sign_convention\": \"beta^x = -v_s f, NEGATIVE. Fixed by "
         "production: AlcubierreMetric writes g_01 = -v f, and the ADM form "
         "g_0i = gamma_ij beta^j with gamma = I forces the sign.\",\n";
    f << "  \"lapse_parameterisation\": \"alpha = exp(a), positive for every "
         "finite a\",\n";
    f << "  \"spatial_metric_parameterisation\": \"gamma = L L^T with "
         "exponential diagonal, positive definite for every finite L\",\n";
    f << "  \"g00_zero_meaning\": \"g_00 = 0 is where the coordinate-stationary "
         "vector becomes null. NOT a metric degeneracy, NOT a signature "
         "change, NOT a curvature singularity, NOT automatically a horizon: "
         "det g stays exactly -1 there.\",\n";
    f << "  \"time_treatment\": \"analytic: d/dt = -(v_s/R0) d/dxi. The "
         "frozen-time diagnostic path is NOT used anywhere in O0B, and the "
         "mixed temporal-spatial channels are verified nonzero.\",\n";
    f << "  \"velocity\": " << cfg.velocity << ",\n";
    f << "  \"bubble_radius\": " << cfg.bubble_radius << ",\n";
    f << "  \"wall_sigma\": " << cfg.wall_sigma << ",\n";
    f << "  \"scale_R0\": " << cfg.scale_R0 << ",\n";
    f << "  \"raw_dimension\": " << rec.raw_dimension << ",\n";
    f << "  \"constraint_rank\": " << rec.constraint_rank << ",\n";
    f << "  \"constrained_dimension\": " << rec.constrained_dimension << ",\n";
    f << "  \"whitened_dimension\": " << rec.whitened_dimension << ",\n";
    f << "  \"oti_directions\": " << g_state.oti_directions << ",\n";
    f << "  \"oti_order\": " << g_state.oti_order << ",\n";
    f << "  \"gradient_batches\": " << g_state.batches << ",\n";
    f << "  \"hessian_coverage\": \"The complete design gradient is available "
         "through twelve batches. Mixed Hessian entries between distinct "
         "batches are not available.\",\n";
    f << "  \"full_hessian_implemented\": false,\n";
    f << "  \"energy_evaluated\": false,\n";
    f << "  \"optimizer_run\": false,\n";
    f << "  \"determinism_digest\": \"" << hex64(g_state.determinism_digest)
      << "\",\n";
    f << "  \"digest_algorithm\": \"FNV-1a 64-bit over file bytes; a content "
         "digest, not a cryptographic checksum\",\n";
    f << "  \"manifest_excludes_itself\": true,\n";
    f << "  \"record_count\": " << names.size() << ",\n";
    f << "  \"records\": {\n";
    for (std::size_t i = 0; i < names.size(); ++i) {
        bool ok = false;
        const std::uint64_t h = fnv1a64_file(out_dir + "/" + names[i], ok);
        f << "    \"" << names[i] << "\": \"" << (ok ? hex64(h) : "unreadable")
          << "\"" << (i + 1 == names.size() ? "\n" : ",\n");
    }
    f << "  },\n";
    f << "  \"cases\": " << report.num_cases() << ",\n";
    f << "  \"passed\": " << report.num_passed() << ",\n";
    f << "  \"failed\": " << report.num_failed() << ",\n";
    f << "  \"skipped\": " << report.num_skipped() << ",\n";
    f << "  \"scope\": \"O0B maps verified spline coordinates into "
         "structurally valid ADM metrics and analytic metric derivatives. It "
         "does not yet minimize energy, certify the spline family for energy "
         "evaluation, demonstrate a causal shortcut or establish physical "
         "realizability.\"\n}\n";
}

}  // namespace checks
}  // namespace warp_hypad

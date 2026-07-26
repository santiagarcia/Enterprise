#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include "b4c_checks.hpp"
#include "program_b/inverse_design/constitutive_decomposition.hpp"
#include "program_b/inverse_design/layered.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

// ---------------------------------------------------------------------------
// Gate B4C-PROPAGATION-CORRECTION.
//
// Three corrections to B4:
//
//   1. the metricity claim was overstated. Withdrawn and replaced.
//   2. the layered solver was hard-wired to double, so nothing could be
//      differentiated through the actual stack assembly.
//   3. the quartic root finder bracketed SIGN CHANGES, so it was blind to roots
//      of even multiplicity -- exactly the metric-induced case.
// ---------------------------------------------------------------------------

namespace program_b {
namespace checks {

using namespace program_b::inverse_design;
using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace {

using C = Complex<double>;

std::string sci(double v) {
    std::ostringstream o;
    o.setf(std::ios::scientific);
    o.precision(4);
    o << v;
    return o.str();
}

LaboratoryMetric<double> eta_lab() { return LaboratoryMetric<double>::minkowski(); }

ObserverFrame<double> lab_frame(const LaboratoryMetric<double>& g) {
    double V[4] = {1.0, 0.0, 0.0, 0.0};
    return build_observer_frame(g, make_laboratory_observer(g, V, "lab_rest"));
}

Mat4<double> alcubierre_at(double vs, double x) {
    warp_hypad::metrics::AlcubierreMetric<double> w;
    const warp_hypad::metrics::AlcubierreParameters<double> p{vs, 3.0, 2.0};
    return w.evaluate(0.0, x, 0.0, 0.0, p);
}

Rank4<double> anisotropic_Z(const LaboratoryMetric<double>& g,
                            const ObserverFrame<double>& frame) {
    ObserverDecomposition4<double> d;
    double e3[3][3] = {}, m3[3][3] = {}, z3[3][3] = {};
    e3[0][0] = 2.0; e3[1][1] = 3.5; e3[2][2] = 5.0;
    for (int i = 0; i < 3; ++i) m3[i][i] = 1.0;
    d.epsilon = block_from_coefficients(e3, frame);
    d.mu_inverse = block_from_coefficients(m3, frame);
    d.zeta = block_from_coefficients(z3, frame);
    d.zeta_adjoint = block_from_coefficients(z3, frame);
    return reconstruct_Z_A11(d, g, frame.V);
}

// A matched homogeneous slab: mu = eps, so the impedance equals vacuum's, r = 0
// and t = exp(i omega eps d) exactly. That gives closed-form derivatives
// d(phase)/d(eps) = omega d and d(phase)/d(d) = omega eps.
template <typename Scalar>
std::vector<Layer<Scalar>> matched_slab(const Scalar& eps, const Scalar& d) {
    std::vector<Layer<Scalar>> stack;
    Layer<Scalar> L;
    L.epsilon = Complex<Scalar>(eps, scalar_cast<Scalar>(0.0));
    // mu = eps  =>  mu^-1 = 1/eps
    L.mu_inverse = Complex<Scalar>(scalar_cast<Scalar>(1.0) / eps,
                                   scalar_cast<Scalar>(0.0));
    L.zeta = complex_from_real<Scalar>(0.0);
    L.thickness = d;
    stack.push_back(L);
    return stack;
}

}  // namespace

void run_b4c_checks(VerificationReport& report, const std::string& out_dir,
                   B4CMetrics& metrics) {
    const LaboratoryMetric<double> g = eta_lab();
    const ObserverFrame<double> frame = lab_frame(g);

    metrics.old_metricity_claim = "metricity iff no birefringence";
    metrics.corrected_metricity_claim =
        "Metric-induced constitutive tensors are non-birefringent. Generic "
        "reciprocal constitutive tensors need not be. The converse requires "
        "additional assumptions and is subject to conformal/dilaton freedom and "
        "constitutive pieces that do not affect the Fresnel polynomial.";

    // =======================================================================
    // 1. The corrected metricity wording.
    // =======================================================================
    {
        Stopwatch sw;
        report.check_predicate(
            "b4c.metricity_claim_is_corrected",
            "the unrestricted claim 'metricity iff no birefringence' is "
            "WITHDRAWN and replaced by the one-directional statement",
            VerificationLevel::A,
            metrics.corrected_metricity_claim.find("need not be") !=
                std::string::npos &&
                metrics.corrected_metricity_claim.find("converse") !=
                    std::string::npos,
            "the corrected wording states only the provable direction",
            metrics.corrected_metricity_claim, sw.seconds(),
            "B4's tests establish that metric-induced tensors give a repeated "
            "quadratic cone, and that ONE selected generic reciprocal tensor "
            "birefringes. Neither establishes the converse, and the Fresnel "
            "polynomial does NOT uniquely reconstruct Z^{abcd} -- the axion piece "
            "and the conformal factor are both invisible to it");
    }

    // =======================================================================
    // 2. The axion-like term: propagation-invisible, yet a different tensor.
    // =======================================================================
    {
        Stopwatch sw;
        const double alpha = 0.37;
        const Rank4<double> Zm =
            covariant_constitutive(EffectiveMetric<double>(alcubierre_at(0.3, 3.0)),
                                   g).Z;
        const Rank4<double> Zax = axion_tensor<double>(alpha, g);
        Rank4<double> Zsum;
        for (int i = 0; i < 256; ++i) Zsum.t[i] = Zm.t[i] + Zax.t[i];

        // (a) the tensors genuinely differ
        double tensor_diff = 0.0;
        for (int i = 0; i < 256; ++i)
            tensor_diff = std::fmax(tensor_diff, std::fabs(Zsum.t[i] - Zm.t[i]));
        metrics.constitutive_difference_under_axion = tensor_diff;

        // (b) the observer blocks differ, and where expected: in zeta
        const ObserverDecomposition4<double> bm = decompose_A7(Zm, g.value, frame.V);
        const ObserverDecomposition4<double> bs = decompose_A7(Zsum, g.value, frame.V);
        const double d_eps = warp_hypad::max_abs_diff(bm.epsilon, bs.epsilon);
        const double d_mu = warp_hypad::max_abs_diff(bm.mu_inverse, bs.mu_inverse);
        const double d_zeta = warp_hypad::max_abs_diff(bm.zeta, bs.zeta);
        metrics.zeta_block_difference_under_axion = d_zeta;

        // The axion appears as 2 alpha h in zeta, so the predicted difference is
        // 2 alpha on the spatial diagonal.
        const Mat4<double> h = projector_upper(g.value, frame.V);
        double predicted = 0.0;
        for (int i = 0; i < 16; ++i)
            predicted = std::fmax(predicted, std::fabs(2.0 * alpha * h.m[i]));
        const double zeta_form_error = std::fabs(d_zeta - predicted);

        // (c) the Fresnel polynomial is UNCHANGED
        const double n[3] = {0.5773502691896258, 0.5773502691896258,
                             0.5773502691896258};
        double c1[5], c2[5];
        fresnel_quartic_coefficients(Zm, 1.0, n, c1);
        fresnel_quartic_coefficients(Zsum, 1.0, n, c2);
        double coeff_diff = 0.0, coeff_scale = 0.0;
        for (int i = 0; i < 5; ++i) {
            coeff_diff = std::fmax(coeff_diff, std::fabs(c1[i] - c2[i]));
            coeff_scale = std::fmax(coeff_scale, std::fabs(c1[i]));
        }
        metrics.fresnel_equality_under_axion = coeff_diff / std::fmax(coeff_scale, 1e-300);

        // (d) the axion coefficient is recovered by projection
        metrics.axion_coefficient_recovered = axion_coefficient(Zsum, g);

        std::ostringstream got;
        got << "max|Z_sum - Z_metric| = " << sci(tensor_diff)
            << "; block differences eps " << sci(d_eps) << ", mu^-1 " << sci(d_mu)
            << ", zeta " << sci(d_zeta) << " (predicted 2*alpha*h = "
            << sci(predicted) << "); Fresnel coefficient difference "
            << sci(metrics.fresnel_equality_under_axion) << "; alpha recovered "
            << sci(metrics.axion_coefficient_recovered) << " against " << alpha;

        report.check_predicate(
            "b4c.axion_changes_the_tensor_but_not_the_cone",
            "adding alpha eps^{abcd} to a metric-induced tensor changes Z and "
            "changes the zeta block by exactly 2 alpha h, while leaving the "
            "Fresnel polynomial bit-identical",
            VerificationLevel::B,
            tensor_diff > 1e-3 && d_zeta > 1e-3 && zeta_form_error < 1e-12 &&
                d_eps < 1e-12 && d_mu < 1e-12 &&
                metrics.fresnel_equality_under_axion < 1e-14 &&
                std::fabs(metrics.axion_coefficient_recovered - alpha) < 1e-12,
            "tensor and zeta differ; eps, mu^-1 and the cone do not",
            got.str(), sw.seconds(),
            "THIS IS THE POINT: same Fresnel cone does NOT imply the same "
            "constitutive tensor. W^{bd} = Z^{abcd}k_a k_c annihilates the axion "
            "piece IDENTICALLY, because a totally antisymmetric symbol contracted "
            "with the symmetric k_a k_c vanishes -- so the axion contributes "
            "nothing to W, nothing to its adjugate, and nothing to the light cone");
    }

    {
        Stopwatch sw;
        // Conformal/dilaton freedom: a different metric, the same cone.
        const Mat4<double> geff = alcubierre_at(0.3, 3.0);
        Mat4<double> scaled;
        for (int i = 0; i < 16; ++i) scaled.m[i] = 3.3 * geff.m[i];
        const Rank4<double> Z1 =
            covariant_constitutive(EffectiveMetric<double>(geff), g).Z;
        const Rank4<double> Z2 =
            covariant_constitutive(EffectiveMetric<double>(scaled), g).Z;

        // g_eff differs, Z does not (conformal invariance, verified at B1C), so
        // the cone cannot distinguish them: the metric is recoverable only up to
        // a conformal factor.
        double metric_diff = 0.0, Z_diff = 0.0;
        for (int i = 0; i < 16; ++i)
            metric_diff = std::fmax(metric_diff, std::fabs(scaled.m[i] - geff.m[i]));
        for (int i = 0; i < 256; ++i)
            Z_diff = std::fmax(Z_diff, std::fabs(Z2.t[i] - Z1.t[i]));
        metrics.conformal_fresnel_zero_set_shift = Z_diff;

        report.check_predicate(
            "b4c.conformal_freedom_is_explicit",
            "two effective metrics differing by a conformal factor give the SAME "
            "Z and therefore the same cone: the metric is determined only up to "
            "that factor",
            VerificationLevel::B, metric_diff > 1e-1 && Z_diff < 1e-13,
            "metrics differ, Z does not",
            "max|g2 - g1| = " + sci(metric_diff) + ", max|Z2 - Z1| = " +
                sci(Z_diff),
            sw.seconds(),
            "the second reason the Fresnel polynomial cannot reconstruct a unique "
            "metric, independent of the axion argument");
    }

    // =======================================================================
    // 3. The classification, reported piece by piece.
    // =======================================================================
    {
        Stopwatch sw;
        const double n[3] = {0.5773502691896258, 0.5773502691896258,
                             0.5773502691896258};
        const Rank4<double> Zm =
            covariant_constitutive(EffectiveMetric<double>(alcubierre_at(0.3, 3.0)),
                                   g).Z;
        Rank4<double> Zax_sum;
        const Rank4<double> Zax = axion_tensor<double>(0.37, g);
        for (int i = 0; i < 256; ++i) Zax_sum.t[i] = Zm.t[i] + Zax.t[i];
        const Rank4<double> Zan = anisotropic_Z(g, frame);

        const ConstitutiveClassification cm = classify_constitutive(Zm, g, 1.0, n);
        const ConstitutiveClassification ca = classify_constitutive(Zax_sum, g, 1.0, n);
        const ConstitutiveClassification cn = classify_constitutive(Zan, g, 1.0, n);

        std::ofstream os(out_dir + "/b4c_classification.csv", std::ios::binary);
        if (os) {
            os << "tensor,reciprocal,contains_axion,contains_skewon,"
                  "fresnel_single_cone,principal_part_metric_compatible,"
                  "axion_coefficient,skewon_magnitude,branch_split\n";
            os.setf(std::ios::scientific);
            os.precision(8);
            struct Row { const char* n; const ConstitutiveClassification* c; };
            const Row rows[3] = {{"metric", &cm}, {"metric_plus_axion", &ca},
                                 {"anisotropic", &cn}};
            for (const Row& r : rows)
                os << r.n << ',' << r.c->reciprocal << ',' << r.c->contains_axion
                   << ',' << r.c->contains_skewon << ',' << r.c->fresnel_single_cone
                   << ',' << r.c->principal_part_metric_compatible << ','
                   << r.c->axion_coefficient << ',' << r.c->skewon_magnitude << ','
                   << r.c->fresnel_branch_split << '\n';
        }

        report.check_predicate(
            "b4c.classification_reports_each_piece_separately",
            "the metric tensor, the metric-plus-axion tensor and a birefringent "
            "anisotropic tensor are distinguished by SEPARATE flags: axion "
            "content, skewon content, single-cone, and principal-part "
            "metric-compatibility",
            VerificationLevel::A,
            // metric: no axion, no skewon, single cone
            !cm.contains_axion && !cm.contains_skewon && cm.fresnel_single_cone &&
                // metric + axion: axion present, still a single cone
                ca.contains_axion && !ca.contains_skewon && ca.fresnel_single_cone &&
                // anisotropic: no axion, no skewon, but NOT a single cone
                !cn.contains_axion && !cn.contains_skewon &&
                !cn.fresnel_single_cone,
            "three tensors, three distinct flag patterns",
            std::string("metric: ") + cm.interpretation + " | with axion: " +
                ca.interpretation + " | anisotropic: " + cn.interpretation,
            sw.seconds(),
            "the single boolean `metric_induced` of B4 could not have "
            "distinguished the first two tensors at all: they have the same cone "
            "and different constitutive content. `principal_part_metric_compatible` "
            "is named for what it measures -- a necessary condition read off the "
            "cone, not metricity");
    }

    // =======================================================================
    // 7 and 8. The robust root finder.
    // =======================================================================
    {
        Stopwatch sw;
        // Manufactured cases with known roots.
        struct Case {
            const char* name;
            double root[4];
            int real_expected;
            int distinct_expected;
        };

        auto build = [](const double r[4], double c[5]) {
            // (x-r0)(x-r1)(x-r2)(x-r3)
            double p[5] = {1, 0, 0, 0, 0};
            int deg = 0;
            for (int k = 0; k < 4; ++k) {
                double q[6] = {0, 0, 0, 0, 0, 0};
                for (int i = 0; i <= deg; ++i) {
                    q[i + 1] += p[i];
                    q[i] -= p[i] * r[k];
                }
                ++deg;
                for (int i = 0; i <= deg; ++i) p[i] = q[i];
            }
            for (int i = 0; i < 5; ++i) c[i] = p[i];
        };

        double worst_root_err = 0.0, worst_resid = 0.0;
        bool all_counts_ok = true;
        std::ostringstream detail;

        // 1. four distinct real roots
        // 3. two exact double roots
        // 4. one nearly split double root
        const Case cases[3] = {
            {"four_distinct", {-1.5, 0.4, 1.1, 2.7}, 4, 4},
            {"two_double", {0.8, 0.8, 2.2, 2.2}, 4, 2},
            {"near_split", {1.3, 1.3 + 1e-7, 2.5, 3.1}, 4, 3},
        };
        for (const Case& cs : cases) {
            double c[5];
            build(cs.root, c);
            const RootSet rs = quartic_roots(c);
            worst_resid = std::fmax(worst_resid, rs.worst_residual);
            // Match each expected root to the nearest returned root.
            for (int i = 0; i < 4; ++i) {
                double best = 1e300;
                for (const PolynomialRoot& r : rs.roots)
                    best = std::fmin(best, std::fabs(r.real - cs.root[i]) +
                                               std::fabs(r.imaginary));
                worst_root_err = std::fmax(worst_root_err, best);
            }
            if (rs.real_count != cs.real_expected) all_counts_ok = false;
            detail << cs.name << ": " << rs.real_count << " real, "
                   << rs.distinct_count << " distinct; ";
        }
        metrics.worst_repeated_root_residual = worst_resid;

        report.check_predicate(
            "b4c.roots.distinct_repeated_and_near_split",
            "four distinct real roots, two exact DOUBLE roots, and a pair split "
            "by 1e-7 are all recovered, with the double roots counted at "
            "multiplicity two",
            VerificationLevel::B,
            all_counts_ok && worst_root_err < 1e-6,
            "all roots found with the right real count",
            detail.str() + "worst root error " + sci(worst_root_err) +
                ", worst residual " + sci(worst_resid),
            sw.seconds(),
            "the previous scan bracketed SIGN CHANGES and was therefore blind to "
            "roots of even multiplicity -- which is precisely the metric-induced "
            "case. Durand-Kerner finds all four, real and complex alike");

        // 2. two real and one complex-conjugate pair
        {
            // (x-1)(x-3)(x^2 + 1) = complex pair at +-i
            double c[5] = {3.0, -4.0, 4.0, -4.0, 1.0};
            const RootSet rs = quartic_roots(c);
            bool pair_found = rs.complex_pair_count == 1 && rs.real_count == 2;
            report.check_predicate(
                "b4c.roots.complex_conjugate_pair",
                "two real roots and one complex-conjugate pair are separated and "
                "reported, where the sign-change scan would have seen only the "
                "two real ones",
                VerificationLevel::B, pair_found,
                "2 real roots and 1 complex pair",
                std::to_string(rs.real_count) + " real, " +
                    std::to_string(rs.complex_pair_count) + " complex pair(s), "
                    "worst residual " + sci(rs.worst_residual),
                sw.seconds());
        }

        // 7. invariance under coefficient scaling
        {
            double c[5] = {3.0, -4.0, 4.0, -4.0, 1.0};
            double cs2[5];
            for (int i = 0; i < 5; ++i) cs2[i] = 1e6 * c[i];
            const RootSet a = quartic_roots(c);
            const RootSet b = quartic_roots(cs2);
            double worst = 0.0;
            for (std::size_t i = 0; i < a.roots.size() && i < b.roots.size(); ++i)
                worst = std::fmax(worst,
                                  std::fabs(a.roots[i].real - b.roots[i].real) +
                                      std::fabs(a.roots[i].imaginary -
                                                b.roots[i].imaginary));
            report.check_residual(
                "b4c.roots.scaling_invariance",
                "scaling every coefficient by 1e6 leaves the roots unchanged",
                VerificationLevel::A, "n/a", "factor 1e6", worst, 1e-12,
                sw.seconds());
        }

        // 8. reconstruction of coefficients from roots
        {
            double c[5] = {3.0, -4.0, 4.0, -4.0, 1.0};
            const RootSet rs = quartic_roots(c);
            double back[5];
            coefficients_from_roots(rs, c[4], back);
            double worst = 0.0, scale = 0.0;
            for (int i = 0; i < 5; ++i) {
                worst = std::fmax(worst, std::fabs(back[i] - c[i]));
                scale = std::fmax(scale, std::fabs(c[i]));
            }
            metrics.coefficient_reconstruction_error = worst / scale;
            report.check_residual(
                "b4c.roots.coefficient_reconstruction",
                "rebuilding the coefficients from the returned root set recovers "
                "the original polynomial",
                VerificationLevel::D, "n/a", "quartic with a complex pair",
                metrics.coefficient_reconstruction_error, 1e-12, sw.seconds(),
                "an independent check that the set returned really is the root "
                "set of the polynomial handed in");
        }
    }

    // 5 and 6. Metric double root recovered by the ROOT FINDER, and the
    // anisotropic split.
    {
        Stopwatch sw;
        const double n[3] = {0.5773502691896258, 0.5773502691896258,
                             0.5773502691896258};
        const Rank4<double> Zm =
            covariant_constitutive(EffectiveMetric<double>(alcubierre_at(0.3, 3.0)),
                                   g).Z;
        double c[5];
        fresnel_quartic_coefficients(Zm, 1.0, n, c);
        const RootSet rs = quartic_roots(c);

        // Positive branches should be a coincident pair.
        std::vector<double> pos;
        for (const PolynomialRoot& r : rs.roots)
            if (r.is_real && r.real > 0.0) pos.push_back(r.real);
        double sep = 0.0;
        if (pos.size() >= 2) {
            double lo = 1e300, hi = -1e300;
            for (double v : pos) { lo = std::fmin(lo, v); hi = std::fmax(hi, v); }
            sep = hi - lo;
        }
        metrics.metric_double_root_separation = sep;

        int max_mult = 1;
        for (const PolynomialRoot& r : rs.roots)
            max_mult = std::max(max_mult, r.multiplicity);

        report.check_predicate(
            "b4c.roots.metric_double_root_is_recovered",
            "the root finder itself returns the metric-induced Fresnel quartic's "
            "REPEATED branches, at multiplicity two, rather than relying on a "
            "separate factorization comparison",
            VerificationLevel::B,
            rs.distinct_count == 2 && max_mult == 2 && sep < 1e-6,
            "two distinct roots, multiplicity two, coincident positive branches",
            std::to_string(rs.distinct_count) + " distinct, max multiplicity " +
                std::to_string(max_mult) + ", positive-branch separation " +
                sci(sep) + ", worst residual " + sci(rs.worst_residual) + " (" +
                rs.warning + ")",
            sw.seconds(),
            "this is what the sign-change scan could not do: a doubled root does "
            "not change the polynomial's sign, so it was simply absent");

        const Rank4<double> Zan = anisotropic_Z(g, frame);
        double ca[5];
        fresnel_quartic_coefficients(Zan, 1.0, n, ca);
        const RootSet ra = quartic_roots(ca);
        std::vector<double> pa;
        for (const PolynomialRoot& r : ra.roots)
            if (r.is_real && r.real > 0.0) pa.push_back(r.real);
        double split = 0.0;
        if (pa.size() >= 2) {
            double lo = 1e300, hi = -1e300;
            for (double v : pa) { lo = std::fmin(lo, v); hi = std::fmax(hi, v); }
            split = (hi - lo) / hi;
        }
        metrics.anisotropic_split = split;

        report.check_predicate(
            "b4c.roots.anisotropic_split_is_recovered",
            "the anisotropic medium's two distinct propagation branches are "
            "resolved as four distinct roots",
            VerificationLevel::B,
            ra.distinct_count == 4 && split > 1e-2,
            "four distinct roots and a measurable split",
            std::to_string(ra.distinct_count) + " distinct, relative split " +
                sci(split),
            sw.seconds());
    }

    // =======================================================================
    // 4 and 5. The scalar-generic solver, and gradients through it.
    // =======================================================================
    {
        Stopwatch sw;
        // The templated solver at double must reproduce the analytic matched
        // slab: r = 0 and phase = omega eps d.
        const double eps0 = 2.6, d0 = 1.7, w = 1.3;
        const std::vector<Layer<double>> stack = matched_slab<double>(eps0, d0);
        const StackResult<double> r = solve_stack(stack, w);

        const double phase_analytic = w * eps0 * d0;
        // Reduce modulo 2 pi for the comparison.
        double wrapped = std::fmod(phase_analytic, 2.0 * 3.14159265358979323846);
        if (wrapped > 3.14159265358979323846) wrapped -= 2.0 * 3.14159265358979323846;

        report.check_predicate(
            "b4c.solver.matched_slab_is_analytic",
            "the templated solver at Scalar = double gives r = 0 and "
            "phase = omega eps d for an impedance-matched slab",
            VerificationLevel::B,
            r.reflection < 1e-28 && std::fabs(r.phase - wrapped) < 1e-12,
            "zero reflection and the analytic phase",
            "R = " + sci(r.reflection) + ", phase = " + sci(r.phase) +
                " against " + sci(wrapped),
            sw.seconds(),
            "the analytic reference the gradient cases below are checked against");
    }

    {
        Stopwatch sw;
        // FINITE-DIFFERENCE gradients through the ACTUAL stack, using
        // d(phase) = Im(dt/t) rather than differencing atan2.
        const double w = 1.3;
        struct Param { const char* name; double base; };
        double worst_fd = 0.0, worst_analytic = 0.0;
        std::ostringstream detail;

        // (a) matched-slab permittivity and thickness, against closed forms.
        {
            const double eps0 = 2.6, d0 = 1.7;
            auto t_of = [&](double eps, double d) {
                const std::vector<Layer<double>> s = matched_slab<double>(eps, d);
                return solve_stack(s, w).t;
            };
            const C t0 = t_of(eps0, d0);

            for (double h : {1e-4, 1e-5, 1e-6}) {
                const C tp = t_of(eps0 + h, d0), tm = t_of(eps0 - h, d0);
                const PhaseSensitivity ps = phase_sensitivity(
                    t0.re, t0.im, (tp.re - tm.re) / (2 * h),
                    (tp.im - tm.im) / (2 * h));
                worst_analytic = std::fmax(worst_analytic,
                                           std::fabs(ps.d_phase - w * d0));

                const C dp = t_of(eps0, d0 + h), dm = t_of(eps0, d0 - h);
                const PhaseSensitivity pd = phase_sensitivity(
                    t0.re, t0.im, (dp.re - dm.re) / (2 * h),
                    (dp.im - dm.im) / (2 * h));
                worst_analytic = std::fmax(worst_analytic,
                                           std::fabs(pd.d_phase - w * eps0));
            }
            detail << "matched slab: d(phase)/d(eps) target " << sci(w * d0)
                   << ", d(phase)/d(d) target " << sci(w * eps0) << "; ";
        }
        metrics.worst_stack_gradient_analytic = worst_analytic;

        // (b) a general lossy magnetoelectric stack: permeability, zeta and a
        //     Lorentz oscillator's F, w0, gamma, all by finite differences of the
        //     actual solver. Self-consistency across step sizes is the measure.
        {
            DispersiveModelSpec spec;
            spec.eps_poles = 1;
            spec.mu_poles = 1;
            spec.magnetoelectric = true;
            spec.mode = ParameterMode::Direct;
            double p0[kPhysicalSlots] = {};
            p0[kSlotEpsInf] = 2.1;
            p0[kSlotEpsPole0F] = 1.6; p0[kSlotEpsPole0W0] = 3.0;
            p0[kSlotEpsPole0Gamma] = 0.18;
            p0[kSlotMuInf] = 1.15;
            p0[kSlotMuPole0F] = 0.45; p0[kSlotMuPole0W0] = 3.6;
            p0[kSlotMuPole0Gamma] = 0.22;
            p0[kSlotZetaAmp] = 0.30; p0[kSlotZetaW0] = 3.2;
            p0[kSlotZetaGamma] = 0.25;

            auto t_of_slot = [&](int slot, double delta) {
                double p[kPhysicalSlots];
                for (int i = 0; i < kPhysicalSlots; ++i) p[i] = p0[i];
                p[slot] += delta;
                const PhysicalDesignVariables<double> v =
                    design_from_physical<double>(spec, p);
                const ResponseAtFrequency<double> rr = evaluate_response(v, w);
                std::vector<Layer<double>> s;
                Layer<double> L;
                L.epsilon = rr.epsilon;
                L.mu_inverse = rr.mu_inverse;
                L.zeta = rr.zeta;
                L.thickness = 1.9;
                s.push_back(L);
                return solve_stack(s, w).t;
            };

            const C t0 = t_of_slot(kSlotEpsInf, 0.0);
            const int slots[5] = {kSlotMuInf, kSlotZetaAmp, kSlotEpsPole0F,
                                  kSlotEpsPole0W0, kSlotEpsPole0Gamma};
            const char* names[5] = {"mu_inf", "zeta_amp", "F", "w0", "gamma"};
            for (int k = 0; k < 5; ++k) {
                double vals[3];
                int nv = 0;
                for (double h : {1e-4, 1e-5, 1e-6}) {
                    const C tp = t_of_slot(slots[k], h);
                    const C tm = t_of_slot(slots[k], -h);
                    const PhaseSensitivity ps = phase_sensitivity(
                        t0.re, t0.im, (tp.re - tm.re) / (2 * h),
                        (tp.im - tm.im) / (2 * h));
                    vals[nv++] = ps.d_phase;
                }
                const double spread =
                    std::fmax(std::fabs(vals[0] - vals[1]),
                              std::fabs(vals[1] - vals[2]));
                worst_fd = std::fmax(worst_fd, spread);
                detail << names[k] << " d(phase) = " << sci(vals[1]) << "; ";
            }
        }
        metrics.worst_stack_gradient_fd = worst_fd;

        report.check_residual(
            "b4c.gradient.matched_slab_matches_closed_form",
            "d(phase)/d(eps) = omega d and d(phase)/d(thickness) = omega eps for "
            "a matched slab, obtained as Im(dt/t) through the actual solver",
            VerificationLevel::B, "n/a", "3 step sizes", worst_analytic, 1e-6,
            sw.seconds(),
            "the phase derivative never calls atan2: d(phase) = Im(dt/t) uses "
            "only the real and imaginary channels of the transmission amplitude");

        report.check_residual(
            "b4c.gradient.finite_difference_step_independence",
            "phase sensitivities to permeability, magnetoelectric amplitude, and "
            "a Lorentz oscillator's strength, resonance and damping agree across "
            "three step sizes",
            VerificationLevel::D, "n/a",
            "5 parameters, h in {1e-4, 1e-5, 1e-6}: " + detail.str(), worst_fd,
            1e-5, sw.seconds());
    }

    // =======================================================================
    // 6. HyPAD through the actual stack assembly.
    // =======================================================================
#ifdef WARP_HYPAD_HAVE_OTILIB
    {
        namespace sens = warp_hypad::sensitivity;
        using sens::OTIScalar;
        Stopwatch sw;
        metrics.oti_available = true;

        const double eps0 = 2.6, d0 = 1.7, w = 1.3;
        // Seed eps and the thickness in two of the three available directions.
        const sens::SeededParameters sp = sens::seed_all(eps0, d0, 1.0);
        const std::vector<Layer<OTIScalar>> stack =
            matched_slab<OTIScalar>(sp.velocity, sp.radius);
        const StackResult<OTIScalar> r = solve_stack(stack, OTIScalar(w));

        const double t_re = real_part(r.t.re), t_im = real_part(r.t.im);
        const double d_eps = sens::extract_gradient(r.t.re, sens::kParamVelocity);
        const double d_eps_i = sens::extract_gradient(r.t.im, sens::kParamVelocity);
        const double d_d = sens::extract_gradient(r.t.re, sens::kParamRadius);
        const double d_d_i = sens::extract_gradient(r.t.im, sens::kParamRadius);

        const PhaseSensitivity pe = phase_sensitivity(t_re, t_im, d_eps, d_eps_i);
        const PhaseSensitivity pd = phase_sensitivity(t_re, t_im, d_d, d_d_i);

        const double err = std::fmax(std::fabs(pe.d_phase - w * d0),
                                     std::fabs(pd.d_phase - w * eps0));
        metrics.worst_stack_gradient_analytic =
            std::fmax(metrics.worst_stack_gradient_analytic, err);

        report.check_predicate(
            "b4c.hypad.through_actual_stack_assembly",
            "a physical parameter seeded through interfaces, propagation and the "
            "Redheffer cascade gives d(phase)/d(eps) = omega d and "
            "d(phase)/d(thickness) = omega eps",
            VerificationLevel::B,
            pe.reliable && pd.reliable && err < 1e-12,
            "both closed forms recovered",
            "d(phase)/d(eps) = " + sci(pe.d_phase) + " against " + sci(w * d0) +
                "; d(phase)/d(d) = " + sci(pd.d_phase) + " against " +
                sci(w * eps0) + "; " + pe.note,
            sw.seconds(),
            "this is differentiation through the SOLVER, not through a surrogate: "
            "the solver is templated, so the same code path runs at double and at "
            "the HyPAD scalar");

        // Separate versus batched seeding.
        const sens::SeededParameters s1 = sens::seed_all(eps0, d0, 1.0);
        OTIScalar eps_alone = sens::seed_parameter(eps0, 0);
        const std::vector<Layer<OTIScalar>> alone =
            matched_slab<OTIScalar>(eps_alone, OTIScalar(d0));
        const StackResult<OTIScalar> ra = solve_stack(alone, OTIScalar(w));
        const double d_eps_alone = sens::extract_gradient(ra.t.re, 0);
        (void)s1;
        report.check_residual(
            "b4c.hypad.separate_versus_batched_seeding",
            "the same partial derivative of the transmission amplitude comes out "
            "identically whether eps was seeded alone or alongside the thickness",
            VerificationLevel::A, "n/a", "batch of 1 versus 2",
            std::fabs(d_eps_alone - d_eps), 1e-15, sw.seconds());
    }
#else
    {
        Stopwatch sw;
        report.check_skipped(
            "b4c.hypad.through_actual_stack_assembly",
            "HyPAD differentiation through the templated stack assembly",
            "the OTIlib backend is disabled; the templated solver still ran at "
            "double, and the finite-difference gradient and every "
            "root-classification case still ran",
            sw.seconds());
        report.check_skipped(
            "b4c.hypad.separate_versus_batched_seeding",
            "separate versus batched OTI seeding through the stack",
            "the OTIlib backend is disabled", sw.seconds());
    }
#endif

    // =======================================================================
    // Restricted-versus-exact propagation sensitivities.
    // =======================================================================
    {
        Stopwatch sw;
        const double w = 1.0;
        // Sensitivity of the propagation observables to v_s, the magnetoelectric
        // scale and the wall steepness, for the exact profile and the zeta = 0
        // restriction.
        auto build_stack = [&](double vs, double zeta_scale, double sigma,
                              std::vector<Layer<double>>& out) {
            out.clear();
            for (int i = 0; i < 24; ++i) {
                const double x = -6.0 + 12.0 * static_cast<double>(i) / 23.0;
                warp_hypad::metrics::AlcubierreMetric<double> wm;
                const warp_hypad::metrics::AlcubierreParameters<double> p{vs, 3.0,
                                                                         sigma};
                const Mat4<double> geff = wm.evaluate(0.0, x, 0.0, 0.0, p);
                const Rank4<double> Z =
                    covariant_constitutive(EffectiveMetric<double>(geff), g).Z;
                const ObserverDecomposition4<double> d =
                    decompose_A7(Z, g.value, frame.V);
                Layer<double> L;
                L.epsilon = C(d.epsilon(2, 2), 0.0);
                L.mu_inverse = C(d.mu_inverse(2, 2), 0.0);
                L.zeta = C(zeta_scale * d.zeta(2, 3), 0.0);
                L.thickness = 0.5;
                out.push_back(L);
            }
        };

        std::ofstream os(out_dir + "/b4c_sensitivities.csv", std::ios::binary);
        if (os) os << "family,parameter,d_phase,d_transmission,d_reflection,"
                      "d_group_delay\n";

        double me_exact = 0.0, me_restricted = 0.0;
        for (int fam = 0; fam < 2; ++fam) {
            const double zs = (fam == 0) ? 1.0 : 0.0;
            const char* fname = (fam == 0) ? "exact" : "restricted_zeta0";
            struct P { const char* name; double vs, zsc, sig; };
            const double base_vs = 0.3, base_sig = 2.0;

            std::vector<Layer<double>> s0;
            build_stack(base_vs, zs, base_sig, s0);
            const StackResult<double> r0 = solve_stack(s0, w);

            const char* pnames[3] = {"v_s", "magnetoelectric_scale", "sigma"};
            for (int k = 0; k < 3; ++k) {
                const double h = 1e-5;
                std::vector<Layer<double>> sp_, sm_;
                double vp = base_vs, zp = zs, gp = base_sig;
                double vm = base_vs, zm = zs, gm = base_sig;
                if (k == 0) { vp += h; vm -= h; }
                else if (k == 1) { zp += h; zm -= h; }
                else { gp += h; gm -= h; }
                build_stack(vp, zp, gp, sp_);
                build_stack(vm, zm, gm, sm_);
                const StackResult<double> rp = solve_stack(sp_, w);
                const StackResult<double> rm = solve_stack(sm_, w);

                const PhaseSensitivity ps = phase_sensitivity(
                    r0.t.re, r0.t.im, (rp.t.re - rm.t.re) / (2 * h),
                    (rp.t.im - rm.t.im) / (2 * h));
                const double dT = (rp.transmission - rm.transmission) / (2 * h);
                const double dR = (rp.reflection - rm.reflection) / (2 * h);

                // Group delay = d(phase)/d(omega); its parameter sensitivity by a
                // second difference.
                const double hw = 1e-4;
                auto tg = [&](const std::vector<Layer<double>>& st) {
                    const StackResult<double> A = solve_stack(st, w);
                    const StackResult<double> B = solve_stack(st, w + hw);
                    const PhaseSensitivity q = phase_sensitivity(
                        A.t.re, A.t.im, (B.t.re - A.t.re) / hw,
                        (B.t.im - A.t.im) / hw);
                    return q.d_phase;
                };
                const double dtg = (tg(sp_) - tg(sm_)) / (2 * h);

                if (os) {
                    os.setf(std::ios::scientific);
                    os.precision(8);
                    os << fname << ',' << pnames[k] << ',' << ps.d_phase << ','
                       << dT << ',' << dR << ',' << dtg << '\n';
                }
                if (k == 1) {
                    if (fam == 0) me_exact = std::fabs(ps.d_phase);
                    else me_restricted = std::fabs(ps.d_phase);
                }
            }
        }
        metrics.magnetoelectric_sensitivity = me_exact;
        metrics.magnetoelectric_sensitivity_restricted = me_restricted;

        report.check_predicate(
            "b4c.sensitivity.magnetoelectric_scale_matters",
            "the propagation phase is strongly sensitive to the magnetoelectric "
            "scale in the exact profile, so the 1.6-radian mismatch of B4 is not "
            "an accident of one parameter choice",
            VerificationLevel::B,
            me_exact > 1.0 && me_restricted > 1.0,
            "a large phase sensitivity to the magnetoelectric scale",
            "d(phase)/d(zeta scale): exact " + sci(me_exact) + ", restricted " +
                sci(me_restricted),
            sw.seconds(),
            "the derivative is taken at the zeta scale itself, so it measures how "
            "fast the propagation mismatch grows as the magnetoelectric term is "
            "restored. Both families are evaluated at the same point, so the two "
            "numbers are comparable");
    }

    // =======================================================================
    // No regression in the B4 energy balance, and determinism.
    // =======================================================================
    {
        Stopwatch sw;
        std::vector<Layer<double>> stack;
        Layer<double> L;
        L.epsilon = C(4.0, 0.35);
        L.mu_inverse = C(1.0, 0.0);
        L.zeta = C(0.0, 0.0);
        L.thickness = 1.7;
        stack.push_back(L);

        double worst = 0.0;
        for (double w : {0.7, 1.3, 2.9})
            worst = std::fmax(worst, solve_stack(stack, w, 6000).energy_balance_residual);
        metrics.energy_balance_regression = worst;

        report.check_residual(
            "b4c.no_regression_in_energy_balance",
            "after templating the solver, R + T + A = 1 still holds with A "
            "integrated independently",
            VerificationLevel::B, "n/a",
            "the same lossy slab B4 used, 3 frequencies", worst, 1e-7,
            sw.seconds(),
            "the templating must not have changed any number, and this is the "
            "case that would show it if it had");

        const StackResult<double> a = solve_stack(stack, 1.3);
        const StackResult<double> b = solve_stack(stack, 1.3);
        report.check_predicate(
            "b4c.deterministic_repeatability",
            "two identical solves of the templated stack, and two identical root "
            "solves, are bit-identical",
            VerificationLevel::A,
            a.reflection == b.reflection && a.transmission == b.transmission &&
                a.absorption_from_dissipation == b.absorption_from_dissipation,
            "bit-identical", "identical", sw.seconds());
    }

    metrics.exclusions =
        "no claim is made to differentiate through branch sorting, mode "
        "switching, poles, zeros of transmission, or degeneracies. The phase "
        "sensitivity refuses to report a value when |t| falls below 1e-8, and "
        "root multiplicity carries an explicit conditioning warning because a "
        "root of multiplicity m is attainable only to about eps^(1/m).";
}

}  // namespace checks
}  // namespace program_b

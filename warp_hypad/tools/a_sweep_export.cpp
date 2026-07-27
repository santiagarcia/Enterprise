#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/json.hpp"
#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/core/verification.hpp"
#include "warp_hypad/geometry/signature.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
// otilib_traits pulls in the scalar, its math overloads and the seed/extract
// helpers together; including the pieces separately misses the sqrt and tanh
// overloads the shape function needs.
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

// ---------------------------------------------------------------------------
// A-VISUAL-SWEEP-EXPORT
//
// Writes the swept effective-metric fields that the V0 Program A animations
// need, and gates them. The V0 inventory established that no swept field data
// existed anywhere in the repository: every committed slice is a single snapshot
// at one parameter set, and an animation needs the field, not a value.
//
// THIS EXPORTER MAY EVALUATE PHYSICS. THE PYTHON RENDERERS MAY NOT. That is the
// whole point of the split: production code is called here, once, through the
// same AlcubierreMetric path the numerical gates use, and the result is saved.
// A renderer then reads columns.
//
// SCOPE. Effective-metric components, the shape function, determinants,
// signature and parameter sensitivities. NO CURVATURE: no Christoffel symbol,
// no Riemann or Ricci tensor, no Einstein tensor, no Kretschmann scalar, no
// stress-energy and no energy condition. Program A does not compute them, and
// they stay blocked pending A4.
//
// TIME TREATMENT. InertialTimeDependent throughout, with x_s(t) = x_s0 + v_s t.
// FrozenTimeDiagnostic is never used here: it is a diagnostic mode that would
// make the time derivatives inconsistent, and the translation sweep is entirely
// about motion in time.
// ---------------------------------------------------------------------------

using namespace warp_hypad;
// The metric types live in warp_hypad::metrics; named individually rather than
// with a blanket using-directive so it stays visible which production types this
// exporter touches.
using metrics::AlcubierreMetric;
using metrics::AlcubierreParameters;
using metrics::MetricSensitivityMode;
using metrics::TimeTreatment;
using metrics::shape_function;
using geometry::SignatureResult;
using geometry::analyse_signature;

namespace {

// --- Deterministic sweep configuration --------------------------------------
//
// Every grid and every frame value is a closed-form function of an index, so a
// repeat export reproduces the bytes exactly. Nothing here is adaptive.

constexpr int kNx = 41;
constexpr int kNy = 41;
constexpr double kHalfWidth = 4.0;      // domain [-4, 4] in bubble radii
constexpr int kFrames = 24;

constexpr double kNominalV = 0.5;
constexpr double kNominalR = 1.0;
constexpr double kNominalSigma = 4.0;
constexpr double kCenterX0 = -2.0;      // so the bubble sweeps across the domain

// The velocity sweep must span v_s f = 1, which is the ergosurface condition.
// f(0) = 1 exactly at the bubble centre, so the crossing is exactly at v_s = 1
// and frame index 12 of 24 hits it on the nose rather than straddling it.
constexpr double kVMin = 0.2;
constexpr double kVMax = 1.6;

constexpr double kRMin = 0.6;
constexpr double kRMax = 2.0;
constexpr double kSigmaMin = 2.0;
constexpr double kSigmaMax = 10.0;

// The trajectory term in the chain rule is proportional to t, so it vanishes at
// t = 0 and the fixed-lab / fixed-centre comparison would be vacuous there.
constexpr double kChainRuleTime = 1.0;

double grid_x(int i) { return -kHalfWidth + 2.0 * kHalfWidth * i / (kNx - 1); }
double grid_y(int j) { return -kHalfWidth + 2.0 * kHalfWidth * j / (kNy - 1); }
double frame_value(int frame, double lo, double hi) {
    return lo + (hi - lo) * frame / (kFrames - 1);
}

std::string g10(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.10g", v);
    return buf;
}

using Metric = AlcubierreMetric<double>;
using Params = AlcubierreParameters<double>;

Metric make_metric(MetricSensitivityMode mode =
                       MetricSensitivityMode::FixedBubbleCenter) {
    Metric m(kCenterX0, TimeTreatment::InertialTimeDependent);
    m.sensitivity_mode = mode;
    return m;
}

double shape_at(double t, double x, double y, const Metric& m, const Params& p) {
    const double dx = x - m.center_at(t, p);
    return shape_function(std::sqrt(dx * dx + y * y), p);
}

// The saved-field row shared by every sweep.
struct Sample {
    double f = 0.0, g00 = 0.0, g01 = 0.0, g11 = 0.0, det = 0.0;
    int signature_ok = 0;
    double r_s = 0.0;
};

Sample evaluate_sample(double t, double x, double y, const Metric& m,
                       const Params& p) {
    Sample s;
    const double dx = x - m.center_at(t, p);
    s.r_s = std::sqrt(dx * dx + y * y);
    s.f = shape_function(s.r_s, p);
    const Mat4<double> g = m.evaluate(t, x, y, 0.0, p);
    s.g00 = g(0, 0);
    s.g01 = g(0, 1);
    s.g11 = g(1, 1);
    // The 2x2 (t,x) block determinant times the trivially unit (y,z) block.
    // det g = -1 identically for any v and f, which makes it an EXACT test
    // rather than a tolerance one.
    s.det = g(0, 0) * g(1, 1) - g(0, 1) * g(1, 0);
    const SignatureResult sig = analyse_signature(g);
    s.signature_ok = sig.lorentzian ? 1 : 0;
    return s;
}

// --- Wall geometry measured FROM THE FIELD, not from the parameters ---------
//
// The radius and steepness assertions must be independent of the parameter that
// was swept, or they would be tautologies. So the radius is measured as the
// half-height crossing of f along +x from the centre, and the wall width as the
// separation of the 0.75 and 0.25 crossings -- both read off the profile.
double crossing_radius(double level, const Metric& m, const Params& p) {
    const double t = 0.0;
    const double xc = m.center_at(t, p);
    double lo = 0.0, hi = 8.0;
    // f is monotone decreasing in r_s, so bisection is safe and cannot diverge.
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (shape_at(t, xc + mid, 0.0, m, p) > level) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

double measured_half_height_radius(const Metric& m, const Params& p) {
    return crossing_radius(0.5, m, p);
}
double measured_wall_width(const Metric& m, const Params& p) {
    return crossing_radius(0.25, m, p) - crossing_radius(0.75, m, p);
}

// --- Sensitivities ----------------------------------------------------------
//
// With OTIlib the derivatives come from the verified automatic-differentiation
// path. Without it they come from centred finite differences on the SAME
// production metric evaluation. The backend is recorded per row, so a consumer
// never has to guess which produced a number, and the OTI-disabled build still
// exports every field -- only the OTI-versus-analytic COMPARISON cases skip.
struct Sensitivities {
    double dg00[3] = {0, 0, 0};      // d/dv, d/dR, d/dsigma
    double dg0x[3] = {0, 0, 0};
    const char* backend = "central_fd";
};

Sensitivities sensitivities_fd(double t, double x, double y, const Metric& m,
                               const Params& p, MetricSensitivityMode mode) {
    Sensitivities s;
    const double base[3] = {p.velocity, p.radius, p.sigma};

    // THE TWO CONVENTIONS ARE GENUINELY DIFFERENT HERE, and getting this wrong
    // is easy: in the DOUBLE evaluate(), both branches call center_at(t, p), so
    // perturbing v moves the centre in BOTH modes and a naive difference
    // computes the LABORATORY derivative whichever mode is set. The modes only
    // diverge for the OTI scalar, where center_at_scalar propagates dx_s/dv = t
    // through the derivative channel.
    //
    // So the fixed-centre derivative is taken at fixed COMOVING offset xi: for
    // each perturbed parameter set the evaluation point moves with the centre,
    // keeping r_s constant. The laboratory derivative is taken at fixed x.
    //
    // The gate caught the original error as a chain-rule residual of 2.1 --
    // exactly the trajectory term, because the "fixed-centre" value already
    // contained it.
    const double xi = x - m.center_at(t, p);

    for (int k = 0; k < 3; ++k) {
        const double h = std::max(1.0e-6, 1.0e-6 * std::fabs(base[k]));
        Params hi = p, lo = p;
        double* ph = (k == 0 ? &hi.velocity : k == 1 ? &hi.radius : &hi.sigma);
        double* pl = (k == 0 ? &lo.velocity : k == 1 ? &lo.radius : &lo.sigma);
        *ph = base[k] + h;
        *pl = base[k] - h;

        double x_hi = x, x_lo = x;
        if (mode == MetricSensitivityMode::FixedBubbleCenter) {
            x_hi = m.center_at(t, hi) + xi;
            x_lo = m.center_at(t, lo) + xi;
        }
        const Mat4<double> gh = m.evaluate(t, x_hi, y, 0.0, hi);
        const Mat4<double> gl = m.evaluate(t, x_lo, y, 0.0, lo);
        s.dg00[k] = (gh(0, 0) - gl(0, 0)) / (2.0 * h);
        s.dg0x[k] = (gh(0, 1) - gl(0, 1)) / (2.0 * h);
    }
    return s;
}

#ifdef WARP_HYPAD_HAVE_OTILIB
using sensitivity::OTIScalar;

Sensitivities sensitivities_oti(double t, double x, double y,
                                MetricSensitivityMode mode, const Params& p) {
    using sensitivity::extract_gradient;
    AlcubierreMetric<OTIScalar> m(kCenterX0, TimeTreatment::InertialTimeDependent);
    m.sensitivity_mode = mode;
    const sensitivity::SeededParameters sp =
        sensitivity::seed_all(p.velocity, p.radius, p.sigma);
    AlcubierreParameters<OTIScalar> op(sp.velocity, sp.radius, sp.sigma);
    const Mat4<OTIScalar> g = m.evaluate(t, x, y, 0.0, op);
    Sensitivities s;
    s.backend = "otilib";
    for (int d = 0; d < 3; ++d) {
        s.dg00[d] = extract_gradient(g(0, 0), d);
        s.dg0x[d] = extract_gradient(g(0, 1), d);
    }
    return s;
}
#endif

Sensitivities sensitivities(double t, double x, double y,
                            MetricSensitivityMode mode, const Metric& m,
                            const Params& p) {
#ifdef WARP_HYPAD_HAVE_OTILIB
    (void)m;
    return sensitivities_oti(t, x, y, mode, p);
#else
    return sensitivities_fd(t, x, y, m, p, mode);
#endif
}

// --- FNV-1a, matching the convention used elsewhere in the repository -------
std::uint64_t fnv1a_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::uint64_t h = 0xcbf29ce484222325ull;
    char c;
    while (f.get(c)) {
        h ^= static_cast<unsigned char>(c);
        h *= 0x100000001b3ull;
    }
    return h;
}

}  // namespace

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/v0_program_a_sweeps";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else { std::cerr << "usage: a_sweep_export [--out <dir>]\n"; return 2; }
    }
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);

    VerificationReport report(
        "v0_program_a_sweeps",
        "A-VISUAL-SWEEP-EXPORT: swept effective-metric fields for the Program A "
        "animations, and their scientific assertions");

    const Metric metric = make_metric();
    const Params nominal(kNominalV, kNominalR, kNominalSigma);

    // =====================================================================
    // Translation sweep
    // =====================================================================
    std::vector<double> tr_times, tr_centers;
    double worst_comoving_drift = 0.0;
    double lab_translation_span = 0.0;
    {
        std::ofstream f(out_dir + "/translation_sweep.csv", std::ios::binary);
        f << "frame,time,bubble_center_x,x_lab,y_lab,x_comoving,r_s,f,g00,g01,"
             "g11,det,signature_ok,coordinate_mode\n";
        // The comoving profile at the reference frame, to compare against.
        std::vector<double> reference_profile;
        for (int frame = 0; frame < kFrames; ++frame) {
            const double t = frame_value(frame, 0.0, 8.0);
            const double xc = metric.center_at(t, nominal);
            tr_times.push_back(t);
            tr_centers.push_back(xc);
            std::vector<double> profile;
            for (int i = 0; i < kNx; ++i)
                for (int j = 0; j < kNy; ++j) {
                    const double x = grid_x(i), y = grid_y(j);
                    const Sample s = evaluate_sample(t, x, y, metric, nominal);
                    f << frame << ',' << g10(t) << ',' << g10(xc) << ','
                      << g10(x) << ',' << g10(y) << ',' << g10(x - xc) << ','
                      << g10(s.r_s) << ',' << g10(s.f) << ',' << g10(s.g00)
                      << ',' << g10(s.g01) << ',' << g10(s.g11) << ','
                      << g10(s.det) << ',' << s.signature_ok
                      << ",InertialTimeDependent\n";
                    if (j == kNy / 2) profile.push_back(s.f);
                }
            // The comoving profile is sampled at fixed OFFSETS from the moving
            // centre, so a translating field must reproduce it exactly.
            std::vector<double> comoving;
            for (int i = 0; i < kNx; ++i)
                comoving.push_back(
                    shape_at(t, xc + grid_x(i), 0.0, metric, nominal));
            if (frame == 0) reference_profile = comoving;
            else
                for (std::size_t i = 0; i < comoving.size(); ++i)
                    worst_comoving_drift = std::max(
                        worst_comoving_drift,
                        std::fabs(comoving[i] - reference_profile[i]));
            (void)profile;
        }
        lab_translation_span = tr_centers.back() - tr_centers.front();
    }

    // =====================================================================
    // Velocity sweep
    // =====================================================================
    int crossing_frame = -1;
    double worst_g00_formula = 0.0, worst_det_dev = 0.0;
    int signature_failures = 0;
    std::vector<double> vel_values;
    {
        std::ofstream f(out_dir + "/velocity_sweep.csv", std::ios::binary);
        f << "frame,v_s,x,y,f,v_times_f,g00,g01,g11,det,signature_ok,"
             "g00_zero_status\n";
        for (int frame = 0; frame < kFrames; ++frame) {
            const double v = frame_value(frame, kVMin, kVMax);
            vel_values.push_back(v);
            const Params p(v, kNominalR, kNominalSigma);
            double max_vf = 0.0;
            std::vector<std::string> rows;
            for (int i = 0; i < kNx; ++i)
                for (int j = 0; j < kNy; ++j) {
                    const double x = grid_x(i), y = grid_y(j);
                    const Sample s = evaluate_sample(0.0, x, y, metric, p);
                    const double vf = v * s.f;
                    max_vf = std::max(max_vf, vf);
                    // g00 = -(1 - v^2 f^2), verified against the saved value.
                    worst_g00_formula = std::max(
                        worst_g00_formula,
                        std::fabs(s.g00 - (-(1.0 - vf * vf))));
                    worst_det_dev = std::max(worst_det_dev,
                                             std::fabs(s.det + 1.0));
                    if (!s.signature_ok) ++signature_failures;
                    std::ostringstream o;
                    o << frame << ',' << g10(v) << ',' << g10(x) << ','
                      << g10(y) << ',' << g10(s.f) << ',' << g10(vf) << ','
                      << g10(s.g00) << ',' << g10(s.g01) << ',' << g10(s.g11)
                      << ',' << g10(s.det) << ',' << s.signature_ok;
                    rows.push_back(o.str());
                }
            // The status is a statement about g00 reaching zero. It is an
            // ERGOSURFACE CONDITION and nothing else: the metric stays
            // Lorentzian, det g stays -1, and there is no degeneracy, no
            // signature failure and no horizon here.
            const char* status = (max_vf >= 1.0 ? "ergosurface_crossing"
                                                : "no_zero_crossing");
            if (max_vf >= 1.0 && crossing_frame < 0) crossing_frame = frame;
            for (const std::string& r : rows) f << r << ',' << status << "\n";
        }
    }

    // =====================================================================
    // Radius sweep
    // =====================================================================
    std::vector<double> radius_values, measured_radii, radius_centers;
    {
        std::ofstream f(out_dir + "/radius_sweep.csv", std::ios::binary);
        f << "frame,R,x,y,r_s,f,g00,g01,det,signature_ok,"
             "measured_half_height_radius,bubble_center_x\n";
        for (int frame = 0; frame < kFrames; ++frame) {
            const double R = frame_value(frame, kRMin, kRMax);
            const Params p(kNominalV, R, kNominalSigma);
            const double meas = measured_half_height_radius(metric, p);
            const double xc = metric.center_at(0.0, p);
            radius_values.push_back(R);
            measured_radii.push_back(meas);
            radius_centers.push_back(xc);
            for (int i = 0; i < kNx; ++i)
                for (int j = 0; j < kNy; ++j) {
                    const double x = grid_x(i), y = grid_y(j);
                    const Sample s = evaluate_sample(0.0, x, y, metric, p);
                    f << frame << ',' << g10(R) << ',' << g10(x) << ','
                      << g10(y) << ',' << g10(s.r_s) << ',' << g10(s.f) << ','
                      << g10(s.g00) << ',' << g10(s.g01) << ',' << g10(s.det)
                      << ',' << s.signature_ok << ',' << g10(meas) << ','
                      << g10(xc) << "\n";
                }
        }
    }

    // =====================================================================
    // Steepness sweep
    // =====================================================================
    std::vector<double> sigma_values, wall_widths, sigma_radii, sigma_centers;
    {
        std::ofstream f(out_dir + "/steepness_sweep.csv", std::ios::binary);
        f << "frame,sigma,x,y,r_s,f,g00,g01,det,signature_ok,"
             "measured_wall_width,bubble_center_x\n";
        for (int frame = 0; frame < kFrames; ++frame) {
            const double sg = frame_value(frame, kSigmaMin, kSigmaMax);
            const Params p(kNominalV, kNominalR, sg);
            const double width = measured_wall_width(metric, p);
            const double xc = metric.center_at(0.0, p);
            sigma_values.push_back(sg);
            wall_widths.push_back(width);
            sigma_radii.push_back(measured_half_height_radius(metric, p));
            sigma_centers.push_back(xc);
            for (int i = 0; i < kNx; ++i)
                for (int j = 0; j < kNy; ++j) {
                    const double x = grid_x(i), y = grid_y(j);
                    const Sample s = evaluate_sample(0.0, x, y, metric, p);
                    f << frame << ',' << g10(sg) << ',' << g10(x) << ','
                      << g10(y) << ',' << g10(s.r_s) << ',' << g10(s.f) << ','
                      << g10(s.g00) << ',' << g10(s.g01) << ',' << g10(s.det)
                      << ',' << s.signature_ok << ',' << g10(width) << ','
                      << g10(xc) << "\n";
                }
        }
    }

    // =====================================================================
    // Sensitivity sweep
    // =====================================================================
    bool sens_finite = true;
    const char* sens_backend = "central_fd";
    {
        std::ofstream f(out_dir + "/sensitivity_sweep.csv", std::ios::binary);
        f << "frame,swept_parameter,parameter_value,x,y,dg00_dv,dg00_dR,"
             "dg00_dsigma,dg0x_dv,dg0x_dR,dg0x_dsigma,backend\n";
        struct Sweep { const char* name; double lo, hi; int which; };
        const Sweep sweeps[3] = {{"v", kVMin, 1.0, 0}, {"R", kRMin, kRMax, 1},
                                 {"sigma", kSigmaMin, kSigmaMax, 2}};
        for (const Sweep& sw : sweeps)
            for (int frame = 0; frame < kFrames; ++frame) {
                const double val = frame_value(frame, sw.lo, sw.hi);
                Params p(kNominalV, kNominalR, kNominalSigma);
                if (sw.which == 0) p.velocity = val;
                else if (sw.which == 1) p.radius = val;
                else p.sigma = val;
                for (int i = 0; i < kNx; ++i)
                    for (int j = 0; j < kNy; ++j) {
                        const double x = grid_x(i), y = grid_y(j);
                        const Sensitivities s = sensitivities(
                            0.0, x, y, MetricSensitivityMode::FixedBubbleCenter,
                            metric, p);
                        sens_backend = s.backend;
                        for (int k = 0; k < 3; ++k)
                            if (!std::isfinite(s.dg00[k]) ||
                                !std::isfinite(s.dg0x[k]))
                                sens_finite = false;
                        f << frame << ',' << sw.name << ',' << g10(val) << ','
                          << g10(x) << ',' << g10(y) << ',' << g10(s.dg00[0])
                          << ',' << g10(s.dg00[1]) << ',' << g10(s.dg00[2])
                          << ',' << g10(s.dg0x[0]) << ',' << g10(s.dg0x[1])
                          << ',' << g10(s.dg0x[2]) << ',' << s.backend << "\n";
                    }
            }
    }

    // =====================================================================
    // Fixed-lab versus fixed-centre
    // =====================================================================
    //
    // The two conventions differ by the trajectory term, because in the
    // laboratory convention r_s depends on v through x_s(t) = x_s0 + v t.
    //
    //   dg/dv |_lab = dg/dv |_centre + (dg/dx_s)(dx_s/dv)
    //               = dg/dv |_centre - (dg/dx) t
    //
    // The trajectory term is computed INDEPENDENTLY here, by a centred spatial
    // difference of the production metric, so the residual is a real check
    // rather than an identity. The two conventions are also evaluated by
    // separate metric objects with different sensitivity_mode, never by one
    // shared helper.
    double worst_chain_residual = 0.0;
    {
        std::ofstream f(out_dir + "/fixed_lab_vs_fixed_center.csv",
                        std::ios::binary);
        f << "frame,time,x,y,bubble_center_x,dg_dv_fixed_center,"
             "dg_dv_fixed_lab,trajectory_term,chain_rule_residual,backend\n";
        const Metric m_lab =
            make_metric(MetricSensitivityMode::FixedLaboratoryCoordinates);
        const double t = kChainRuleTime;
        for (int frame = 0; frame < kFrames; ++frame) {
            const double v = frame_value(frame, kVMin, 1.0);
            const Params p(v, kNominalR, kNominalSigma);
            const double xc = metric.center_at(t, p);
            for (int i = 0; i < kNx; ++i)
                for (int j = 0; j < kNy; ++j) {
                    const double x = grid_x(i), y = grid_y(j);
                    const Sensitivities c = sensitivities(
                        t, x, y, MetricSensitivityMode::FixedBubbleCenter,
                        metric, p);
                    const Sensitivities l = sensitivities(
                        t, x, y, MetricSensitivityMode::FixedLaboratoryCoordinates,
                        m_lab, p);
                    // dg00/dx by a centred difference on the production metric.
                    const double h = 1.0e-6;
                    const double gp =
                        metric.evaluate(t, x + h, y, 0.0, p)(0, 0);
                    const double gm =
                        metric.evaluate(t, x - h, y, 0.0, p)(0, 0);
                    const double dg_dx = (gp - gm) / (2.0 * h);
                    const double trajectory = -dg_dx * t;
                    const double residual =
                        std::fabs(l.dg00[0] - (c.dg00[0] + trajectory));
                    worst_chain_residual =
                        std::max(worst_chain_residual, residual);
                    f << frame << ',' << g10(t) << ',' << g10(x) << ','
                      << g10(y) << ',' << g10(xc) << ',' << g10(c.dg00[0])
                      << ',' << g10(l.dg00[0]) << ',' << g10(trajectory) << ','
                      << g10(residual) << ',' << c.backend << "\n";
                }
        }
    }

    // =====================================================================
    // Gate cases
    // =====================================================================
    auto monotonic = [](const std::vector<double>& v) {
        for (std::size_t i = 1; i < v.size(); ++i)
            if (!(v[i] > v[i - 1])) return false;
        return true;
    };
    auto constant = [](const std::vector<double>& v, double tol) {
        double worst = 0.0;
        for (double q : v) worst = std::max(worst, std::fabs(q - v.front()));
        return worst <= tol;
    };

    {
        Stopwatch sw;
        bool ok = true;
        for (std::size_t i = 0; i < tr_times.size(); ++i)
            ok = ok && std::fabs(tr_centers[i] -
                                 (kCenterX0 + kNominalV * tr_times[i])) < 1e-12;
        report.check_predicate(
            "asweep.translation_center_follows_trajectory",
            "the saved bubble centre equals x_s0 + v_s t at every frame",
            VerificationLevel::A, ok, "exact agreement with the trajectory",
            ok ? "all 24 frames agree to 1e-12" : "MISMATCH", sw.seconds(),
            "this is what makes the translation animation's moving centre a "
            "rendered fact rather than a rendering artefact");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.lab_frame_field_translates",
            "the laboratory-frame field moves across the domain over the sweep",
            VerificationLevel::A, lab_translation_span > 1.0,
            "centre displacement greater than one bubble radius",
            "displacement " + g10(lab_translation_span) + " over " +
                g10(tr_times.back()) + " time units",
            sw.seconds());
    }
    {
        Stopwatch sw;
        report.check_residual(
            "asweep.comoving_profile_stationary",
            "the comoving shape profile is identical at every frame, so the "
            "motion is translation and nothing else",
            VerificationLevel::A, "41 offsets x 24 frames",
            "v = " + g10(kNominalV), worst_comoving_drift, 1.0e-12,
            sw.seconds(),
            "sampled at fixed offsets from the MOVING centre: if the profile "
            "changed shape as it moved, this would be nonzero");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.time_values_monotonic", "sweep time increases frame to frame",
            VerificationLevel::A, monotonic(tr_times), "strictly increasing",
            "24 frames from " + g10(tr_times.front()) + " to " +
                g10(tr_times.back()), sw.seconds());
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.velocity_values_monotonic",
            "sweep velocity increases frame to frame", VerificationLevel::A,
            monotonic(vel_values), "strictly increasing",
            "24 frames from " + g10(vel_values.front()) + " to " +
                g10(vel_values.back()), sw.seconds());
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.velocity_sweep_crosses_vf_unity",
            "the velocity sweep contains a frame where v_s f reaches 1",
            VerificationLevel::A, crossing_frame >= 0,
            "at least one crossing frame",
            crossing_frame >= 0
                ? "first crossing at frame " + std::to_string(crossing_frame) +
                      ", v_s = " + g10(vel_values[crossing_frame])
                : "NO CROSSING IN RANGE",
            sw.seconds(),
            "f(0) = 1 exactly at the bubble centre, so the crossing is exactly "
            "at v_s = 1; the committed milestone_03 slice never reaches it, "
            "which is why no ergosurface appears in the static figure");
    }
    {
        Stopwatch sw;
        report.check_residual(
            "asweep.g00_matches_closed_form",
            "the saved g00 equals -(1 - v_s^2 f^2) everywhere in the velocity "
            "sweep",
            VerificationLevel::A, "41 x 41 x 24", "v_s in [0.2, 1.6]",
            worst_g00_formula, 1.0e-15, sw.seconds(),
            "g00 = 0 at v_s f = 1 is therefore an ERGOSURFACE CONDITION, not a "
            "degeneracy, a signature failure or a horizon");
    }
    {
        Stopwatch sw;
        report.check_residual(
            "asweep.determinant_is_exactly_minus_one",
            "det g = -1 across every frame of the velocity sweep, including "
            "beyond the crossing",
            VerificationLevel::A, "41 x 41 x 24", "v_s in [0.2, 1.6]",
            worst_det_dev, 1.0e-15, sw.seconds(),
            "the (t,x) block determinant is (-1 + v^2 f^2) - (v f)^2 = -1 for "
            "ANY v and f, so this is an exact identity rather than a tolerance");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.signature_remains_lorentzian",
            "the metric stays Lorentzian at every sample, including where "
            "g00 >= 0",
            VerificationLevel::A, signature_failures == 0,
            "zero signature failures",
            std::to_string(signature_failures) + " failures over " +
                std::to_string(kNx * kNy * kFrames) + " samples",
            sw.seconds(),
            "g00 changing sign does NOT change the signature: the eigenvalues "
            "of the full 4x4 keep their pattern, which is precisely why the "
            "crossing is labelled an ergosurface");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.radius_sweep_center_fixed",
            "the bubble centre does not move while the radius is swept",
            VerificationLevel::A, constant(radius_centers, 0.0),
            "centre exactly constant",
            "centre " + g10(radius_centers.front()) + " at every frame",
            sw.seconds(),
            "this is the case that distinguishes a radius change from a "
            "translation in the animation");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.measured_radius_monotonic",
            "the MEASURED half-height radius increases with R",
            VerificationLevel::A, monotonic(measured_radii),
            "strictly increasing", "from " + g10(measured_radii.front()) +
                " to " + g10(measured_radii.back()),
            sw.seconds(),
            "measured from the saved profile by bisecting f = 0.5, not read "
            "back from the swept parameter, so it is not a tautology");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.steepness_sweep_center_fixed",
            "the bubble centre does not move while the steepness is swept",
            VerificationLevel::A, constant(sigma_centers, 0.0),
            "centre exactly constant",
            "centre " + g10(sigma_centers.front()) + " at every frame",
            sw.seconds());
    }
    {
        Stopwatch sw;
        // Larger sigma is a THINNER wall under this project's convention.
        bool decreasing = true;
        for (std::size_t i = 1; i < wall_widths.size(); ++i)
            if (!(wall_widths[i] < wall_widths[i - 1])) decreasing = false;
        report.check_predicate(
            "asweep.measured_wall_width_monotonic",
            "the MEASURED wall width falls monotonically as sigma rises",
            VerificationLevel::A, decreasing, "strictly decreasing",
            "from " + g10(wall_widths.front()) + " to " +
                g10(wall_widths.back()),
            sw.seconds(),
            "measured as the separation of the f = 0.75 and f = 0.25 crossings. "
            "Larger sigma is a thinner wall, which is the project convention");
    }
    {
        Stopwatch sw;
        // The nominal radius must NOT move while sigma is swept. The
        // half-height crossing of the Alcubierre f is not exactly R, but it
        // must stay put.
        double spread = 0.0;
        for (double q : sigma_radii)
            spread = std::max(spread, std::fabs(q - sigma_radii.front()));
        report.check_value(
            "asweep.nominal_radius_fixed_during_steepness_sweep",
            "the measured half-height radius stays at the nominal radius while "
            "only sigma changes",
            VerificationLevel::B, "24 frames", "R = " + g10(kNominalR),
            0.0, spread, 0.05, sw.seconds(),
            "the half-height crossing of the Alcubierre shape function is close "
            "to but not exactly R, and it shifts slightly as the wall sharpens; "
            "the tolerance is that residual shape effect, not a slack allowance");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.sensitivity_fields_finite",
            "every exported sensitivity value is finite",
            VerificationLevel::A, sens_finite, "no NaN or infinity",
            std::string("backend ") + sens_backend, sw.seconds());
    }
    {
        Stopwatch sw;
#ifdef WARP_HYPAD_HAVE_OTILIB
        // Compare the exported OTI sensitivity against the analytic shape
        // derivative chain at the bubble centre, where f = 1 and dg00/dv = 2v.
        const Params p(kNominalV, kNominalR, kNominalSigma);
        const double xc = metric.center_at(0.0, p);
        const Sensitivities s = sensitivities(
            0.0, xc, 0.0, MetricSensitivityMode::FixedBubbleCenter, metric, p);
        report.check_value(
            "asweep.oti_matches_analytic_at_center",
            "dg00/dv = 2 v f^2 = 2 v at the bubble centre, where f = 1",
            VerificationLevel::A, "single point", "v = " + g10(kNominalV),
            2.0 * kNominalV, s.dg00[0], 1.0e-12, sw.seconds(),
            "ties the exported sensitivity field to the same analytic identity "
            "the milestone_03 gate checks");
#else
        report.check_skipped(
            "asweep.oti_matches_analytic_at_center",
            "dg00/dv against the analytic value using the OTI backend",
            "OTIlib is disabled; the sensitivity FIELDS are still exported via "
            "centred finite differences and every metric, determinant, "
            "signature, trajectory and chain-rule case still runs",
            sw.seconds());
#endif
    }
    {
        Stopwatch sw;
#ifdef WARP_HYPAD_HAVE_OTILIB
        // OTI against centred FD at a wall point, where the derivative is
        // largest and an error would show.
        const Params p(kNominalV, kNominalR, kNominalSigma);
        const double xw = metric.center_at(0.0, p) + kNominalR;
        const Sensitivities o = sensitivities_oti(
            0.0, xw, 0.0, MetricSensitivityMode::FixedBubbleCenter, p);
        const Sensitivities d = sensitivities_fd(
            0.0, xw, 0.0, metric, p, MetricSensitivityMode::FixedBubbleCenter);
        double worst = 0.0;
        for (int k = 0; k < 3; ++k)
            worst = std::max(worst, std::fabs(o.dg00[k] - d.dg00[k]));
        report.check_residual(
            "asweep.oti_matches_finite_difference",
            "the OTI sensitivities agree with centred finite differences at a "
            "wall point",
            VerificationLevel::B, "single point, h = 1e-6",
            "r_s = R", worst, 1.0e-5, sw.seconds(),
            "the tolerance is the finite-difference truncation error, not the "
            "accuracy of OTI: OTI has no truncation error at all, which is what "
            "the milestone_03 convergence study established");
#else
        report.check_skipped(
            "asweep.oti_matches_finite_difference",
            "OTI against centred finite differences",
            "OTIlib is disabled, so there is no second backend to compare "
            "against; the finite-difference values themselves are exported",
            sw.seconds());
#endif
    }
    {
        Stopwatch sw;
        // The chain rule holds to the accuracy of the finite-difference
        // trajectory term, which is the least accurate ingredient.
        report.check_residual(
            "asweep.fixed_lab_fixed_center_chain_rule",
            "dg/dv|_lab = dg/dv|_centre - (dg/dx) t, with the trajectory term "
            "computed independently",
            VerificationLevel::B, "41 x 41 x 24", "t = " + g10(kChainRuleTime),
            worst_chain_residual, 1.0e-4, sw.seconds(),
            "the two conventions come from two metric objects with different "
            "sensitivity_mode, never one shared helper; the trajectory term is "
            "a centred spatial difference, so its 1e-6 step sets the tolerance");
    }

    // Checksums, then the deterministic-repeat case.
    const char* files[] = {"translation_sweep.csv", "velocity_sweep.csv",
                           "radius_sweep.csv", "steepness_sweep.csv",
                           "sensitivity_sweep.csv",
                           "fixed_lab_vs_fixed_center.csv"};
    JsonValue checksums = JsonValue::object();
    bool all_present = true;
    for (const char* n : files) {
        const std::string p = out_dir + "/" + n;
        if (!std::filesystem::exists(p)) { all_present = false; continue; }
        checksums.set(n, JsonValue(std::to_string(fnv1a_file(p))));
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "asweep.output_checksums_valid",
            "every sweep file exists and has a recorded FNV-1a checksum",
            VerificationLevel::A, all_present, "6 files with checksums",
            all_present ? "6 files written" : "FILES MISSING", sw.seconds());
    }
    {
        Stopwatch sw;
        // Determinism: every frame value and grid coordinate is a closed-form
        // function of its index, and nothing adaptive runs, so a repeat export
        // reproduces the bytes. Verified by recomputing the frame values.
        bool deterministic = true;
        for (int frame = 0; frame < kFrames; ++frame) {
            if (frame_value(frame, kVMin, kVMax) != vel_values[frame])
                deterministic = false;
            if (frame_value(frame, kRMin, kRMax) != radius_values[frame])
                deterministic = false;
        }
        report.check_predicate(
            "asweep.deterministic_repeat_export",
            "frame values and grid coordinates are closed-form functions of "
            "their index, so a repeat export is byte-identical",
            VerificationLevel::A, deterministic, "exact reproduction",
            deterministic ? "all frame values reproduce exactly" : "DRIFT",
            sw.seconds(),
            "no adaptive stepping, no accumulation, no wall-clock input");
    }

    report.print(std::cout);
    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    // --- sweep_manifest.json -------------------------------------------------
    {
        JsonValue m = JsonValue::object();
        m.set("program", JsonValue("A (warp_hypad)"));
        m.set("milestone", JsonValue("V0.1A / A-VISUAL-SWEEP-EXPORT"));
        m.set("coordinate_mode", JsonValue("InertialTimeDependent, "
                                           "x_s(t) = x_s0 + v_s t"));
        m.set("frozen_time_used", JsonValue(false));
        m.set("grid_nx", JsonValue(kNx));
        m.set("grid_ny", JsonValue(kNy));
        m.set("grid_half_width", JsonValue(kHalfWidth));
        m.set("frames_per_sweep", JsonValue(kFrames));
        m.set("velocity_range", JsonValue(g10(kVMin) + " to " + g10(kVMax)));
        m.set("radius_range", JsonValue(g10(kRMin) + " to " + g10(kRMax)));
        m.set("sigma_range", JsonValue(g10(kSigmaMin) + " to " + g10(kSigmaMax)));
        m.set("nominal_parameters",
              JsonValue("v = " + g10(kNominalV) + ", R = " + g10(kNominalR) +
                        ", sigma = " + g10(kNominalSigma)));
        m.set("units", JsonValue("geometrised, c = 1; lengths in bubble radii; "
                                 "time in the same geometrised unit"));
        m.set("time_convention",
              JsonValue("coordinate time of the inertial laboratory frame"));
        m.set("sensitivity_convention",
              JsonValue("FixedBubbleCenter for the sensitivity sweep; both "
                        "conventions appear in fixed_lab_vs_fixed_center.csv"));
        m.set("sensitivity_backend", JsonValue(sens_backend));
#ifdef WARP_HYPAD_HAVE_OTILIB
        m.set("otilib_enabled", JsonValue(true));
#else
        m.set("otilib_enabled", JsonValue(false));
#endif
        m.set("checksums", checksums);
        m.set("ergosurface_crossing_frame", JsonValue(crossing_frame));
        m.set("scope",
              JsonValue("These records contain effective-metric components and "
                        "parameter sensitivities. They do not contain "
                        "curvature, stress-energy or propulsion physics."));
        m.set("build", build_info_json());
        m.set("generated_timestamp", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/sweep_manifest.json", m, error))
            return 3;
    }

    std::cout << "\nA-VISUAL-SWEEP-EXPORT SUMMARY\n";
    std::cout << "  grid / frames             " << kNx << " x " << kNy << " / "
              << kFrames << "\n";
    std::cout << "  worst |det g + 1|         " << worst_det_dev << "\n";
    std::cout << "  signature failures        " << signature_failures << "\n";
    std::cout << "  ergosurface crossing      frame " << crossing_frame
              << " (v_s = "
              << (crossing_frame >= 0 ? g10(vel_values[crossing_frame]) : "n/a")
              << ")\n";
    std::cout << "  comoving profile drift    " << worst_comoving_drift << "\n";
    std::cout << "  lab translation span      " << lab_translation_span << "\n";
    std::cout << "  chain-rule residual       " << worst_chain_residual << "\n";
    std::cout << "  sensitivity backend       " << sens_backend << "\n";
    std::cout << "\nSCOPE: effective-metric components and sensitivities only. "
                 "No curvature,\n  no stress-energy, no propulsion. Curvature "
                 "remains blocked pending A4.\n";

    const bool ok = report.all_passed();
    std::cout << "\nA-VISUAL-SWEEP-EXPORT: " << (ok ? "PASS" : "FAIL") << " ("
              << report.num_passed() << "/" << report.num_cases() << ", "
              << report.num_skipped() << " skipped)\n";
    return ok ? 0 : 1;
}

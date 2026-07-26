#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include "otilib_checks.hpp"
#include "warp_hypad/sensitivity/otilib_traits.hpp"

// ---------------------------------------------------------------------------
// OTI versus centered finite differences.
//
// This is a Level D check: an independently implemented numerical method. It is
// weaker than the analytic comparisons, because centered FD carries its own
// truncation and cancellation error -- which is precisely the point. The
// expected picture is a V-shaped error curve in the step size h: truncation
// error O(h^2) falling as h shrinks, then roundoff error O(eps/h) rising again.
// The best achievable agreement is around 1e-11 near h ~ 1e-5, NOT machine
// precision.
//
// If OTI agreed with FD to 1e-16 at every step size, that would indicate the FD
// reference was accidentally computed from the same code path, not that both
// were right. So the check asserts a floor on the best agreement and does not
// demand more.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using sensitivity::extract_gradient;
using sensitivity::extract_hessian;
using sensitivity::OTIScalar;

constexpr double kA = 1.30;
constexpr double kB = 0.70;
constexpr double kC = 0.45;

template <typename S>
S composite(const S& a, const S& b, const S& c) {
    using std::exp;
    using std::sin;
    const S a2 = a * a;
    const S sb = sin(b);
    const S t1 = a2 * sb;
    const S ac = a * c;
    const S t2 = exp(ac);
    return t1 + t2;
}

double eval_real(double a, double b, double c) { return composite<double>(a, b, c); }

double fd_first(int direction, double h) {
    double p[3] = {kA, kB, kC};
    double lo[3] = {kA, kB, kC};
    double hi[3] = {kA, kB, kC};
    (void)p;
    lo[direction] -= h;
    hi[direction] += h;
    return (eval_real(hi[0], hi[1], hi[2]) - eval_real(lo[0], lo[1], lo[2])) /
           (2.0 * h);
}

double fd_second(int direction, double h) {
    double lo[3] = {kA, kB, kC};
    double mid[3] = {kA, kB, kC};
    double hi[3] = {kA, kB, kC};
    lo[direction] -= h;
    hi[direction] += h;
    return (eval_real(hi[0], hi[1], hi[2]) - 2.0 * eval_real(mid[0], mid[1], mid[2]) +
            eval_real(lo[0], lo[1], lo[2])) /
           (h * h);
}

}  // namespace

void run_otilib_fd_checks(VerificationReport& report, const std::string& out_dir) {
    const sensitivity::SeededParameters p = sensitivity::seed_all(kA, kB, kC);
    const OTIScalar f = composite<OTIScalar>(p.velocity, p.radius, p.sigma);

    const double steps[7] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8};
    const char* names[3] = {"velocity_slot(a)", "radius_slot(b)", "sigma_slot(c)"};

    std::ofstream os(out_dir + "/finite_difference_comparison.csv", std::ios::binary);
    if (os) {
        os << "parameter,derivative_order,step_size,oti,finite_difference,"
              "absolute_error,relative_error\n";
        os.setf(std::ios::scientific);
        os.precision(12);
    }

    double best_first = 1e300;
    double best_second = 1e300;

    for (int d = 0; d < 3; ++d) {
        const double oti_g = extract_gradient(f, d);
        const double oti_h = extract_hessian(f, d, d);

        for (double h : steps) {
            const double g = fd_first(d, h);
            const double hh = fd_second(d, h);

            const double ge = std::fabs(g - oti_g);
            const double he = std::fabs(hh - oti_h);
            const double grel = ge / std::fmax(1e-300, std::fabs(oti_g));
            const double hrel = he / std::fmax(1e-300, std::fabs(oti_h));

            best_first = std::fmin(best_first, grel);
            best_second = std::fmin(best_second, hrel);

            if (os) {
                os << names[d] << ",1," << h << ',' << oti_g << ',' << g << ','
                   << ge << ',' << grel << '\n';
                os << names[d] << ",2," << h << ',' << oti_h << ',' << hh << ','
                   << he << ',' << hrel << '\n';
            }
        }
    }

    {
        Stopwatch sw;
        std::ostringstream note;
        note << "7 step sizes from 1e-1 to 1e-8; best relative agreement "
             << best_first;
        report.check_residual(
            "oti.vs_finite_difference.first",
            "the best centered-FD first derivative agrees with OTI",
            VerificationLevel::D, "n/a", "a=1.3, b=0.7, c=0.45", best_first,
            1e-9, sw.seconds(), note.str());
    }
    {
        Stopwatch sw;
        std::ostringstream note;
        note << "second differences are far noisier (O(eps/h^2) roundoff); best "
                "relative agreement "
             << best_second;
        report.check_residual(
            "oti.vs_finite_difference.second",
            "the best centered-FD second derivative agrees with OTI",
            VerificationLevel::D, "n/a", "a=1.3, b=0.7, c=0.45", best_second,
            1e-5, sw.seconds(), note.str());
    }

    // The V-shape itself is evidence the two methods are independent. If the
    // error did not grow again at very small h, FD would not be behaving like
    // FD, and the comparison would prove nothing.
    {
        Stopwatch sw;
        const double coarse = std::fabs(fd_first(0, 1e-1) - extract_gradient(f, 0));
        const double sweet = std::fabs(fd_first(0, 1e-5) - extract_gradient(f, 0));
        const double tiny = std::fabs(fd_first(0, 1e-12) - extract_gradient(f, 0));
        const bool v_shaped = (sweet < coarse) && (sweet < tiny);
        std::ostringstream got;
        got.setf(std::ios::scientific);
        got.precision(3);
        got << "h=1e-1: " << coarse << ", h=1e-5: " << sweet << ", h=1e-12: " << tiny;
        report.check_predicate(
            "oti.fd_error_is_v_shaped",
            "FD error falls then rises as h shrinks (truncation then roundoff), "
            "confirming FD is a genuinely independent method",
            VerificationLevel::D, v_shaped, "error at h=1e-5 below both ends",
            got.str(), sw.seconds(),
            "monotone agreement would suggest the reference shares a code path "
            "with OTI");
    }
}

}  // namespace checks
}  // namespace warp_hypad

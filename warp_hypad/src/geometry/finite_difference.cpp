#include "warp_hypad/geometry/finite_difference.hpp"

#include <sstream>

namespace warp_hypad {
namespace geometry {

std::string validate_for_derivatives(const Grid& g, bool need_second) {
    const int minimum = need_second ? kMinPointsSecond : kMinPointsFirst;
    const char* what = need_second ? "second" : "first";

    const int n[3] = {g.nx(), g.ny(), g.nz()};
    const char* names[3] = {"nx", "ny", "nz"};

    for (int a = 0; a < 3; ++a) {
        if (n[a] < minimum) {
            std::ostringstream os;
            os << "grid." << names[a] << " = " << n[a] << " is too small for "
               << what << " derivatives at second order; the one-sided boundary "
               << "stencil needs " << minimum << " points along each "
               << "differentiated axis";
            return os.str();
        }
    }

    const double h[3] = {g.dx(), g.dy(), g.dz()};
    for (int a = 0; a < 3; ++a) {
        if (!(h[a] > 0.0)) {
            std::ostringstream os;
            os << "spacing along " << names[a]
               << " is zero; a degenerate axis cannot be differentiated";
            return os.str();
        }
    }
    return std::string();
}

}  // namespace geometry
}  // namespace warp_hypad

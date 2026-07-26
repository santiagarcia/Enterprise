#include "warp_hypad/sensitivity/otilib_seed.hpp"

namespace warp_hypad {
namespace sensitivity {

const char* parameter_name(int direction) {
    switch (direction) {
        case kParamVelocity: return "velocity";
        case kParamRadius: return "radius";
        case kParamSigma: return "sigma";
        default: return "unknown";
    }
}

}  // namespace sensitivity
}  // namespace warp_hypad

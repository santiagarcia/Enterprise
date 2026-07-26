/* Scalar-only translation unit for OTIlib's onumm3n2 static type.
 *
 * OTIlib is a unity-build project: the files under src/c/static/onumm3n2/ have
 * no includes of their own and are meant to be #included into one translation
 * unit. Upstream's src/c/static/onumm3n2.c does that for the scalar AND array
 * paths; the array path pulls in oarrm3n2_* symbols that live in other modules
 * and would drag in most of the library.
 *
 * warp_hypad has its own Mat4/Rank3/Rank4 containers and needs only the OTI
 * SCALAR type, so this TU includes the scalar half only. That is what makes the
 * dependency small enough to build natively with MSVC.
 *
 * Deliberately excluded:
 *   - onumm3n2/array/*        (needs oarrm3n2_*, dense array module)
 *   - onumm3n2/scalar/gauss.c (FEM Gauss-point support, unused here)
 *
 * Upstream states OTIlib on Windows works "only under WSL". That is true of the
 * full CMake/conda build, which also runs `make gendata` to produce precomputed
 * sparse-direction tables. The static scalar types carry no such dependency --
 * verified by grep for dHelp/precompute across the static sources -- which is
 * why this subset compiles and links with cl.exe unmodified. See
 * docs/otilib_discovery_report.md.
 *
 * Licence: OTIlib is GPL-3.0. See external/otilib/LICENSE and the licence
 * section of the discovery report before distributing anything linked to this.
 */

#include "oti/static/onumm3n2.h"

#include "static/onumm3n2/scalar/base.c"
#include "static/onumm3n2/scalar/algebra.c"
#include "static/onumm3n2/scalar/algebra_to.c"
#include "static/onumm3n2/scalar/functions.c"
#include "static/onumm3n2/scalar/functions_to.c"
#include "static/onumm3n2/scalar/utils.c"

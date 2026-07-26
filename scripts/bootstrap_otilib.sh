#!/usr/bin/env bash
#
# Fetch the pinned OTIlib checkout.
#
# OTIlib is a PINNED EXTERNAL DEPENDENCY, not vendored source. GPL-3.0 is
# accepted for this research prototype; keeping the checkout external is an
# architectural choice, not an unresolved licensing question. It keeps a large
# upstream history out of this repository while the backend stays optional and
# isolated behind include/warp_hypad/sensitivity/.
#
# This script does NOT modify or patch upstream source. It clones, checks out the
# recorded revision, verifies the resulting SHA, and exits nonzero on any
# mismatch. CMake performs the same SHA check at configure time, so a drifted
# checkout cannot silently produce a green sensitivity gate.
#
#   ./scripts/bootstrap_otilib.sh
#   cmake -S warp_hypad -B warp_hypad/build -DWARP_HYPAD_WITH_OTILIB=ON

set -euo pipefail

OTILIB_URL="https://github.com/mauriaristi/otilib.git"
OTILIB_PIN="de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
dest="${repo_root}/warp_hypad/external/otilib"

fail() { printf 'bootstrap_otilib: %s\n' "$1" >&2; exit 1; }

command -v git >/dev/null 2>&1 || fail "git not found on PATH"

if [ -e "${dest}/.git" ]; then
    printf 'bootstrap_otilib: existing checkout at %s\n' "${dest}"
    # Never discard local state silently. If the working tree is dirty, say so
    # and stop: an unexpected modification to a pinned dependency is exactly the
    # thing this script exists to surface.
    if [ -n "$(git -C "${dest}" status --porcelain)" ]; then
        fail "the existing checkout has local modifications. Inspect it, then
       remove ${dest} and rerun to fetch a clean copy."
    fi
    git -C "${dest}" fetch --depth 1 origin "${OTILIB_PIN}" 2>/dev/null || \
        git -C "${dest}" fetch origin
    git -C "${dest}" checkout --detach "${OTILIB_PIN}"
else
    printf 'bootstrap_otilib: cloning %s\n' "${OTILIB_URL}"
    mkdir -p "$(dirname "${dest}")"
    # Shallow, then deepen only if the pinned commit is not reachable.
    git clone --no-checkout "${OTILIB_URL}" "${dest}"
    git -C "${dest}" checkout --detach "${OTILIB_PIN}"
fi

actual="$(git -C "${dest}" rev-parse HEAD)"
if [ "${actual}" != "${OTILIB_PIN}" ]; then
    fail "SHA mismatch after checkout.
       expected ${OTILIB_PIN}
       found    ${actual}"
fi

# The one file the build actually needs. If upstream ever moves it, fail loudly
# here rather than at link time.
header="${dest}/include/oti/static/onumm3n2.h"
[ -f "${header}" ] || fail "checkout verified at ${OTILIB_PIN} but ${header} is missing"

printf 'bootstrap_otilib: OK\n'
printf '  path %s\n' "${dest}"
printf '  pin  %s\n' "${actual}"
printf '\nEnable the backend with:\n'
printf '  cmake -S warp_hypad -B warp_hypad/build -DWARP_HYPAD_WITH_OTILIB=ON\n'

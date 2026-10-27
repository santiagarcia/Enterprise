<#
.SYNOPSIS
    Fetch the pinned OTIlib checkout.

.DESCRIPTION
    OTIlib is GPL-3.0, accepted for this research prototype. Since commit
    bd262fb the working tree at the pinned revision IS VENDORED in
    warp_hypad/external/otilib/ (see VENDORING_NOTE.md there), so a normal clone
    does not need this script. It remains for fetching a fresh upstream checkout
    with git history. The backend stays optional and isolated behind
    include/warp_hypad/sensitivity/.

    This script does NOT modify or patch upstream source. It clones, checks out
    the recorded revision, verifies the resulting SHA, and exits nonzero on any
    mismatch. CMake performs the same SHA check at configure time, so a drifted
    checkout cannot silently produce a green sensitivity gate.

.EXAMPLE
    .\scripts\bootstrap_otilib.ps1
    cmake -S warp_hypad -B warp_hypad/build -DWARP_HYPAD_WITH_OTILIB=ON
#>

$ErrorActionPreference = 'Stop'

$OtilibUrl = 'https://github.com/mauriaristi/otilib.git'
$OtilibPin = 'de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54'

function Fail($message) {
    Write-Error "bootstrap_otilib: $message"
    exit 1
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot  = Split-Path -Parent $scriptDir
$dest      = Join-Path $repoRoot 'warp_hypad\external\otilib'

if ($null -eq (Get-Command git -ErrorAction SilentlyContinue)) {
    Fail 'git not found on PATH'
}

if (Test-Path (Join-Path $dest '.git')) {
    Write-Host "bootstrap_otilib: existing checkout at $dest"

    # Never discard local state silently. An unexpected modification to a pinned
    # dependency is exactly the thing this script exists to surface.
    $dirty = & git -C $dest status --porcelain
    if ($dirty) {
        Fail "the existing checkout has local modifications. Inspect it, then remove $dest and rerun to fetch a clean copy."
    }

    & git -C $dest fetch --depth 1 origin $OtilibPin 2>$null
    if (-not $?) { & git -C $dest fetch origin }
    & git -C $dest checkout --detach $OtilibPin
    if (-not $?) { Fail "could not check out $OtilibPin" }
}
else {
    Write-Host "bootstrap_otilib: cloning $OtilibUrl"
    $parent = Split-Path -Parent $dest
    if (-not (Test-Path $parent)) { New-Item -ItemType Directory -Force $parent | Out-Null }

    & git clone --no-checkout $OtilibUrl $dest
    if (-not $?) { Fail 'clone failed' }
    & git -C $dest checkout --detach $OtilibPin
    if (-not $?) { Fail "could not check out $OtilibPin" }
}

$actual = (& git -C $dest rev-parse HEAD).Trim()
if ($actual -ne $OtilibPin) {
    Fail "SHA mismatch after checkout.`n       expected $OtilibPin`n       found    $actual"
}

# The one file the build actually needs. If upstream ever moves it, fail loudly
# here rather than at link time.
$header = Join-Path $dest 'include\oti\static\onumm3n2.h'
if (-not (Test-Path $header)) {
    Fail "checkout verified at $OtilibPin but $header is missing"
}

Write-Host 'bootstrap_otilib: OK'
Write-Host "  path $dest"
Write-Host "  pin  $actual"
Write-Host ''
Write-Host 'Enable the backend with:'
Write-Host '  cmake -S warp_hypad -B warp_hypad/build -DWARP_HYPAD_WITH_OTILIB=ON'

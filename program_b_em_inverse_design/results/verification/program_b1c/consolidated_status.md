# B1C consolidated status

B1C is covered by three gates. This one holds only the cross-cutting cases; the other two are separate binaries and separate CTest tests, and this gate DEPENDS on them rather than repeating their assertions.

| gate | binary | CTest test | record |
|---|---|---|---|
| B1C-S23 | `b1c_s23_verify` | `programb1c_s23.gate` | `results/verification/program_b1c_s23/` |
| B1C-S456 | `b1c_s456_verify` | `programb1c_s456.gate` | `results/verification/program_b1c_s456/` |
| B1C aggregate | `b1c_verify` | `programb1c.gate` | `results/verification/program_b1c/` |

## This gate

- cases: 53
- passed: 53
- failed: 0
- skipped: 0

## Velocity sensitivity

- worst OTI versus analytic: 2.220446e-16
- worst analytic versus centred finite differences: 6.869283e-11
- worst chain-rule residual: 1.110223e-16

## Withdrawn assertions

See `docs/b1c_aggregate_gate_audit.md` for the per-case record of what was withdrawn, why, and what supersedes it. No active gate includes `legacy_unverified.hpp`.

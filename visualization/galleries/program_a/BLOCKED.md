# Blocked V0-A products

`ricci_scalar_slice`, `curvature_norm_slice` and
`curvature_wall_evolution` are NOT generated.

Program A contains no curvature code. `warp_hypad` computes
metric components, determinants, signature masks and parameter
sensitivities; the word 'curvature' appears in it only as a
guard forbidding the frozen-time mode from being fed to a
Christoffel routine that does not yet exist.

Computing a Ricci scalar inside this plotting script would be
exactly the failure the visualization contract prevents: a
field invented in Python, rendered as a figure, and thereafter
believed. It would look identical to a verified one.

Unblocking requires a Program A scientific milestone:
Christoffel symbols from the verified finite-difference
operators, Riemann and Ricci assembly, and verification
against known answers -- flat space exactly zero, Schwarzschild
R = 0 with a nonzero Kretschmann invariant, and a convergence
study. See visualization/V0_ARCHITECTURE.md.

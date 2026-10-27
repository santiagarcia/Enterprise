# Enterprise

Research prototypes for warp-spacetime analysis and related programs:

- [`warp_hypad/`](warp_hypad/README.md): C++ metric, curvature and parameter
  sensitivity pipeline with hypercomplex automatic differentiation.
- [`program_b_em_inverse_design/`](program_b_em_inverse_design/README.md):
  electromagnetic inverse design.
- `program_c_astrogation/`: relativistic astrogation using Gaia DR3 snapshots.
- `visualization/`, `scripts/`: figures, animations and helper scripts.

## License

- **First-party code** in this repository is licensed under the GNU General
  Public License, version 3 only (SPDX: `GPL-3.0-only`). The full text is in
  [`LICENSE`](LICENSE).
- **`warp_hypad/external/otilib/`** is a vendored copy of upstream
  [OTIlib](https://github.com/mauriaristi/otilib) at commit
  `de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54`, distributed under its own
  upstream GPL-3.0 licence (`warp_hypad/external/otilib/LICENSE`). Differences
  from upstream are listed in `warp_hypad/external/otilib/VENDORING_NOTE.md`.
  Some upstream build files carry their own notices (for example the Kitware
  CMake modules, Apache License 2.0).
- **Third-party data** are **not** covered by the code licence: the Gaia DR3
  snapshots in `program_c_astrogation/snapshots/` (with SIMBAD identifiers and
  NASA Exoplanet Archive records), and the outputs derived from them under
  `program_c_astrogation/results/verification/` (`program_c0_real/`,
  `program_c1/`, and the C0 Gaia demonstration file). Gaia data are
  CC BY-NC 3.0 IGO. See
  [`program_c_astrogation/snapshots/DATA_LICENSES.md`](program_c_astrogation/snapshots/DATA_LICENSES.md)
  for the terms, the full list of affected files and the required
  acknowledgements.

## AI assistance

Some commits carry a `Co-Authored-By` trailer that names an AI
assistant. All 45 commits up to `477141d` (2026-08-15) carry one naming Claude
(Anthropic): 38 name "Claude Opus 5" and 7 name "Claude Opus 5 (1M context)".
The trailers record that an AI coding assistant contributed to those commits;
they do not identify which lines it wrote or how the changes were checked.

> **[AUTHOR STATEMENT REQUIRED: the author's own statement of how the
> AI-assisted work was reviewed and verified, and of responsibility for this
> repository's contents. To be written by the author; not filled in by the AI
> assistant.]**

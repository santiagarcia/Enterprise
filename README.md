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
- **Data snapshots** in `program_c_astrogation/snapshots/` (Gaia DR3, with
  SIMBAD identifiers and NASA Exoplanet Archive records) are third-party data
  under their own terms and are **not** covered by the code licence. Gaia data
  are CC BY-NC 3.0 IGO. See
  [`program_c_astrogation/snapshots/DATA_LICENSES.md`](program_c_astrogation/snapshots/DATA_LICENSES.md)
  for the terms and the required acknowledgements.

## AI assistance

All 45 commits up to `477141d` (2026-08-15) carry a `Co-Authored-By` trailer
naming Claude (Anthropic): 38 name "Claude Opus 5" and 7 name "Claude Opus 5
(1M context)". The trailers record that an AI coding assistant contributed to
those commits; they do not identify which lines it wrote. The human author
reviewed the changes and is responsible for the contents of this repository.

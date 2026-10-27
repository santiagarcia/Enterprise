# Vendored copy of OTIlib: provenance and differences from upstream

This directory is a copy of the working tree of upstream OTIlib
(https://github.com/mauriaristi/otilib) at commit
`de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54` (2025-08-25). It was committed to
this repository in `bd262fb` (2026-07-28). The upstream git history is not
included.

OTIlib is licensed under the GNU General Public License, version 3. The licence
text is in `LICENSE` in this directory. Upstream copyright and licence notices
are kept unchanged. A few build files carry their own notices (for example
`cmake/UseCython.cmake` and `cmake/FindCython.cmake`: Copyright 2011 Kitware,
Inc., Apache License 2.0).

`PINNED_SHA.txt` holds only the pinned commit. `warp_hypad/CMakeLists.txt`
reads that file and compares its whole content with the pin, so this note is
kept in a separate file.

## Differences from the upstream tree at the pin

Compared on 2026-10-27, path by path and blob by blob, against the upstream
tree `8f0728c1dd65714054c7bf4dc0aa1838ec1689ac` (the tree of `de3f7fa`):

- **Files added here:** `PINNED_SHA.txt` and this note (`VENDORING_NOTE.md`).
- **Files omitted (14 of 2,737):** 6 Jupyter checkpoint copies and 8 compiled
  build artefacts present in the upstream tree. No upstream source file is
  omitted.
  - `examples/notebooks/.ipynb_checkpoints/core_tutorial_01-checkpoint.ipynb`
  - `examples/notebooks/.ipynb_checkpoints/sparse_tutorial_01-checkpoint.ipynb`
  - `examples/notebooks/.ipynb_checkpoints/sparse_tutorial_02-checkpoint.ipynb`
  - `examples/notebooks/.ipynb_checkpoints/sparse_tutorial_03-checkpoint.ipynb`
  - `examples/notebooks/.ipynb_checkpoints/static_tutorial_01-checkpoint.ipynb`
  - `examples/notebooks/.ipynb_checkpoints/static_tutorial_02-checkpoint.ipynb`
  - `src/cpp/a.o`
  - `src/fortran/core/libmaster_parameters.a`
  - `src/fortran/core/libreal_utils.a`
  - `src/fortran/core/master_parameters.o`
  - `src/fortran/core/real_utils.o`
  - `src/fortran/sparse/libsparse_oti.a`
  - `src/fortran/sparse/sparse_oti.o`
  - `src/fortran/tests/test.exe`
- **Line endings (41 files):** these upstream files use CRLF line endings. This
  repository stores them with LF (`.gitattributes`: `* text=auto eol=lf`).
  After removing carriage returns, all 41 are byte-identical to upstream
  (22 `.f90`, 5 `.py`, 4 `.rst`, 3 `.c`, 2 `.txt`, 2 `Makefile`, 1 `.pxi`,
  1 `.h`, 1 `.bat`).
- **All other 2,682 files** are byte-identical to upstream.

No other change was made to the upstream content.

Reproduce (from the repository root, with an upstream clone in `UP`):

```bash
git -C UP ls-tree -r de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54
git ls-tree -r HEAD:warp_hypad/external/otilib
# compare the blob ids per path; for differing paths compare
# `git cat-file blob <id> | tr -d '\r'` on both sides.
```

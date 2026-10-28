# study3: files kept out of the repository

These build products were left out because they are large compiled binaries
(object files and Windows executables, 3-11 MB each) that are rebuilt from the
committed driver sources by driver/build.sh. Their SHA-256 hashes are listed so
a rebuilt binary can be compared with the one that produced the evidence.

Every other file of the study (sources, scripts, logs, CSV/JSON results,
assessments, final fields, figures, reports
and reviews) is committed byte-identical; see .gitattributes.

| File | Size (bytes) | SHA-256 |
|---|---:|---|
| driver/o0d.o | 3129987 | 8279994b804bf16557d1302ba36c7f82681a9ea93611e4848c83db6fe5be7c49 |
| driver/s3.exe | 9302402 | 366e9c721619e9255b411e81b2704d5a99a40de6df1946e1e6358a1c2c5881ec |
| driver/s3.o | 4438126 | f9a1288abb42f3f88aaf9fd35d058fb1b826e22f0d1a29aee97e31f8f00d6ba9 |

Also omitted: Python bytecode caches (__pycache__/) and the ruff lint cache
(.ruff_cache/), which are tool caches rather than evidence.

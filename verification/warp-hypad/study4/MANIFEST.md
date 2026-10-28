# study4: files kept out of the repository

These build products were left out because they are large compiled binaries
(object files and Windows executables, 3-11 MB each) that are rebuilt from the
committed driver sources by driver/build.sh. Their SHA-256 hashes are listed so
a rebuilt binary can be compared with the one that produced the evidence.

Every other file of the study (sources, scripts, logs, CSV/JSON results,
assessments, Hessian summaries and Hessian .bin/.out files, figures, reports
and reviews) is committed byte-identical; see .gitattributes.

| File | Size (bytes) | SHA-256 |
|---|---:|---|
| driver/o0d.o | 3129987 | 049ae9ee8e33b77696859b2040765a7845316e6b06c7203e88fc77442d93a638 |
| driver/s4.exe | 10733542 | b6c0ad75c6fc4b32d88c8ab63f05edabec77e399ec520b2c0e137a0f38dd52ff |
| driver/s4.o | 6572365 | 4f961d063a4fd296f4c00d3a1a179015685fe6bc2c6a820328335c786bc28951 |
| driver/v2_s3harness.exe | 9412523 | 01f863e4d420794519e8b82268b4276d0c24d46a5391ea4dd58b693628cbfb3a |
| driver/v2h.o | 4605685 | e8f1c96ad8e0dcbb710ea43baa704e9620c28d2dd1ce74a2d30e18ca634d1651 |

Also omitted: Python bytecode caches (__pycache__/) and the ruff lint cache
(.ruff_cache/), which are tool caches rather than evidence.

Also omitted: 37 empty (0-byte) files in assess/snap/ whose names end in a
carriage return (shown on Windows as y_<run>_it<k><U+F00D>.txt). They were left
by a first run of driver/assess.sh with CRLF line endings; each has a correctly
named, non-empty twin (y_<run>_it<k>.txt) that is committed and that
assess/iters.csv references.

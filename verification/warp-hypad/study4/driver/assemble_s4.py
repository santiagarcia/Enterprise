"""Assemble s4.cpp from verbatim line ranges of study3/driver/s3.cpp and the study-4 pieces.

Usage: python -I assemble_s4.py WARP_HYPAD_ROOT
"""
import hashlib, sys
from pathlib import Path

root = Path(sys.argv[1])
s3 = (root / "study3/driver/s3.cpp").read_text().splitlines(keepends=True)
pc = root / "study4/driver/src_pieces"


def lines(a, b):  # 1-based inclusive
    return "".join(s3[a - 1:b])


parts = [
    (pc / "p1_header.cpp").read_text(),
    "// ==== verbatim s3.cpp lines 32-190 ====\n", lines(32, 190),
    (pc / "p2_sums.cpp").read_text(),
    "// ==== verbatim s3.cpp lines 216-324 ====\n", lines(216, 324),
    (pc / "p3_refine.cpp").read_text(),
    "// ==== verbatim s3.cpp lines 377-437 ====\n", lines(377, 437),
    (pc / "p4_run.cpp").read_text(),
    (pc / "p5_selftest.cpp").read_text(),
    "// ==== verbatim s3.cpp lines 693-1087 ====\n", lines(693, 1087),
    (pc / "p6_main.cpp").read_text(),
]
out = root / "study4/driver/s4.cpp"
out.write_text("".join(parts))
print("s3.cpp sha256", hashlib.sha256((root / "study3/driver/s3.cpp").read_bytes()).hexdigest())
print("s4.cpp sha256", hashlib.sha256(out.read_bytes()).hexdigest())

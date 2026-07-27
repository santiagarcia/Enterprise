"""Shared argument handling. NO SCRIPT DEPENDS ON THE CURRENT WORKING DIRECTORY.

Every script takes an explicit input path and an explicit output directory, both
resolved to absolute paths before anything is read or written. Importing the
visualization package from an arbitrary directory therefore works, and so does
running a script from anywhere.
"""

from __future__ import annotations

import argparse
import os
import sys


def bootstrap_path() -> None:
    """Put the visualization package root on sys.path.

    Scripts live in static/ and animation/ and are run as files, not modules, so
    the package root is two levels up from the script. Derived from __file__,
    never from the working directory.
    """
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if root not in sys.path:
        sys.path.insert(0, root)


def parser(description: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=description)
    p.add_argument("input", help="explicit path to the saved solver record")
    p.add_argument("--metadata", required=True,
                   help="path to c1_metadata.json, the provenance every figure "
                        "inherits")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--theme", default="scientific",
                   choices=["scientific", "scientific_print", "cinematic"])
    return p

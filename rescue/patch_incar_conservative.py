#!/usr/bin/env python3
"""
Patch INCAR files for conservative VASP restart/rescue runs.

Default patch:
  ALGO  = Normal
  POTIM = 0.1
  AMIX  = 0.1
  BMIX  = 0.0001
  AMIN  = 0.01

Example:
  python patch_incar_conservative.py "level*_*-ads" --apply
"""

from __future__ import annotations

import argparse
import glob
import re
import shutil
import time
from pathlib import Path


PATCH = {
    "ALGO": "Normal",
    "POTIM": "0.1",
    "AMIX": "0.1",
    "BMIX": "0.0001",
    "AMIN": "0.01",
}


def patch_text(text: str) -> str:
    lines = text.splitlines()
    keys_seen = set()

    for idx, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        m = re.match(r"^([A-Za-z][A-Za-z0-9_]*)\s*=", stripped)
        if not m:
            continue
        key = m.group(1).upper()
        if key in PATCH:
            lines[idx] = f"{key} = {PATCH[key]}"
            keys_seen.add(key)

    missing = [k for k in PATCH if k not in keys_seen]
    if missing:
        lines.append("")
        lines.append("# Conservative rescue settings")
        for k in missing:
            lines.append(f"{k} = {PATCH[k]}")

    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("glob_pattern", nargs="?", default="level*_*-ads")
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    dirs = [Path(p) for p in sorted(glob.glob(args.glob_pattern)) if Path(p).is_dir()]
    for d in dirs:
        incar = d / "INCAR"
        if not incar.exists():
            print(f"{d}: MISSING INCAR")
            continue

        old = incar.read_text(errors="ignore")
        new = patch_text(old)
        if old == new:
            print(f"{d}: OK already patched")
            continue

        if args.apply:
            stamp = time.strftime("%Y%m%d_%H%M%S")
            shutil.copy2(incar, d / f"INCAR.before_conservative_patch_{stamp}")
            incar.write_text(new)
            print(f"{d}: PATCHED")
        else:
            print(f"{d}: WOULD_PATCH")


if __name__ == "__main__":
    main()

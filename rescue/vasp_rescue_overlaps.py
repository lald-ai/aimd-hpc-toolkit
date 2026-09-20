#!/usr/bin/env python3
"""
Find VASP POSCARs with atom overlaps and replace them using safe OLD_* geometries.

If no safe geometry is found, writes DO_NOT_SUBMIT_OVERLAP.

Examples:
  python vasp_rescue_overlaps.py "level*_*-ads" --min-dist 0.65
  python vasp_rescue_overlaps.py "level*_*-ads" --min-dist 0.65 --apply
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import time
from pathlib import Path
from typing import Optional

import numpy as np
from ase.io import read


def min_dist(path: Path) -> Optional[tuple[float, int, str]]:
    try:
        atoms = read(str(path))
        dmat = atoms.get_all_distances(mic=True)
        np.fill_diagonal(dmat, 999.0)
        i, j = np.unravel_index(np.argmin(dmat), dmat.shape)
        return float(dmat[i, j]), len(atoms), f"{i+1}-{atoms[i].symbol}/{j+1}-{atoms[j].symbol}"
    except Exception:
        return None


def candidates(d: Path) -> list[Path]:
    out: list[Path] = []
    for name in ["POSCAR", "CONTCAR"]:
        p = d / name
        if p.exists() and p.stat().st_size > 0:
            out.append(p)
    for p in sorted(d.glob("POSCAR.before*")):
        if p.exists() and p.stat().st_size > 0:
            out.append(p)
    for old in sorted(d.glob("OLD_*"), reverse=True):
        for name in ["CONTCAR", "POSCAR"]:
            p = old / name
            if p.exists() and p.stat().st_size > 0:
                out.append(p)

    seen: set[str] = set()
    unique: list[Path] = []
    for p in out:
        rp = str(p.resolve())
        if rp not in seen:
            unique.append(p)
            seen.add(rp)
    return unique


def rescue_overlap(d: Path, cutoff: float, apply: bool) -> str:
    current = min_dist(d / "POSCAR")
    if current is None:
        return "BAD_PARSE"
    cur_min, cur_n, cur_pair = current
    if cur_min >= cutoff:
        return "OK"

    good: list[tuple[float, int, str, Path]] = []
    for p in candidates(d):
        r = min_dist(p)
        if r is None:
            continue
        md, nat, pair = r
        if md >= cutoff:
            good.append((md, nat, pair, p))

    if not good:
        if apply:
            (d / "DO_NOT_SUBMIT_OVERLAP").write_text(
                f"No replacement with min_dist >= {cutoff}; current={cur_min:.3f} pair={cur_pair}\n"
            )
        return f"NO_REPLACEMENT current={cur_min:.3f} pair={cur_pair}"

    best_min, best_n, best_pair, best_path = sorted(good, reverse=True)[0]
    if not apply:
        return f"WOULD_FIX current={cur_min:.3f} pair={cur_pair} using={best_path} new_min={best_min:.3f}"

    stamp = time.strftime("%Y%m%d_%H%M%S")
    shutil.copy2(d / "POSCAR", d / f"POSCAR.before_overlap_rescue_{stamp}")
    if (d / "CONTCAR").exists():
        shutil.copy2(d / "CONTCAR", d / f"CONTCAR.before_overlap_rescue_{stamp}")
    shutil.copy2(best_path, d / "POSCAR")
    shutil.copy2(best_path, d / "CONTCAR")
    return f"FIXED current={cur_min:.3f} pair={cur_pair} using={best_path} new_min={best_min:.3f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("glob_pattern", nargs="?", default="level*_*-ads")
    parser.add_argument("--min-dist", type=float, default=0.65)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    dirs = [Path(p) for p in sorted(glob.glob(args.glob_pattern)) if Path(p).is_dir()]
    for d in dirs:
        msg = rescue_overlap(d, args.min_dist, args.apply)
        if msg != "OK":
            print(f"{d}: {msg}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Repair malformed/truncated POSCAR files in VASP calculation directories.

Strategy:
  species order  <- POTCAR TITEL order
  atom counts    <- OUTCAR ions per type, POSCAR-like files, or OLD_* fallbacks
  cell           <- POSCAR/CONTCAR/XDATCAR/OLD_*
  geometry       <- current POSCAR, CONTCAR, XDATCAR last frame, OUTCAR last ionic step, OLD_*

Requires ASE for final validation.

Examples:
  python vasp_rescue_poscar.py "level*_*-ads"
  python vasp_rescue_poscar.py "level*_*-ads" --apply --also-contcar
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import time
from pathlib import Path
from typing import Optional

from ase.io import read, write


def ase_ok(path: Path) -> bool:
    try:
        read(str(path))
        return True
    except Exception:
        return False


def is_int_tokens(tokens: list[str]) -> bool:
    if not tokens:
        return False
    try:
        [int(x) for x in tokens]
        return True
    except Exception:
        return False


def is_coord_mode(line: str) -> bool:
    s = line.strip().lower()
    return s == "direct" or s == "cartesian" or s.startswith("direct ") or s.startswith("cartesian ")


def is_selective(line: str) -> bool:
    return line.strip().lower().startswith("selective")


def potcar_species(path: Path) -> list[str]:
    if not path.exists():
        return []
    out: list[str] = []
    for line in path.read_text(errors="ignore").splitlines():
        if "TITEL" in line:
            parts = line.split()
            if len(parts) >= 4:
                out.append(parts[3].split("_")[0])
    return out


def counts_from_outcar(path: Path) -> list[int]:
    if not path.exists():
        return []
    last: list[int] = []
    for line in path.read_text(errors="ignore").splitlines():
        if "ions per type" in line:
            nums = re.findall(r"\d+", line.split("=", 1)[-1])
            if nums:
                last = [int(x) for x in nums]
    return last


def counts_from_poscar_like(path: Path, nspecies: int) -> list[int]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    lines = path.read_text(errors="ignore").splitlines()
    counts: list[int] = []
    for line in lines[5:50]:
        if is_selective(line) or is_coord_mode(line):
            break
        toks = line.split()
        if is_int_tokens(toks):
            counts.extend([int(x) for x in toks])
            if len(counts) >= nspecies:
                return counts[:nspecies]
    return []


def cell_from(path: Path) -> Optional[tuple[str, list[str]]]:
    if not path.exists() or path.stat().st_size == 0:
        return None
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < 5:
        return None
    try:
        scale = lines[1].strip()
        float(scale.split()[0])
        cell = lines[2:5]
        for row in cell:
            vals = row.split()[:3]
            if len(vals) < 3:
                return None
            [float(v) for v in vals]
        return scale, cell
    except Exception:
        return None


def natoms_from_counts(counts: list[int]) -> int:
    return sum(counts)


def write_poscar_from_atoms(source: Path, target: Path) -> bool:
    try:
        atoms = read(str(source))
        write(str(target), atoms, format="vasp", direct=True, vasp5=True, sort=False)
        return True
    except Exception:
        return False


def xdatcar_to_atoms(source: Path, target: Path) -> bool:
    try:
        atoms = read(str(source), index=-1)
        write(str(target), atoms, format="vasp", direct=True, vasp5=True, sort=False)
        return True
    except Exception:
        return False


def outcar_to_atoms(source: Path, target: Path) -> bool:
    try:
        atoms = read(str(source), index=-1)
        write(str(target), atoms, format="vasp", direct=True, vasp5=True, sort=False)
        return True
    except Exception:
        return False


def candidate_sources(d: Path) -> list[Path]:
    out: list[Path] = []
    for name in ["POSCAR", "CONTCAR", "XDATCAR", "OUTCAR"]:
        p = d / name
        if p.exists() and p.stat().st_size > 0:
            out.append(p)
    for old in sorted(d.glob("OLD_*"), reverse=True):
        for name in ["CONTCAR", "POSCAR", "XDATCAR", "OUTCAR"]:
            p = old / name
            if p.exists() and p.stat().st_size > 0:
                out.append(p)
    return out


def rescue_dir(d: Path, apply: bool, also_contcar: bool) -> str:
    poscar = d / "POSCAR"
    if ase_ok(poscar):
        return "OK already readable"

    species = potcar_species(d / "POTCAR")
    if not species:
        species = potcar_species(Path.cwd() / "POTCAR")
    if not species:
        return "BAD no POTCAR species"

    counts: list[int] = []
    for p in [d / "OUTCAR", d / "POSCAR", d / "CONTCAR"] + [q for old in sorted(d.glob("OLD_*"), reverse=True) for q in [old/"OUTCAR", old/"POSCAR", old/"CONTCAR"]]:
        if p.name == "OUTCAR":
            counts = counts_from_outcar(p)
        else:
            counts = counts_from_poscar_like(p, len(species))
        if counts:
            break

    if not counts:
        return "BAD no atom counts found"
    if len(counts) != len(species):
        return f"BAD species/count mismatch species={len(species)} counts={len(counts)}"

    # Prefer ASE-readable sources for geometry because that preserves a sane cell/coords.
    source = None
    tmp = d / f".tmp_rescue_{int(time.time())}.POSCAR"
    for cand in candidate_sources(d):
        ok = False
        if cand.name in ("POSCAR", "CONTCAR"):
            ok = write_poscar_from_atoms(cand, tmp)
        elif cand.name == "XDATCAR":
            ok = xdatcar_to_atoms(cand, tmp)
        elif cand.name == "OUTCAR":
            ok = outcar_to_atoms(cand, tmp)

        if ok and ase_ok(tmp):
            source = cand
            break

    if source is None:
        if tmp.exists():
            tmp.unlink()
        return "BAD no readable geometry source"

    if not apply:
        tmp.unlink(missing_ok=True)
        return f"WOULD_RESCUE using {source}"

    stamp = time.strftime("%Y%m%d_%H%M%S")
    if poscar.exists():
        shutil.copy2(poscar, d / f"POSCAR.before_rescue_{stamp}")
    shutil.copy2(tmp, poscar)
    if also_contcar:
        if (d / "CONTCAR").exists():
            shutil.copy2(d / "CONTCAR", d / f"CONTCAR.before_rescue_{stamp}")
        shutil.copy2(tmp, d / "CONTCAR")
    tmp.unlink(missing_ok=True)

    if ase_ok(poscar):
        return f"OK rescued using {source}"
    return "BAD wrote rescued POSCAR but ASE still fails"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("glob_pattern", nargs="?", default="level*_*-ads")
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--also-contcar", action="store_true")
    args = parser.parse_args()

    dirs = [Path(p) for p in sorted(glob.glob(args.glob_pattern)) if Path(p).is_dir()]
    for d in dirs:
        msg = rescue_dir(d, args.apply, args.also_contcar)
        if not msg.startswith("OK already"):
            print(f"{d}: {msg}")


if __name__ == "__main__":
    main()

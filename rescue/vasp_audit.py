#!/usr/bin/env python3
"""
Audit VASP calculation folders for POSCAR/POTCAR consistency, atom overlaps,
completion, fatal errors, and key input-file presence.

Requires ASE for POSCAR reading and distance checks.

Example:
    python vasp_audit.py "level*_*-ads" --min-dist 0.65 --out audit.csv
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from ase.io import read


FATAL_PATTERNS = re.compile(
    r"ZHEGV|EDDDAV|LAPACK|segmentation fault|segfault|out of memory|OOM",
    re.IGNORECASE,
)


@dataclass
class AuditResult:
    directory: str
    poscar_ok: bool
    poscar_error: str
    natoms: Optional[int]
    min_dist: Optional[float]
    min_pair: str
    poscar_species: str
    potcar_species: str
    species_match: bool
    has_incar: bool
    has_kpoints: bool
    has_potcar: bool
    has_submit: bool
    completed: bool
    fatal_history: bool
    algo: str
    potim: str
    amix: str
    bmix: str
    amin: str
    status: str


def potcar_species(path: Path) -> list[str]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    species: list[str] = []
    for line in path.read_text(errors="ignore").splitlines():
        if "TITEL" in line:
            parts = line.split()
            if len(parts) >= 4:
                species.append(parts[3].split("_")[0])
    return species


def poscar_species(path: Path) -> list[str]:
    try:
        lines = path.read_text(errors="ignore").splitlines()
        if len(lines) >= 6:
            return lines[5].split()
    except Exception:
        pass
    return []


def grep_last_key(path: Path, key: str) -> str:
    if not path.exists():
        return ""
    last = ""
    for line in path.read_text(errors="ignore").splitlines():
        if re.match(rf"^\s*{re.escape(key)}\b", line, re.IGNORECASE):
            last = line.strip()
    return last


def file_contains(path: Path, pattern: str) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    return pattern in path.read_text(errors="ignore")


def has_fatal(d: Path) -> bool:
    candidates = [d / "OUTCAR"] + list(d.glob("slurm*.out")) + list(d.glob("slurm*.err"))
    for p in candidates:
        if not p.exists() or p.stat().st_size == 0:
            continue
        if FATAL_PATTERNS.search(p.read_text(errors="ignore")):
            return True
    return False


def min_distance(poscar: Path) -> tuple[int, float, str]:
    atoms = read(str(poscar))
    dist = atoms.get_all_distances(mic=True)
    np.fill_diagonal(dist, 999.0)
    i, j = np.unravel_index(np.argmin(dist), dist.shape)
    return len(atoms), float(dist[i, j]), f"{i+1}-{atoms[i].symbol}/{j+1}-{atoms[j].symbol}"


def audit_dir(d: Path, min_dist_cutoff: float) -> AuditResult:
    poscar = d / "POSCAR"
    potcar = d / "POTCAR"
    incar = d / "INCAR"

    pos_ok = False
    pos_err = ""
    natoms: Optional[int] = None
    min_dist: Optional[float] = None
    min_pair = ""

    try:
        natoms, min_dist, min_pair = min_distance(poscar)
        pos_ok = True
    except Exception as exc:
        pos_err = str(exc)

    psp = poscar_species(poscar)
    potsp = potcar_species(potcar)
    species_match = bool(psp and potsp and psp == potsp)

    completed = file_contains(d / "OUTCAR", "reached required accuracy")
    fatal = has_fatal(d)

    issues: list[str] = []
    if not pos_ok:
        issues.append("BAD_POSCAR")
    if not species_match:
        issues.append("SPECIES_MISMATCH")
    if min_dist is not None and min_dist < min_dist_cutoff:
        issues.append("OVERLAP")
    for name in ["INCAR", "KPOINTS", "POTCAR"]:
        if not (d / name).exists():
            issues.append(f"MISSING_{name}")
    if not (d / "first_vasp.run").exists() and not (d / "vasp.run").exists():
        issues.append("MISSING_SUBMIT")
    if fatal:
        issues.append("FATAL_HISTORY")

    return AuditResult(
        directory=str(d),
        poscar_ok=pos_ok,
        poscar_error=pos_err,
        natoms=natoms,
        min_dist=min_dist,
        min_pair=min_pair,
        poscar_species=" ".join(psp),
        potcar_species=" ".join(potsp),
        species_match=species_match,
        has_incar=(d / "INCAR").exists(),
        has_kpoints=(d / "KPOINTS").exists(),
        has_potcar=potcar.exists(),
        has_submit=(d / "first_vasp.run").exists() or (d / "vasp.run").exists(),
        completed=completed,
        fatal_history=fatal,
        algo=grep_last_key(incar, "ALGO"),
        potim=grep_last_key(incar, "POTIM"),
        amix=grep_last_key(incar, "AMIX"),
        bmix=grep_last_key(incar, "BMIX"),
        amin=grep_last_key(incar, "AMIN"),
        status="OK" if not issues else ";".join(issues),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("glob_pattern", nargs="?", default="level*_*-ads")
    parser.add_argument("--min-dist", type=float, default=0.65)
    parser.add_argument("--out", default="vasp_audit.csv")
    args = parser.parse_args()

    dirs = [Path(p) for p in sorted(glob.glob(args.glob_pattern)) if Path(p).is_dir()]
    results = [audit_dir(d, args.min_dist) for d in dirs]

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(AuditResult.__dataclass_fields__.keys()))
        writer.writeheader()
        for r in results:
            writer.writerow(r.__dict__)

    bad = [r for r in results if r.status != "OK"]
    print(f"Scanned: {len(results)} directories")
    print(f"Issues:  {len(bad)}")
    print(f"Wrote:   {args.out}")
    for r in bad[:50]:
        md = "" if r.min_dist is None else f" min={r.min_dist:.3f} pair={r.min_pair}"
        print(f"{r.directory:35s} {r.status}{md}")
    if len(bad) > 50:
        print(f"... plus {len(bad)-50} more")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summarize constrained-MD/PMF window progress without modifying calculations."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from ase.io import read


def count_oszicar_steps(path: Path) -> int:
    if not path.exists():
        return 0
    n = 0
    for line in path.read_text(errors="replace").splitlines():
        if " T=" in line or line.lstrip().startswith("T="):
            n += 1
    return n


def count_blue_moon(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for line in path.read_text(errors="replace").splitlines() if line.strip().startswith("b_m>"))


def distance(path: Path, a: int, b: int) -> float | None:
    if not path.exists():
        return None
    try:
        atoms = read(path, format="vasp")
        return float(atoms.get_distance(a - 1, b - 1, mic=True))
    except Exception:
        return None


def target_from_name(name: str) -> float | None:
    matches = re.findall(r"\d+(?:\.\d+)?", name)
    return float(matches[-1]) if matches else None


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", nargs="?", default=".")
    p.add_argument("--glob", default="**/*pmf*", help="directory glob under root")
    p.add_argument("--carbon-index", type=int, help="1-based C index for final distance check")
    p.add_argument("--dopant-index", type=int, help="1-based dopant/surface index")
    p.add_argument("--csv", default="pmf_progress.csv")
    args = p.parse_args()

    root = Path(args.root)
    dirs = sorted(p for p in root.glob(args.glob) if p.is_dir())
    rows = []
    for d in dirs:
        outcar = d / "OUTCAR"
        contcar = d / "CONTCAR"
        report = d / "REPORT"
        row = {
            "window": str(d),
            "target_from_folder": target_from_name(d.name),
            "oszicar_steps": count_oszicar_steps(d / "OSZICAR"),
            "blue_moon_samples": count_blue_moon(report),
            "has_OUTCAR": outcar.exists(),
            "has_CONTCAR": contcar.exists() and contcar.stat().st_size > 0,
            "has_REPORT": report.exists() and report.stat().st_size > 0,
            "final_C_dopant": "",
        }
        if args.carbon_index and args.dopant_index:
            val = distance(contcar, args.carbon_index, args.dopant_index)
            row["final_C_dopant"] = "" if val is None else f"{val:.6f}"
        rows.append(row)

    fields = ["window", "target_from_folder", "oszicar_steps", "blue_moon_samples",
              "has_OUTCAR", "has_CONTCAR", "has_REPORT", "final_C_dopant"]
    with open(args.csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"{'window':55s} {'steps':>7s} {'BM':>7s} {'CONTCAR':>8s} {'C-dop':>10s}")
    for r in rows:
        print(f"{r['window'][-55:]:55s} {r['oszicar_steps']:7d} {r['blue_moon_samples']:7d} "
              f"{str(r['has_CONTCAR']):>8s} {str(r['final_C_dopant']):>10s}")
    print(f"wrote {args.csv}")


if __name__ == "__main__":
    main()

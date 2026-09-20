#!/usr/bin/env python3
"""Extract Blue-Moon reaction coordinates and free-energy gradients from VASP REPORT files.

VASP writes constrained-coordinate records as ``cc>`` and Blue-Moon records as
``b_m>`` when LBLUEOUT is enabled.  This parser pairs rows by order within each
MD step and can therefore select one constrained coordinate from a run that has
multiple constraints (for example C--dopant plus C--O).

Outputs:
  * raw CSV: one row per parsed MD sample
  * summary CSV: mean coordinate, mean gradient, SD, SEM, and sample count per window
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
STEP_RE = re.compile(r"MD\s+step\s+No\.\s*(\d+)", re.I)


def fnum(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def parse_report(path: Path, constraint_index: int) -> list[dict]:
    samples: list[dict] = []
    step = None
    cc_rows: list[float] = []
    bm_rows: list[float] = []

    def flush() -> None:
        nonlocal cc_rows, bm_rows
        k = constraint_index - 1
        if step is not None and k < len(cc_rows) and k < len(bm_rows):
            samples.append({"step": step, "coordinate": cc_rows[k], "gradient": bm_rows[k]})
        cc_rows, bm_rows = [], []

    for line in path.read_text(errors="replace").splitlines():
        m = STEP_RE.search(line)
        if m:
            flush()
            step = int(m.group(1))
            continue
        s = line.strip()
        if s.startswith("cc>"):
            nums = re.findall(FLOAT, s.replace("D", "E"))
            # For distance coordinates VASP prints e.g. cc> R value target error.
            # The first numeric field after the coordinate label is the actual coordinate.
            if nums:
                cc_rows.append(fnum(nums[0]))
        elif s.startswith("b_m>"):
            nums = re.findall(FLOAT, s.replace("D", "E"))
            # Columns: lambda, |z|^-1/2, GkT, |z|^-1/2*(lambda+GkT)
            if len(nums) >= 4:
                bm_rows.append(fnum(nums[3]))
    flush()
    return samples


def infer_window(path: Path) -> str:
    for part in reversed(path.parts):
        m = re.search(r"pmf[_-](\d+(?:\.\d+)?)", part, re.I)
        if m:
            return m.group(1)
    for part in reversed(path.parts):
        m = re.search(r"(\d+\.\d+)", part)
        if m:
            return m.group(1)
    return path.parent.name


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("paths", nargs="+", help="REPORT files or directories searched recursively for REPORT/report.*")
    p.add_argument("--constraint-index", type=int, default=1, help="1-based constrained-coordinate row to analyze")
    p.add_argument("--discard", type=int, default=0, help="discard this many parsed samples from the start of each REPORT")
    p.add_argument("--raw", default="blue_moon_raw.csv")
    p.add_argument("--summary", default="blue_moon_summary.csv")
    args = p.parse_args()
    if args.constraint_index < 1:
        p.error("--constraint-index must be >= 1")

    reports: list[Path] = []
    for text in args.paths:
        path = Path(text)
        if path.is_file():
            reports.append(path)
        elif path.is_dir():
            reports.extend(sorted(path.rglob("REPORT")))
            reports.extend(sorted(path.rglob("report.*")))
    # preserve order while deduplicating
    reports = list(dict.fromkeys(p.resolve() for p in reports))
    if not reports:
        raise SystemExit("No REPORT/report.* files found")

    all_rows: list[dict] = []
    for report in reports:
        rows = parse_report(report, args.constraint_index)
        rows = rows[args.discard:]
        window = infer_window(report)
        for row in rows:
            all_rows.append({"window": window, "report": str(report), **row})

    if not all_rows:
        raise SystemExit("No paired cc>/b_m> samples parsed. Check REPORT format and constraint index.")

    with open(args.raw, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["window", "report", "step", "coordinate", "gradient"])
        writer.writeheader()
        writer.writerows(all_rows)

    grouped: dict[str, list[dict]] = {}
    for row in all_rows:
        grouped.setdefault(row["window"], []).append(row)

    summary = []
    for window, rows in grouped.items():
        coords = np.array([r["coordinate"] for r in rows], dtype=float)
        grads = np.array([r["gradient"] for r in rows], dtype=float)
        n = len(grads)
        sd = float(np.std(grads, ddof=1)) if n > 1 else float("nan")
        sem = sd / math.sqrt(n) if n > 1 else float("nan")
        summary.append({
            "window": window,
            "coordinate_mean": float(np.mean(coords)),
            "coordinate_sd": float(np.std(coords, ddof=1)) if n > 1 else float("nan"),
            "gradient_mean": float(np.mean(grads)),
            "gradient_sd": sd,
            "gradient_sem": sem,
            "n": n,
        })
    summary.sort(key=lambda r: r["coordinate_mean"])

    with open(args.summary, "w", newline="", encoding="utf-8") as fh:
        fields = ["window", "coordinate_mean", "coordinate_sd", "gradient_mean", "gradient_sd", "gradient_sem", "n"]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)

    print(f"parsed {len(all_rows)} samples from {len(reports)} REPORT file(s)")
    print(f"wrote {args.raw}")
    print(f"wrote {args.summary}")


if __name__ == "__main__":
    main()

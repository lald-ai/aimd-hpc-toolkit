#!/usr/bin/env python3
"""Analyze convergence/equilibration of chunked VASP AIMD from OSZICAR files.

The script stitches numeric ``chunk_*`` directories, extracts T, F, E0 and EK,
plots raw/rolling/cumulative behavior, estimates tail drift, and writes block
statistics.  It is intended as a diagnostic; a numerical threshold is only
applied when the user explicitly supplies one.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PATTERNS = {
    "temperature": re.compile(r"\bT=\s*([-+0-9.Ee]+)"),
    "free_energy_F": re.compile(r"\bF=\s*([-+0-9.Ee]+)"),
    "energy_E0": re.compile(r"\bE0=\s*([-+0-9.Ee]+)"),
    "kinetic_EK": re.compile(r"\bEK=\s*([-+0-9.Ee]+)"),
}


def chunk_key(path: Path):
    m = re.search(r"chunk_(\d+)", path.name)
    return int(m.group(1)) if m else 10**12


def parse_oszicar(path: Path):
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        if "T=" not in line:
            continue
        row = {}
        for name, pat in PATTERNS.items():
            m = pat.search(line)
            row[name] = float(m.group(1)) if m else math.nan
        rows.append(row)
    return rows


def rolling_mean(y: np.ndarray, n: int) -> np.ndarray:
    out = np.full(len(y), np.nan)
    if n <= 1:
        return y.copy()
    if len(y) >= n:
        kernel = np.ones(n) / n
        out[n - 1:] = np.convolve(y, kernel, mode="valid")
    return out


def block_stats(y: np.ndarray, block: int):
    result = []
    for start in range(0, len(y), block):
        z = y[start:start + block]
        z = z[np.isfinite(z)]
        if len(z):
            result.append((start, start + len(z) - 1, len(z), float(np.mean(z)), float(np.std(z, ddof=1)) if len(z) > 1 else math.nan))
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", nargs="?", default=".")
    p.add_argument("--timestep-fs", type=float, default=1.0)
    p.add_argument("--rolling-ps", type=float, nargs="+", default=[0.5, 1.0, 2.5])
    p.add_argument("--block-ps", type=float, default=2.5)
    p.add_argument("--tail-ps", type=float, default=5.0, help="tail interval used for mean/SD/drift diagnostics")
    p.add_argument("--drift-threshold", type=float, help="optional absolute drift threshold in eV/ps")
    p.add_argument("--prefix", default="aimd_convergence")
    args = p.parse_args()

    root = Path(args.root)
    osz = sorted(root.glob("chunk_*/OSZICAR"), key=lambda p: chunk_key(p.parent))
    if not osz and (root / "OSZICAR").exists():
        osz = [root / "OSZICAR"]
    if not osz:
        raise SystemExit("No OSZICAR found in root or chunk_*/OSZICAR")

    rows = []
    for f in osz:
        rows.extend(parse_oszicar(f))
    if not rows:
        raise SystemExit("No ionic MD lines containing T= were parsed")

    names = list(PATTERNS)
    data = {name: np.array([r[name] for r in rows], dtype=float) for name in names}
    step = np.arange(1, len(rows) + 1)
    time_ps = step * args.timestep_fs / 1000.0

    block_n = max(1, round(args.block_ps * 1000 / args.timestep_fs))
    tail_n = max(2, round(args.tail_ps * 1000 / args.timestep_fs))
    tail_start = max(0, len(step) - tail_n)

    # Summary and block statistics use E0 when available, otherwise F.
    energy_name = "energy_E0" if np.isfinite(data["energy_E0"]).any() else "free_energy_F"
    energy = data[energy_name]
    mask = np.isfinite(energy[tail_start:])
    x_tail = time_ps[tail_start:][mask]
    y_tail = energy[tail_start:][mask]
    slope = float(np.polyfit(x_tail, y_tail, 1)[0]) if len(y_tail) >= 2 else math.nan
    tail_mean = float(np.mean(y_tail)) if len(y_tail) else math.nan
    tail_sd = float(np.std(y_tail, ddof=1)) if len(y_tail) > 1 else math.nan

    summary_path = f"{args.prefix}_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(f"files: {len(osz)}\n")
        fh.write(f"MD samples: {len(rows)}\n")
        fh.write(f"timestep_fs: {args.timestep_fs}\n")
        fh.write(f"energy_metric: {energy_name}\n")
        fh.write(f"tail_ps: {args.tail_ps}\n")
        fh.write(f"tail_mean_eV: {tail_mean:.10g}\n")
        fh.write(f"tail_sd_eV: {tail_sd:.10g}\n")
        fh.write(f"tail_linear_drift_eV_per_ps: {slope:.10g}\n")
        if args.drift_threshold is not None and math.isfinite(slope):
            fh.write(f"drift_threshold_eV_per_ps: {args.drift_threshold}\n")
            fh.write(f"within_drift_threshold: {abs(slope) <= args.drift_threshold}\n")

    blocks = block_stats(energy, block_n)
    with open(f"{args.prefix}_blocks.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["start_step", "end_step", "n", "mean_eV", "sd_eV"])
        writer.writerows(blocks)

    # Energy convergence plot.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(time_ps, energy, linewidth=0.7, alpha=0.45, label=energy_name)
    for ps in args.rolling_ps:
        n = max(1, round(ps * 1000 / args.timestep_fs))
        ax.plot(time_ps, rolling_mean(energy, n), linewidth=1.5, label=f"{ps:g} ps rolling mean")
    finite = np.isfinite(energy)
    if finite.any():
        cumulative = np.cumsum(np.where(finite, energy, 0.0)) / np.maximum(1, np.cumsum(finite))
        ax.plot(time_ps, cumulative, linewidth=1.2, label="cumulative mean")
    ax.axvspan(time_ps[tail_start], time_ps[-1], alpha=0.08, label=f"last {args.tail_ps:g} ps")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(f"{energy_name} (eV)")
    ax.set_title("AIMD energy convergence")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.prefix}_energy.png", dpi=200)

    # Temperature plot.
    temp = data["temperature"]
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.plot(time_ps, temp, linewidth=0.7, alpha=0.5, label="T")
    for ps in args.rolling_ps:
        n = max(1, round(ps * 1000 / args.timestep_fs))
        ax.plot(time_ps, rolling_mean(temp, n), linewidth=1.4, label=f"{ps:g} ps rolling mean")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("AIMD temperature convergence")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.prefix}_temperature.png", dpi=200)

    print(f"wrote {summary_path}")
    print(f"wrote {args.prefix}_blocks.csv")
    print(f"wrote {args.prefix}_energy.png")
    print(f"wrote {args.prefix}_temperature.png")


if __name__ == "__main__":
    main()

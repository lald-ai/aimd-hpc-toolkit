#!/usr/bin/env python3
"""Integrate window-averaged Blue-Moon free-energy gradients into a PMF.

Input is normally ``blue_moon_summary.csv`` from extract_blue_moon.py.  The PMF
is obtained by trapezoidal integration of dA/dxi over xi.  The sign convention
is intentionally user-controlled because reaction-coordinate definitions differ.
Uncertainties are propagated from independent per-window SEM values.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_csv(path: str):
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    if not rows:
        raise ValueError("input CSV is empty")
    x = np.array([float(r["coordinate_mean"]) for r in rows])
    g = np.array([float(r["gradient_mean"]) for r in rows])
    sem = np.array([float(r.get("gradient_sem", "nan")) for r in rows])
    order = np.argsort(x)
    return [rows[i] for i in order], x[order], g[order], sem[order]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="Blue-Moon window summary CSV")
    p.add_argument("--sign", type=float, default=1.0, help="multiply gradients by this sign before integration")
    p.add_argument("--reference", choices=["first", "last", "minimum"], default="first",
                   help="zero of free energy after integration")
    p.add_argument("--output", default="pmf_profile.csv")
    p.add_argument("--plot", default="pmf_profile.png")
    p.add_argument("--xlabel", default="Reaction coordinate (Angstrom)")
    args = p.parse_args()

    rows, x, g, sem = load_csv(args.input)
    g = args.sign * g
    n = len(x)
    if n < 2:
        raise SystemExit("Need at least two windows to integrate a PMF")

    A = np.zeros(n)
    var = np.zeros(n)
    for i in range(1, n):
        dx = x[i] - x[i - 1]
        A[i] = A[i - 1] + 0.5 * dx * (g[i - 1] + g[i])
        if np.isfinite(sem[i - 1]) and np.isfinite(sem[i]):
            var[i] = var[i - 1] + (0.5 * dx) ** 2 * (sem[i - 1] ** 2 + sem[i] ** 2)
        else:
            var[i] = np.nan

    if args.reference == "first":
        shift = A[0]
    elif args.reference == "last":
        shift = A[-1]
    else:
        shift = np.nanmin(A)
    A = A - shift
    sigma = np.sqrt(var)

    fields = ["window", "coordinate", "gradient", "gradient_sem", "pmf", "pmf_sigma"]
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for src, xi, gi, sei, ai, si in zip(rows, x, g, sem, A, sigma):
            writer.writerow({
                "window": src.get("window", ""),
                "coordinate": xi,
                "gradient": gi,
                "gradient_sem": sei,
                "pmf": ai,
                "pmf_sigma": si,
            })

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(x, A, marker="o", linewidth=1.6)
    if np.all(np.isfinite(sigma)):
        ax.fill_between(x, A - sigma, A + sigma, alpha=0.2)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel("Relative free energy / PMF (eV)")
    ax.set_title("Blue-Moon potential of mean force")
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)
    print(f"wrote {args.output}")
    print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()

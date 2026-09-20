#!/usr/bin/env python3
"""Stage a constrained CO--surface distance window while preserving the C--O bond.

This is a reconstructed utility based on the retained production workflow behavior.
It is NOT a byte-for-byte copy of the historical Gautschi script.

The script reads a starting POSCAR (normally the *starting* POSCAR of the previous
window), identifies the bonded oxygen as the O atom nearest the selected carbon,
and translates C and O together so that the C--dopant minimum-image distance
matches the requested target.  Optionally it writes a two-distance ICONST file:
C--dopant and C--O.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
from ase.io import read, write


def one_based_to_zero(idx: int, n: int, label: str) -> int:
    if idx < 1 or idx > n:
        raise ValueError(f"{label} index {idx} outside 1..{n}")
    return idx - 1


def unique_symbol_index(atoms, symbol: str) -> int:
    matches = [i for i, atom in enumerate(atoms) if atom.symbol == symbol]
    if len(matches) != 1:
        raise ValueError(
            f"Auto-detection requires exactly one {symbol}; found {len(matches)}. "
            f"Pass an explicit atom index."
        )
    return matches[0]


def nearest_symbol(atoms, center: int, symbol: str) -> tuple[int, float]:
    candidates = [i for i, atom in enumerate(atoms) if atom.symbol == symbol and i != center]
    if not candidates:
        raise ValueError(f"No {symbol} atoms found")
    distances = [(i, atoms.get_distance(center, i, mic=True)) for i in candidates]
    return min(distances, key=lambda item: item[1])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="source POSCAR-like structure")
    p.add_argument("output", help="output POSCAR")
    p.add_argument("--target", type=float, required=True, help="target C--dopant distance in Angstrom")
    p.add_argument("--dopant-index", type=int, required=True, help="1-based surface/dopant atom index")
    p.add_argument("--carbon-index", type=int, help="1-based C index; auto-detects only if exactly one C exists")
    p.add_argument("--co-min", type=float, default=0.95, help="minimum plausible C--O bond length")
    p.add_argument("--co-max", type=float, default=1.40, help="maximum plausible C--O bond length")
    p.add_argument("--tolerance", type=float, default=2e-3, help="distance verification tolerance in Angstrom")
    p.add_argument("--iconst", help="optional ICONST output with C--dopant and C--O hard constraints")
    p.add_argument("--backup", action="store_true", help="copy input to <output>.source_before_distance_adjust")
    args = p.parse_args()

    atoms = read(args.input, format="vasp")
    n = len(atoms)
    c = one_based_to_zero(args.carbon_index, n, "carbon") if args.carbon_index else unique_symbol_index(atoms, "C")
    dop = one_based_to_zero(args.dopant_index, n, "dopant")
    if c == dop:
        raise ValueError("Carbon and dopant indices must be different")

    o, co_before = nearest_symbol(atoms, c, "O")
    if not (args.co_min <= co_before <= args.co_max):
        raise ValueError(
            f"Nearest O gives C--O={co_before:.4f} A, outside accepted "
            f"range [{args.co_min}, {args.co_max}] A. Refusing to move atoms."
        )

    cd_before = atoms.get_distance(c, dop, mic=True)
    vec_c_to_dop = atoms.get_distance(c, dop, mic=True, vector=True)
    norm = np.linalg.norm(vec_c_to_dop)
    if norm < 1e-12:
        raise ValueError("C and dopant are coincident")

    displacement = vec_c_to_dop / norm * (cd_before - args.target)
    atoms.positions[c] += displacement
    atoms.positions[o] += displacement
    atoms.wrap()

    cd_after = atoms.get_distance(c, dop, mic=True)
    co_after = atoms.get_distance(c, o, mic=True)
    if abs(cd_after - args.target) > args.tolerance:
        raise RuntimeError(
            f"Final C--dopant distance {cd_after:.6f} A differs from target "
            f"{args.target:.6f} A by more than {args.tolerance:g} A"
        )
    if abs(co_after - co_before) > args.tolerance:
        raise RuntimeError(
            f"Rigid CO check failed: C--O changed from {co_before:.6f} to {co_after:.6f} A"
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if args.backup:
        shutil.copy2(args.input, str(out) + ".source_before_distance_adjust")
    write(out, atoms, format="vasp", direct=True, vasp5=True, sort=False)

    if args.iconst:
        # ICONST atom numbers are 1-based. STATUS=0 means constrained coordinate.
        Path(args.iconst).write_text(
            f"R {c + 1} {dop + 1} 0\n"
            f"R {c + 1} {o + 1} 0\n",
            encoding="utf-8",
        )

    print(f"C index              : {c + 1}")
    print(f"bonded O index       : {o + 1}")
    print(f"dopant/surface index : {dop + 1}")
    print(f"C--dopant             : {cd_before:.6f} -> {cd_after:.6f} A")
    print(f"C--O                  : {co_before:.6f} -> {co_after:.6f} A")
    print(f"wrote                 : {out}")
    if args.iconst:
        print(f"wrote                 : {args.iconst}")


if __name__ == "__main__":
    main()

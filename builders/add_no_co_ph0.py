#!/usr/bin/env python3
import numpy as np
from ase.io import read, write
from ase import Atoms

# -------- settings --------
out_poscar = "POSCAR_solvent_NO_CO_pH0"
min_dist = 1.25  # slightly softer so proton placement is possible in dense water
rng = np.random.default_rng(7)

# diatomic bond lengths (A)
bond_NO = 1.15
bond_CO = 1.13

# hydronium extra O-H bond tries (A): start physical, then relax a bit if crowded
oh_try = [0.98, 1.02, 1.06, 1.10, 1.15, 1.20]

# candidate insertion points in fractional coords (spread out)
frac_points = [
    (0.20, 0.20, 0.20),
    (0.80, 0.80, 0.80),
    (0.20, 0.80, 0.50),
    (0.80, 0.20, 0.50),
    (0.50, 0.50, 0.20),
    (0.50, 0.50, 0.80),
]

def wrap_to_cell(cart, cell):
    frac = np.linalg.solve(cell.T, cart)
    frac = frac - np.floor(frac)
    return frac @ cell

def ok_distance(newpos, positions, cell):
    frac_new = np.linalg.solve(cell.T, newpos)
    frac_all = np.linalg.solve(cell.T, positions.T).T
    dfrac = frac_all - frac_new
    dfrac -= np.round(dfrac)               # minimum-image in fractional
    dcart = dfrac @ cell
    d2 = np.sum(dcart**2, axis=1)
    return np.all(d2 >= min_dist**2)

def min_dist_to_all(newpos, positions, cell):
    frac_new = np.linalg.solve(cell.T, newpos)
    frac_all = np.linalg.solve(cell.T, positions.T).T
    dfrac = frac_all - frac_new
    dfrac -= np.round(dfrac)
    dcart = dfrac @ cell
    d2 = np.sum(dcart**2, axis=1)
    return float(np.sqrt(np.min(d2)))

def place_diatomic(sym_a, sym_b, bond, cell, positions):
    # try candidate points then random search
    for base in frac_points:
        cart0 = np.array(base) @ cell
        v = rng.normal(size=3); v /= np.linalg.norm(v)
        cart1 = wrap_to_cell(cart0 + v * bond, cell)
        if ok_distance(cart0, positions, cell) and ok_distance(cart1, positions, cell):
            return (sym_a, cart0), (sym_b, cart1)

    for _ in range(8000):
        cart0 = rng.random(3) @ cell
        v = rng.normal(size=3); v /= np.linalg.norm(v)
        cart1 = wrap_to_cell(cart0 + v * bond, cell)
        if ok_distance(cart0, positions, cell) and ok_distance(cart1, positions, cell):
            return (sym_a, cart0), (sym_b, cart1)

    raise RuntimeError("Failed to place diatomic without overlaps. Consider larger box or lower min_dist.")

def place_atom(sym, cell, positions):
    for base in frac_points:
        cart = np.array(base) @ cell
        if ok_distance(cart, positions, cell):
            return sym, cart
    for _ in range(8000):
        cart = rng.random(3) @ cell
        if ok_distance(cart, positions, cell):
            return sym, cart
    raise RuntimeError("Failed to place atom without overlaps. Consider larger box or lower min_dist.")

def pick_best_water_oxygen(atoms):
    # choose oxygen with maximum clearance (best for forming H3O+)
    syms = atoms.get_chemical_symbols()
    O_idx = [i for i,s in enumerate(syms) if s == "O"]
    if not O_idx:
        raise RuntimeError("No oxygen found in POSCAR.")

    cell = atoms.cell.array
    best_i = O_idx[0]
    best_clear = -1.0

    # evaluate clearance ignoring the oxygen itself
    pos = atoms.positions
    for i in O_idx:
        # clearance = min distance from this O to any other atom
        # (exclude itself)
        others = np.delete(pos, i, axis=0)
        clear = min_dist_to_all(pos[i], others, cell)
        if clear > best_clear:
            best_clear = clear
            best_i = i

    return best_i, best_clear

def add_extra_H_near_O(atoms, O_idx):
    cell = atoms.cell.array
    Opos = atoms.positions[O_idx].copy()

    # try many directions + several bond lengths
    for r in oh_try:
        for _ in range(400):
            v = rng.normal(size=3); v /= np.linalg.norm(v)
            Hpos = wrap_to_cell(Opos + v * r, cell)
            if ok_distance(Hpos, atoms.positions, cell):
                atoms += Atoms("H", positions=[Hpos], cell=atoms.cell, pbc=True)
                return True

    return False

def add_extra_H_in_void(atoms):
    # fallback: place H in a safe void; relaxation will “find” hydronium-like environment
    cell = atoms.cell.array
    for _ in range(15000):
        cart = rng.random(3) @ cell
        if ok_distance(cart, atoms.positions, cell):
            atoms += Atoms("H", positions=[cart], cell=atoms.cell, pbc=True)
            return True
    return False

def main():
    atoms = read("POSCAR")
    atoms.set_pbc(True)
    cell = atoms.cell.array

    # Add NO
    a0, a1 = place_diatomic("N", "O", bond_NO, cell, atoms.positions)
    atoms += Atoms([a0[0], a1[0]], positions=[a0[1], a1[1]], cell=atoms.cell, pbc=True)

    # Add CO
    b0, b1 = place_diatomic("C", "O", bond_CO, cell, atoms.positions)
    atoms += Atoms([b0[0], b1[0]], positions=[b0[1], b1[1]], cell=atoms.cell, pbc=True)

    # Add Cl-
    c = place_atom("Cl", cell, atoms.positions)
    atoms += Atoms(c[0], positions=[c[1]], cell=atoms.cell, pbc=True)

    # Add “excess proton” robustly
    O_idx, clear = pick_best_water_oxygen(atoms)
    ok = add_extra_H_near_O(atoms, O_idx)
    if not ok:
        # fallback: void placement
        ok2 = add_extra_H_in_void(atoms)
        if not ok2:
            raise RuntimeError("Failed to place extra H anywhere without overlaps (box too dense).")

    # Reorder species for clean POTCAR order: H O N C Cl
    order = {"H": 0, "O": 1, "N": 2, "C": 3, "Cl": 4}
    idx = sorted(range(len(atoms)), key=lambda i: order.get(atoms[i].symbol, 99))
    atoms = atoms[idx]

    write(out_poscar, atoms, format="vasp", vasp5=True, direct=False)
    counts = {s: atoms.get_chemical_symbols().count(s) for s in ["H","O","N","C","Cl"]}
    print(f"Wrote: {out_poscar}")
    print("Counts:", counts)

if __name__ == "__main__":
    main()

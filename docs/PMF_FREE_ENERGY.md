# PMF and free-energy utilities

The repository now distinguishes two quantities that were both called "free energy" in parts of the historical workflow.

## 1. `OSZICAR` `F=` during AIMD

`analysis/aimd_convergence.py` stitches `chunk_*/OSZICAR`, extracts temperature plus `F=`, `E0=`, and `EK=`, and produces rolling/cumulative convergence diagnostics. This is an equilibration diagnostic. It is **not** the Blue-Moon PMF.

## 2. Blue-Moon PMF from constrained MD

For VASP constrained MD, define constrained coordinates in `ICONST` and enable `LBLUEOUT = .TRUE.`. VASP writes constrained-coordinate data (`cc>`) and Blue-Moon data (`b_m>`) to `REPORT`. The fourth numeric `b_m>` field is the Blue-Moon free-energy-gradient quantity shown by the VASP tutorial.

Workflow:

```bash
# 1. Create the next distance window while translating CO rigidly.
python pmf/stage_distance_window.py previous/POSCAR next/POSCAR \
  --target 3.4 --carbon-index 46 --dopant-index 1 --iconst next/ICONST --backup

# 2. Inspect progress across window folders.
python pmf/check_pmf_progress.py . --glob '**/02_pmf_*' \
  --carbon-index 46 --dopant-index 1

# 3. Extract the reaction-coordinate gradient from REPORT files.
#    For dual constraints, choose the C--dopant row with --constraint-index.
python pmf/extract_blue_moon.py . \
  --constraint-index 1 --discard 1000

# 4. Integrate <dA/dxi> over xi and plot the PMF.
python pmf/integrate_pmf.py blue_moon_summary.csv --reference minimum
```

## Important scientific checks

- Atom indices must be verified for each structure; the scripts intentionally do not guess the dopant atom.
- The staging script identifies the CO oxygen as the O atom nearest the selected carbon and refuses to proceed if its C--O distance is outside a configurable physical range.
- C and O are translated together, preserving the C--O bond while changing the C--surface distance.
- If two hard constraints are present, confirm which `cc>`/`b_m>` row is the reaction coordinate before integrating.
- Choose the sign of the integrated gradient consistently with the reaction-coordinate definition. `integrate_pmf.py --sign -1` is available rather than silently imposing a convention.
- Discard equilibration samples separately for every window before averaging gradients.
- The propagated uncertainty in `integrate_pmf.py` treats window SEM values as independent; correlated trajectories require a more careful statistical treatment.

## Provenance

These four PMF utilities and `analysis/aimd_convergence.py` are **new reconstructed tools** written from the retained workflow behavior and official VASP constrained-MD output conventions. They are not claimed to be byte-for-byte copies of the five unavailable historical `pmf_co__111__*.py` files.

Official VASP references:

- https://vasp.at/wiki/Constrained_molecular_dynamics
- https://vasp.at/tutorials/latest/transition_states/part3/

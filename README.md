# AIMD-HPC Toolkit

[![quality](https://github.com/lald-ai/aimd-hpc-toolkit/actions/workflows/quality.yml/badge.svg)](https://github.com/lald-ai/aimd-hpc-toolkit/actions/workflows/quality.yml)

A practical research-software toolkit for **VASP ab initio molecular dynamics (AIMD)** on HPC systems, with utilities for convergence analysis, constrained Blue-Moon sampling, potential-of-mean-force reconstruction, structure preparation, and calculation recovery.

This repository grew out of computational-catalysis workflows where long AIMD simulations are split across many jobs and must remain restartable, auditable, and scientifically interpretable.

## What this project demonstrates

- Production-oriented **VASP AIMD** analysis across chunked trajectories
- Temperature, energy, drift, rolling-average, and block-convergence diagnostics
- **Blue-Moon constrained molecular dynamics** analysis
- Mean-force extraction and **PMF / free-energy integration**
- Safe staging of adjacent reaction-coordinate windows
- Structure-building utilities for catalytic interfaces
- Recovery tools for malformed or failed VASP calculations
- HPC-aware workflow design and reproducibility checks

## Workflow

```text
Structure / restart geometry
          ↓
      VASP AIMD
          ↓
Convergence diagnostics
          ↓
Constrained MD windows
          ↓
Blue-Moon mean forces
          ↓
PMF / free-energy profile
          ↓
Audit / restart / recovery
```

## Repository layout

| Directory | Purpose |
| --- | --- |
| `analysis/` | AIMD temperature/energy diagnostics and convergence analysis |
| `pmf/` | Blue-Moon window staging, progress checks, gradient extraction, and PMF integration |
| `builders/` | Structure-generation utilities for catalytic-interface workflows |
| `vasp_input/` | ASE/VASP input-generation example |
| `rescue/` | VASP auditing, POSCAR recovery, overlap repair, conservative INCAR patching, and Slurm helpers |
| `docs/` | Workflow notes and scientific-method documentation |
| `scripts/` | Lightweight repository validation |

## Selected tools

### AIMD convergence

`analysis/aimd_convergence.py` stitches chunked `OSZICAR` output and evaluates temperature and energy behavior using rolling means, cumulative averages, block statistics, tail fluctuations, and linear drift.

```bash
python analysis/aimd_convergence.py . \
  --timestep-fs 1.0 \
  --rolling-ps 0.5 1.0 2.5 \
  --block-ps 2.5 \
  --tail-ps 5
```

### Blue-Moon / PMF workflow

Extract constrained mean-force information from VASP `REPORT` files:

```bash
python pmf/extract_blue_moon.py . \
  --constraint-index 1 \
  --discard 1000
```

Integrate the resulting mean force into a free-energy profile:

```bash
python pmf/integrate_pmf.py blue_moon_summary.csv \
  --reference minimum
```

Stage the next constrained distance window while preserving the molecular geometry of the adsorbate:

```bash
python pmf/stage_distance_window.py previous/POSCAR next/POSCAR \
  --target 3.4 \
  --carbon-index 46 \
  --dopant-index 1 \
  --iconst next/ICONST \
  --backup
```

See [`docs/PMF_FREE_ENERGY.md`](docs/PMF_FREE_ENERGY.md) for the PMF workflow and conventions.

### VASP rescue / audit

Audit many VASP calculations without modifying them:

```bash
python rescue/vasp_audit.py "level*_*-ads" --min-dist 0.65 --out audit.csv
```

Preview candidate POSCAR recovery operations:

```bash
python rescue/vasp_rescue_poscar.py "level*_*-ads"
```

## Installation

```bash
python -m pip install -r requirements.txt
```

Some examples depend on external research software or system-specific VASP/Slurm environments. Those dependencies are documented in [`docs/EXTERNAL_INPUTS.md`](docs/EXTERNAL_INPUTS.md).

## Validation

A lightweight GitHub Actions workflow checks Python and shell syntax on each push and pull request. The same checks can be run locally with:

```bash
bash scripts/verify_repo.sh
```

## Research context

The toolkit is aimed at atomistic simulations of catalytic interfaces, especially workflows where long AIMD trajectories, constrained reaction coordinates, and HPC restarts must be managed together. Example paths and account-specific information have been sanitized; large simulation outputs, datasets, and unpublished research results are intentionally not distributed.

Some PMF/convergence utilities were reconstructed from documented workflow behavior when the original cluster sources were unavailable; this is noted in the accompanying documentation rather than presenting them as byte-for-byte historical originals.

## Author

**Dhruv Lal**  
PhD Researcher in Chemical Engineering  
University of Wisconsin–Madison

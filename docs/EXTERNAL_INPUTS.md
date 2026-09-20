# External inputs and environment

These are research utilities, not a turnkey VASP distribution. Typical inputs are existing VASP calculation folders and standard VASP files (`POSCAR`, `CONTCAR`, `POTCAR`, `OUTCAR`, `OSZICAR`, `XDATCAR`, `INCAR`, submission scripts).

## Python dependencies

- Python 3.9+
- NumPy
- ASE
- Matplotlib

`builders/heat_adsorption_setup.py` additionally expects an external `heat` package used by the original research environment. Set `HEAT_APPS` if that package needs an extra import path, plus `HEAT_CLUSTER` and `HEAT_EMAIL` as appropriate for your environment.

`vasp_input/generate_aimd_inputs_ase.py` uses ASE's VASP calculator and therefore assumes a separately configured VASP/ASE environment. Its parameter choices are a historical example, not universal recommendations.

## Safety before batch operations

The rescue utilities default to inspection/dry-run behavior where possible. Review generated reports and backups before applying changes across a large calculation tree. Cluster submission/cancellation helpers should be tested on a small path first.

# Workflow map

This recovered snapshot covers the parts of the original AIMD workflow for which complete source survived.

| Area | Included utilities | Intended use |
|---|---|---|
| Analysis | `chunk_diagnostics.py`, `chunk_temperature.py`, `diagnostics_dat_plot.py`, `oszicar_progress.py` | Parse chunked `OSZICAR` data and visualize temperature/energy progress. |
| Builders | `add_no_co_ph0.py` | Add NO/CO/protonated-water species to a periodic solvent structure while screening short contacts. |
| Builders | `heat_adsorption_setup.py` | Example driver for an external HEAT surface-adsorption package. |
| VASP input | `generate_aimd_inputs_ase.py` | Historical ASE/VASP input-generation example for an AIMD heating segment. |
| Rescue | `vasp_audit.py` | Bulk audit of VASP folders for parseability, species/count consistency, overlaps, completion, and fatal errors. |
| Rescue | `vasp_rescue_poscar.py` | Rebuild malformed/truncated POSCAR geometry from POTCAR/OUTCAR/XDATCAR/OLD_* evidence. |
| Rescue | `vasp_rescue_overlaps.py` | Find short-contact structures and replace them with safe historical geometries. |
| Rescue | `patch_incar_conservative.py` | Apply conservative restart settings to INCAR files. |
| Queue | `slurm_workdir_tools.sh` | Query/cancel/submit Slurm jobs by work directory rather than fragile job-name matching. |

The original repository also contained trajectory stitching/wrapping, PMF staging and checking, production-health checks, queue mapping, additional builders, and system-specific legacy examples. Their filenames are preserved in `ORIGINAL_TARGET_MANIFEST.txt`, but their source bodies are not published unless they were fully recovered.

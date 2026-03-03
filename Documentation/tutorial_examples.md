# EntDetect Tutorial (Examples with Included Test Assets)

This tutorial walks through end-to-end example runs using the small test datasets shipped in this repo under `assets/`.

All commands below assume you are running from the EntDetect repo root.

## Setup (Linux or macOS via Miniconda)

Create the conda environment and install EntDetect in editable mode:

```bash
conda env create -f environment.yml
conda activate entdetect
```

macOS notes:
- If you are on macOS and `conda env create -f environment.yml` fails, try:
  - `conda env create -f environment-mac.yml`
- On Apple Silicon, some scientific packages may not have native osx-arm64 builds for your exact Python version. If you hit solver/build issues, a common workaround is to use the x86_64 conda subdir:

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment-mac.yml
conda activate entdetect
```

## Example 1: Native entanglements (all-atom PDB; heavy-atom contacts)

This uses the all-atom PDB shipped in the repo.

```bash
python scripts/run_nativeNCLE.py \
  --struct assets/test_structures/1zmr.pdb \
  --outdir Testing/1zmr_aapdb/ \
  --organism Human \
  --Accession P0A799 \
  --model EXP \
  --ent_detection_method 3
```

What this produces (under `Testing/1zmr_aapdb/`):
- `Native_GE/` native entanglements
- `Native_HQ_GE/` high-quality subset
- `Native_clustered_HQ_GE/` clustered HQ entanglements
- `Native_clustered_HQ_GE_features/` feature tables

## Example 2: Native entanglements (C-alpha coarse-grained PDB; C-alpha contacts)

This uses the C-alpha-only PDB shipped in the repo. For CG structures, use C-alpha contacts:

```bash
python scripts/run_nativeNCLE.py \
  --struct assets/test_structures/1zmr_model_clean_ca.pdb \
  --outdir Testing/1zmr_cg/ \
  --organism Human \
  --Accession P0A799 \
  --model EXP \
  --ent_detection_method 3 \
  --resolution cg \
  --contacts calpha
```

## Example 3: Trajectory OPs (coarse-grained trajectory; C-alpha contacts)

This example is fully runnable with the bundled CHARMM-style C-alpha topology (`.psf`/`.cor`) and the CG trajectory (`.dcd`).

```bash
python scripts/run_OP_on_simulation_traj.py \
  --Traj 420 \
  --PSF assets/test_simulations/1zmr_model_clean_ca.psf \
  --DCD assets/test_simulations/420_prod.dcd \
  --ID 1zmr_ca_t420 \
  --COR assets/test_simulations/1zmr_model_clean_ca.cor \
  --sec_elements assets/test_simulations/secondary_struc_defs.txt \
  --domain assets/test_simulations/domain_def.dat \
  --start -10 \
  --outdir Testing/1zmr_ca_t420/ \
  --ent_detection_method 3 \
  --resolution cg \
  --contacts calpha
```

Notes:
- `--resolution cg --contacts calpha` matches the typical CG workflow.
- As of v1.1.3, the per-frame trajectory entanglement CSV output (`Traj_GE/*_GE.csv`) is filtered so it only reports contact pairs that are also present in the same-run reference contact set.

Legacy flags are still supported in `run_nativeNCLE.py` (`--cg` and `--Calpha/--calpha`), but the tutorial uses `--resolution/--contacts` for consistency with the trajectory script.

## Example 4: Trajectory OPs (all-atom trajectory; heavy-atom contacts)

The repo includes an example all-atom trajectory file:
- `assets/test_simulations/420_prod_aa.dcd`

However, `scripts/run_OP_on_simulation_traj.py` requires matching topology/reference inputs in CHARMM `.psf`/`.cor` format (plus secondary structure and domain definition files).

If you have an all-atom simulation with matching files, the workflow looks like:

```bash
python scripts/run_OP_on_simulation_traj.py \
  --Traj 420 \
  --PSF /path/to/all_atom.psf \
  --DCD /path/to/all_atom.dcd \
  --ID my_protein_aa \
  --COR /path/to/all_atom.cor \
  --sec_elements /path/to/secondary_struc_defs.txt \
  --domain /path/to/domain_def.dat \
  --start 0 \
  --outdir results/my_protein_aa/ \
  --ent_detection_method 3 \
  --resolution aa \
  --contacts heavy
```

For all-atom runs, `--resolution aa --contacts heavy` ensures native contacts are defined using heavy atoms.

## Other scripts (examples / templates)

The repo includes additional workflow scripts under `scripts/`. Some require large external datasets that are not shipped with this repo.

### Convert CHARMM CG topology to PDB (`convert_cor_psf_to_pdb.py`)

```bash
python scripts/convert_cor_psf_to_pdb.py \
  --psf assets/test_simulations/1zmr_model_clean_ca.psf \
  --cor assets/test_simulations/1zmr_model_clean_ca.cor \
  --outpdb Testing/converted/1zmr_model_clean_ca.pdb
```

### Resolution conversion (`run_change_resolution.py`)

Template (requires domain definitions for your system):

```bash
python scripts/run_change_resolution.py \
  --outdir results/change_resolution/ \
  --pdbfile /path/to/structure.pdb \
  --domain_file /path/to/domain_def.dat \
  --ID my_protein
```

### Non-native entanglement clustering (`run_nonnative_entanglement_clustering.py`)

Template (requires combined GE pickles across trajectories):

```bash
python scripts/run_nonnative_entanglement_clustering.py \
  --outdir results/nonnative_clustering/ \
  --pkl_file_path /path/to/Combined_GE_pickles/ \
  --trajnum2pklfile_path /path/to/trajnum_to_pickle_mapping.txt \
  --traj_dir_prefix /path/to/trajectory_dirs/
```

### MSM building (`run_MSM.py`)

Template (requires OP outputs across trajectories):

```bash
python scripts/run_MSM.py \
  --outdir results/MSM/ \
  --ID my_protein_msm \
  --OPpath /path/to/OP_outputs/ \
  --lagtime 20 \
  --n_large_states 10
```

### Simulation vs experiment comparison (`run_compare_sim2exp.py`)

Template (requires MSM outputs + experimental datasets):

```bash
python scripts/run_compare_sim2exp.py --help
```

### Population modeling / Monte Carlo / folding pathway

These scripts require external datasets. See help for required inputs:

```bash
python scripts/run_population_modeling.py --help
python scripts/run_montecarlo.py --help
python scripts/run_Foldingpathway.py --help
```

# Workflow 3: Compare Structural Ensembles to High-Throughput Experimental Data

← [Back to Master Index](index.md)

---

## Goal

Test whether specific metastable states from simulation are statistically consistent with experimental conformational signals (LiP-MS and XL-MS). Identify and extract representative structures from the best-supported metastable states.

---

## Table of Contents

- [Step 1. Activate your environment and set paths](#step-1-activate-your-environment-and-set-paths)
- [Step 2. Prepare experimental inputs](#step-2-prepare-experimental-inputs)
- [Step 3. Back-map coarse-grained trajectories to all-atom structures](#step-3-back-map-coarse-grained-trajectories-to-all-atom-structures)
- [Step 4. Collect per-trajectory SASA and XP outputs into consolidated arrays](#step-4-collect-per-trajectory-sasa-and-xp-outputs-into-consolidated-arrays)
- [Step 5. Run the LiP-MS / XL-MS consistency test](#step-5-run-the-lip-ms--xl-ms-consistency-test)
- [Step 6. Visualize representative structures](#step-6-visualize-representative-structures)
- [Running the consistency test as a single script](#running-the-consistency-test-as-a-single-script)

---

## Typical runtime

| Step | Runtime |
|------|---------|
| Back-mapping (per frame) | 5–30 minutes per frame |
| Collecting SASA / Jwalk arrays | Minutes |
| LiP-MS / XL-MS consistency test | Minutes to hours |

> **Note:** Per-trajectory SASA and XP (Jwalk) calculations are performed in **Workflow 2**. Step 4 in this workflow collects those per-trajectory files into the consolidated `.npy` arrays required by the consistency test.

> **Strategy for this tutorial:** Pre-collected `SASA.npy` and `Jwalk.npy` are already available in `$OUTDIR`. Steps 5 and 6 are fully runnable from those pre-computed inputs.

---

## Pre-computed reference outputs

All upstream inputs from Workflow 2 are in the DATASTORE:

```
$DATASTORE/outputs/workflow2/
├── MSM/
│   ├── 1ZMR_prod_MSMmapping.csv
│   └── 1ZMR_prod_meta_dist.npy
└── nonnative_clustering/
    └── cluster_data_topoly_linking_number.npz
```

This workflow writes its outputs to:

```
$DATASTORE/outputs/workflow3/
├── SASA.npy                                          # collected in Step 4
├── Jwalk.npy                                         # collected in Step 4
└── MassSpec_ConsistencyTest/
    ├── LiPMS_XLMS_consist_pvalues_metastates_*.xlsx
    ├── consist_signal_struct_data.npz
    └── Consistent_structures_*.xlsx
```

Where:
```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow3
```

> **Strategy for this tutorial:** Pre-collected `SASA.npy` and `Jwalk.npy` are already available in `$OUTDIR`. Steps 5 and 6 are fully runnable from those pre-computed inputs.

---

## Required input files

| File | Path | Notes |
|------|------|-------|
| Native all-atom PDB | `$REFSTRUCT/1zmr_model_clean.pdb` | All-atom reference |
| Cα PSF topology | `$REFSTRUCT/1zmr_model_clean_ca.psf` | Cα topology for back-mapping |
| CG trajectory DCDs | `$CG_TRAJ_DIR/{N}_prod.dcd` (N=1–1000) | For back-mapping frame extraction |
| AA trajectory DCDs | `$AA_TRAJ_DIR/{N}_prod_aa.dcd` (N=1–1000) | All-atom MD trajectories |
| Per-traj SASA files | `$DATASTORE/outputs/workflow2/OP_AA/SASA/{ID}_Traj{N}.SASA` | From Workflow 2 Step 3 |
| Per-traj XP files | `$DATASTORE/outputs/workflow2/OP_AA/XP/{ID}_Traj{N}.XP` | From Workflow 2 Step 3 |
| SASA array (collected) | `$OUTDIR/SASA.npy` | Built in Step 4 |
| Jwalk array (collected) | `$OUTDIR/Jwalk.npy` | Built in Step 4 |
| LiP-MS experimental file | `$DATASTORE/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx` | |
| XL-MS experimental file | `$DATASTORE/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx` | |
| MSM mapping | `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv` | From Workflow 2 |
| MSM meta distribution | `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy` | From Workflow 2 |
| Non-native clustering data | `$DATASTORE/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz` | From Workflow 2 |
| OP directory | `$DATASTORE/outputs/workflow2/OP_AA/` | From Workflow 2 Step 3 |

```bash
REFSTRUCT=$DATASTORE/user_input/reference_structures
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
```

---

## Step 1. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
OUTDIR=$DATASTORE/outputs/workflow3

mkdir -p $OUTDIR/BackMapping \
         $OUTDIR/MassSpec_ConsistencyTest \
         $OUTDIR/logs
```

---

## Step 2. Prepare experimental inputs

The processed experimental files are already available in `$DATASTORE/user_input/experimental_data/`. For your own system, you would:

1. Run your LiP-MS statistical analysis pipeline to identify significantly changed peptides.
2. Map significant peptides back to residue-level assignments.
3. Format the output as an Excel file matching the expected column structure.
4. Repeat for XL-MS cross-link data.

> **Critical:** There are many valid ways to define statistical significance in LiP-MS data. The choice of significance threshold and mapping strategy can strongly affect downstream conclusions. Be consistent and document your decisions.

To verify the experimental input files are accessible:

```python
import pandas as pd

DATASTORE = "/scratch/ims86/EntDetect_Datastore"

lipms = pd.read_excel(f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx")
xlms  = pd.read_excel(f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx")

print("LiP-MS shape:", lipms.shape)
print(lipms.head())
print("\nXL-MS shape:", xlms.shape)
print(xlms.head())
```

---

## Step 3. Back-map coarse-grained trajectories to all-atom structures

> **This step is computationally expensive and has already been completed** for the 1ZMR example system. The all-atom DCD trajectories in `$AA_TRAJ_DIR` are the result. This section describes the procedure for reference and for running on your own system.

Back-mapping reconstructs full all-atom models from Cα-only (coarse-grained) trajectory frames using the native all-atom structure as a template.

### 3a. Save individual Cα PDB frames from a trajectory

The frame extraction step is system-specific and typically done with MDAnalysis or VMD. Example with MDAnalysis:

```python
import MDAnalysis as mda

DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
CG_TRAJ_DIR = f"{DATASTORE}/user_input/cg_trajectories"
OUTDIR      = f"{DATASTORE}/outputs/workflow3"

import os
os.makedirs(f"{OUTDIR}/cg_frames", exist_ok=True)

psf = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean_ca.psf"
dcd = f"{CG_TRAJ_DIR}/420_prod.dcd"

u = mda.Universe(psf, dcd)
protein = u.select_atoms("protein")

for i, ts in enumerate(u.trajectory[-67:]):   # last 67 frames
    out_pdb = f"{OUTDIR}/cg_frames/frame_{i:04d}.pdb"
    protein.write(out_pdb)
```

### 3b. Back-map each Cα frame to all-atom

```python
from EntDetect.change_resolution import BackMapping

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow3"

# cg_pdb : a single CG (Cα-only) PDB frame extracted above
cg_pdb  = f"{OUTDIR}/cg_frames/frame_0000.pdb"
aa_pdb  = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
ID      = "1ZMR"

backMapper = BackMapping(outdir=f"{OUTDIR}/BackMapping")
backMapper.backmap(cg_pdb=cg_pdb, aa_pdb=aa_pdb, ID=ID)
```

### What back-mapping produces

Depending on configuration, outputs include:

- reconstructed all-atom structures (`.pdb`);
- intermediate PD2 / Pulchra outputs;
- OpenMM energy-minimization logs and energy-minimized final structures.

> **After back-mapping:** Inspect a representative subset of reconstructed structures in VMD before proceeding. Collate the validated per-frame PDBs into an all-atom DCD for downstream use.

---

## Step 4. Collect per-trajectory SASA and XP outputs into consolidated arrays

Per-trajectory SASA and XP (Jwalk) data were computed in **Workflow 2** and are stored as one file per trajectory:

```
$DATASTORE/outputs/workflow2/OP_AA/SASA/{ID}_Traj{N}.SASA   ← CSV, long-form, units nm²
$DATASTORE/outputs/workflow2/OP_AA/XP/{ID}_Traj{N}.XP       ← TSV, per-frame Jwalk + XP scores
```

This step reads all of those files and writes two consolidated arrays:

| Output file | Shape | Description |
|-------------|-------|-------------|
| `SASA.npy` | `(n_traj, n_frames, prot_len)` float64 | Per-residue SASA in Å² |
| `Jwalk.npy` | `(n_traj, n_frames)` object | Per-frame dicts of `{pair_key: {'Euclidean': float, 'Jwalk': float}}` |

The pair key format matches the reference: `'RESNUM|CHAIN-RESNUM|CHAIN'` (e.g. `'1|A-124|A'`).

Trajectories whose files are absent (not yet completed) are filled with NaN (SASA) or `None` (Jwalk) and are automatically skipped by `MassSpec` via its existing NaN-filtering logic.

### 4a. Collect SASA

```python
from EntDetect.order_params import CollectOP

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow3"

collector = CollectOP(
    sasa_dir = f"{DATASTORE}/outputs/workflow2/OP_AA/SASA",
    xp_dir   = f"{DATASTORE}/outputs/workflow2/OP_AA/XP",
    outdir   = OUTDIR,
    ID       = "1ZMR",
    n_traj   = 1000,
    n_frames = 335,    # frames stored per trajectory (match Workflow 2 settings)
    prot_len = 387,    # number of residues
)

sasa_npy = collector.collect_SASA()   # writes $OUTDIR/SASA.npy
print("SASA saved to:", sasa_npy)
```

Unit conversion from nm² (CalculateOP output) to Å² (MassSpec input) is applied automatically (×100).

### 4b. Collect Jwalk / XP

```python
jwalk_npy = collector.collect_Jwalk()  # writes $OUTDIR/Jwalk.npy
print("Jwalk saved to:", jwalk_npy)
```

The `SASD` column from each `.XP` file becomes the `'Jwalk'` value in the dict; `Euclidean Distance` becomes `'Euclidean'`.

### Verify the output shapes

```python
import numpy as np

sasa  = np.load(sasa_npy)
jwalk = np.load(jwalk_npy, allow_pickle=True)

print("SASA  shape:", sasa.shape)   # expected (1000, 335, 387)
print("Jwalk shape:", jwalk.shape)  # expected (1000, 335)
print("Example pair keys:", list(jwalk[0, 0].keys())[:3])
```

Compare to pre-collected reference arrays in `$DATASTORE/outputs/workflow3/`.

---

## Step 5. Run the LiP-MS / XL-MS consistency test

This is the key integrative step. It tests whether specific metastable states identified in Workflow 2 are statistically consistent with the experimental signals.

### 5a. Initialize `MassSpec`

```python
from EntDetect.compare_sim2exp import MassSpec

DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
AA_TRAJ_DIR = f"{DATASTORE}/user_input/aa_trajectories"
OUTDIR      = f"{DATASTORE}/outputs/workflow3"

msm_data_file    = f"{DATASTORE}/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv"
meta_dist_file   = f"{DATASTORE}/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy"
LiPMS_exp_file   = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx"
sasa_data_file   = f"{OUTDIR}/SASA.npy"       # collected in Step 4
XLMS_exp_file    = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx"
dist_data_file   = f"{OUTDIR}/Jwalk.npy"      # collected in Step 4
cluster_data_file = f"{DATASTORE}/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz"
OPpath           = f"{DATASTORE}/outputs/workflow2/OP_AA/"
AAdcd_dir        = AA_TRAJ_DIR
native_AA_pdb    = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
outdir           = f"{OUTDIR}/MassSpec_ConsistencyTest"

# Protocol-specific parameters for the 1ZMR / ecPGK example system
state_idx_list   = [4, 6, 8]   # metastable state indices to test
prot_len         = 387          # protein length in residues
last_num_frames  = 335          # number of frames per trajectory analyzed
rm_traj_list     = []           # trajectories excluded in Workflow 2 Step 10
native_state_idx = 9            # index of the native/folded metastable state
start            = 6600         # first frame index used in OP calculations
ID               = "1ZMR"

MS = MassSpec(
    msm_data_file=msm_data_file,
    meta_dist_file=meta_dist_file,
    LiPMS_exp_file=LiPMS_exp_file,
    sasa_data_file=sasa_data_file,
    XLMS_exp_file=XLMS_exp_file,
    dist_data_file=dist_data_file,
    cluster_data_file=cluster_data_file,
    OPpath=OPpath,
    AAdcd_dir=AAdcd_dir,
    native_AA_pdb=native_AA_pdb,
    state_idx_list=state_idx_list,
    prot_len=prot_len,
    last_num_frames=last_num_frames,
    rm_traj_list=rm_traj_list,
    native_state_idx=native_state_idx,
    outdir=outdir,
    ID=ID,
    start=start
)
```

### 5b. Run the consistency test

```python
consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
print("Consistency data:", consist_data_file)
print("Consistency results:", consist_result_file)
```

### 5c. Select representative structures

```python
MS.select_rep_structs(
    consist_data_file,
    consist_result_file,
    total_traj_num_frames=335,   # total frames per trajectory in the full run
    last_num_frames=67            # frames from which representative structs are selected
)
```

### What the consistency test does

For each metastable state in `state_idx_list`, the test evaluates whether:

- LiP-MS signals (elevated protease cleavage) are concentrated in residues that are **more solvent-exposed** in that state than in the native state.
- XL-MS signals (cross-links formed / lost) are consistent with the **inter-residue distances** in that state.

A p-value is computed for each state/experiment combination.

### Expected outputs

| File | Contents |
|------|----------|
| `LiPMS_XLMS_consist_pvalues_metastates_*.xlsx` | per-state p-values for LiP-MS and XL-MS consistency |
| `consist_signal_struct_data.npz` | raw per-residue consistency arrays |
| `Consistent_structures_v8.xlsx` | selected representative structures from consistent states |

Compare to pre-computed reference files in `$DATASTORE/outputs/workflow3/MassSpec_ConsistencyTest/`.

---

## Step 6. Visualize representative structures

Using the structure identifiers in `Consistent_structures_v8.xlsx`, load the corresponding all-atom trajectory frames into VMD:

```tcl
# In the VMD TkConsole — load the all-atom trajectory for the representative trajectory
set traj_num 420
mol new /scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb
mol addfile /scratch/ims86/EntDetect_Datastore/user_input/aa_trajectories/420_prod_aa.dcd
# Navigate to the frame index from the results file
animate goto <frame_index>
```

---

## Running the consistency test as a single script

Steps 4 and 5 (CollectOP + MassSpec consistency test) are automated by `scripts/run_compare_sim2exp.py`. The script accepts `--sasa_dir` / `--xp_dir` to trigger collection, or `--sasa_data_file` / `--dist_data_file` to skip collection and use pre-built arrays directly.

### Collect + run (standard use)

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow3

mkdir -p $OUTDIR/MassSpec_ConsistencyTest $OUTDIR/logs

python scripts/run_compare_sim2exp.py \
    --sasa_dir        $DATASTORE/outputs/workflow2/OP_AA/SASA \
    --xp_dir          $DATASTORE/outputs/workflow2/OP_AA/XP \
    --n_traj          1000 \
    --n_frames        335 \
    --msm_data_file   $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv \
    --meta_dist_file  $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy \
    --LiPMS_exp_file  $DATASTORE/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx \
    --XLMS_exp_file   $DATASTORE/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx \
    --cluster_data_file $DATASTORE/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz \
    --OPpath          $DATASTORE/outputs/workflow2/OP_AA \
    --AAdcd_dir       $DATASTORE/user_input/aa_trajectories \
    --native_AA_pdb   $DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb \
    --state_idx_list  4 6 8 \
    --prot_len        387 \
    --last_num_frames 335 \
    --rm_traj_list    65 75 155 162 199 231 264 286 296 314 354 417 448 472 473 474 577 579 591 703 704 732 758 812 833 870 876 944 967 \
    --native_state_idx 9 \
    --outdir          $OUTDIR/MassSpec_ConsistencyTest \
    --ID              1ZMR \
    --start           6600 \
    --end             -1 \
    --stride          1 \
    --num_perm        1000 \
    --n_boot          100 \
    --lag_frame       20 \
    --nproc           10
```

### Skip collection (pre-built arrays already exist)

If `SASA.npy` and `Jwalk.npy` are already in `$OUTDIR`, pass them directly and omit `--sasa_dir`/`--xp_dir`:

```bash
python scripts/run_compare_sim2exp.py \
    --sasa_data_file  $OUTDIR/SASA.npy \
    --dist_data_file  $OUTDIR/Jwalk.npy \
    --msm_data_file   $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv \
    ...
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `CollectOP` skips many trajectories with "Missing" warnings | Workflow 2 SASA/XP jobs not yet complete | Check job status; re-run missing trajectories before collecting |
| `SASA.npy` shape mismatch (`n_frames` axis wrong) | `n_frames` argument doesn't match how many frames CalculateOP wrote | Inspect one `.SASA` file: `wc -l` ÷ `prot_len` = actual `n_frames` |
| `FileNotFoundError` for `SASA.npy` or `Jwalk.npy` | Step 4 (CollectOP) not yet run, or `$OUTDIR` path wrong | Run Step 4 first or verify `OUTDIR` |
| `MassSpec` raises shape mismatch | `prot_len`, `last_num_frames`, or `start` inconsistent with collected data | Confirm values match those used in Workflow 2 Step 3 |
| Consistency test returns no significant states | `state_idx_list` out of range | Check valid state indices in `1ZMR_prod_MSMmapping.csv` |
| Back-mapping produces clashing structures | Pulchra or PD2 reconstruction failure | Inspect CG input frame; remove frames with very distorted Cα geometry |

---

← [Workflow 2](workflow2_trajectory_analysis.md) | [Back to Master Index](index.md) | Next → [Workflow 4: Population Analysis](workflow4_population.md)

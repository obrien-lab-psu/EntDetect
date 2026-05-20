# Workflow 2: Detect Changes in Entanglement Across Simulation Trajectories

← [Back to Master Index](index.md)

---

## Goal

Compute order parameters (Q, G, K, SASA, XP) across simulation trajectories, cluster non-native entanglement changes, build a Markov state model (MSM), and analyze metastable state behavior.

---

## Typical runtime

| Step | Runtime |
|------|---------|
| Q, G, K for one CG trajectory (nproc=10) | 12–20 hours |
| SASA and XP for one AA trajectory (nproc=10) | 2–6 hours |
| Q, G, K for all 1000 trajectories (cluster) | Hours to days |
| Non-native clustering | Hours (memory-intensive) |
| MSM construction | Minutes |
| MSM statistics / folding pathways | Minutes |

---

## Pre-computed outputs

All 1000 CG and all-atom trajectories have already been analyzed and their outputs are stored in the DATASTORE:

```
$DATASTORE/outputs/workflow2/
├── OP/                        # CG order parameters (all 1000 trajectories)
│   ├── Q/1ZMR_Traj{N}.Q
│   ├── G/1ZMR_Traj{N}.G
│   │   └── Combined_GE/            # Per-trajectory entanglement pkl files
│   │       └── 1ZMR_traj{N}_GE.pkl
│   └── K/K_{N}_prod.dat
└── OP_AA/                     # All-atom order parameters (all 1000 trajectories)
    ├── SASA/1ZMR.SASA
    └── XP/
        ├── Jwalk_results/
        ├── XLresidue_pairs_Full.csv
        └── XLresidue_pairs.txt
```

Where:
```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow2
```

---

## Required input files

| File | Path | Notes |
|------|------|-------|
| Cα PSF topology | `$REFSTRUCT/1zmr_model_clean_ca.psf` | CG Q/G/K |
| Cα COR reference | `$REFSTRUCT/1zmr_model_clean_ca.cor` | CG Q/G/K |
| CG trajectories (all) | `$CG_TRAJ_DIR/{N}_prod.dcd` (N=1–1000) | Full production run |
| All-atom PDB topology | `$REFSTRUCT/1zmr_model_clean.pdb` | SASA/XP |
| AA trajectories (all) | `$AA_TRAJ_DIR/{N}_prod_aa.dcd` (N=1–1000) | Full AA production run |
| Secondary structure defs | `$REFSTRUCT/secondary_struc_defs.txt` | **Required for Q/G/K** |
| Domain boundary file | `$REFSTRUCT/domain_def.dat` | **Required for Q/G/K** |

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
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
REFSTRUCT=$DATASTORE/user_input/reference_structures
OUTDIR=$DATASTORE/outputs/workflow2

mkdir -p $OUTDIR/OP/G $OUTDIR/OP/Q $OUTDIR/OP/K \
         $OUTDIR/OP_AA/SASA $OUTDIR/OP_AA/XP \
         $OUTDIR/nonnative_clustering $OUTDIR/MSM \
         $OUTDIR/MSM_StateProbabilityStats
```

---

## Step 2. Compute order parameters Q, G, K on the CG trajectory

Compute the three canonical order parameters for the coarse-grained trajectory using the parameters from `assets/slurm/scripts/run_OP_traj420.slurm`.

### 2a. Initialize `CalculateOP` for the CG trajectory

```python
from EntDetect.order_params import CalculateOP

DATASTORE    = "/scratch/ims86/EntDetect_Datastore"
REFSTRUCT    = f"{DATASTORE}/user_input/reference_structures"
CG_TRAJ_DIR  = f"{DATASTORE}/user_input/cg_trajectories"
OUTDIR       = f"{DATASTORE}/outputs/workflow2"

Traj         = 420
PSF          = f"{REFSTRUCT}/1zmr_model_clean_ca.psf"
COR          = f"{REFSTRUCT}/1zmr_model_clean_ca.cor"
DCD          = f"{CG_TRAJ_DIR}/420_prod.dcd"
ID           = "1ZMR"
sec_elements = f"{REFSTRUCT}/secondary_struc_defs.txt"
domain       = f"{REFSTRUCT}/domain_def.dat"
outdir       = f"{OUTDIR}/OP"

# start: first frame to include (0-indexed).
# Adjust to skip early equilibration frames for production runs.
start = 0

CalcOP = CalculateOP(
    outdir=outdir,
    Traj=Traj,
    ID=ID,
    psf=PSF,
    cor=COR,
    sec_elements=sec_elements,
    dcd=DCD,
    domain=domain,
    start=start,
    ent_detection_method=1,   # 1 = GLN-only; matches production run settings
)
```

### 2b. Compute Q — fraction of native contacts

```python
Qdata_dict = CalcOP.Q()
```

`Q` measures how many of the native residue–residue contacts present in the reference structure are also present in each trajectory frame. A value near 1.0 indicates a native-like conformation.

**Output:** A `.Q` file in `$OUTDIR/OP/Q/`

### 2c. Compute G — entanglement order parameter

```python
Gdata_dict = CalcOP.G(topoly=False, Calpha=True, CG=True, nproc=10, chunk_frames=100, chunk_suffix='_chunk')
```

`G` captures the fraction of native contacts that exhibit a **change in entanglement state** relative to the native structure.

| Argument | Value | Meaning |
|----------|-------|----------|
| `topoly` | `False` | Use GLN-only workflow (no Topoly linking numbers) |
| `Calpha` | `True` | Use Cα-defined contacts (appropriate for CG trajectories) |
| `CG` | `True` | Input trajectory is coarse-grained |
| `nproc` | `10` | Number of CPU cores to use |
| `chunk_frames` | `100` | Write intermediate results every N frames (reduces memory usage on large trajectories) |
| `chunk_suffix` | `'_chunk'` | Suffix for chunked output files |

**Output:** A `.G` file and per-frame entanglement metadata `.pkl` in `$OUTDIR/OP/G/`

> **Runtime note:** `G` is the most expensive order parameter to compute. Expect 12–20 hours per CG trajectory at ~6700 frames on 10 cores. Submit via the cluster for a full 1000-trajectory run (see [Running the full script](#running-the-full-workflow-as-a-single-script)).

### 2d. Compute K — mirror-symmetry order parameter

```python
Kdata_dict = CalcOP.K()
```

`K` detects frames where the protein has adopted a **mirror-image conformation** relative to the native structure. These frames are artifacts that must be removed before clustering.

**Output:** A `K_*.dat` file in `$OUTDIR/OP/K/`

---

## Step 3. Compute order parameters SASA and XP on the all-atom trajectory

> **Why separate from Q, G, K?** The CG Cα-only representation does not carry enough atomic detail for accurate SASA (Shrake-Rupley requires explicit side-chain/backbone atoms) or cross-link SASD (Jwalk uses solvent-accessible surface geometry). The AA trajectories are produced by back-mapping the CG frames. Q, G, K are computed only on CG trajectories (Step 2).

### 3a. Initialize `CalculateOP` for the all-atom trajectory

```python
from EntDetect.order_params import CalculateOP

DATASTORE    = "/scratch/ims86/EntDetect_Datastore"
REFSTRUCT    = f"{DATASTORE}/user_input/reference_structures"
AA_TRAJ_DIR  = f"{DATASTORE}/user_input/aa_trajectories"
OUTDIR       = f"{DATASTORE}/outputs/workflow2"

Traj     = 420
AA_PDB   = f"{REFSTRUCT}/1zmr_model_clean.pdb"    # all-atom topology / reference
AA_DCD   = f"{AA_TRAJ_DIR}/420_prod_aa.dcd"
ID       = "1ZMR"

CalcOP = CalculateOP(
    outdir=f"{OUTDIR}/OP_AA",
    Traj=Traj,
    ID=ID,
    psf=AA_PDB,
    cor=AA_PDB,
    dcd=AA_DCD,
    sec_elements=f"{REFSTRUCT}/secondary_struc_defs.txt",
    domain=f"{REFSTRUCT}/domain_def.dat",
    start=0,
)
```

### 3b. Compute SASA — solvent-accessible surface area per residue

```python
SASAdata_dict = CalcOP.SASA()
```

Uses the Shrake-Rupley algorithm (mdtraj, `probe_radius=0.14 nm`, 1000 sphere points) to compute per-residue SASA for every frame. Used downstream (Workflow 3) to test LiP-MS signals: high SASA residues are more accessible to the protease.

| Key | Contents |
|-----|----------|
| `outfile` | Path to `{ID}.SASA` CSV written to `$OUTDIR/OP_AA/SASA/` |
| `result`  | DataFrame with columns `Time(ns)`, `Frame`, `resid`, `SASA(nm^2)` |

### 3c. Compute XP — Jwalk cross-link probability

```python
XPdata_dict = CalcOP.XP(pdb=AA_PDB)
```

`XP` computes the **solvent-accessible surface distance (SASD)** between all pairs of cross-linkable residue types (K, S, T, Y, M) and converts each SASD to a cross-link probability score using a Gaussian parameterised on K–K linker geometry. Used downstream to test XL-MS signals.

**Output:** Per-residue-pair XP scores in `$OUTDIR/OP_AA/XP/Jwalk_results/{stem}_crosslink_list.txt`

---

## Running the order parameter analysis as a single script

Order parameters are computed by `scripts/run_OP_on_simulation_traj.py`. Q, G, K are computed for CG trajectories; SASA and XP are computed for AA trajectories. This is the script used for the full 1000-trajectory production run.

### CG run — Q, G, K

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

python scripts/run_OP_on_simulation_traj.py \
    --Traj 420 \
    --PSF  $REFSTRUCT/1zmr_model_clean_ca.psf \
    --COR  $REFSTRUCT/1zmr_model_clean_ca.cor \
    --DCD  $DATASTORE/user_input/cg_trajectories/420_prod.dcd \
    --resolution cg \
    --contacts calpha \
    --ID   1ZMR \
    --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
    --domain       $REFSTRUCT/domain_def.dat \
    --outdir       $DATASTORE/outputs/workflow2/OP \
    --logdir       $DATASTORE/outputs/workflow2/OP/logs \
    --start        0 \
    --ent_detection_method 1 \
    --nproc        10 \
    --ops Q G K \
    --no_topoly \
    --chunk_frames 100 \
    --chunk_suffix _chunk
```

### AA run — SASA, XP

```bash
python scripts/run_OP_on_simulation_traj.py \
    --Traj   420 \
    --ID     1ZMR \
    --PSF    $REFSTRUCT/1zmr_model_clean.pdb \
    --COR    $REFSTRUCT/1zmr_model_clean.pdb \
    --DCD    $DATASTORE/user_input/aa_trajectories/420_prod_aa.dcd \
    --resolution aa \
    --contacts calpha \
    --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
    --domain       $REFSTRUCT/domain_def.dat \
    --outdir $DATASTORE/outputs/workflow2/OP_AA \
    --logdir $DATASTORE/outputs/workflow2/OP_AA/logs \
    --start  0 \
    --nproc        10 \
    --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
    --ops SASA XP
```
---

## Step 4. Identify and remove artificial mirror conformations

Before clustering entanglement changes, remove trajectories or frames that are mirror-image artifacts.

### Recommended procedure

1. Load Q and K time series for all trajectories.
2. Identify trajectories where K is persistently low (<=0.6) and Q is (>=0.2).
3. Visually inspect flagged frames in VMD for mirror image conformation. 
4. Record the trajectory numbers that are confirmed mirrors in a list (used as `rm_traj_list` in [Workflow 3: Sim-to-Experiment](workflow3_sim2exp.md)).

In this tutorial we identified the mirror artifact trajectories as: 65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967   

> **Important:** The cutoff values for flagging mirror conformations must be tuned for your system by examining the Q and K distributions. 

---

## Step 5. Cluster non-native entanglement changes

This step identifies non-redundant changes in entanglement topology across all trajectories. The pkl file paths are specified in the trajectory-to-pkl mapping CSV (`trajnum2pklfile_path`), which serves as the single source of truth for which pkl files to analyze.

### `trajnum2file.txt` — maps trajectory numbers to pkl files

This file maps trajectory numbers to their **G-order-parameter pkl files** (not DCD files). The format is comma-separated:

```
trajnum,pklfile
<trajectory_number>,<path_to_pkl_file>
```

The pre-populated copy at `$DATASTORE/user_input/metadata/trajnum2file.txt` maps all 1000 trajectories to the Combined_GE pkl files in the DATASTORE. It can be regenerated if needed:

```bash
echo "trajnum,pklfile" > $DATASTORE/user_input/metadata/trajnum2file.txt
for pkl in $DATASTORE/outputs/workflow2/OP/G/Combined_GE/*.pkl; do
    num=$(basename $pkl | sed 's/1ZMR_traj\([0-9]*\)_GE.pkl/\1/')
    echo "$num,$pkl"
done >> $DATASTORE/user_input/metadata/trajnum2file.txt
```
Keep in mind that if you used the `--chunk_frames` argument when calculating G that the .pkl files will be chunked and you will need to specify which chunk to use for each traj. EntDetect currently do not support using multiple chunks per trajectory number. 

> **Memory warning:** This step can require tens of gigabytes of RAM. Run on a high-memory node or reduce the number of frames per trajectory in the clustering pool. 

### 5a. Build the clustering with Python

```python
from EntDetect.clustering import ClusterNonNativeEntanglements

DATASTORE            = "/scratch/ims86/EntDetect_Datastore"
OUTDIR               = f"{DATASTORE}/outputs/workflow2"
CG_TRAJ_DIR          = f"{DATASTORE}/user_input/cg_trajectories"

trajnum2pklfile_path = f"{DATASTORE}/user_input/metadata/trajnum2file.txt"
traj_dir_prefix      = CG_TRAJ_DIR
outdir               = f"{OUTDIR}/nonnative_clustering"

clustering_NNents = ClusterNonNativeEntanglements(
    trajnum2pklfile_path=trajnum2pklfile_path,
    traj_dir_prefix=traj_dir_prefix,
    outdir=outdir,
)

# start_frame, end_frame: frame range to analyze (0-indexed)
# Default is to use all frames; here we use the last 67 frames for a faster demo
clustering_NNents.cluster(start_frame=6600, end_frame=6667)
```

### 5b. Using the command-line interface

For convenience, use the `scripts/run_nonnative_entanglement_clustering.py` script directly (see "Running non-native clustering as a single script" below for full details).

### 5c. Inspecting clustering results

```python
import pandas as pd
import numpy as np

DATASTORE  = "/scratch/ims86/EntDetect_Datastore"
clust_dir  = f"{DATASTORE}/outputs/workflow2/nonnative_clustering"

# Representative entanglement changes
rep_df = pd.read_csv(f"{clust_dir}/rep_chg_ent_topoly_linking_number.csv")
print(f"Number of representative entanglement changes: {len(rep_df)}")
print(rep_df.head())

# Per-frame structural assignments
chg_df = pd.read_csv(f"{clust_dir}/chg_ent_struct_topoly_linking_number.csv")
print(chg_df.head())

# Cluster data array
cluster_data = np.load(f"{clust_dir}/cluster_data_topoly_linking_number.npz", allow_pickle=True)
print("Available arrays:", list(cluster_data.keys()))
```

### 5d. Expected outputs

| File | Contents |
|------|----------|
| `rep_chg_ent_topoly_linking_number.csv` | Representative entanglement changes (one per cluster) |
| `chg_ent_struct_topoly_linking_number.csv` | Per-frame cluster assignment |
| `cluster_data_topoly_linking_number.npz` | Compressed cluster data array |
| `cluster_tree_topoly_linking_number.dat` | Text representation of the clustering hierarchy |
| `rep_chg_ent_list_topoly_linking_number.pkl` | List of representative entanglement objects |
| `chg_ent_topoly_linking_number_distribution.pdf` | Distribution plot of loop/crossing residues |

---

## Running non-native clustering as a single script

The `scripts/run_nonnative_entanglement_clustering.py` script automates the clustering workflow. For production runs with all 1000 trajectories, submit via SLURM using `assets/slurm/scripts/nonNativeClustering.slurm`.

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore

mkdir -p $DATASTORE/outputs/workflow2/nonnative_clustering
mkdir -p $DATASTORE/outputs/workflow2/nonnative_clustering/logs

python scripts/run_nonnative_entanglement_clustering.py \
    --outdir               $DATASTORE/outputs/workflow2/nonnative_clustering \
    --trajnum2pklfile_path $DATASTORE/user_input/metadata/trajnum2file.txt \
    --traj_dir_prefix      $DATASTORE/user_input/cg_trajectories \
    --start_frame          6600 \
    --end_frame            6667 \
    --logdir               $DATASTORE/outputs/workflow2/nonnative_clustering/logs \
    --nproc                4 \
    --log_level            DEBUG
```

For full production clustering (all frames, all trajectories):

```bash
sbatch assets/slurm/scripts/nonNativeClustering.slurm
```

---

## Step 6. Build a Markov state model (MSM)

Organize simulation frames into microstates and metastable states using the full-trajectory order-parameter data using parameters from `assets/slurm/scripts/run_MSM.slurm`.

### 6a. Build the MSM with Python

```python
from EntDetect.clustering import MSMNonNativeEntanglementClustering

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

outdir         = f"{OUTDIR}/MSM"
ID             = "1ZMR_prod"

# OPpath must point to a directory containing Q/, G/, K/ subdirectories
# with per-trajectory OP files for ALL trajectories.
OPpath         = f"{OUTDIR}/OP/"
n_large_states = 10   # number of metastable macro-states requested
lagtime        = 20   # lag time in frames

MSM = MSMNonNativeEntanglementClustering(
    outdir=outdir,
    ID=ID,
    OPpath=OPpath,
    start=0,
    n_large_states=n_large_states,
    lagtime=lagtime
)
MSM.run()
```

### 6b. Using the command-line interface

For convenience, use the `scripts/run_MSM.py` script directly (see "Running MSM construction as a single script" below for full details).

### 6c. Critical notes

- Try **multiple values** of `n_large_states` (e.g., 5, 10, 15). The protocol notes ≤15 often works well.
- The final number of states may be lower than requested if empty states are discarded.
- The default lag time of 1 frame is appropriate for visualization and exploratory grouping. For kinetic interpretation, test for Markovian behavior explicitly and choose lag time accordingly.

### 6d. Expected outputs

| File | Contents |
|------|----------|
| `1ZMR_prod_MSMmapping.csv` | Per-frame microstate and metastable-state assignments |
| `1ZMR_prod_meta_set.csv` | Metastable-state summary |
| `1ZMR_prod_meta_dist.npy` | Metastable-state probability distribution |
| `1ZMR_prod_StateAndFEplot.png` | Order-parameter landscape and state assignments |

---

## Running MSM construction as a single script

The `scripts/run_MSM.py` script automates MSM construction. For production runs, submit via SLURM using `assets/slurm/scripts/run_MSM.slurm`.

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore

mkdir -p $DATASTORE/outputs/workflow2/MSM
mkdir -p $DATASTORE/outputs/workflow2/MSM/logs

python scripts/run_MSM.py \
    --outdir         $DATASTORE/outputs/workflow2/MSM \
    --ID             1ZMR_prod \
    --OPpath         $DATASTORE/outputs/workflow2/OP/ \
    --start          0 \
    --n_large_states 10 \
    --lagtime        20 \
    --logdir         $DATASTORE/outputs/workflow2/MSM/logs
```

For full production MSM:

```bash
sbatch assets/slurm/scripts/run_MSM.slurm
```

---

## Step 7. Label MSM Data — Define Your Analysis Cases

Before computing statistics or folding pathways, each trajectory needs a **type label** (e.g. `A` or `B`) representing the biological comparison of interest. Any column with consistent string labels per trajectory can serve as the type column downstream.

We demonstrate two contrasting cases:

| Case | Labelling rule | Expected signal |
|------|---------------|----------------|
| **Case 1 — Biologically-informed** | A = max Q ≥ 0.80 **and** max G ≤ 0.05 (native-like, non-entangled); B = all others | High JS divergence |
| **Case 2 — Random (negative control)** | A/B assigned randomly per trajectory (seed=42) | JS divergence ≈ 0 |

> **Define your own cases:** Replace the labelling rule below with whatever biological comparison makes sense for your data — e.g. different temperature conditions, mutation variants, or folding outcomes. The only requirement is a column with consistent string labels per trajectory.

### 7a. Generate annotated MSM mapping files

```python
import pandas as pd
import numpy as np

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

msm_mapping = pd.read_csv(f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping.csv")
```

**Case 1 — Biologically-informed split (Q ≥ 0.8 and G ≤ 0.05):**

```python
# A = trajectories whose max Q >= 0.80 AND max G <= 0.05 (native-like and non-entangled)
# B = all others (misfolded or with significant entanglement changes)
traj_stats = msm_mapping.groupby('traj').agg(max_Q=('Q', 'max'), max_G=('G', 'max'))
native_trajs = traj_stats[(traj_stats['max_Q'] >= 0.80) & (traj_stats['max_G'] <= 0.05)].index
msm_mapping['traj_type_QG_native'] = msm_mapping['traj'].isin(native_trajs).map({True: 'A', False: 'B'})

annotated_file_case1 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_QG_native.csv"
msm_mapping[['traj', 'frame', 'microstate', 'metastablestate', 'Q', 'G', 'traj_type_QG_native']].to_csv(
    annotated_file_case1, index=False)
print("Case 1 — trajectory type distribution:")
print(msm_mapping.groupby('traj')['traj_type_QG_native'].first().value_counts())
```

> **Note:** The Q and G thresholds are system-specific. Adjust to match the folding and entanglement distributions of your system.

**Case 2 — Random split (negative control):**

```python
# Randomly assign A/B labels to trajectories (fixed seed for reproducibility)
rng = np.random.default_rng(seed=42)
all_trajs = msm_mapping['traj'].unique()
random_labels = dict(zip(all_trajs, rng.choice(['A', 'B'], size=len(all_trajs))))
msm_mapping['traj_type_random'] = msm_mapping['traj'].map(random_labels)

annotated_file_case2 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_random.csv"
msm_mapping[['traj', 'frame', 'microstate', 'metastablestate', 'Q', 'G', 'traj_type_random']].to_csv(
    annotated_file_case2, index=False)
print("Case 2 — random type distribution:")
print(msm_mapping.groupby('traj')['traj_type_random'].first().value_counts())
```

### 7b. Visualize representative metastable structures

Identify the representative frame for each metastable state from `1ZMR_prod_MSMmapping.csv`, extract the corresponding structure from the trajectory, and inspect it in VMD.

---

## Step 8. Metastable-State Probability Evolution (MSMStats)

Using `MSMStats`, compute how each trajectory population (A and B) distributes across metastable states over simulation time. Run separately for each labelled case to compare the biological signal against the random baseline.

**Input:** Annotated MSM mapping CSV files from Step 7  
**Output:** `$OUTDIR/MSM_StateProbabilityStats_{case}/` — probability plots and summary tables

### 8a. Compute state probability statistics

```python
from EntDetect.statistics import MSMStats
import os

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

meta_set_file  = f"{OUTDIR}/MSM/1ZMR_prod_meta_set.csv"
traj_type_list = ['A', 'B']
rm_traj_list   = []
```

**Case 1 — Biologically-informed split:**

```python
outdir_stats_case1 = f"{OUTDIR}/MSM_StateProbabilityStats_QG_native"
os.makedirs(outdir_stats_case1, exist_ok=True)

MS1 = MSMStats(
    outdir=outdir_stats_case1,
    msm_data_file=f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_QG_native.csv",
    meta_set_file=meta_set_file,
    tarj_type_col='traj_type_QG_native',
    rm_traj_list=rm_traj_list,
    traj_type_list=traj_type_list,
)

df1 = MS1.StateProbabilityStats()
MS1.Plot_StateProbabilityStats(df=df1)
```

**Case 2 — Random split (negative control):**

```python
outdir_stats_case2 = f"{OUTDIR}/MSM_StateProbabilityStats_random"
os.makedirs(outdir_stats_case2, exist_ok=True)

MS2 = MSMStats(
    outdir=outdir_stats_case2,
    msm_data_file=f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_random.csv",
    meta_set_file=meta_set_file,
    tarj_type_col='traj_type_random',
    rm_traj_list=rm_traj_list,
    traj_type_list=traj_type_list,
)

df2 = MS2.StateProbabilityStats()
MS2.Plot_StateProbabilityStats(df=df2)
```

**Output:** Time series plots showing the population (probability) of each metastable state over simulation time, one plot per case. Useful for identifying which states are transiently vs persistently populated by each trajectory subpopulation.

---

## Step 9. Folding Pathways and Jensen-Shannon Divergence

Analyse how trajectories transition between metastable states and quantify the divergence between the two populations using Jensen-Shannon (JS) divergence.

We demonstrate **two contrasting cases** to illustrate what the JS signal looks like under biologically meaningful vs. meaningless partitioning.

| Case | Labelling rule | Expected JS signal |
|------|---------------|-------------------|
| **Case 1 — Biologically-informed split** | A = correctly folded (max Q ≥ 0.8 **and** max G ≤ 0.05); B = misfolded/entangled | High divergence — A and B follow distinct metastable-state progressions |
| **Case 2 — Random split (negative control)** | A/B assigned randomly per trajectory (fixed seed) | Near-zero divergence — populations are statistically identical by construction |

Running both cases back-to-back makes the biological signal immediately recognisable against baseline noise.

### 9a. Compute folding pathways and Jensen-Shannon divergence

Run `FoldingPathwayStats` for each case. `post_trans()` traces state-to-state transitions for each trajectory, removing loops to yield the minimal directed pathway. `JS_divergence()` computes a windowed Jensen-Shannon divergence between the two populations over simulation time.

**Case 1 — Biologically-informed split:**

```python
from EntDetect.statistics import FoldingPathwayStats

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

msm_data_case1 = pd.read_csv(f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_QG_native.csv")
meta_set_file  = f"{OUTDIR}/MSM/1ZMR_prod_meta_set.csv"
outdir_case1   = f"{OUTDIR}/FoldingPathway_QG_native"

FP1 = FoldingPathwayStats(
    msm_data=msm_data_case1,
    meta_set_file=meta_set_file,
    tarj_type_col='traj_type_QG_native',
    traj_type_list=['A', 'B'],
    outdir=outdir_case1,
    rm_traj_list=[],
)

folding_pathways_case1 = FP1.post_trans()
FP1.JS_divergence()
```

**Case 2 — Random split:**

```python
msm_data_case2 = pd.read_csv(f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_random.csv")
outdir_case2   = f"{OUTDIR}/FoldingPathway_random"

FP2 = FoldingPathwayStats(
    msm_data=msm_data_case2,
    meta_set_file=meta_set_file,
    tarj_type_col='traj_type_random',
    traj_type_list=['A', 'B'],
    outdir=outdir_case2,
    rm_traj_list=[],
)

folding_pathways_case2 = FP2.post_trans()
FP2.JS_divergence()
```

| JS divergence | Interpretation |
|---------------|----------------|
| Near 0 | A and B explore similar state distributions |
| Near 1 | A and B have divergent state usage |

**Expected outcome:**
- Case 1 should show **elevated JS divergence** — the native-like (A) and misfolded (B) populations traverse metastable states differently.
- Case 2 should show **JS divergence near 0** throughout — random labels produce no systematic separation.

### 9b. Using the command-line interface

Use `scripts/run_Foldingpathway.py` for each case (see "Running folding pathway analysis as a single script" below):

```bash
# Case 1 — biologically-informed split
python scripts/run_Foldingpathway.py \
    --msm_data_file $OUTDIR/MSM/1ZMR_prod_MSMmapping_QG_native.csv \
    --meta_set_file $OUTDIR/MSM/1ZMR_prod_meta_set.csv \
    --traj_type_col traj_type_QG_native \
    --traj_type_list A B \
    --outdir        $OUTDIR/FoldingPathway_QG_native

# Case 2 — random split (negative control)
python scripts/run_Foldingpathway.py \
    --msm_data_file $OUTDIR/MSM/1ZMR_prod_MSMmapping_random.csv \
    --meta_set_file $OUTDIR/MSM/1ZMR_prod_meta_set.csv \
    --traj_type_col traj_type_random \
    --traj_type_list A B \
    --outdir        $OUTDIR/FoldingPathway_random
```

### 9c. Expected outputs

Each case produces the same file set in its respective output directory:

| File | Contents |
|------|----------|
| `FoldingPathways_metastablestate_A-B.csv` | Per-type folding pathway probabilities |
| `JS_div_metastablestate_A-B.dat` | Windowed JS divergence time series |

---

## Running folding pathway analysis as a single script

The `scripts/run_Foldingpathway.py` script computes both folding pathways and JS divergence in a single call. It requires the MSM mapping CSV to already contain the trajectory-type column (see Step 7). Run it separately for each case.

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow2

# Case 1 — biologically-informed split (Q >= 0.8 and G <= 0.05)
mkdir -p $OUTDIR/FoldingPathway_QG_native/logs

python scripts/run_Foldingpathway.py \
    --msm_data_file $OUTDIR/MSM/1ZMR_prod_MSMmapping_QG_native.csv \
    --meta_set_file $OUTDIR/MSM/1ZMR_prod_meta_set.csv \
    --traj_type_col traj_type_QG_native \
    --traj_type_list A B \
    --outdir        $OUTDIR/FoldingPathway_QG_native \
    --logdir        $OUTDIR/FoldingPathway_QG_native/logs \
    --log_level     INFO

# Case 2 — random split (negative control)
mkdir -p $OUTDIR/FoldingPathway_random/logs

python scripts/run_Foldingpathway.py \
    --msm_data_file $OUTDIR/MSM/1ZMR_prod_MSMmapping_random.csv \
    --meta_set_file $OUTDIR/MSM/1ZMR_prod_meta_set.csv \
    --traj_type_col traj_type_random \
    --traj_type_list A B \
    --outdir        $OUTDIR/FoldingPathway_random \
    --logdir        $OUTDIR/FoldingPathway_random/logs \
    --log_level     INFO
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `G()` hangs or is very slow | Normal for long trajectories; CG + Topoly is expensive | Reduce `nproc` if memory-limited; use a shorter demo trajectory |
| `SASA()` returns NaN for some frames | Bad backmapped AA structure or topology mismatch | Inspect suspect frames in VMD; these frames are filtered in Workflow 3 |
| Jwalk error: PDB not found | `xp_pdb` path incorrect or file missing | Verify `1zmr_model_clean.pdb` is in `$DATASTORE/user_input/reference_structures/` |
| Jwalk error: `freesasa` not found | Package not installed in env | `pip install freesasa` inside the `entdetect` conda env |
| `Combined_GE/` is empty | Full OP run not yet complete | Use pre-computed results from `$DATASTORE/outputs/workflow2/OP/G/Combined_GE/` |
| MSM produces fewer states than `n_large_states` | Empty states discarded | Normal; try a higher `n_large_states` |
| `secondary_struc_defs.txt` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |
| `domain_def.dat` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |

---

← [Workflow 1](workflow1_native_ncle.md) | [Back to Master Index](index.md) | Next → [Workflow 3: Sim-to-Experiment](workflow3_sim2exp.md)

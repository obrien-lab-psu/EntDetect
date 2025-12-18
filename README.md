# EntDetect

A comprehensive Python package for studying non-covalent lasso entanglements in protein folding through molecular dynamics simulations and experimental data analysis.

## Overview

EntDetect provides a complete toolkit for analyzing protein entanglements across multiple scales - from individual structures to large-scale proteomic datasets. The package enables researchers to:

- **Identify and characterize** native entanglements in protein structures
- **Calculate order parameters** for simulation trajectories (Q, G, K, SASA, Jwalk, XP)
- **Build Markov State Models** from coarse-grained simulation ensembles
- **Compare simulations to experiments** using LiP-MS and XL-MS data
- **Perform population-level analysis** across heterogeneous protein datasets
- **Coarse-grain and back-map** between atomic resolutions

## Key Features

- **Multi-scale Analysis**: From single proteins to proteome-wide studies
- **Experimental Integration**: Direct comparison with mass spectrometry data
- **Advanced Clustering**: Non-native entanglement clustering and MSM construction
- **Statistical Methods**: Monte Carlo simulations and logistic regression modeling
- **Flexible Resolution**: Seamless conversion between all-atom and coarse-grained representations

## Installation

Create a new conda environment using the provided environment file:

```bash
conda env create -f environment.yml --name EntDetect_env
conda activate EntDetect_env
```

Install the package in the EntDetect directory:

```bash
pip install .
```

## Quick Start

### Basic Order Parameter Calculation

Calculate fraction of native contacts (Q) and entanglement changes (G) for a simulation trajectory:

```bash
python scripts/run_OP_on_simulation_traj.py \
  --Traj 1 \
  --PSF /path/to/structure.psf \
  --DCD /path/to/trajectory.dcd \
  --ID protein_name \
  --COR /path/to/structure.cor \
  --sec_elements /path/to/secondary_structure.txt \
  --domain /path/to/domain_def.dat \
  --outdir results/OP_analysis/ \
  --start 0
```

### Native Entanglement Analysis

Identify and cluster native entanglements in a protein structure:

```bash
python scripts/run_nativeNCLE.py \
  --struct /path/to/structure.pdb \
  --outdir results/native_entanglements/ \
  --ID protein_name \
  --organism Ecoli
```

## Complete Workflow Scripts
## Complete Workflow Scripts

### 1. Order Parameter Analysis

#### Simulation Trajectory Analysis (`run_OP_on_simulation_traj.py`)
Calculate Q, G, and K order parameters for coarse-grained simulation trajectories:

```bash
python scripts/run_OP_on_simulation_traj.py \
  --Traj 1 \
  --PSF /path/to/structure.psf \
  --DCD /path/to/trajectory.dcd \
  --ID 1ZMR \
  --COR /path/to/structure.cor \
  --sec_elements /path/to/secondary_structure.txt \
  --domain /path/to/domain_def.dat \
  --outdir results/OP_analysis/ \
  --start 0
```

**Outputs**: Q (native contacts), G (entanglement changes), K (mirror symmetry)

### 2. Entanglement Analysis

#### Native Entanglement Identification (`run_nativeNCLE.py`)
Identify, cluster, and characterize native entanglements in protein structures:

```bash
python scripts/run_nativeNCLE.py \
  --struct /path/to/structure.pdb \
  --outdir results/native_analysis/ \
  --ID protein_name \
  --organism Ecoli
```

**Outputs**: Native entanglements, high-quality filtered entanglements, clustered entanglements, entanglement features

#### Non-Native Entanglement Clustering (`run_nonnative_entanglement_clustering.py`)
Cluster non-native entanglements across multiple simulation trajectories:

```bash
python scripts/run_nonnative_entanglement_clustering.py \
  --outdir results/nonnative_clustering/ \
  --pkl_file_path /path/to/pkl/files/ \
  --trajnum2pklfile_path /path/to/mapping.txt \
  --traj_dir_prefix /path/to/trajectories/
```

### 3. Resolution Conversion

#### Coarse-Graining and Back-Mapping (`run_change_resolution.py`)
Convert between all-atom and coarse-grained representations:

```bash
python scripts/run_change_resolution.py \
  --outdir results/coarse_graining/ \
  --pdbfile /path/to/structure.pdb \
  --nscal 2 \
  --domain_file /path/to/domain_def.dat \
  --ID protein_name
```

**Outputs**: Coarse-grained structure, force field files, back-mapped all-atom structure

### 4. Markov State Model Construction

#### MSM Building (`run_MSM.py`)
Build Markov State Models from ensemble simulation data:

```bash
python scripts/run_MSM.py \
  --outdir results/MSM/ \
  --ID protein_msm \
  --OPpath /path/to/OP/data/ \
  --start 0 \
  --n_large_states 10 \
  --lagtime 20 \
  --rm_traj_list 65 75 155
```

**Outputs**: MSM states, transition matrices, representative structures

### 5. Experimental Comparison

#### Simulation-Experiment Consistency (`run_compare_sim2exp.py`)
Compare MSM states with LiP-MS and XL-MS experimental data:

```bash
python scripts/run_compare_sim2exp.py \
  --msm_data_file /path/to/msm_data.csv \
  --meta_dist_file /path/to/meta_distances.csv \
  --LiPMS_exp_file /path/to/LiPMS_data.csv \
  --XLMS_exp_file /path/to/XLMS_data.csv \
  --sasa_data_file /path/to/sasa_data.csv \
  --dist_data_file /path/to/distances.csv \
  --cluster_data_file /path/to/clusters.csv \
  --OPpath /path/to/OP/ \
  --AAdcd_dir /path/to/trajectories/ \
  --native_AA_pdb /path/to/native.pdb \
  --state_idx_list 1 2 3 4 5 \
  --prot_len 390 \
  --last_num_frames 100 \
  --rm_traj_list 65 75 155 \
  --native_state_idx 0 \
  --outdir results/consistency/ \
  --ID protein_comparison \
  --start 0 \
  --end 1000 \
  --stride 1 \
  --num_perm 1000 \
  --n_boot 100 \
  --lag_frame 20 \
  --nproc 10
```

### 6. Population-Level Analysis

#### Proteome Logistic Regression (`run_population_modeling.py`)
Analyze entanglement associations across heterogeneous protein datasets:

```bash
python scripts/run_population_modeling.py \
  --dataframe_files /path/to/proteome/data/ \
  --outdir results/population_analysis/ \
  --gene_list /path/to/gene_list.txt \
  --tag proteome_study \
  --reg_formula "cut_C_Rall ~ AA + region"
```

#### Monte Carlo Analysis (`run_montecarlo.py`)
Identify candidate proteins with poor per-protein statistics using Monte Carlo methods:

```bash
python scripts/run_montecarlo.py \
  --dataframe_files /path/to/proteome/data/ \
  --outpath results/monte_carlo/ \
  --gene_list /path/to/gene_list.txt \
  --tag mc_analysis \
  --steps 100000 \
  --n_groups 4 \
  --C1 1.0 \
  --C2 2.5 \
  --beta 0.05
```

### 7. Folding Pathway Analysis

#### Pathway Statistics (`run_Foldingpathway.py`)
Analyze post-transitional folding pathways from temperature quenching simulations:

```bash
python scripts/run_Foldingpathway.py \
  --msm_data_file /path/to/msm_meta_set.csv \
  --meta_set_file /path/to/meta_set.csv \
  --traj_type_col traj_type_column \
  --outdir results/folding_pathways/ \
  --rm_traj_list 65 75 155
```

## Package Structure

```
EntDetect/
├── EntDetect/                    # Main package
│   ├── __init__.py
│   ├── gaussian_entanglement.py  # Core entanglement calculations
│   ├── clustering.py             # Entanglement clustering methods
│   ├── order_params.py          # Order parameter calculations
│   ├── compare_sim2exp.py       # Simulation-experiment comparison
│   ├── statistics.py            # Statistical analysis methods
│   ├── entanglement_features.py # Feature generation
│   ├── change_resolution.py     # Resolution conversion
│   └── utilities.py             # Helper functions
├── scripts/                     # Example workflow scripts
├── Documentation/               # Detailed module documentation
└── TestingGrounds/             # Test data and examples
```

## Core Modules

- **`gaussian_entanglement`**: Calculate Gaussian linking numbers and identify entanglements
- **`clustering`**: Cluster native and non-native entanglements, build MSMs
- **`order_params`**: Compute Q, G, K, SASA, Jwalk, and cross-linking propensity
- **`compare_sim2exp`**: Integrate LiP-MS and XL-MS experimental data
- **`statistics`**: Population modeling and Monte Carlo analysis
- **`entanglement_features`**: Generate structural features for entanglements
- **`change_resolution`**: Convert between all-atom and coarse-grained representations

## Documentation

Detailed documentation for each module is available:

- [Gaussian Entanglement](Documentation/gaussian_entanglement.md)
- [Clustering](Documentation/clustering.md)
- [Order Parameters](Documentation/order_params.md)
- [Simulation-Experiment Comparison](Documentation/compare_sim2exp.md)
- [Statistical Analysis](Documentation/statistics.md)
- [Entanglement Features](Documentation/entanglement_features.md)
- [Resolution Conversion](Documentation/change_resolution.md)
- [Utilities](Documentation/utilities.md)

## Requirements

- Python 3.8+
- NumPy
- Pandas
- SciPy
- MDAnalysis
- OpenMM (for force field operations)
- Matplotlib (for visualization)
- See `environment.yml` for complete dependencies

## Citation

If you use EntDetect in your research, please cite:

```bibtex
@software{entdetect2024,
  title={EntDetect: A Python Package for Protein Entanglement Analysis},
  author={Your Name},
  year={2024},
  url={https://github.com/SparkyDaBear/EntDetect}
}
```

## Contributing

Contributions are welcome! Please see our contributing guidelines and submit pull requests for any improvements.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

For questions and support, please open an issue on GitHub or contact the developers. 
  


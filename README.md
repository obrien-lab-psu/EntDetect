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

Create a new conda environment and install EntDetect (from this repo checkout):

```bash
conda env create -f environment.yml
conda activate entdetect
```

Notes:
- The provided conda environment targets Python 3.11 for best compatibility with the scientific stack.
- Run `conda env create` from the EntDetect repo root (the environment file uses `pip -e .`).
- If you prefer installing into an existing env, use `pip install -e .` from the repo root.

macOS (Miniconda) notes:
- The default `environment.yml` is intended to work on both Linux and macOS via conda-forge.
- If you hit solver/build issues on macOS, try:

```bash
conda env create -f environment-mac.yml
conda activate entdetect
```

- On Apple Silicon, if a dependency is missing for `osx-arm64`, a common workaround is:

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment-mac.yml
conda activate entdetect
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

## Tutorial

For step-by-step, runnable examples (including both all-atom vs C-alpha structure examples, and a bundled CG trajectory example), see:

- [Documentation/tutorial_examples.md](Documentation/tutorial_examples.md)

For all scripts, `--help` is the most up-to-date reference for required inputs:

```bash
python scripts/run_nativeNCLE.py --help
python scripts/run_OP_on_simulation_traj.py --help
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

- [Tutorial (examples with included assets)](Documentation/tutorial_examples.md)
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
  


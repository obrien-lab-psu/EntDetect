# Gaussian Entanglement Module

The `gaussian_entanglement` module provides the core functionality for calculating and analyzing non-covalent lasso entanglements in protein structures using Gaussian entanglement theory.

## Overview

This module implements methods to:
- Calculate Gaussian linking numbers between protein loops
- Identify native entanglements in protein structures  
- Filter high-quality entanglements based on structural criteria
- Support both all-atom and coarse-grained representations

## Classes

### GaussianEntanglement

Main class for entanglement calculations and analysis.

#### Initialization
```python
from EntDetect.gaussian_entanglement import GaussianEntanglement

ge = GaussianEntanglement(g_threshold=0.6, density=0.0, Calpha=False, CG=False)
```

**Parameters:**
- `g_threshold` (float): Threshold for Gaussian entanglement score (default: 0.6)
- `density` (float): Density parameter for contact filtering (default: 0.0)  
- `Calpha` (bool): Use Cα atoms for native contacts vs 4.5Å heavy atom cutoff (default: False)
- `CG` (bool): Enable coarse-grained mode (default: False)

#### Key Methods

##### calculate_native_entanglements()
Identifies all native entanglements in a protein structure.

```python
result = ge.calculate_native_entanglements(
    pdb_file="/path/to/structure.pdb",
    outdir="results/native_entanglements/", 
    ID="protein_name"
)
```

**Parameters:**
- `pdb_file` (str): Path to cleaned PDB structure file
- `outdir` (str): Output directory for results
- `ID` (str): Identifier for the analysis

**Returns:**
- Dictionary with output file paths and entanglement data

##### select_high_quality_entanglements()
Filters entanglements based on structural quality criteria.

```python
hq_result = ge.select_high_quality_entanglements(
    entanglement_file="native_entanglements.csv",
    pdb_file="/path/to/structure.pdb",
    outdir="results/hq_entanglements/",
    ID="protein_name", 
    model="EXP"
)
```

**Parameters:**
- `entanglement_file` (str): Path to entanglement data file
- `pdb_file` (str): Path to structure file
- `outdir` (str): Output directory
- `ID` (str): Analysis identifier  
- `model` (str): Model type ("EXP" for experimental, "AF" for AlphaFold)

## Usage Examples

### Basic Native Entanglement Analysis

```python
from EntDetect.gaussian_entanglement import GaussianEntanglement

# Initialize for experimental structure analysis
ge = GaussianEntanglement(g_threshold=0.6, Calpha=False, CG=False)

# Calculate native entanglements
native_ents = ge.calculate_native_entanglements(
    pdb_file="structure.pdb",
    outdir="results/native/", 
    ID="1ZMR"
)

# Filter for high quality entanglements
hq_ents = ge.select_high_quality_entanglements(
    entanglement_file=native_ents['outfile'],
    pdb_file="structure.pdb",
    outdir="results/hq/",
    ID="1ZMR",
    model="EXP"
)
```

### Coarse-Grained Analysis

```python
# Initialize for coarse-grained trajectory analysis
ge_cg = GaussianEntanglement(g_threshold=0.6, Calpha=True, CG=True)

# Process coarse-grained structure
cg_result = ge_cg.calculate_native_entanglements(
    pdb_file="cg_structure.pdb",
    outdir="results/cg_native/",
    ID="1ZMR_CG"
)
```

## Integration with Scripts

This module is used extensively in the standardized workflow scripts:

### run_nativeNCLE.py
Complete pipeline for native entanglement analysis:

```bash
python scripts/run_nativeNCLE.py \
  --struct structure.pdb \
  --outdir results/native_analysis/ \
  --ID protein_name \
  --organism Ecoli
```

### run_OP_on_simulation_traj.py  
Uses entanglement calculations for G parameter computation:

```bash
python scripts/run_OP_on_simulation_traj.py \
  --Traj 1 \
  --PSF structure.psf \
  --DCD trajectory.dcd \
  --ID protein_name \
  --COR structure.cor \
  --sec_elements secondary.txt \
  --domain domain_def.dat \
  --outdir results/OP/ \
  --start 0
```

## Output Files

### Native Entanglement Files
- `*_native_entanglements.csv`: Raw entanglement data with GLN values
- `*_hq_entanglements.csv`: High-quality filtered entanglements  
- `*_clustered_entanglements.csv`: Clustered entanglements (via clustering module)

### File Format
CSV files contain columns:
- `loop1_start`, `loop1_end`: First loop boundaries
- `loop2_start`, `loop2_end`: Second loop boundaries  
- `GLN`: Gaussian linking number
- `entanglement_type`: Classification of entanglement
- `quality_score`: Structural quality metrics

## Theory Background

The Gaussian entanglement method calculates linking numbers between protein loops using:

1. **Loop Identification**: Protein backbone segments between secondary structures
2. **GLN Calculation**: Gaussian linking number between loop pairs  
3. **Threshold Application**: Entanglements identified above g_threshold
4. **Quality Filtering**: Structural validation based on contact patterns

## Performance Notes

- **Memory Usage**: Scales quadratically with protein size
- **Computation Time**: ~1-10 minutes for typical protein structures
- **Parallelization**: Thread-safe for multiple structure analysis
- **File I/O**: Optimized CSV output for large datasets

## Dependencies

- NumPy: Numerical computations
- Pandas: Data manipulation  
- MDAnalysis: Structure file parsing
- SciPy: Statistical functions
select_high_quality_entanglements(entanglement_file, pdb_file, outdir, ID, model, mapping)
```
Filters and selects high-quality entanglements based on additional criteria.
- `entanglement_file` (str): Path to the file containing raw entanglement data.
- `pdb_file` (str): Path to the cleaned PDB file.
- `outdir` (str): Output directory for filtered results.
- `ID` (str): Identifier for the current analysis.
- `model` (str): Model type (e.g., 'EXP' or 'AF').
- `mapping` (str): Mapping file or 'None'.

**Returns:**
A dictionary containing output file paths and filtered entanglement data.

---

## Usage Example
```python
from EntDetect.gaussian_entanglement import GaussianEntanglement

ge = GaussianEntanglement(g_threshold=0.6, density=0.0, Calpha=False, CG=False)
native_ent = ge.calculate_native_entanglements("protein_clean.pdb", outdir="results/Native_GE", ID="protein1")
hq_ent = ge.select_high_quality_entanglements(native_ent['outfile'], "protein_clean.pdb", outdir="results/Native_HQ_GE", ID="protein1", model="AF", mapping="None")
```

**Notes:**
- The class is designed to be used in protein structure analysis pipelines, especially for large-scale or automated workflows.
- Output directories will be created if they do not exist.
- For best results, input PDB files should be pre-processed and cleaned.

For further details, refer to the source code in `EntDetect/gaussian_entanglement.py` or the package documentation.

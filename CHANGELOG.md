# Changelog

All notable changes to EntDetect will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-01-26

### Added
- **Configurable ENT Detection Methods** - Three distinct methods for determining entanglement status:
  - Method 1: GLN-based (any nonzero Gaussian Linking Number)
  - Method 2: TLN-based (any nonzero Topological Linking Number) - default
  - Method 3: Consensus (both GLN and TLN must agree on same terminus)
  - New `--ent_detection_method` command-line argument in `run_nativeNCLE.py`
  - New `determine_ent_status()` method in `GaussianEntanglement` class

- **Memory Optimization for AlphaFold Models**
  - Auto-enable C-alpha mode for AF structures to reduce memory usage
  - Prevents out-of-memory errors on large structures

- **Version Tracking**
  - Added `__version__` to `EntDetect/__init__.py`
  - Added CHANGELOG.md for release documentation

### Fixed
- **Crossing Residue Validation Bug** - Critical fix in `find_crossing()` method:
  - Original code used alternating pair comparison which missed consecutive crossings
  - Replaced with greedy algorithm ensuring all retained crossings are ≥10 residues apart
  - Applied fix to both N-terminal and C-terminal crossing validation

### Removed
- **Placeholder '?' Logic** - Removed no-longer-needed sentinel value system:
  - Removed `mark_absent_crossings()` method from `GaussianEntanglement`
  - Removed `replace('?', '-100000')` from clustering pipeline
  - Crossing strings now only contain actual identified residues or remain empty

### Changed
- Default `ent_detection_method` in `run_nativeNCLE.py` changed to 3 (consensus method)
- Crossing validation now uses greedy filtering instead of alternating pair comparison

### Technical Details
- **Files Modified**:
  - `EntDetect/gaussian_entanglement.py` - Added ENT detection methods, fixed crossing validation
  - `EntDetect/clustering.py` - Removed placeholder logic
  - `scripts/run_nativeNCLE.py` - Added CLI argument, AF model auto-optimization
  - `EntDetect/__init__.py` - Added version and package metadata
  - `setup.py` - Updated version to 0.2.0

## [0.1.0] - Initial Release

### Features
- Native Gaussian entanglement detection in protein structures
- Support for both experimental (PDB) and AlphaFold (AF) structures
- Coarse-grained and all-atom representation options
- Entanglement clustering to remove degeneracies
- High-quality entanglement filtering
- Entanglement feature generation
- Disulfide bond detection
- Support for C-alpha (CA) and coarse-grained (CG) models

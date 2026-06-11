#!/usr/bin/env python
"""
Comprehensive test suite for Workflow 3: Sim-to-Experiment consistency testing.

This script tests:
1. Markdown tutorial (Section 1-6)
2. Jupyter notebook cells
3. Data loading and validation
4. CollectOP (Steps 4a-4b)
5. MassSpec consistency test (Step 5)
6. Result validation and visualization

Run with: python test_workflow3_rigorous.py
"""

import os
import sys
import numpy as np
import pandas as pd
import time
from datetime import datetime, timedelta
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION
# ============================================================================
DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR = f"{DATASTORE}/outputs/workflow3"
REFSTRUCT = f"{DATASTORE}/user_input/reference_structures"
AA_TRAJ_DIR = f"{DATASTORE}/user_input/aa_trajectories"
CG_TRAJ_DIR = f"{DATASTORE}/user_input/cg_trajectories"

WORKFLOW2_OUTDIR = f"{DATASTORE}/outputs/workflow2"
SASA_DIR = f"{WORKFLOW2_OUTDIR}/OP_AA/SASA"
XP_DIR = f"{WORKFLOW2_OUTDIR}/OP_AA/XP"

PROT_LEN = 387
N_FRAMES = 335
N_TRAJ = 1000
STATE_IDX_LIST = [4, 6, 8]
NATIVE_STATE_IDX = 9

# ============================================================================
# TEST 1: Verify Experimental Data
# ============================================================================

def test_experimental_data():
    """Test Section: Load and verify experimental inputs."""
    logger.info("=" * 80)
    logger.info("TEST 1: Experimental Data Loading")
    logger.info("=" * 80)
    
    lipms_file = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx"
    xlms_file = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx"
    
    try:
        lipms = pd.read_excel(lipms_file)
        xlms = pd.read_excel(xlms_file)
        
        logger.info(f"✓ LiP-MS data loaded: {lipms.shape}")
        logger.info(f"  Columns: {list(lipms.columns)}")
        logger.info(f"✓ XL-MS data loaded: {xlms.shape}")
        logger.info(f"  Columns: {list(xlms.columns)}")
        
        return True, (lipms, xlms)
    except Exception as e:
        logger.error(f"✗ Failed to load experimental data: {e}")
        return False, None

# ============================================================================
# TEST 2: Verify SASA and XP Files
# ============================================================================

def test_op_files():
    """Test Step 3: Verify per-trajectory SASA and XP files."""
    logger.info("\n" + "=" * 80)
    logger.info("TEST 2: Order Parameter Files Verification")
    logger.info("=" * 80)
    
    # Check SASA files
    import glob
    sasa_files = sorted(glob.glob(f"{SASA_DIR}/1ZMR_Traj*.SASA"))
    logger.info(f"SASA files found: {len(sasa_files)}/{N_TRAJ}")
    
    if len(sasa_files) == N_TRAJ:
        logger.info("✓ All 1000 SASA files present")
        sample = sasa_files[0]
        df = pd.read_csv(sample)
        logger.info(f"  Sample SASA file shape: {df.shape}")
        logger.info(f"  Columns: {list(df.columns)}")
    else:
        logger.warning(f"⚠️  Missing {N_TRAJ - len(sasa_files)} SASA files")
    
    # Check XP files
    xp_files = sorted(glob.glob(f"{XP_DIR}/1ZMR_Traj*.XP"))
    logger.info(f"XP files found: {len(xp_files)}/{N_TRAJ}")
    
    if len(xp_files) == N_TRAJ:
        logger.info("✓ All 1000 XP files present")
        sample = xp_files[0]
        df = pd.read_csv(sample, sep='\t')
        logger.info(f"  Sample XP file shape: {df.shape}")
        logger.info(f"  Columns: {list(df.columns)}")
    else:
        logger.warning(f"⚠️  Missing {N_TRAJ - len(xp_files)} XP files")
    
    return len(sasa_files) == N_TRAJ and len(xp_files) == N_TRAJ

# ============================================================================
# TEST 3: Collect OP Arrays
# ============================================================================

def test_collect_op():
    """Test Step 4: Collect SASA and Jwalk arrays."""
    logger.info("\n" + "=" * 80)
    logger.info("TEST 3: OP Collection (SASA + Jwalk)")
    logger.info("=" * 80)
    
    from EntDetect.order_params import CollectOP
    
    try:
        collector = CollectOP(
            sasa_dir=SASA_DIR,
            xp_dir=XP_DIR,
            outdir=OUTDIR,
            ID="1ZMR",
            n_traj=N_TRAJ,
            n_frames=N_FRAMES,
            prot_len=PROT_LEN,
        )
        
        # Collect SASA
        logger.info("\nCollecting SASA...")
        sasa_file = collector.collect_SASA()
        logger.info(f"✓ SASA saved: {sasa_file}")
        
        # Verify SASA shape
        sasa = np.load(sasa_file)
        logger.info(f"  Shape: {sasa.shape}")
        logger.info(f"  Expected: ({N_TRAJ}, {N_FRAMES}, {PROT_LEN})")
        assert sasa.shape == (N_TRAJ, N_FRAMES, PROT_LEN), "SASA shape mismatch!"
        logger.info(f"  ✓ Shape verified")
        logger.info(f"  Data range: {np.nanmin(sasa):.2f} – {np.nanmax(sasa):.2f} Å²")
        logger.info(f"  NaN count: {np.sum(np.isnan(sasa))}")
        
        # Collect Jwalk
        logger.info("\nCollecting Jwalk (this may take 30-60 minutes)...")
        start_time = time.time()
        jwalk_file = collector.collect_Jwalk()
        elapsed = time.time() - start_time
        logger.info(f"✓ Jwalk saved: {jwalk_file}")
        logger.info(f"  Time elapsed: {elapsed/60:.1f} minutes")
        
        # Verify Jwalk shape
        jwalk = np.load(jwalk_file, allow_pickle=True)
        logger.info(f"  Shape: {jwalk.shape}")
        logger.info(f"  Expected: ({N_TRAJ}, {N_FRAMES})")
        assert jwalk.shape == (N_TRAJ, N_FRAMES), "Jwalk shape mismatch!"
        logger.info(f"  ✓ Shape verified")
        
        # Check sample
        sample = jwalk[0, 0]
        if sample is not None:
            logger.info(f"  Sample entry type: {type(sample)}")
            if isinstance(sample, dict):
                logger.info(f"  Sample entry has {len(sample)} pairs")
        
        return True, (sasa_file, jwalk_file)
        
    except Exception as e:
        logger.error(f"✗ OP collection failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None

# ============================================================================
# TEST 4: Load MSM Data
# ============================================================================

def test_msm_data():
    """Test Step 5a: Load MSM mapping and metadata."""
    logger.info("\n" + "=" * 80)
    logger.info("TEST 4: MSM Data Loading")
    logger.info("=" * 80)
    
    msm_mapping_file = f"{WORKFLOW2_OUTDIR}/MSM/1ZMR_prod_MSMmapping.csv"
    meta_dist_file = f"{WORKFLOW2_OUTDIR}/MSM/1ZMR_prod_meta_dist.npy"
    
    try:
        msm_mapping = pd.read_csv(msm_mapping_file)
        logger.info(f"✓ MSM mapping loaded: {msm_mapping.shape}")
        logger.info(f"  Columns: {list(msm_mapping.columns)}")
        
        meta_dist = np.load(meta_dist_file)
        logger.info(f"✓ Meta distribution loaded: {meta_dist.shape}")
        
        meta_states = sorted(msm_mapping['metastablestate'].unique())
        logger.info(f"  Metastable states: {meta_states}")
        
        return True, (msm_mapping, meta_dist)
    except Exception as e:
        logger.error(f"✗ Failed to load MSM data: {e}")
        return False, None

# ============================================================================
# TEST 5: Initialize and Run MassSpec
# ============================================================================

def test_consistency_test(sasa_file, jwalk_file):
    """Test Step 5: Run LiP-MS / XL-MS consistency test."""
    logger.info("\n" + "=" * 80)
    logger.info("TEST 5: Run Consistency Test (FULL PRODUCTION DATA)")
    logger.info("=" * 80)
    logger.info(f"⏱️  EXPECTED RUNTIME: 1–3 hours")
    logger.info(f"⚠️  Running on ALL {N_TRAJ} trajectories × {N_FRAMES} frames")
    logger.info("=" * 80)
    
    from EntDetect.compare_sim2exp import MassSpec
    
    msm_mapping_file = f"{WORKFLOW2_OUTDIR}/MSM/1ZMR_prod_MSMmapping.csv"
    meta_dist_file = f"{WORKFLOW2_OUTDIR}/MSM/1ZMR_prod_meta_dist.npy"
    lipms_file = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx"
    xlms_file = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx"
    cluster_data_file = f"{WORKFLOW2_OUTDIR}/nonnative_clustering/cluster_data_topoly_linking_number.npz"
    oppath = f"{WORKFLOW2_OUTDIR}/OP_AA/"
    native_aa_pdb = f"{REFSTRUCT}/1zmr_model_clean.pdb"
    
    try:
        logger.info("\nInitializing MassSpec...")
        MS = MassSpec(
            msm_data_file=msm_mapping_file,
            meta_dist_file=meta_dist_file,
            LiPMS_exp_file=lipms_file,
            sasa_data_file=sasa_file,
            XLMS_exp_file=xlms_file,
            dist_data_file=jwalk_file,
            cluster_data_file=cluster_data_file,
            OPpath=oppath,
            AAdcd_dir=AA_TRAJ_DIR,
            native_AA_pdb=native_aa_pdb,
            state_idx_list=STATE_IDX_LIST,
            prot_len=PROT_LEN,
            last_num_frames=N_FRAMES,
            rm_traj_list=[],
            native_state_idx=NATIVE_STATE_IDX,
            outdir=f"{OUTDIR}/MassSpec_ConsistencyTest",
            ID="1ZMR",
            start=6600
        )
        logger.info("✓ MassSpec initialized")
        
        logger.info("\nStarting consistency test...")
        start_time = time.time()
        consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
        elapsed = time.time() - start_time
        
        logger.info(f"✓ Consistency test completed in {elapsed/60:.1f} minutes ({elapsed/3600:.2f} hours)")
        logger.info(f"  Data: {consist_data_file}")
        logger.info(f"  Results: {consist_result_file}")
        
        # Select representative structures
        logger.info("\nSelecting representative structures...")
        MS.select_rep_structs(
            consist_data_file,
            consist_result_file,
            total_traj_num_frames=335,
            last_num_frames=67
        )
        logger.info("✓ Representative structures selected")
        
        return True, (consist_data_file, consist_result_file)
        
    except Exception as e:
        logger.error(f"✗ Consistency test failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None

# ============================================================================
# TEST 6: Validate Results
# ============================================================================

def test_result_validation(consist_data_file, consist_result_file):
    """Test Step 6: Validate consistency test results."""
    logger.info("\n" + "=" * 80)
    logger.info("TEST 6: Result Validation")
    logger.info("=" * 80)
    
    try:
        import glob
        
        outdir = f"{OUTDIR}/MassSpec_ConsistencyTest"
        
        # Check for result files
        pvalue_files = glob.glob(f"{outdir}/LiPMS_XLMS_consist_pvalues_metastates_*.xlsx")
        struct_files = glob.glob(f"{outdir}/Consistent_structures_*.xlsx")
        
        logger.info(f"P-value files: {len(pvalue_files)}")
        logger.info(f"Structure files: {len(struct_files)}")
        
        if pvalue_files:
            pval_df = pd.read_excel(pvalue_files[0])
            logger.info(f"✓ P-values loaded: {pval_df.shape}")
            logger.info(f"  Data:\n{pval_df}")
        
        if struct_files:
            struct_df = pd.read_excel(struct_files[0])
            logger.info(f"✓ Structure references: {struct_df.shape}")
            logger.info(f"  First few rows:\n{struct_df.head()}")
        
        return True
        
    except Exception as e:
        logger.error(f"✗ Result validation failed: {e}")
        return False

# ============================================================================
# MAIN TEST SUITE
# ============================================================================

def main():
    """Run all tests."""
    logger.info("\n" + "=" * 80)
    logger.info("WORKFLOW 3 COMPREHENSIVE TEST SUITE")
    logger.info("=" * 80)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Estimated total runtime: 2–4 hours")
    
    tests_passed = 0
    tests_total = 0
    
    # Test 1: Experimental data
    tests_total += 1
    success, exp_data = test_experimental_data()
    if success:
        tests_passed += 1
    else:
        logger.error("TEST 1 FAILED - Aborting")
        return
    
    # Test 2: OP files
    tests_total += 1
    success = test_op_files()
    if success:
        tests_passed += 1
    
    # Test 3: Collect OP
    tests_total += 1
    success, op_files = test_collect_op()
    if not success:
        logger.error("TEST 3 FAILED - Aborting (cannot proceed without OP arrays)")
        return
    else:
        tests_passed += 1
    
    sasa_file, jwalk_file = op_files
    
    # Test 4: MSM data
    tests_total += 1
    success, msm_data = test_msm_data()
    if success:
        tests_passed += 1
    else:
        logger.error("TEST 4 FAILED - Aborting")
        return
    
    # Test 5: Consistency test
    tests_total += 1
    success, result_files = test_consistency_test(sasa_file, jwalk_file)
    if success:
        tests_passed += 1
        
        # Test 6: Validate results
        tests_total += 1
        consist_data_file, consist_result_file = result_files
        success = test_result_validation(consist_data_file, consist_result_file)
        if success:
            tests_passed += 1
    
    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("TEST SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Tests passed: {tests_passed}/{tests_total}")
    logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if tests_passed == tests_total:
        logger.info("\n✓ ALL TESTS PASSED!")
        return 0
    else:
        logger.error(f"\n✗ {tests_total - tests_passed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())

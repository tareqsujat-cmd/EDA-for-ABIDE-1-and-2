"""
Time Series Exploratory Data Analysis (EDA) for ABIDE 1 and ABIDE 2
Analyzes resting-state fMRI time series data, signal quality, and functional connectivity.
"""

import os
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import nibabel as nib
from pathlib import Path
from tqdm import tqdm
from scipy.stats import pearsonr

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIG
# ============================================================================
SCRIPT_DIR = Path(__file__).parent

CONFIG = {
    'abide1_ts_dir': SCRIPT_DIR / 'ABIDE_1/preprocessed/rois_cc200',
    'abide2_ts_dir': SCRIPT_DIR / 'ABIDE_2/preprocessed/dcan_gordon2014',
    'output_dir': SCRIPT_DIR / 'plots/timeseries',
    'figsize': (12, 6),
    'dpi': 150,
    'style': 'whitegrid',
    'context': 'talk',
    'n_sample_subjects': 10,
    'n_sample_rois': 5,
}

# ============================================================================
# SETUP
# ============================================================================
sns.set_theme(style=CONFIG['style'], context=CONFIG['context'])
CONFIG['output_dir'].mkdir(parents=True, exist_ok=True)


# ============================================================================
# LOAD TIME SERIES DATA
# ============================================================================
def load_abide1_timeseries():
    """Load ABIDE1 CC200 time series (.1D files)."""
    print("Loading ABIDE1 time series (CC200)...")
    
    ts_data = {}
    n_rois = None
    
    if not CONFIG['abide1_ts_dir'].exists():
        print(f"  WARNING: ABIDE1 time series directory not found: {CONFIG['abide1_ts_dir']}")
        return ts_data
    
    files = sorted(CONFIG['abide1_ts_dir'].glob('*.1D'))
    print(f"  Found {len(files)} .1D files")
    
    for filepath in tqdm(files, desc="Loading ABIDE1 .1D files"):
        try:
            # Load with np.loadtxt, skipping comments
            data = np.loadtxt(filepath, comments='#')
            
            # Handle 1D arrays (single ROI)
            if data.ndim == 1:
                data = data.reshape(-1, 1)
            
            subject_id = filepath.stem
            ts_data[subject_id] = {
                'data': data,
                'n_tr': data.shape[0],
                'n_roi': data.shape[1],
                'dataset': 'ABIDE1',
                'filepath': filepath,
            }
            
            if n_rois is None:
                n_rois = data.shape[1]
        except Exception as e:
            print(f"    Error loading {filepath}: {e}")
    
    print(f"  Loaded {len(ts_data)} subjects. Expected N_ROI = {n_rois}")
    return ts_data


def load_abide2_timeseries():
    """Load ABIDE2 Gordon2014 time series (.ptseries.nii files)."""
    print("Loading ABIDE2 time series (Gordon2014)...")
    
    ts_data = {}
    n_rois = None
    
    if not CONFIG['abide2_ts_dir'].exists():
        print(f"  WARNING: ABIDE2 time series directory not found: {CONFIG['abide2_ts_dir']}")
        return ts_data
    
    files = sorted(CONFIG['abide2_ts_dir'].glob('*.ptseries.nii'))
    print(f"  Found {len(files)} .ptseries.nii files")
    
    for filepath in tqdm(files, desc="Loading ABIDE2 .ptseries.nii files"):
        try:
            # Load with nibabel
            img = nib.load(filepath)
            data = np.asarray(img.dataobj)
            
            # Time series should be (T, n_roi) or have time in last dim
            if data.ndim > 1:
                # Usually (n_roi, T) or (T, n_roi), try to infer
                if data.shape[0] > data.shape[1]:
                    # Likely (T, n_roi)
                    pass
                else:
                    # Likely (n_roi, T), transpose
                    data = data.T
            else:
                data = data.reshape(-1, 1)
            
            subject_id = filepath.stem.replace('.ptseries.nii', '').replace('.ptseries', '')
            ts_data[subject_id] = {
                'data': data,
                'n_tr': data.shape[0],
                'n_roi': data.shape[1] if data.ndim > 1 else 1,
                'dataset': 'ABIDE2',
                'filepath': filepath,
            }
            
            if n_rois is None:
                n_rois = data.shape[1] if data.ndim > 1 else 1
        except Exception as e:
            print(f"    Error loading {filepath}: {e}")
    
    print(f"  Loaded {len(ts_data)} subjects. Expected N_ROI = {n_rois}")
    return ts_data


# ============================================================================
# TASK 1: TIME SERIES PROPERTIES SUMMARY
# ============================================================================
def create_timeseries_summary(ts_data_1, ts_data_2):
    """Create summary of time series properties."""
    print("\n1. Creating time series properties summary...")
    
    summary_data = []
    
    for dataset_name, ts_data in [('ABIDE1', ts_data_1), ('ABIDE2', ts_data_2)]:
        if not ts_data:
            continue
        
        n_files = len(ts_data)
        n_trs = [v['n_tr'] for v in ts_data.values()]
        n_rois = [v['n_roi'] for v in ts_data.values()]
        
        summary_data.append({
            'Dataset': dataset_name,
            'N Files': n_files,
            'Mean TRs': np.mean(n_trs),
            'Min TRs': np.min(n_trs),
            'Max TRs': np.max(n_trs),
            'N ROIs': np.mean(n_rois),
        })
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    print()


# ============================================================================
# TASK 2: SIGNAL QUALITY CHECKS
# ============================================================================
def check_signal_quality(ts_data_1, ts_data_2):
    """Check signal quality for sample subjects."""
    print("\n2. Checking signal quality (sample of 10 subjects per dataset)...")
    
    for dataset_name, ts_data in [('ABIDE1', ts_data_1), ('ABIDE2', ts_data_2)]:
        if not ts_data:
            print(f"  No {dataset_name} data available")
            continue
        
        # Sample subjects
        sample_subjects = list(np.random.choice(list(ts_data.keys()), 
                                         min(CONFIG['n_sample_subjects'], len(ts_data)), 
                                         replace=False))
        
        dead_rois = []
        
        for subject_id in sample_subjects:
            data = ts_data[subject_id]['data']
            
            # Check for dead ROIs (all zeros or very low variance)
            for roi_idx in range(data.shape[1]):
                roi_ts = data[:, roi_idx]
                if np.allclose(roi_ts, 0) or np.std(roi_ts) < 1e-6:
                    dead_rois.append((subject_id, roi_idx))
        
        if dead_rois:
            print(f"  {dataset_name} - Found {len(dead_rois)} dead ROIs:")
            for subj, roi in dead_rois[:5]:
                print(f"    {subj}, ROI {roi}")
            if len(dead_rois) > 5:
                print(f"    ... and {len(dead_rois) - 5} more")
        else:
            print(f"  {dataset_name} - No dead ROIs detected in sample")
        
        # Plot raw time series for 5 random ROIs of 1 subject
        if sample_subjects:
            subject_id = sample_subjects[0]
            data = ts_data[subject_id]['data']
            n_rois_to_plot = min(CONFIG['n_sample_rois'], data.shape[1])
            roi_indices = np.random.choice(data.shape[1], n_rois_to_plot, replace=False)
            
            fig, axes = plt.subplots(n_rois_to_plot, 1, figsize=(14, 2*n_rois_to_plot), dpi=CONFIG['dpi'])
            if n_rois_to_plot == 1:
                axes = [axes]
            
            for idx, roi in enumerate(roi_indices):
                ts = data[:, roi]
                axes[idx].plot(ts, linewidth=0.5)
                axes[idx].set_ylabel(f'ROI {roi}')
                axes[idx].set_xlim([0, len(ts)])
                axes[idx].grid(True, alpha=0.3)
            
            axes[-1].set_xlabel('TR')
            fig.suptitle(f'{dataset_name} - Raw Time Series (Subject: {subject_id})', fontsize=14)
            plt.tight_layout()
            plt.savefig(CONFIG['output_dir'] / f'raw_timeseries_{dataset_name}_{subject_id[:20]}.png', dpi=CONFIG['dpi'])
            plt.close()
            
            # Plot mean signal ± SD band
            fig, ax = plt.subplots(figsize=CONFIG['figsize'], dpi=CONFIG['dpi'])
            mean_ts = np.mean(data, axis=1)
            std_ts = np.std(data, axis=1)
            tr_indices = np.arange(len(mean_ts))
            
            ax.plot(tr_indices, mean_ts, 'b-', linewidth=2, label='Mean')
            ax.fill_between(tr_indices, mean_ts - std_ts, mean_ts + std_ts, alpha=0.3, label='±1 SD')
            ax.set_xlabel('TR')
            ax.set_ylabel('Signal')
            ax.set_title(f'{dataset_name} - Mean Signal ± SD (Subject: {subject_id})')
            ax.legend()
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(CONFIG['output_dir'] / f'mean_signal_{dataset_name}_{subject_id[:20]}.png', dpi=CONFIG['dpi'])
            plt.close()


# ============================================================================
# TASK 3: FUNCTIONAL CONNECTIVITY MATRICES
# ============================================================================
def compute_functional_connectivity(ts_data, phenotypic_file=None):
    """Compute functional connectivity matrices and plot."""
    
    if not ts_data:
        print("  No time series data available")
        return
    
    dataset_name = list(ts_data.values())[0]['dataset']
    print(f"  Computing FC matrices for {dataset_name}...")
    
    # Load phenotypic data to split by diagnosis if available
    pheno = None
    if phenotypic_file and Path(phenotypic_file).exists():
        try:
            pheno = pd.read_csv(phenotypic_file)
        except:
            pheno = None
    
    # Compute average FC across all subjects
    fc_matrices = []
    fc_asd = []
    fc_tc = []
    
    for subject_id, ts_info in tqdm(ts_data.items(), desc=f"Computing FC for {dataset_name}"):
        data = ts_info['data']
        
        # Handle NaNs and standardize
        data = np.nan_to_num(data)
        if data.shape[0] < 2:
            continue
        
        # Standardize each ROI
        data_z = (data - np.mean(data, axis=0)) / (np.std(data, axis=0) + 1e-8)
        
        # Compute correlation matrix
        fc = np.corrcoef(data_z.T)
        fc = np.nan_to_num(fc)
        fc_matrices.append(fc)
        
        # Try to get diagnosis from phenotypic data
        if pheno is not None:
            # Try to match subject ID
            try:
                matching_rows = pheno[pheno['SUB_ID'].astype(str).str.contains(subject_id.split('_')[-1])]
                if not matching_rows.empty:
                    dx = matching_rows.iloc[0].get('DX_GROUP', None)
                    if dx == 1:
                        fc_asd.append(fc)
                    elif dx == 2:
                        fc_tc.append(fc)
            except:
                pass
    
    if not fc_matrices:
        print(f"  No valid FC matrices computed for {dataset_name}")
        return
    
    # Handle variable FC matrix sizes (filter to most common size)
    if len(fc_matrices) > 0:
        sizes = [fc.shape[0] for fc in fc_matrices]
        from collections import Counter
        most_common_size = Counter(sizes).most_common(1)[0][0]
        
        # Filter to most common size
        fc_matrices_filtered = [fc for fc in fc_matrices if fc.shape[0] == most_common_size]
        fc_asd_filtered = [fc for fc in fc_asd if fc.shape[0] == most_common_size]
        fc_tc_filtered = [fc for fc in fc_tc if fc.shape[0] == most_common_size]
        
        if len(fc_matrices_filtered) < len(fc_matrices):
            print(f"    Filtered FC matrices to most common size ({most_common_size}×{most_common_size}): {len(fc_matrices_filtered)}/{len(fc_matrices)} subjects")
        
        if not fc_matrices_filtered:
            print(f"  No valid FC matrices after size filtering")
            return
        
        fc_matrices = fc_matrices_filtered
        fc_asd = fc_asd_filtered
        fc_tc = fc_tc_filtered
    
    # Average FC
    fc_avg = np.mean(fc_matrices, axis=0)
    
    # Plot average FC heatmap
    fig, ax = plt.subplots(figsize=(10, 9), dpi=CONFIG['dpi'])
    im = ax.imshow(fc_avg, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax.set_xlabel('ROI')
    ax.set_ylabel('ROI')
    ax.set_title(f'{dataset_name} - Average Functional Connectivity')
    plt.colorbar(im, ax=ax, label='Pearson Correlation')
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / f'fc_average_{dataset_name}.png', dpi=CONFIG['dpi'])
    plt.close()
    
    # Plot ASD vs TC difference if available
    if fc_asd and fc_tc:
        fc_asd_avg = np.mean(fc_asd, axis=0)
        fc_tc_avg = np.mean(fc_tc, axis=0)
        fc_diff = fc_asd_avg - fc_tc_avg
        
        fig, ax = plt.subplots(figsize=(10, 9), dpi=CONFIG['dpi'])
        vmax = np.max(np.abs(fc_diff))
        im = ax.imshow(fc_diff, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
        ax.set_xlabel('ROI')
        ax.set_ylabel('ROI')
        ax.set_title(f'{dataset_name} - FC Difference (ASD - TC)')
        plt.colorbar(im, ax=ax, label='Correlation Difference')
        plt.tight_layout()
        plt.savefig(CONFIG['output_dir'] / f'fc_diff_asd_tc_{dataset_name}.png', dpi=CONFIG['dpi'])
        plt.close()
        
        # Find top 10 ROI pairs with largest difference
        # Get upper triangle indices
        triu_indices = np.triu_indices(fc_diff.shape[0], k=1)
        diffs = fc_diff[triu_indices]
        top_indices = np.argsort(np.abs(diffs))[-10:]
        
        print(f"\n  Top 10 ROI pairs with largest FC difference (ASD - TC) in {dataset_name}:")
        for idx in reversed(top_indices):
            roi1, roi2 = triu_indices[0][idx], triu_indices[1][idx]
            diff_val = diffs[idx]
            print(f"    ROI {roi1} - ROI {roi2}: {diff_val:.4f}")
    
    print(f"  FC plots saved for {dataset_name}")


# ============================================================================
# TASK 4: ROI VARIANCE ANALYSIS
# ============================================================================
def analyze_roi_variance(ts_data):
    """Analyze variance of each ROI."""
    
    if not ts_data:
        print("  No time series data available")
        return
    
    dataset_name = list(ts_data.values())[0]['dataset']
    print(f"\n4. Analyzing ROI variance for {dataset_name}...")
    
    # Compute variance per ROI per subject
    roi_variances = []
    
    for subject_id, ts_info in ts_data.items():
        data = ts_info['data']
        data = np.nan_to_num(data)
        
        variances = np.var(data, axis=0)
        roi_variances.append(variances)
    
    if not roi_variances:
        print("  No variance data computed")
        return
    
    # Handle variable ROI counts (filter to most common size)
    if len(roi_variances) > 0:
        sizes = [len(var_array) for var_array in roi_variances]
        from collections import Counter
        most_common_size = Counter(sizes).most_common(1)[0][0]
        
        # Filter to most common size
        roi_variances_filtered = [var for var in roi_variances if len(var) == most_common_size]
        
        if len(roi_variances_filtered) < len(roi_variances):
            print(f"    Filtered ROI variances to most common size ({most_common_size} ROIs): {len(roi_variances_filtered)}/{len(roi_variances)} subjects")
        
        if not roi_variances_filtered:
            print(f"  No valid ROI variances after size filtering")
            return
        
        roi_variances = roi_variances_filtered
    
    # Average variance across subjects
    mean_roi_var = np.mean(roi_variances, axis=0)
    
    # Sort and identify top/bottom
    sorted_indices = np.argsort(mean_roi_var)[::-1]
    
    print(f"\n  Top 10 highest-variance ROIs in {dataset_name}:")
    for rank, idx in enumerate(sorted_indices[:10], 1):
        print(f"    {rank}. ROI {idx}: variance = {mean_roi_var[idx]:.6f}")
    
    print(f"\n  Top 10 lowest-variance ROIs in {dataset_name}:")
    for rank, idx in enumerate(sorted_indices[-10:], 1):
        print(f"    {rank}. ROI {idx}: variance = {mean_roi_var[idx]:.6f}")
    
    # Plot bar chart of ROI variance
    fig, ax = plt.subplots(figsize=(16, 6), dpi=CONFIG['dpi'])
    roi_indices = np.arange(len(mean_roi_var))
    sorted_order = np.argsort(mean_roi_var)[::-1]
    
    ax.bar(roi_indices, mean_roi_var[sorted_order], color='steelblue')
    ax.set_xlabel('ROI (sorted by variance)')
    ax.set_ylabel('Mean Variance')
    ax.set_title(f'{dataset_name} - Mean ROI Variance (sorted)')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / f'roi_variance_{dataset_name}.png', dpi=CONFIG['dpi'])
    plt.close()
    
    print(f"  Variance plot saved: roi_variance_{dataset_name}.png")


# ============================================================================
# TASK 5: MOTION SCRUBBING PREVIEW
# ============================================================================
def motion_scrubbing_preview(ts_data):
    """Preview motion scrubbing (note: FD not available in .1D/.ptseries.nii)."""
    
    if not ts_data:
        return
    
    dataset_name = list(ts_data.values())[0]['dataset']
    print(f"\n5. Motion scrubbing preview for {dataset_name}...")
    print("  Note: Motion parameters not available in .1D or .ptseries.nii files")
    print("  Motion scrubbing would typically use 6-parameter motion estimates from preprocessing.")
    print("  Consider loading motion parameters from FSF/FMRIPREP outputs separately.")


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("TIME SERIES EDA - ABIDE 1 and 2")
    print("=" * 70)
    
    # Load time series data
    ts_data_1 = load_abide1_timeseries()
    ts_data_2 = load_abide2_timeseries()
    
    # Task 1: Summary
    create_timeseries_summary(ts_data_1, ts_data_2)
    
    # Task 2: Signal quality
    check_signal_quality(ts_data_1, ts_data_2)
    
    # Task 3: Functional connectivity
    print("\n3. Computing functional connectivity matrices...")
    compute_functional_connectivity(ts_data_1, phenotypic_file='ABIDE_1/raw/Phenotypic_V1_0b_preprocessed1.csv')
    compute_functional_connectivity(ts_data_2, phenotypic_file='ABIDE_2/raw/ABIDEII_Composite_Phenotypic.csv')
    
    # Task 4: ROI variance
    analyze_roi_variance(ts_data_1)
    analyze_roi_variance(ts_data_2)
    
    # Task 5: Motion scrubbing preview
    motion_scrubbing_preview(ts_data_1)
    motion_scrubbing_preview(ts_data_2)
    
    print("\n" + "=" * 70)
    print("Time series EDA complete! All plots saved to:", CONFIG['output_dir'])
    print("=" * 70)


if __name__ == '__main__':
    main()

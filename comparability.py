"""
Cross-Dataset Comparability Analysis
Examines site harmonization, atlas differences, and diagnosis balance between ABIDE 1 and 2.
"""

import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIG
# ============================================================================
SCRIPT_DIR = Path(__file__).parent

CONFIG = {
    'abide1_pheno_path': SCRIPT_DIR / 'ABIDE_1/raw/Phenotypic_V1_0b_preprocessed1.csv',
    'abide2_pheno_path': SCRIPT_DIR / 'ABIDE_2/raw/ABIDEII_Composite_Phenotypic.csv',
    'output_dir': SCRIPT_DIR / 'plots/phenotypic',
    'figsize': (12, 6),
    'dpi': 150,
    'style': 'whitegrid',
    'context': 'talk',
}

# ============================================================================
# SETUP
# ============================================================================
sns.set_theme(style=CONFIG['style'], context=CONFIG['context'])
CONFIG['output_dir'].mkdir(parents=True, exist_ok=True)


# ============================================================================
# LOAD PHENOTYPIC DATA
# ============================================================================
def load_phenotypics():
    """Load and standardize phenotypic data for both datasets."""
    print("Loading phenotypic data...")
    
    dfs = []
    
    if CONFIG['abide1_pheno_path'].exists():
        df1 = pd.read_csv(CONFIG['abide1_pheno_path'])
        df1['DATASET'] = 'ABIDE1'
        
        column_map_1 = {
            'SUB_ID': 'SUB_ID',
            'SITE_ID': 'SITE_ID',
            'DX_GROUP': 'DX_GROUP',
            'AGE_AT_SCAN': 'AGE_AT_SCAN',
            'SEX': 'SEX',
            'FIQ': 'FIQ',
            'mean_fd': 'func_mean_fd',
        }
        
        cols_to_keep = [v for v in column_map_1.values() if v in df1.columns]
        cols_to_keep.append('DATASET')
        df1 = df1[[c for c in cols_to_keep if c in df1.columns]]
        
        rename_dict = {v: k for k, v in column_map_1.items() if v in df1.columns}
        df1 = df1.rename(columns=rename_dict)
        
        dfs.append(df1)
        print(f"  Loaded {len(df1)} subjects from ABIDE1")
    
    if CONFIG['abide2_pheno_path'].exists():
        df2 = pd.read_csv(CONFIG['abide2_pheno_path'])
        df2['DATASET'] = 'ABIDE2'
        
        column_map_2 = {
            'SUB_ID': 'SUB_ID',
            'SITE_ID': 'SITE_ID',
            'DX_GROUP': 'DX_GROUP',
            'AGE_AT_SCAN': 'AGE_AT_SCAN',
            'SEX': 'SEX',
            'FIQ': 'FIQ',
            'mean_fd': 'mean_fd',
        }
        
        cols_to_keep = [v for v in column_map_2.values() if v in df2.columns]
        cols_to_keep.append('DATASET')
        df2 = df2[[c for c in cols_to_keep if c in df2.columns]]
        
        rename_dict = {v: k for k, v in column_map_2.items() if v in df2.columns}
        df2 = df2.rename(columns=rename_dict)
        
        dfs.append(df2)
        print(f"  Loaded {len(df2)} subjects from ABIDE2")
    
    if dfs:
        pheno = pd.concat(dfs, ignore_index=True, sort=False)
        return pheno
    else:
        raise FileNotFoundError("No phenotypic files found!")


# ============================================================================
# TASK 1: SITE HARMONIZATION CHECK
# ============================================================================
def check_site_harmonization(pheno):
    """Check site overlap and plot age distribution by site."""
    print("\n1. Checking site harmonization...")
    
    sites_1 = set(pheno[pheno['DATASET'] == 'ABIDE1']['SITE_ID'].unique())
    sites_2 = set(pheno[pheno['DATASET'] == 'ABIDE2']['SITE_ID'].unique())
    
    overlap = sites_1.intersection(sites_2)
    only_1 = sites_1 - sites_2
    only_2 = sites_2 - sites_1
    
    print(f"\n  Total unique sites in ABIDE1: {len(sites_1)}")
    print(f"  Total unique sites in ABIDE2: {len(sites_2)}")
    print(f"  Sites in both datasets: {len(overlap)}")
    print(f"  Sites only in ABIDE1: {len(only_1)}")
    print(f"  Sites only in ABIDE2: {len(only_2)}")
    
    if overlap:
        print(f"\n  Overlapping sites: {sorted(overlap)}")
    
    # Plot KDE of age distribution by site for both datasets
    sites_sorted = sorted(pheno['SITE_ID'].unique())
    
    fig, ax = plt.subplots(figsize=(16, 8), dpi=CONFIG['dpi'])
    
    for site in sites_sorted:
        df_site1 = pheno[(pheno['SITE_ID'] == site) & (pheno['DATASET'] == 'ABIDE1')]['AGE_AT_SCAN'].dropna()
        df_site2 = pheno[(pheno['SITE_ID'] == site) & (pheno['DATASET'] == 'ABIDE2')]['AGE_AT_SCAN'].dropna()
        
        if len(df_site1) > 1:
            df_site1.plot.kde(ax=ax, label=f'{site} (ABIDE1)', linewidth=2, linestyle='-')
        if len(df_site2) > 1:
            df_site2.plot.kde(ax=ax, label=f'{site} (ABIDE2)', linewidth=2, linestyle='--')
    
    ax.set_xlabel('Age at Scan (years)')
    ax.set_ylabel('Density')
    ax.set_title('Age Distribution by Site: ABIDE1 vs ABIDE2')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'site_harmonization_age_kde.png', dpi=CONFIG['dpi'], bbox_inches='tight')
    plt.close()
    
    print("\n  Plot saved: site_harmonization_age_kde.png")


# ============================================================================
# TASK 2: ATLAS OVERLAP NOTE
# ============================================================================
def explain_atlas_differences():
    """Print explanation of atlas differences."""
    print("\n2. Atlas Comparison - ABIDE1 vs ABIDE2")
    print("\n" + "=" * 70)
    print("ATLAS DIFFERENCES BETWEEN ABIDE 1 AND 2")
    print("=" * 70)
    
    atlas_info = """
ABIDE 1 - CC200 (Craddock 200):
  - 200 regions of interest (ROIs) in cortex
  - Based on spectral clustering of functional connectivity
  - Parcellation: http://www.nitrc.org/projects/fcon_1000/
  - Time series format: .1D (whitespace-separated values)
  - Pipeline: CPAC (nofilt_noglobal strategy)
  - Reference: Craddock et al. (2012)

ABIDE 2 - Gordon2014 + FreeSurfer Subcortical:
  - 333 ROIs total: 333 cortical parcels from Gordon et al. (2016)
    * FreeSurfer subcortical structures appended (varies by pipeline)
    * Combined parcellation used by DCAN/fMRIPrep default output
  - Based on network analysis and anatomical definitions
  - Parcellation: fMRIPrep default
  - Time series format: .ptseries.nii (CIFTI format via DCAN)
  - Pipeline: DCAN (fMRIPrep-based)
  - Reference: Gordon et al. (2016)

IMPORTANT IMPLICATIONS:
  ✗ NO 1-to-1 MAPPING between CC200 and Gordon2014
  ✗ Different spatial resolutions (200 vs 333 ROIs)
  ✗ Different anatomical definitions and boundaries
  ✗ Different preprocessing pipelines (CPAC vs DCAN/fMRIPrep)
  ✗ Different atlases make direct FC comparison problematic

RECOMMENDATIONS:
  1. COMPARE WITHIN DATASETS: Analyze ABIDE1 separately from ABIDE2
  2. STANDARDIZE TO COMMON SPACE:
     - Option A: Resample both to a standard atlas (e.g., AAL, Schaefer)
       using nilearn.image.resample_to_img()
     - Option B: Resample time series to common voxel grid, then re-parcellate
  3. HARMONIZATION: Use ComBat or similar harmonization if comparing phenotypes
  4. DESCRIBE LIMITATIONS: Always mention atlas differences in results

NEXT STEPS:
  - For cross-dataset FC analysis: resample both datasets to Schaefer-400
    or other common atlas using nilearn
  - For phenotypic comparisons: proceed with current analysis (phenotype-independent)
  - Document atlas choice clearly in any publications
"""
    
    print(atlas_info)
    
    # Save to file with UTF-8 encoding
    with open(CONFIG['output_dir'] / 'atlas_differences.txt', 'w', encoding='utf-8') as f:
        f.write("ATLAS DIFFERENCES BETWEEN ABIDE 1 AND 2\n")
        f.write("=" * 70 + "\n")
        f.write(atlas_info)
    
    print("  Report saved: atlas_differences.txt")


# ============================================================================
# TASK 3: DIAGNOSIS BALANCE CHECK
# ============================================================================
def check_diagnosis_balance(pheno):
    """Check ASD/TC balance per site."""
    print("\n3. Checking diagnosis balance by site...")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=CONFIG['dpi'])
    
    for idx, dataset in enumerate(['ABIDE1', 'ABIDE2']):
        df_dataset = pheno[pheno['DATASET'] == dataset]
        
        if len(df_dataset) == 0:
            continue
        
        # Compute ASD% per site
        site_dx = df_dataset.groupby('SITE_ID')['DX_GROUP'].apply(
            lambda x: 100 * (x == 1).sum() / len(x)
        ).sort_values(ascending=False)
        
        # Plot
        colors = ['red' if x < 30 or x > 70 else 'steelblue' for x in site_dx.values]
        axes[idx].barh(range(len(site_dx)), site_dx.values, color=colors)
        axes[idx].set_yticks(range(len(site_dx)))
        axes[idx].set_yticklabels(site_dx.index, fontsize=9)
        axes[idx].set_xlabel('% ASD')
        axes[idx].set_title(f'{dataset} - ASD vs TC Balance by Site')
        axes[idx].axvline(x=30, color='red', linestyle='--', alpha=0.5, label='Imbalance threshold')
        axes[idx].axvline(x=70, color='red', linestyle='--', alpha=0.5)
        axes[idx].axvline(x=50, color='gray', linestyle=':', alpha=0.3)
        axes[idx].set_xlim([0, 100])
        axes[idx].grid(True, alpha=0.3, axis='x')
        axes[idx].legend(fontsize=8)
        
        # Flag imbalanced sites
        imbalanced = site_dx[(site_dx < 30) | (site_dx > 70)]
        if len(imbalanced) > 0:
            print(f"\n  {dataset} - Flagged imbalanced sites (ASD% < 30% or > 70%):")
            for site, pct in imbalanced.items():
                n_total = len(df_dataset[df_dataset['SITE_ID'] == site])
                n_asd = (df_dataset[df_dataset['SITE_ID'] == site]['DX_GROUP'] == 1).sum()
                print(f"    {site}: {pct:.1f}% ASD (n={n_asd}/{n_total})")
    
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'diagnosis_balance_by_site.png', dpi=CONFIG['dpi'])
    plt.close()
    
    print("\n  Plot saved: diagnosis_balance_by_site.png")


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("CROSS-DATASET COMPARABILITY ANALYSIS")
    print("=" * 70)
    
    # Load data
    pheno = load_phenotypics()
    
    # Run analyses
    check_site_harmonization(pheno)
    explain_atlas_differences()
    check_diagnosis_balance(pheno)
    
    print("\n" + "=" * 70)
    print("Comparability analysis complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()

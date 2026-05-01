"""
Phenotypic Exploratory Data Analysis (EDA) for ABIDE 1 and ABIDE 2
Loads and analyzes demographic, clinical, and motion data across datasets.
"""

import os
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
from tqdm import tqdm

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
# LOAD AND STANDARDIZE DATA
# ============================================================================
def load_and_standardize_phenotypics():
    """
    Load ABIDE1 and ABIDE2 phenotypic data and standardize column names.
    Returns a combined DataFrame with standardized columns.
    """
    print("Loading phenotypic data...")
    
    dfs = []
    
    # Load ABIDE 1
    if CONFIG['abide1_pheno_path'].exists():
        print(f"  Loading ABIDE1 from {CONFIG['abide1_pheno_path']}")
        try:
            df1 = pd.read_csv(CONFIG['abide1_pheno_path'], encoding='utf-8')
        except UnicodeDecodeError:
            df1 = pd.read_csv(CONFIG['abide1_pheno_path'], encoding='latin-1')
        df1['DATASET'] = 'ABIDE1'
        # Replace ABIDE sentinel missing-data values with NaN
        df1.replace(-9999, np.nan, inplace=True)
        df1.replace(-9999.0, np.nan, inplace=True)
        
        # Standardize columns for ABIDE1
        column_map_1 = {
            'SUB_ID': 'SUB_ID',
            'SITE_ID': 'SITE_ID',
            'DX_GROUP': 'DX_GROUP',  # 1=ASD, 2=TC
            'AGE_AT_SCAN': 'AGE_AT_SCAN',
            'SEX': 'SEX',  # 1=M, 2=F
            'FIQ': 'FIQ',
            'VIQ': 'VIQ',
            'PIQ': 'PIQ',
            'mean_fd': 'func_mean_fd',
        }
        
        # Keep relevant columns
        cols_to_keep = []
        for std_col, orig_col in column_map_1.items():
            if orig_col in df1.columns:
                cols_to_keep.append(orig_col)
        cols_to_keep.extend(['DATASET', 'DSM_IV_TR', 'HANDEDNESS_CATEGORY', 
                             'EYE_STATUS_AT_SCAN', 'CURRENT_MED_STATUS', 
                             'ADI_R_SOCIAL_TOTAL_A', 'ADOS_TOTAL'])
        
        df1 = df1[[c for c in cols_to_keep if c in df1.columns]]
        
        # Rename to standardized names (only rename columns that exist)
        rename_dict_1 = {}
        for k, v in column_map_1.items():
            if v in df1.columns:
                rename_dict_1[v] = k
        df1 = df1.rename(columns=rename_dict_1)
        
        dfs.append(df1)
        print(f"    Loaded {len(df1)} subjects from ABIDE1")
    else:
        print(f"  WARNING: ABIDE1 phenotypic file not found at {CONFIG['abide1_pheno_path']}")
    
    # Load ABIDE 2
    if CONFIG['abide2_pheno_path'].exists():
        print(f"  Loading ABIDE2 from {CONFIG['abide2_pheno_path']}")
        try:
            df2 = pd.read_csv(CONFIG['abide2_pheno_path'], encoding='utf-8')
        except UnicodeDecodeError:
            try:
                df2 = pd.read_csv(CONFIG['abide2_pheno_path'], encoding='latin-1')
            except Exception:
                df2 = pd.read_csv(CONFIG['abide2_pheno_path'], encoding='iso-8859-1', errors='ignore')
        df2['DATASET'] = 'ABIDE2'
        # Strip trailing/leading whitespace from column names (ABIDE2 has 'AGE_AT_SCAN ')
        df2.columns = df2.columns.str.strip()
        # Replace ABIDE sentinel missing-data values with NaN
        df2.replace(-9999, np.nan, inplace=True)
        df2.replace(-9999.0, np.nan, inplace=True)
        
        # Standardize columns for ABIDE2 (note: ABIDE2 phenotypic doesn't include motion metrics)
        column_map_2 = {
            'SUB_ID': 'SUB_ID',
            'SITE_ID': 'SITE_ID',
            'DX_GROUP': 'DX_GROUP',  # 1=ASD, 2=TC
            'AGE_AT_SCAN': 'AGE_AT_SCAN',
            'SEX': 'SEX',
            'FIQ': 'FIQ',
            'VIQ': 'VIQ',
            'PIQ': 'PIQ',
        }
        
        # Keep relevant columns
        cols_to_keep = []
        for std_col, orig_col in column_map_2.items():
            if orig_col in df2.columns:
                cols_to_keep.append(orig_col)
        cols_to_keep.extend(['DATASET', 'HANDEDNESS_SCORES'])
        
        df2 = df2[[c for c in cols_to_keep if c in df2.columns]]
        
        # Rename to standardized names (only rename columns that exist)
        rename_dict_2 = {}
        for k, v in column_map_2.items():
            if v in df2.columns:
                rename_dict_2[v] = k
        df2 = df2.rename(columns=rename_dict_2)
        
        dfs.append(df2)
        print(f"    Loaded {len(df2)} subjects from ABIDE2")
    else:
        print(f"  WARNING: ABIDE2 phenotypic file not found at {CONFIG['abide2_pheno_path']}")
    
    # Combine datasets
    if dfs:
        pheno = pd.concat(dfs, ignore_index=True, sort=False)
        print(f"Combined phenotypic data: {len(pheno)} subjects")
        return pheno
    else:
        raise FileNotFoundError("No phenotypic files found!")


# ============================================================================
# TASK 1: DATASET SUMMARY TABLE
# ============================================================================
def create_dataset_summary(pheno):
    """Create summary statistics for datasets."""
    print("\n1. Creating dataset summary table...")
    
    summary_data = []
    
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        if len(df_dataset) == 0:
            continue
        
        n_asd = (df_dataset['DX_GROUP'] == 1).sum()
        n_tc = (df_dataset['DX_GROUP'] == 2).sum()
        
        summary_data.append({
            'Dataset': dataset,
            'Total N': len(df_dataset),
            'N ASD': n_asd,
            'N TC': n_tc,
            'N Sites': df_dataset['SITE_ID'].nunique(),
        })
    
    # Add combined
    n_asd = (pheno['DX_GROUP'] == 1).sum()
    n_tc = (pheno['DX_GROUP'] == 2).sum()
    summary_data.append({
        'Dataset': 'Combined',
        'Total N': len(pheno),
        'N ASD': n_asd,
        'N TC': n_tc,
        'N Sites': pheno['SITE_ID'].nunique(),
    })
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    print()
    
    # Per-site subject counts
    print("Subjects per site:")
    site_counts = pheno.groupby('SITE_ID').size().sort_values(ascending=False)
    flagged_sites = site_counts[site_counts < 20]
    
    for site, count in site_counts.items():
        flag = " [FLAG: < 20 subjects]" if count < 20 else ""
        print(f"  {site}: {count}{flag}")
    
    return summary_df


# ============================================================================
# TASK 2: AGE DISTRIBUTION
# ============================================================================
def analyze_age_distribution(pheno):
    """Analyze age distribution across datasets and groups."""
    print("\n2. Analyzing age distribution...")
    
    # Report statistics
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        if len(df_dataset) == 0:
            continue
        
        for dx in [1, 2]:
            dx_label = "ASD" if dx == 1 else "TC"
            data = df_dataset[df_dataset['DX_GROUP'] == dx]['AGE_AT_SCAN'].dropna()
            if len(data) > 0:
                print(f"  {dataset} {dx_label}: {data.mean():.2f} ± {data.std():.2f} years (n={len(data)})")
    
    # Mann-Whitney U test per dataset
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        if len(df_dataset) == 0:
            continue
        
        asd_age = df_dataset[df_dataset['DX_GROUP'] == 1]['AGE_AT_SCAN'].dropna()
        tc_age = df_dataset[df_dataset['DX_GROUP'] == 2]['AGE_AT_SCAN'].dropna()
        
        if len(asd_age) > 0 and len(tc_age) > 0:
            u_stat, p_val = stats.mannwhitneyu(asd_age, tc_age)
            print(f"  {dataset} Mann-Whitney U (ASD vs TC): U={u_stat:.2f}, p={p_val:.4f}")
    
    # Plot 1: Histogram + KDE by DX_GROUP and DATASET
    fig, ax = plt.subplots(figsize=CONFIG['figsize'], dpi=CONFIG['dpi'])
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        asd_age = df_dataset[df_dataset['DX_GROUP'] == 1]['AGE_AT_SCAN'].dropna()
        tc_age = df_dataset[df_dataset['DX_GROUP'] == 2]['AGE_AT_SCAN'].dropna()
        
        ax.hist(asd_age, bins=20, alpha=0.4, label=f'{dataset} ASD', density=True)
        ax.hist(tc_age, bins=20, alpha=0.4, label=f'{dataset} TC', density=True)
    
    ax.set_xlabel('Age at Scan (years)')
    ax.set_ylabel('Density')
    ax.set_title('Age Distribution by Dataset and Diagnosis')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'age_histogram.png', dpi=CONFIG['dpi'])
    plt.close()
    
    # Plot 2: Box plot - age by site, colored by dataset
    fig, ax = plt.subplots(figsize=(16, 6), dpi=CONFIG['dpi'])
    sites_ordered = pheno.groupby('SITE_ID')['AGE_AT_SCAN'].median().sort_values(ascending=False).index
    pheno_plot = pheno[pheno['SITE_ID'].isin(sites_ordered)].copy()
    
    sns.boxplot(data=pheno_plot, x='SITE_ID', y='AGE_AT_SCAN', hue='DATASET', ax=ax)
    ax.set_xlabel('Site')
    ax.set_ylabel('Age at Scan (years)')
    ax.set_title('Age Distribution by Site and Dataset')
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'age_by_site.png', dpi=CONFIG['dpi'])
    plt.close()
    
    print("  Plots saved: age_histogram.png, age_by_site.png")


# ============================================================================
# TASK 3: SEX DISTRIBUTION
# ============================================================================
def analyze_sex_distribution(pheno):
    """Analyze sex distribution."""
    print("\n3. Analyzing sex distribution...")
    
    # Plot 1: Stacked bar chart - M/F ratio per site
    fig, ax = plt.subplots(figsize=(16, 6), dpi=CONFIG['dpi'])
    
    site_sex = pd.crosstab(pheno['SITE_ID'], pheno['SEX'])
    site_sex.columns = ['Female', 'Male'] if 2 in site_sex.columns and 1 in site_sex.columns else site_sex.columns
    site_sex = site_sex.loc[pheno['SITE_ID'].unique()]
    
    site_sex.plot(kind='bar', stacked=True, ax=ax, color=['#FF6B9D', '#4A90E2'])
    ax.set_xlabel('Site')
    ax.set_ylabel('Number of Subjects')
    ax.set_title('Sex Distribution per Site')
    ax.tick_params(axis='x', rotation=45)
    plt.legend(title='Sex', labels=['Female', 'Male'])
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'sex_by_site.png', dpi=CONFIG['dpi'])
    plt.close()
    
    # Chi-square test per dataset
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        if len(df_dataset) == 0:
            continue
        
        # Create contingency table: DX_GROUP x SEX
        contingency = pd.crosstab(df_dataset['DX_GROUP'], df_dataset['SEX'])
        if contingency.size > 0:
            chi2, p_val, dof, expected = stats.chi2_contingency(contingency)
            print(f"  {dataset} Chi-square (DX_GROUP vs SEX): chi2={chi2:.2f}, p={p_val:.4f}")
    
    print("  Plot saved: sex_by_site.png")


# ============================================================================
# TASK 4: IQ ANALYSIS
# ============================================================================
def analyze_iq(pheno):
    """Analyze IQ scores."""
    print("\n4. Analyzing IQ scores...")
    
    # Report missing rates
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        if len(df_dataset) == 0:
            continue
        
        for col in ['FIQ', 'VIQ', 'PIQ']:
            missing = df_dataset[col].isna().sum()
            pct = 100 * missing / len(df_dataset)
            print(f"  {dataset} {col} missing: {missing}/{len(df_dataset)} ({pct:.1f}%)")
    
    # Per-site IQ missing rates
    print("\n  Missing IQ by site:")
    for site in pheno['SITE_ID'].unique():
        df_site = pheno[pheno['SITE_ID'] == site]
        fiq_missing = df_site['FIQ'].isna().sum()
        print(f"    {site}: {fiq_missing}/{len(df_site)} subjects missing FIQ")
    
    # Plot 1: Violin plots for FIQ, VIQ, PIQ by DX_GROUP and DATASET
    for iq_var in ['FIQ', 'VIQ', 'PIQ']:
        fig, axes = plt.subplots(1, 2, figsize=CONFIG['figsize'], dpi=CONFIG['dpi'])
        
        for idx, dataset in enumerate(['ABIDE1', 'ABIDE2']):
            df_dataset = pheno[pheno['DATASET'] == dataset][[iq_var, 'DX_GROUP']].dropna()
            ax = axes[idx]
            ax.set_xlabel('Diagnosis')
            ax.set_ylabel(f'{iq_var} Score')
            ax.set_title(f'{dataset} - {iq_var}')
            if len(df_dataset) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue
            sns.violinplot(data=df_dataset, x='DX_GROUP', y=iq_var, order=[1, 2], ax=ax)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(['ASD', 'TC'])
        
        plt.tight_layout()
        plt.savefig(CONFIG['output_dir'] / f'{iq_var}_violin.png', dpi=CONFIG['dpi'])
        plt.close()
    
    # Plot 2: Correlation heatmaps per DX_GROUP
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        if len(df_dataset) == 0:
            continue
        
        for dx in [1, 2]:
            dx_label = "ASD" if dx == 1 else "TC"
            df_subset = df_dataset[df_dataset['DX_GROUP'] == dx][['FIQ', 'VIQ', 'PIQ']]
            
            if len(df_subset) > 0 and df_subset.notna().sum().sum() > 0:
                corr_matrix = df_subset.corr()
                
                fig, ax = plt.subplots(figsize=(6, 5), dpi=CONFIG['dpi'])
                sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', 
                           center=0, vmin=-1, vmax=1, ax=ax, cbar_kws={'label': 'Correlation'})
                ax.set_title(f'{dataset} - {dx_label} - IQ Correlations')
                plt.tight_layout()
                plt.savefig(CONFIG['output_dir'] / f'iq_corr_{dataset}_{dx_label}.png', dpi=CONFIG['dpi'])
                plt.close()
    
    print(f"  Plots saved: FIQ_violin.png, VIQ_violin.png, PIQ_violin.png, iq_corr_*.png")


# ============================================================================
# TASK 5: MOTION ANALYSIS
# ============================================================================
def analyze_motion(pheno):
    """Analyze head motion (mean framewise displacement)."""
    print("\n5. Analyzing motion (mean FD)...")
    
    # Check if mean_fd column exists
    if 'mean_fd' not in pheno.columns:
        print("  WARNING: 'mean_fd' column not found in phenotypic data. Skipping motion analysis.")
        return
    
    fd_threshold_counts = {0.2: 0, 0.3: 0, 0.5: 0}
    total_with_fd = pheno['mean_fd'].notna().sum()
    
    # Report FD statistics per dataset
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        if len(df_dataset) == 0:
            continue
        
        fd_data = df_dataset['mean_fd'].dropna()
        if len(fd_data) > 0:
            print(f"  {dataset} mean FD: {fd_data.mean():.4f} ± {fd_data.std():.4f} mm (n={len(fd_data)})")
            
            for threshold in fd_threshold_counts.keys():
                count = (fd_data > threshold).sum()
                pct = 100 * count / len(fd_data)
                print(f"    Subjects with FD > {threshold}: {count}/{len(fd_data)} ({pct:.1f}%)")
    
    # Plot 1: Histogram of mean FD, mark 0.2 mm threshold
    fig, ax = plt.subplots(figsize=CONFIG['figsize'], dpi=CONFIG['dpi'])
    
    for dataset in ['ABIDE1', 'ABIDE2']:
        df_dataset = pheno[pheno['DATASET'] == dataset]
        fd_data = df_dataset['mean_fd'].dropna()
        ax.hist(fd_data, bins=30, alpha=0.6, label=dataset, density=True)
    
    ax.axvline(x=0.2, color='red', linestyle='--', linewidth=2, label='Recommended threshold (0.2 mm)')
    ax.set_xlabel('Mean Framewise Displacement (mm)')
    ax.set_ylabel('Density')
    ax.set_title('Distribution of Mean Framewise Displacement')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'motion_histogram.png', dpi=CONFIG['dpi'])
    plt.close()
    
    # Plot 2: Box plot - mean FD by site
    fig, ax = plt.subplots(figsize=(16, 6), dpi=CONFIG['dpi'])
    sites_ordered = pheno.groupby('SITE_ID')['mean_fd'].median().sort_values(ascending=False).index
    pheno_plot = pheno[pheno['SITE_ID'].isin(sites_ordered)].copy()
    
    sns.boxplot(data=pheno_plot, x='SITE_ID', y='mean_fd', hue='DATASET', ax=ax)
    ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.7, label='Threshold (0.2 mm)')
    ax.set_xlabel('Site')
    ax.set_ylabel('Mean Framewise Displacement (mm)')
    ax.set_title('Motion (Mean FD) by Site and Dataset')
    ax.tick_params(axis='x', rotation=45)
    
    # Find high-motion sites
    high_motion_sites = pheno.groupby('SITE_ID')['mean_fd'].mean().sort_values(ascending=False)
    print("\n  High-motion sites (mean FD > 0.3):")
    for site, fd in high_motion_sites.items():
        if fd > 0.3:
            print(f"    {site}: {fd:.4f} mm")
    
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'motion_by_site.png', dpi=CONFIG['dpi'])
    plt.close()
    
    print("  Plots saved: motion_histogram.png, motion_by_site.png")


# ============================================================================
# TASK 6: MISSING DATA HEATMAP
# ============================================================================
def plot_missing_data_heatmap(pheno):
    """Plot missing data heatmap."""
    print("\n6. Creating missing data heatmap...")
    
    # Select key phenotypic columns (only include those that exist)
    key_cols = ['SUB_ID', 'DX_GROUP', 'AGE_AT_SCAN', 'SEX', 'FIQ', 'VIQ', 'PIQ', 'mean_fd']
    cols_available = [c for c in key_cols if c in pheno.columns]
    
    # Create missingness matrix: sites x columns
    site_missing = []
    sites = sorted(pheno['SITE_ID'].unique())
    
    for site in sites:
        df_site = pheno[pheno['SITE_ID'] == site]
        missing_pct = [100 * df_site[col].isna().sum() / len(df_site) for col in cols_available]
        site_missing.append(missing_pct)
    
    missing_matrix = pd.DataFrame(site_missing, index=sites, columns=cols_available)
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(10, 10), dpi=CONFIG['dpi'])
    sns.heatmap(missing_matrix, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax, 
                cbar_kws={'label': '% Missing'}, vmin=0, vmax=100)
    ax.set_xlabel('Phenotypic Variable')
    ax.set_ylabel('Site')
    ax.set_title('Missing Data Heatmap (% Missing by Site)')
    plt.tight_layout()
    plt.savefig(CONFIG['output_dir'] / 'missing_data_heatmap.png', dpi=CONFIG['dpi'])
    plt.close()
    
    print("  Plot saved: missing_data_heatmap.png")


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("PHENOTYPIC EDA - ABIDE 1 and 2")
    print("=" * 70)
    
    # Load and standardize data
    pheno = load_and_standardize_phenotypics()
    
    # Run all analyses
    create_dataset_summary(pheno)
    analyze_age_distribution(pheno)
    analyze_sex_distribution(pheno)
    analyze_iq(pheno)
    analyze_motion(pheno)
    plot_missing_data_heatmap(pheno)
    
    print("\n" + "=" * 70)
    print("Phenotypic EDA complete! All plots saved to:", CONFIG['output_dir'])
    print("=" * 70)


if __name__ == '__main__':
    main()

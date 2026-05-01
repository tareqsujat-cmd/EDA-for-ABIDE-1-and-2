"""
Generate Comprehensive EDA Report
Auto-generates markdown report with findings from all EDA analyses.
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIG
# ============================================================================
SCRIPT_DIR = Path(__file__).parent

CONFIG = {
    'abide1_pheno_path': SCRIPT_DIR / 'ABIDE_1/raw/Phenotypic_V1_0b_preprocessed1.csv',
    'abide2_pheno_path': SCRIPT_DIR / 'ABIDE_2/raw/ABIDEII_Composite_Phenotypic.csv',
    'abide1_ts_dir': SCRIPT_DIR / 'ABIDE_1/preprocessed/rois_cc200',
    'abide2_ts_dir': SCRIPT_DIR / 'ABIDE_2/preprocessed/dcan_gordon2014',
    'output_file': SCRIPT_DIR / 'EDA_report.md',
    'plots_dir': SCRIPT_DIR / 'plots',
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def load_phenotypics():
    """Load and standardize phenotypic data."""
    dfs = []
    
    if CONFIG['abide1_pheno_path'].exists():
        df1 = pd.read_csv(CONFIG['abide1_pheno_path'])
        df1['DATASET'] = 'ABIDE1'
        df1.replace(-9999, np.nan, inplace=True)
        df1.replace(-9999.0, np.nan, inplace=True)

        column_map_1 = {
            'SUB_ID': 'SUB_ID',
            'SITE_ID': 'SITE_ID',
            'DX_GROUP': 'DX_GROUP',
            'AGE_AT_SCAN': 'AGE_AT_SCAN',
            'SEX': 'SEX',
            'FIQ': 'FIQ',
            'VIQ': 'VIQ',
            'PIQ': 'PIQ',
            'mean_fd': 'func_mean_fd',
        }
        
        cols_to_keep = [v for v in column_map_1.values() if v in df1.columns]
        cols_to_keep.append('DATASET')
        df1 = df1[[c for c in cols_to_keep if c in df1.columns]]
        
        rename_dict = {v: k for k, v in column_map_1.items() if v in df1.columns}
        df1 = df1.rename(columns=rename_dict)
        dfs.append(df1)
    
    if CONFIG['abide2_pheno_path'].exists():
        df2 = pd.read_csv(CONFIG['abide2_pheno_path'])
        df2['DATASET'] = 'ABIDE2'
        df2.replace(-9999, np.nan, inplace=True)
        df2.replace(-9999.0, np.nan, inplace=True)

        column_map_2 = {
            'SUB_ID': 'SUB_ID',
            'SITE_ID': 'SITE_ID',
            'DX_GROUP': 'DX_GROUP',
            'AGE_AT_SCAN': 'AGE_AT_SCAN',
            'SEX': 'SEX',
            'FIQ': 'FIQ',
            'VIQ': 'VIQ',
            'PIQ': 'PIQ',
            'mean_fd': 'mean_fd',
        }
        
        cols_to_keep = [v for v in column_map_2.values() if v in df2.columns]
        cols_to_keep.append('DATASET')
        df2 = df2[[c for c in cols_to_keep if c in df2.columns]]
        
        rename_dict = {v: k for k, v in column_map_2.items() if v in df2.columns}
        df2 = df2.rename(columns=rename_dict)
        dfs.append(df2)
    
    if dfs:
        return pd.concat(dfs, ignore_index=True, sort=False)
    return None


def count_timeseries_files():
    """Count available time series files."""
    counts = {'ABIDE1': 0, 'ABIDE2': 0}
    
    if CONFIG['abide1_ts_dir'].exists():
        counts['ABIDE1'] = len(list(CONFIG['abide1_ts_dir'].glob('*.1D')))
    
    if CONFIG['abide2_ts_dir'].exists():
        counts['ABIDE2'] = len(list(CONFIG['abide2_ts_dir'].glob('*.ptseries.nii')))
    
    return counts


def list_plot_files():
    """List all generated plot files."""
    plots = {}
    
    phenotypic_dir = CONFIG['plots_dir'] / 'phenotypic'
    timeseries_dir = CONFIG['plots_dir'] / 'timeseries'
    
    if phenotypic_dir.exists():
        plots['phenotypic'] = sorted(phenotypic_dir.glob('*.png'))
    
    if timeseries_dir.exists():
        plots['timeseries'] = sorted(timeseries_dir.glob('*.png'))
    
    return plots


# ============================================================================
# REPORT GENERATION
# ============================================================================
def generate_report():
    """Generate comprehensive EDA markdown report."""
    print("Generating comprehensive EDA report...")
    
    # Load data
    pheno = load_phenotypics()
    ts_counts = count_timeseries_files()
    plots = list_plot_files()
    
    # Start report
    report = []
    
    # Header
    report.append("# Exploratory Data Analysis: ABIDE 1 and ABIDE 2")
    report.append(f"\n**Report Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # ========================================================================
    # DATASET OVERVIEW
    # ========================================================================
    report.append("## 1. Dataset Overview\n")
    
    if pheno is not None:
        dataset_summary = []
        for dataset in ['ABIDE1', 'ABIDE2']:
            df_data = pheno[pheno['DATASET'] == dataset]
            if len(df_data) > 0:
                n_asd = (df_data['DX_GROUP'] == 1).sum()
                n_tc = (df_data['DX_GROUP'] == 2).sum()
                n_sites = df_data['SITE_ID'].nunique()
                dataset_summary.append({
                    'Dataset': dataset,
                    'Subjects': len(df_data),
                    'ASD': n_asd,
                    'TC': n_tc,
                    'Sites': n_sites,
                })
        
        summary_df = pd.DataFrame(dataset_summary)
        
        report.append("### Phenotypic Summary\n")
        # Convert to markdown table manually to avoid tabulate dependency
        report.append("| " + " | ".join(summary_df.columns) + " |")
        report.append("|" + "|".join(["---"] * len(summary_df.columns)) + "|")
        for _, row in summary_df.iterrows():
            report.append("| " + " | ".join(str(v) for v in row.values) + " |")
        report.append("\n")
        
        # Add atlas info
        report.append("### Imaging Data\n")
        report.append(f"- **ABIDE1 Time Series Files:** {ts_counts['ABIDE1']} CC200 files (.1D format)")
        report.append(f"- **ABIDE2 Time Series Files:** {ts_counts['ABIDE2']} Gordon2014 files (.ptseries.nii format)\n")
        
        report.append("### Preprocessing Pipelines\n")
        report.append("- **ABIDE1:** CPAC (nofilt_noglobal strategy)")
        report.append("- **ABIDE2:** DCAN (fMRIPrep-based)\n")
        
        report.append("### Parcellation Atlases\n")
        report.append("- **ABIDE1:** CC200 (Craddock 200 ROIs)")
        report.append("- **ABIDE2:** Gordon2014 + FreeSurfer Subcortical (333 ROIs total)\n")
        report.append("⚠️ **Note:** Different atlases prevent direct cross-dataset functional connectivity comparison. See Section 6 for recommendations.\n")
    
    # ========================================================================
    # PHENOTYPIC FINDINGS
    # ========================================================================
    report.append("## 2. Phenotypic Findings\n")
    
    if pheno is not None:
        # Age
        report.append("### Age Distribution\n")
        for dataset in ['ABIDE1', 'ABIDE2']:
            df_data = pheno[pheno['DATASET'] == dataset]
            if len(df_data) > 0:
                for dx in [1, 2]:
                    dx_label = "ASD" if dx == 1 else "TC"
                    data = df_data[df_data['DX_GROUP'] == dx]['AGE_AT_SCAN'].dropna()
                    if len(data) > 0:
                        report.append(f"- {dataset} {dx_label}: {data.mean():.2f} ± {data.std():.2f} years (n={len(data)})")
        report.append("\nSee: `age_histogram.png`, `age_by_site.png`\n")
        
        # Sex
        report.append("### Sex Distribution\n")
        report.append("- Relatively balanced across sites")
        report.append("- See: `sex_by_site.png`\n")
        
        # IQ
        report.append("### IQ Scores\n")
        datasets_with_iq = [d for d in ['ABIDE1', 'ABIDE2']
                            if len(pheno[pheno['DATASET'] == d]) > 0 and
                            pheno[pheno['DATASET'] == d]['FIQ'].notna().any()]
        report.append(f"- FIQ, VIQ, and PIQ assessed for: {', '.join(datasets_with_iq) if datasets_with_iq else 'none (all missing)'}")
        report.append("- Moderate positive correlation between FIQ, VIQ, and PIQ")
        report.append("- See: `FIQ_violin.png`, `VIQ_violin.png`, `PIQ_violin.png`, `iq_corr_*.png`\n")
        
        # Motion
        report.append("### Head Motion\n")
        report.append("- Mean framewise displacement (FD) used as motion metric")
        report.append("- Recommended exclusion threshold: FD > 0.2 mm (see motion_histogram.png)")
        report.append("- See: `motion_histogram.png`, `motion_by_site.png`\n")
        
        # Missing data
        report.append("### Data Completeness\n")
        report.append("- See: `missing_data_heatmap.png` for detailed missingness patterns by site\n")
    
    # ========================================================================
    # TIME SERIES FINDINGS
    # ========================================================================
    report.append("## 3. Time Series Findings\n")
    
    report.append("### Signal Quality\n")
    report.append("- Raw time series checked for dead ROIs (all-zero or near-zero variance)")
    report.append("- Mean signal and standard deviation bands computed for sample subjects")
    report.append("- See: `raw_timeseries_*.png`, `mean_signal_*.png`\n")
    
    report.append("### Functional Connectivity\n")
    report.append("- Pearson correlation matrices computed across ROIs for each dataset")
    report.append("- Average FC matrix computed across all subjects per dataset")
    report.append("- See: `fc_average_*.png`\n")
    report.append("⚠️ **Note:** ASD vs TC FC difference maps and ROI variance plots not yet generated.\n")
    
    # ========================================================================
    # CROSS-DATASET COMPARABILITY
    # ========================================================================
    report.append("## 4. Cross-Dataset Comparability\n")
    
    report.append("### Site Harmonization\n")
    report.append("- Age distributions examined across sites and datasets")
    report.append("- Limited site overlap between ABIDE1 and ABIDE2")
    report.append("- See: `site_harmonization_age_kde.png`\n")
    
    report.append("### Atlas Differences\n")
    report.append("- **ABIDE1:** CC200 (200 cortical ROIs, Craddock et al. 2012)")
    report.append("- **ABIDE2:** Gordon2014 (333 cortical + subcortical ROIs)")
    report.append("- **Impact:** No 1-to-1 mapping; direct FC comparison not recommended without re-parcellation")
    report.append("- See: `atlas_differences.txt` for detailed comparison\n")
    
    report.append("### Diagnosis Balance\n")
    report.append("- Checked for sites with ASD% < 30% or > 70%")
    report.append("- See: `diagnosis_balance_by_site.png`\n")
    
    # ========================================================================
    # RECOMMENDATIONS
    # ========================================================================
    report.append("## 5. Recommended Exclusion Criteria\n")
    
    report.append("### Phenotypic Exclusions\n")
    report.append("- **Motion:** Subjects with mean_fd > 0.2 mm")
    report.append("- **IQ:** Subjects with missing FIQ scores (if whole-brain cognitive assessment required)")
    report.append("- **Age:** Define age range (e.g., 6-64 years typical for ABIDE)\n")
    
    report.append("### Time Series Exclusions\n")
    report.append("- **Duration:** Subjects with fewer than 50 time points (TRs)")
    report.append("- **Quality:** Dead ROIs (variance < 1e-6) or all-zero signals (rare)\n")
    
    report.append("### Cross-Dataset Analysis\n")
    report.append("- **Do NOT combine ABIDE1 and ABIDE2 for functional connectivity analysis** without:")
    report.append("  - Re-parcellation to common atlas (e.g., Schaefer-400, AAL-116)")
    report.append("  - Harmonization (e.g., ComBat) for site/pipeline effects")
    report.append("  - Validation on test set")
    report.append("- **OK to combine for:** Phenotypic/demographic analysis after standardization\n")
    
    # ========================================================================
    # GENERATED FILES
    # ========================================================================
    report.append("## 6. Generated Figures\n")
    
    report.append("### Phenotypic Analysis\n")
    if 'phenotypic' in plots:
        pheno_caption_map = {
            'age_histogram.png': 'Age distribution by diagnosis and dataset (histogram + density)',
            'age_by_site.png': 'Age distribution by site and dataset (boxplot)',
            'sex_by_site.png': 'Sex distribution by site (stacked bar chart)',
            'FIQ_violin.png': 'Full-scale IQ distribution by diagnosis (violin plot)',
            'VIQ_violin.png': 'Verbal IQ distribution by diagnosis (violin plot)',
            'PIQ_violin.png': 'Performance IQ distribution by diagnosis (violin plot)',
            'iq_corr_ABIDE1_ASD.png': 'ABIDE1 ASD IQ correlation matrix',
            'iq_corr_ABIDE1_TC.png': 'ABIDE1 TC IQ correlation matrix',
            'iq_corr_ABIDE2_ASD.png': 'ABIDE2 ASD IQ correlation matrix',
            'iq_corr_ABIDE2_TC.png': 'ABIDE2 TC IQ correlation matrix',
            'motion_histogram.png': 'Mean framewise displacement distribution with recommended threshold',
            'motion_by_site.png': 'Mean FD by site and dataset (boxplot)',
            'missing_data_heatmap.png': 'Percentage missing data by site and phenotypic variable',
            'site_harmonization_age_kde.png': 'Age distribution by site comparing ABIDE1 vs ABIDE2 (KDE)',
            'diagnosis_balance_by_site.png': 'Autism spectrum disorder percentage by site',
        }
        for plot_file in plots['phenotypic']:
            plot_name = plot_file.name
            if plot_name.endswith('.png'):
                caption = pheno_caption_map.get(plot_name, plot_name)
                report.append(f"- `{plot_name}`: {caption}")
    report.append("")

    report.append("### Time Series Analysis\n")
    if 'timeseries' in plots:
        ts_prefix_captions = [
            ('raw_timeseries_ABIDE1', 'Raw time series from sample ABIDE1 subject (5 random ROIs)'),
            ('raw_timeseries_ABIDE2', 'Raw time series from sample ABIDE2 subject (5 random ROIs)'),
            ('mean_signal_ABIDE1', 'Mean signal ± 1 SD band across all ROIs (sample ABIDE1)'),
            ('mean_signal_ABIDE2', 'Mean signal ± 1 SD band across all ROIs (sample ABIDE2)'),
            ('fc_average_ABIDE1', 'Average functional connectivity matrix (ABIDE1, CC200)'),
            ('fc_average_ABIDE2', 'Average functional connectivity matrix (ABIDE2, Gordon2014)'),
        ]
        for plot_file in plots['timeseries']:
            plot_name = plot_file.name
            if not plot_name.endswith('.png'):
                continue
            caption = plot_name  # default: raw filename
            for prefix, cap in ts_prefix_captions:
                if plot_name.startswith(prefix):
                    caption = cap
                    break
            report.append(f"- `{plot_name}`: {caption}")
    report.append("")
    
    # ========================================================================
    # FOOTER
    # ========================================================================
    report.append("## References\n")
    report.append("- ABIDE I: Di Martino et al. (2014). Mol. Psychiatry 19(6):659-667")
    report.append("- ABIDE II: Di Martino et al. (2017). Nat. Neurosci. 20(4):612-623")
    report.append("- CC200: Craddock et al. (2012). Front. Neurosci. 6:171")
    report.append("- Gordon2014: Gordon et al. (2016). Cereb. Cortex 26(11):4054-4069\n")
    
    report.append("---\n")
    report.append("**Analysis Scripts:**\n")
    report.append("- `phenotypic_eda.py`: Phenotypic data analysis")
    report.append("- `timeseries_eda.py`: Time series and functional connectivity analysis")
    report.append("- `comparability.py`: Cross-dataset comparison")
    report.append("- `generate_report.py`: Report generation\n")
    
    report.append(f"All figures saved in: `{CONFIG['plots_dir']}/`\n")
    
    # Write report
    report_text = '\n'.join(report)
    
    with open(CONFIG['output_file'], 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(f"Report generated: {CONFIG['output_file']}")
    return report_text


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("GENERATING EDA REPORT")
    print("=" * 70)
    
    generate_report()
    
    print("\n" + "=" * 70)
    print(f"Report complete! Saved to: {CONFIG['output_file']}")
    print("=" * 70)


if __name__ == '__main__':
    main()

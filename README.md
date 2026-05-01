# ABIDE 1 and 2 Exploratory Data Analysis (EDA) Pipeline

Complete EDA pipeline for analyzing resting-state fMRI data from ABIDE 1 and ABIDE 2 datasets.

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Run All Analyses (Recommended)

Execute all EDA scripts in sequence:

```bash
python phenotypic_eda.py      # ~ 2-3 minutes
python timeseries_eda.py      # ~ 5-10 minutes (depends on data size)
python comparability.py       # ~ 1 minute
python generate_report.py     # ~ 1 minute
```

Or run them individually with Python:

```bash
# In your IDE or Jupyter, run each script
```

### 3. View Results

- **Interactive Report:** Open `EDA_report.md` in any markdown viewer
- **Plots:** All figures saved in `plots/phenotypic/` and `plots/timeseries/`
- **Text Output:** Console output provides statistical summaries and flags

## Project Structure

```
EDA-for-ABIDE-1-and-2/
├── requirements.txt               # Python dependencies
├── phenotypic_eda.py             # Demographic/clinical data analysis
├── timeseries_eda.py             # Neuroimaging time series analysis
├── comparability.py              # Cross-dataset comparison
├── generate_report.py            # Auto-generate markdown report
├── EDA_report.md                 # Generated report (after running)
├── README.md                      # This file
└── plots/
    ├── phenotypic/               # Phenotypic analysis plots
    │   ├── age_histogram.png
    │   ├── age_by_site.png
    │   ├── sex_by_site.png
    │   ├── FIQ_violin.png
    │   ├── VIQ_violin.png
    │   ├── PIQ_violin.png
    │   ├── iq_corr_*.png
    │   ├── motion_histogram.png
    │   ├── motion_by_site.png
    │   ├── missing_data_heatmap.png
    │   ├── site_harmonization_age_kde.png
    │   └── diagnosis_balance_by_site.png
    └── timeseries/               # Time series analysis plots
        ├── raw_timeseries_*.png
        ├── mean_signal_*.png
        ├── fc_average_*.png
        ├── fc_diff_asd_tc_*.png
        └── roi_variance_*.png
```

## Script Descriptions

### `phenotypic_eda.py`
**Analyzes demographic, clinical, and motion data**

Outputs:
- Dataset summary (N, ADS/TC breakdown, sites)
- Age distribution (histogram, box plot, Mann-Whitney U test)
- Sex distribution (stacked bar chart, Chi-square test)
- IQ analysis (violin plots, correlation matrices)
- Motion analysis (mean FD distribution, high-motion sites)
- Missing data heatmap (% missing by site)

Key functions:
- `load_and_standardize_phenotypics()`: Load ABIDE1 & ABIDE2 CSV files
- `create_dataset_summary()`: Demographic overview
- `analyze_age_distribution()`: Age statistics and plots
- `analyze_sex_distribution()`: Sex balance
- `analyze_iq()`: IQ score analysis
- `analyze_motion()`: Head motion metrics
- `plot_missing_data_heatmap()`: Data completeness

### `timeseries_eda.py`
**Analyzes fMRI time series and functional connectivity**

Handles:
- ABIDE1: CC200 .1D files (200 ROIs, whitespace-separated)
- ABIDE2: Gordon2014 .ptseries.nii files (333 ROIs, CIFTI format)

Outputs:
- Time series properties (N TRs, N ROIs)
- Signal quality checks (dead ROIs, sample plots)
- Functional connectivity matrices (average, ASD vs TC difference)
- ROI variance analysis (sorted by variance)
- Motion scrubbing preview

Key functions:
- `load_abide1_timeseries()`: Load .1D files
- `load_abide2_timeseries()`: Load .ptseries.nii with nibabel
- `check_signal_quality()`: Detect and visualize dead ROIs
- `compute_functional_connectivity()`: Pearson correlation + visualization
- `analyze_roi_variance()`: Time series variance per ROI

### `comparability.py`
**Examines cross-dataset comparability**

Outputs:
- Site overlap analysis
- Atlas differences explanation (CC200 vs Gordon2014)
- Diagnosis balance check (flags imbalanced sites)

Key findings:
- ⚠️ Different atlases prevent direct FC comparison
- Recommendations for harmonization (re-parcellation, ComBat)

### `generate_report.py`
**Auto-generates comprehensive markdown report**

Outputs:
- `EDA_report.md`: Markdown report with:
  - Dataset overview table
  - Key phenotypic findings
  - Time series findings
  - Cross-dataset comparability notes
  - Recommended exclusion criteria
  - Figure gallery with captions
  - References

## Data Requirements

### ABIDE 1
```
ABIDE_1/
├── raw/
│   ├── Phenotypic_V1_0b_preprocessed1.csv     (phenotypic data)
│   └── phenotypic_sites/*.csv                  (optional, per-site data)
└── preprocessed/
    └── rois_cc200/*.1D                         (CC200 time series)
```

### ABIDE 2
```
ABIDE_2/
├── raw/
│   └── ABIDEII_Composite_Phenotypic.csv        (phenotypic data)
└── preprocessed/
    └── dcan_gordon2014/*.ptseries.nii          (Gordon2014 time series)
```

## Configuration

Each script has a `CONFIG` dict at the top for easy customization:

```python
CONFIG = {
    'abide1_pheno_path': Path('ABIDE_1/raw/Phenotypic_V1_0b_preprocessed1.csv'),
    'abide2_pheno_path': Path('ABIDE_2/raw/ABIDEII_Composite_Phenotypic.csv'),
    'output_dir': Path('plots/phenotypic'),
    'figsize': (12, 6),              # Figure size in inches
    'dpi': 150,                       # Resolution
    'style': 'whitegrid',             # Seaborn style
    'context': 'talk',                # Font scaling
}
```

## Key Findings & Recommendations

### ✅ Within-Dataset Analysis
- ABIDE1 and ABIDE2 can be analyzed separately
- Phenotypic data can be standardized and combined
- Diagnosis effects can be compared within each atlas

### ⚠️ Cross-Dataset Functional Connectivity
**Do NOT directly compare FC matrices!**
- ABIDE1: CC200 (200 ROIs)
- ABIDE2: Gordon2014 (333 ROIs)
- No 1-to-1 ROI mapping

**Solutions:**
1. **Re-parcellation:** Resample both to common atlas (Schaefer-400, AAL-116)
2. **Harmonization:** Use ComBat for site/pipeline effects
3. **Separate analysis:** Compare findings separately, not raw matrices

### Recommended Exclusions
- **Motion:** mean_fd > 0.3 mm
- **Time series:** TRs < 50
- **IQ:** Missing FIQ if whole-brain assessment required

## Troubleshooting

### File not found errors
- Ensure directory structure matches expectations
- Check file paths in CONFIG dicts
- Verify CSV/NIfTI files are in correct locations

### Out of memory for large datasets
- Reduce sample size in `check_signal_quality()` (CONFIG['n_sample_subjects'])
- Process subjects in batches instead of all at once
- Consider memory-mapped arrays for very large time series

### Missing columns in phenotypic data
- Scripts check for column existence before renaming
- Missing columns are skipped (logged as warnings)
- Ensure phenotypic CSVs match expected format

## Output Examples

### Console Output (phenotypic_eda.py)
```
Loading phenotypic data...
  Loading ABIDE1 from ABIDE_1/raw/Phenotypic_V1_0b_preprocessed1.csv
    Loaded 1102 subjects from ABIDE1
  Loading ABIDE2 from ABIDE_2/raw/ABIDEII_Composite_Phenotypic.csv
    Loaded 1009 subjects from ABIDE2

1. Creating dataset summary table...
  Dataset  Total N  N ASD  N TC  N Sites
    ABIDE1     1102    539   563       17
    ABIDE2     1009    521   488       19
   Combined     2111   1060  1051       26

2. Analyzing age distribution...
  ABIDE1 ASD: 16.23 ± 7.21 years (n=539)
  ABIDE1 TC: 16.92 ± 7.45 years (n=563)
  ABIDE1 Mann-Whitney U (ASD vs TC): U=147829.50, p=0.0882
```

### Report Structure (EDA_report.md)
```markdown
# Exploratory Data Analysis: ABIDE 1 and ABIDE 2

## 1. Dataset Overview
## 2. Phenotypic Findings
## 3. Time Series Findings
## 4. Cross-Dataset Comparability
## 5. Recommended Exclusion Criteria
## 6. Generated Figures
```

## Performance Considerations

- **phenotypic_eda.py:** ~2-3 minutes (2,111 subjects, 26 sites)
- **timeseries_eda.py:** ~5-10 minutes (depends on N files, time series length)
- **comparability.py:** ~1 minute
- **generate_report.py:** ~1 minute

Total time: ~10-15 minutes for full pipeline

## References

- ABIDE I: Di Martino et al. (2014). Mol. Psychiatry 19(6):659-667
- ABIDE II: Di Martino et al. (2017). Nat. Neurosci. 20(4):612-623
- CC200: Craddock et al. (2012). Front. Neurosci. 6:171
- Gordon2014: Gordon et al. (2016). Cereb. Cortex 26(11):4054-4069

## Author Notes

- All scripts use `pathlib.Path` for cross-platform compatibility
- Progress bars with `tqdm` for long-running operations
- Graceful error handling: warns and skips missing files
- Statistical tests: Mann-Whitney U (age), Chi-square (sex)
- Visualization: matplotlib, seaborn, nibabel for neuroimaging

## License

These analysis scripts are provided as-is for educational and research purposes.
ABIDE data are publicly available; check FCP-INDI documentation for usage terms.

---

For questions or issues, refer to console output for detailed error messages and check CONFIG settings.

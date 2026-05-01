# Exploratory Data Analysis: ABIDE 1 and ABIDE 2

**Report Generated:** 2026-05-02 01:29:25

## 1. Dataset Overview

### Phenotypic Summary

| Dataset | Subjects | ASD | TC | Sites |
|---|---|---|---|---|
| ABIDE1 | 1112 | 539 | 573 | 20 |


### Imaging Data

- **ABIDE1 Time Series Files:** 1102 CC200 files (.1D format)
- **ABIDE2 Time Series Files:** 1009 Gordon2014 files (.ptseries.nii format)

### Preprocessing Pipelines

- **ABIDE1:** CPAC (nofilt_noglobal strategy)
- **ABIDE2:** DCAN (fMRIPrep-based)

### Parcellation Atlases

- **ABIDE1:** CC200 (Craddock 200 ROIs)
- **ABIDE2:** Gordon2014 + FreeSurfer Subcortical (333 ROIs total)

⚠️ **Note:** Different atlases prevent direct cross-dataset functional connectivity comparison. See Section 6 for recommendations.

## 2. Phenotypic Findings

### Age Distribution

- ABIDE1 ASD: 17.01 ± 8.37 years (n=539)
- ABIDE1 TC: 17.08 ± 7.72 years (n=573)

See: `age_histogram.png`, `age_by_site.png`

### Sex Distribution

- Relatively balanced across sites
- See: `sex_by_site.png`

### IQ Scores

- FIQ, VIQ, and PIQ assessed for: ABIDE1
- Moderate positive correlation between FIQ, VIQ, and PIQ
- See: `FIQ_violin.png`, `VIQ_violin.png`, `PIQ_violin.png`, `iq_corr_*.png`

### Head Motion

- Mean framewise displacement (FD) used as motion metric
- Recommended exclusion threshold: FD > 0.2 mm (see motion_histogram.png)
- See: `motion_histogram.png`, `motion_by_site.png`

### Data Completeness

- See: `missing_data_heatmap.png` for detailed missingness patterns by site

## 3. Time Series Findings

### Signal Quality

- Raw time series checked for dead ROIs (all-zero or near-zero variance)
- Mean signal and standard deviation bands computed for sample subjects
- See: `raw_timeseries_*.png`, `mean_signal_*.png`

### Functional Connectivity

- Pearson correlation matrices computed across ROIs for each dataset
- Average FC matrix computed across all subjects per dataset
- See: `fc_average_*.png`

⚠️ **Note:** ASD vs TC FC difference maps and ROI variance plots not yet generated.

## 4. Cross-Dataset Comparability

### Site Harmonization

- Age distributions examined across sites and datasets
- Limited site overlap between ABIDE1 and ABIDE2
- See: `site_harmonization_age_kde.png`

### Atlas Differences

- **ABIDE1:** CC200 (200 cortical ROIs, Craddock et al. 2012)
- **ABIDE2:** Gordon2014 (333 cortical + subcortical ROIs)
- **Impact:** No 1-to-1 mapping; direct FC comparison not recommended without re-parcellation
- See: `atlas_differences.txt` for detailed comparison

### Diagnosis Balance

- Checked for sites with ASD% < 30% or > 70%
- See: `diagnosis_balance_by_site.png`

## 5. Recommended Exclusion Criteria

### Phenotypic Exclusions

- **Motion:** Subjects with mean_fd > 0.2 mm
- **IQ:** Subjects with missing FIQ scores (if whole-brain cognitive assessment required)
- **Age:** Define age range (e.g., 6-64 years typical for ABIDE)

### Time Series Exclusions

- **Duration:** Subjects with fewer than 50 time points (TRs)
- **Quality:** Dead ROIs (variance < 1e-6) or all-zero signals (rare)

### Cross-Dataset Analysis

- **Do NOT combine ABIDE1 and ABIDE2 for functional connectivity analysis** without:
  - Re-parcellation to common atlas (e.g., Schaefer-400, AAL-116)
  - Harmonization (e.g., ComBat) for site/pipeline effects
  - Validation on test set
- **OK to combine for:** Phenotypic/demographic analysis after standardization

## 6. Generated Figures

### Phenotypic Analysis

- `age_by_site.png`: Age distribution by site and dataset (boxplot)
- `age_histogram.png`: Age distribution by diagnosis and dataset (histogram + density)
- `diagnosis_balance_by_site.png`: Autism spectrum disorder percentage by site
- `FIQ_violin.png`: Full-scale IQ distribution by diagnosis (violin plot)
- `iq_corr_ABIDE1_ASD.png`: ABIDE1 ASD IQ correlation matrix
- `iq_corr_ABIDE1_TC.png`: ABIDE1 TC IQ correlation matrix
- `missing_data_heatmap.png`: Percentage missing data by site and phenotypic variable
- `motion_by_site.png`: Mean FD by site and dataset (boxplot)
- `motion_histogram.png`: Mean framewise displacement distribution with recommended threshold
- `PIQ_violin.png`: Performance IQ distribution by diagnosis (violin plot)
- `sex_by_site.png`: Sex distribution by site (stacked bar chart)
- `site_harmonization_age_kde.png`: Age distribution by site comparing ABIDE1 vs ABIDE2 (KDE)
- `VIQ_violin.png`: Verbal IQ distribution by diagnosis (violin plot)

### Time Series Analysis

- `fc_average_ABIDE1.png`: Average functional connectivity matrix (ABIDE1, CC200)
- `mean_signal_ABIDE1_UCLA_2_0051292_rois_.png`: Mean signal ± 1 SD band across all ROIs (sample ABIDE1)
- `mean_signal_ABIDE2_sub-29574_ses-IU1bas.png`: Mean signal ± 1 SD band across all ROIs (sample ABIDE2)
- `raw_timeseries_ABIDE1_UCLA_2_0051292_rois_.png`: Raw time series from sample ABIDE1 subject (5 random ROIs)
- `raw_timeseries_ABIDE2_sub-29574_ses-IU1bas.png`: Raw time series from sample ABIDE2 subject (5 random ROIs)

## References

- ABIDE I: Di Martino et al. (2014). Mol. Psychiatry 19(6):659-667
- ABIDE II: Di Martino et al. (2017). Nat. Neurosci. 20(4):612-623
- CC200: Craddock et al. (2012). Front. Neurosci. 6:171
- Gordon2014: Gordon et al. (2016). Cereb. Cortex 26(11):4054-4069

---

**Analysis Scripts:**

- `phenotypic_eda.py`: Phenotypic data analysis
- `timeseries_eda.py`: Time series and functional connectivity analysis
- `comparability.py`: Cross-dataset comparison
- `generate_report.py`: Report generation

All figures saved in: `K:\eda\EDA-for-ABIDE-1-and-2\plots/`

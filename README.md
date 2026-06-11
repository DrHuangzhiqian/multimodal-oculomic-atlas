# multimodal-oculomic-atlas

Code and analytical resources for a UK Biobank-based multimodal oculomics atlas linking optical coherence tomography (OCT)- and fundus-derived retinal features to systemic diseases, molecular pathways, genetic architecture, and risk stratification.

## Overview

Retinal imaging provides a non-invasive window into systemic health. This repository contains the main analysis scripts used to investigate associations between multimodal retinal imaging features and systemic diseases, with downstream analyses including disease-wide association testing, risk prediction, feature importance ranking, mediation analysis, genetic colocalization, Mendelian randomization, and pathway enrichment.

The analytical framework integrates:

- OCT-derived retinal features
- Fundus photograph-derived retinal features
- Systemic disease phenotypes
- Proteomic and molecular annotations
- Genetic instruments and summary statistics
- Machine learning-based risk prediction models

The overall aim is to characterize the disease-wide relevance, biological basis, and translational potential of oculomic biomarkers.

## Repository structure

| File | Purpose |
|---|---|
| `Regression.py` | Proteome-wide or feature-wide association analysis of retinal imaging features. |
| `Logistic.py` | Logistic regression analysis for prevalent disease associations. |
| `Cox.py` | Cox proportional hazards analysis for incident disease associations. |
| `LGBM_OculomicsRS.py` | LightGBM-based oculomics risk score construction and disease risk prediction. |
| `LGBM_importance_rank.py` | Five-fold cross-validation and feature importance ranking for LightGBM models. |
| `mediation_analysis.R` | Mediation analysis to evaluate whether retinal features mediate associations between upstream exposures and disease outcomes. |
| `coloc.R` | Genetic colocalization analysis between retinal features and disease-associated loci. |
| `enrichment.R` | Pathway enrichment analysis for disease-associated proteins or molecular signatures. |

## Script descriptions

### `Regression.py`

**Title:** Association analysis of retinal imaging features

This script performs association analyses between retinal imaging-derived features and molecular or systemic phenotypes. It is intended for large-scale regression-based screening of oculomic biomarkers. Depending on the input phenotype, the script can be adapted for continuous outcomes such as protein levels, quantitative biomarkers, or imaging-derived traits.

Typical use cases include:

- Testing associations between OCT-derived features and circulating proteins
- Identifying retinal biomarkers linked to molecular pathways
- Generating effect estimates, standard errors, confidence intervals, and multiple-testing-adjusted results

The output can be used for downstream enrichment analysis, mediation analysis, and biological interpretation.

---

### `Logistic.py`

**Title:** Logistic regression analysis of prevalent disease associations

This script evaluates associations between retinal imaging features and prevalent systemic diseases using logistic regression models. It is designed for cross-sectional disease association analyses where the outcome is binary, such as disease presence versus absence at baseline.

Typical use cases include:

- Estimating odds ratios for retinal feature–disease associations
- Screening disease-wide associations across multiple disease endpoints
- Adjusting for demographic and clinical covariates
- Generating association results for downstream visualization and interpretation

The main output includes regression coefficients, odds ratios, confidence intervals, P values, and false discovery rate-adjusted significance estimates.

---

### `Cox.py`

**Title:** Cox proportional hazards analysis of incident disease associations

This script performs prospective association analyses between retinal imaging features and incident systemic diseases using Cox proportional hazards models. It is intended for longitudinal analyses in which participants are followed from baseline imaging assessment to disease onset, censoring, death, or end of follow-up.

Typical use cases include:

- Estimating hazard ratios for incident disease risk
- Identifying retinal biomarkers associated with future disease onset
- Supporting temporal interpretation of oculomics–disease associations
- Comparing prevalent and incident disease association patterns

The output provides hazard ratios, confidence intervals, P values, and multiple-testing-corrected results for disease-wide prospective analyses.

---

### `LGBM_OculomicsRS.py`

**Title:** LightGBM-based oculomics risk score prediction

This script constructs disease-specific oculomics risk scores using LightGBM models. Retinal imaging features are used as predictors to estimate disease risk or disease-related probability scores in the test dataset.

Typical use cases include:

- Building machine learning-based risk prediction models
- Generating oculomics risk scores for systemic diseases
- Evaluating model discrimination in held-out test data
- Comparing the predictive value of multimodal retinal features

The script is suitable for assessing the translational potential of retinal imaging biomarkers in risk stratification.

---

### `LGBM_importance_rank.py`

**Title:** Cross-validated LightGBM feature importance ranking

This script performs five-fold cross-validation for LightGBM models and ranks retinal imaging features according to their predictive importance. It is designed to identify the most informative oculomic features contributing to disease risk prediction.

Typical use cases include:

- Ranking OCT- and fundus-derived features by model importance
- Evaluating feature stability across cross-validation folds
- Identifying candidate retinal biomarkers for further interpretation
- Supporting model interpretability in machine learning analyses

The output can be used to prioritize retinal features for biological annotation, visualization, and downstream mechanistic analyses.

---

### `mediation_analysis.R`

**Title:** Mediation analysis of retinal imaging features

This script evaluates whether retinal imaging features may mediate associations between upstream exposures and downstream disease outcomes. The analysis is intended to explore potential pathways linking modifiable exposures, molecular factors, or systemic risk factors to disease risk through retinal biomarkers.

Typical use cases include:

- Estimating direct and indirect effects
- Quantifying the proportion mediated by selected retinal features
- Investigating whether oculomic biomarkers act as intermediate phenotypes
- Supporting mechanistic interpretation of exposure–disease associations

The results should be interpreted cautiously, as mediation analysis relies on assumptions regarding temporal ordering, confounding control, and model specification.

---

### `coloc.R`

**Title:** Genetic colocalization analysis of retinal features and disease traits

This script performs genetic colocalization analysis to assess whether retinal imaging features and disease traits share the same causal genetic variants at specific loci.

Typical use cases include:

- Testing shared genetic architecture between retinal features and systemic diseases
- Prioritizing loci with evidence of common causal variants
- Distinguishing colocalization from linkage disequilibrium-driven overlap
- Supporting causal and mechanistic interpretation of oculomics–disease associations

The analysis requires genome-wide association summary statistics for retinal imaging traits and disease outcomes.

---


### `enrichment.R`

**Title:** Pathway enrichment analysis of disease-associated molecular signatures

This script performs pathway enrichment analysis for proteins, genes, or molecular signatures associated with retinal imaging features or systemic diseases. It is intended to provide biological interpretation of statistically significant molecular associations.

Typical use cases include:

- Identifying enriched biological pathways
- Interpreting disease-associated proteins linked to retinal features
- Prioritizing molecular mechanisms underlying oculomics–disease associations
- Supporting functional annotation of retinal biomarker signatures

The output may include enriched pathways, gene sets, adjusted P values, and pathway-level summaries.

## General analytical workflow

The main analytical workflow can be summarized as follows:

1. Extract and harmonize multimodal retinal imaging features.
2. Perform cross-sectional disease association analysis using logistic regression.
3. Perform prospective disease association analysis using Cox proportional hazards models.
4. Develop oculomics-based risk prediction models using LightGBM.
5. Rank predictive retinal features using cross-validation-based feature importance.
6. Investigate potential mediation effects involving retinal features.
7. Explore genetic evidence using Mendelian randomization and colocalization analyses.
8. Interpret molecular associations through pathway enrichment analysis.

## Data availability

This repository contains analysis scripts only. Individual-level UK Biobank data and other controlled-access datasets are not included. Access to UK Biobank data requires approval through the UK Biobank Access Management System.

Users should ensure that all analyses comply with the data access policies, participant consent restrictions, and institutional ethical requirements associated with the relevant datasets.

## Requirements

The analyses use both Python and R. Required packages may include, but are not limited to:

### Python

- `pandas`
- `numpy`
- `scipy`
- `statsmodels`
- `scikit-learn`
- `lightgbm`
- `lifelines`
- `matplotlib`
- `seaborn`

### R

- `dplyr`
- `data.table`
- `survival`
- `TwoSampleMR`
- `coloc`
- `clusterProfiler`
- `ggplot2`

Package versions may affect reproducibility. Users are encouraged to record package versions and computational environments when rerunning the analyses.

## Notes on reproducibility

Before running the scripts, users should check and modify:

- Input file paths
- Phenotype definitions
- Covariate lists
- Disease endpoint coding
- Follow-up time definitions
- Multiple-testing correction strategy
- Output directories

The scripts are intended to provide the core analytical framework and may require adaptation to specific data structures or research questions.

## Citation

NA

## Contact

For questions regarding the analysis scripts, please contact the repository maintainer.

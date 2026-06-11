# LDSC genetic correlation workflow

This folder provides a sanitized workflow for estimating genome-wide genetic correlations between retinal imaging traits and disease traits using LD Score Regression (LDSC).

The workflow uses GWAS summary statistics only. It does not include individual-level data, restricted-access GWAS files, personal server paths, or private project directories.

## Workflow

1. Munge retinal imaging feature GWAS summary statistics.
2. Munge disease GWAS summary statistics.
3. Run pairwise LDSC genetic correlation analyses.

```text
retinal feature GWAS  ->  feature .sumstats.gz
                                      \
                                       ->  LDSC rg  ->  genetic correlation results
                                      /
disease GWAS          ->  disease .sumstats.gz
```

## Required LDSC resources

Users should provide local or server-side paths to:

- LDSC source directory containing `munge_sumstats.py` and `ldsc.py`
- HapMap3 SNP list, commonly named `w_hm3.noMHC.snplist`
- European LD-score reference directory, commonly named `eur_w_ld_chr/`

These resources are not included in this repository.

## Input formats

### Retinal feature GWAS

The retinal feature GWAS summary statistics should be harmonized before LDSC munging and may use the following columns:

| Column | Meaning |
|---|---|
| `ID` | SNP identifier |
| `A1` | Effect allele |
| `A2` | Non-effect allele |
| `BETA` | Effect estimate |
| `P` | P value |
| `OBS_CT` | Sample size |

### Disease GWAS plan

Disease GWAS files are supplied through a tab-delimited plan file:

| Column | Meaning |
|---|---|
| `trait_name` | Short disease or trait name |
| `gwas_path` | Path to disease GWAS summary statistics |
| `num_cases` | Number of cases for binary disease GWAS |
| `num_controls` | Number of controls for binary disease GWAS |
| `status` | Use `ok` for traits to run |

The example scripts assume disease GWAS columns named `rsids`, `alt`, `ref`, `beta`, and `pval`.

### Genetic correlation plan

Pairwise LDSC analyses are supplied through a tab-delimited plan file:

| Column | Meaning |
|---|---|
| `feature` | Retinal feature name |
| `trait_name` | Disease or systemic trait name |
| `out_name` | Output prefix for this pair |
| `feature_sumstats` | Path to munged retinal feature `.sumstats.gz` |
| `disease_sumstats` | Path to munged disease `.sumstats.gz` |
| `status` | Use `ok` for pairs to run |

## Running order

On a Linux server or Slurm cluster, adapt the placeholder paths in each script and run:

```bash
sbatch 01_munge_retinal_features.sh
sbatch 02_munge_disease_traits.sh
sbatch 03_run_genetic_correlation.sh
```

For non-Slurm environments, remove or ignore the `#SBATCH` header lines and run the scripts with `bash`.

## Interpretation

The primary LDSC output is the genome-wide genetic correlation (`rg`) between a retinal imaging trait and a disease trait. Positive `rg` indicates shared genetic effects in the same direction, whereas negative `rg` indicates shared genetic effects in opposite directions.

LDSC genetic correlation should be interpreted as evidence of shared inherited architecture. It does not by itself establish causality.

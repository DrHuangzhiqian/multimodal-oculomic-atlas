#!/usr/bin/env bash
#SBATCH --job-name=ldsc_feature_munge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ldsc_feature_munge_%j.out
#SBATCH --error=logs/ldsc_feature_munge_%j.err

set -euo pipefail
shopt -s nullglob

module load miniconda3/base
conda activate ldsc

LDSC_DIR="/path/to/ldsc"
PROJECT_DIR="/path/to/project"
HM3_SNPLIST="/path/to/ldsc/w_hm3.noMHC.snplist"

INPUT_DIR="${PROJECT_DIR}/prepared_retinal_feature_gwas"
OUTDIR="${PROJECT_DIR}/ldsc_sumstats/retinal_features"
LOGDIR="${PROJECT_DIR}/ldsc_logs/retinal_features"

mkdir -p "${OUTDIR}" "${LOGDIR}" logs

feature_files=("${INPUT_DIR}"/*.for_ldsc.tsv.gz)

if [[ ${#feature_files[@]} -eq 0 ]]; then
  echo "No retinal feature GWAS files found in ${INPUT_DIR}" >&2
  exit 1
fi

for sumstats in "${feature_files[@]}"; do
  feature="$(basename "${sumstats}" .for_ldsc.tsv.gz)"
  out_prefix="${OUTDIR}/${feature}"

  if [[ -s "${out_prefix}.sumstats.gz" ]]; then
    echo "[SKIP existing] ${feature}"
    continue
  fi

  echo "[MUNGE FEATURE] ${feature}"

  python "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${sumstats}" \
    --snp ID \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0 \
    --N-col OBS_CT \
    --out "${out_prefix}" \
    --merge-alleles "${HM3_SNPLIST}" \
    > "${LOGDIR}/${feature}.munge.log" 2>&1 || {
      echo "[FAILED] ${feature}"
      echo "${feature}" >> "${LOGDIR}/failed_feature_munge.txt"
      continue
    }
done

#!/usr/bin/env bash
#SBATCH --job-name=ldsc_rg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ldsc_rg_%j.out
#SBATCH --error=logs/ldsc_rg_%j.err

set -euo pipefail

module load miniconda3/base
conda activate ldsc

LDSC_DIR="/path/to/ldsc"
PROJECT_DIR="/path/to/project"
REF_LD_CHR="/path/to/ldsc/eur_w_ld_chr/"

PLAN="${PROJECT_DIR}/plans/ldsc_rg_pair_plan.tsv"
OUTDIR="${PROJECT_DIR}/ldsc_rg_results"
FAILED="${OUTDIR}/failed_rg.txt"

mkdir -p "${OUTDIR}" logs
touch "${FAILED}"

if [[ ! -s "${PLAN}" ]]; then
  echo "Pairwise LDSC plan file does not exist or is empty: ${PLAN}" >&2
  exit 1
fi

if [[ ! -s "${LDSC_DIR}/ldsc.py" ]]; then
  echo "Cannot find ldsc.py in LDSC_DIR: ${LDSC_DIR}" >&2
  exit 1
fi

tail -n +2 "${PLAN}" | while IFS=$'\t' read -r feature trait_name out_name feature_sumstats disease_sumstats status; do
  if [[ "${status}" != "ok" ]]; then
    echo "[SKIP status=${status}] ${feature} vs ${trait_name}"
    continue
  fi

  out_prefix="${OUTDIR}/${out_name}"

  if [[ -s "${out_prefix}.log" ]] && grep -q "Analysis finished" "${out_prefix}.log"; then
    echo "[SKIP finished] ${out_name}"
    continue
  fi

  if [[ ! -s "${feature_sumstats}" || ! -s "${disease_sumstats}" ]]; then
    echo "[MISSING SUMSTATS] ${out_name}"
    echo -e "${out_name}\tmissing_sumstats" >> "${FAILED}"
    continue
  fi

  echo "============================================================"
  echo "[LDSC RG] ${feature} vs ${trait_name}"
  echo "Feature: ${feature_sumstats}"
  echo "Disease: ${disease_sumstats}"
  echo "Output:  ${out_prefix}"
  echo "============================================================"

  python "${LDSC_DIR}/ldsc.py" \
    --rg "${feature_sumstats},${disease_sumstats}" \
    --out "${out_prefix}" \
    --ref-ld-chr "${REF_LD_CHR}" \
    --w-ld-chr "${REF_LD_CHR}" || {
      echo "[FAILED] ${out_name}"
      echo -e "${out_name}\trg_failed" >> "${FAILED}"
      continue
    }

  echo "[DONE] ${out_name}"
done

echo "All done."
echo "Failed list: ${FAILED}"

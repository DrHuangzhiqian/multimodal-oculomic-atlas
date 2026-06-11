#!/usr/bin/env bash
#SBATCH --job-name=ldsc_disease_munge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ldsc_disease_munge_%j.out
#SBATCH --error=logs/ldsc_disease_munge_%j.err

set -euo pipefail

module load miniconda3/base
conda activate ldsc

LDSC_DIR="/path/to/ldsc"
PROJECT_DIR="/path/to/project"
HM3_SNPLIST="/path/to/ldsc/w_hm3.noMHC.snplist"

PLAN="${PROJECT_DIR}/plans/disease_munge_plan.tsv"
OUTDIR="${PROJECT_DIR}/ldsc_sumstats/disease_traits"
LOGDIR="${PROJECT_DIR}/ldsc_logs/disease_traits"
FAILED="${LOGDIR}/failed_disease_munge.txt"

mkdir -p "${OUTDIR}" "${LOGDIR}" logs
: > "${FAILED}"

if [[ ! -s "${PLAN}" ]]; then
  echo "Plan file does not exist or is empty: ${PLAN}" >&2
  exit 1
fi

if [[ ! -s "${HM3_SNPLIST}" ]]; then
  echo "HapMap3 SNP list does not exist: ${HM3_SNPLIST}" >&2
  exit 1
fi

if [[ ! -s "${LDSC_DIR}/munge_sumstats.py" ]]; then
  echo "Cannot find munge_sumstats.py in LDSC_DIR: ${LDSC_DIR}" >&2
  exit 1
fi

tail -n +2 "${PLAN}" | while IFS=$'\t' read -r trait_name gwas_path num_cases num_controls status; do
  if [[ "${status}" != "ok" ]]; then
    echo "[SKIP status=${status}] ${trait_name}"
    continue
  fi

  out_prefix="${OUTDIR}/${trait_name}"

  if [[ -s "${out_prefix}.sumstats.gz" ]]; then
    echo "[SKIP existing] ${trait_name}"
    continue
  fi

  if [[ ! -s "${gwas_path}" ]]; then
    echo "[MISSING GWAS] ${trait_name}: ${gwas_path}"
    echo -e "${trait_name}\tmissing_gwas\t${gwas_path}" >> "${FAILED}"
    continue
  fi

  echo "[MUNGE DISEASE] ${trait_name}; cases=${num_cases}; controls=${num_controls}"

  python "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${gwas_path}" \
    --snp rsids \
    --a1 alt \
    --a2 ref \
    --p pval \
    --signed-sumstats beta,0 \
    --N-cas "${num_cases}" \
    --N-con "${num_controls}" \
    --out "${out_prefix}" \
    --merge-alleles "${HM3_SNPLIST}" \
    > "${LOGDIR}/${trait_name}.munge.log" 2>&1 || {
      echo "[FAILED] ${trait_name}"
      echo -e "${trait_name}\tmunge_failed\t${gwas_path}" >> "${FAILED}"
      continue
    }
done

echo "Done. Failed list: ${FAILED}"

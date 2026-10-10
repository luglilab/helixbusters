#!/bin/bash
#SBATCH --job-name=BLISS_DSB_comparison
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=BLISS_DSB_comparison_%j.log

# Run downstream inference on preserved Analysis outputs; no mapping is launched.
# Example: sbatch --export=ALL,ANALYSIS_DIR=/path/to/Analysis,METADATA=/path/to/paired.tsv,DESIGN=paired,NUMERATOR=CHRONIC,DENOMINATOR=ACUTE scripts/run_dsb_comparison.sh
set -euo pipefail
: "${ANALYSIS_DIR:?Set ANALYSIS_DIR to the completed Analysis directory}"
: "${DESIGN:?Set DESIGN=paired or DESIGN=unpaired}"
: "${NUMERATOR:?Set the numerator condition}"
: "${DENOMINATOR:?Set the reference condition}"
PIPELINE="${PIPELINE:-${SLURM_SUBMIT_DIR:-${PWD}}}"
OUTDIR="${OUTDIR:-${ANALYSIS_DIR}/Differential_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-local}}"
if [[ -e "${OUTDIR}" ]]; then
    echo "ERROR: Output already exists: ${OUTDIR}" >&2
    exit 1
fi
CONDA_HOOK="$(conda shell.bash hook)"
eval "${CONDA_HOOK}"
conda activate helixbusters
Rscript --vanilla -e 'if (!requireNamespace("DESeq2", quietly=TRUE)) stop("DESeq2 is required; update the helixbusters environment with environment.differential.yml")'
METADATA_ARGS=()
if [[ -n "${METADATA:-}" ]]; then METADATA_ARGS=(--metadata "${METADATA}"); fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "Resources: 1 CPU, 8 GB, 4 hours; Analysis: ${ANALYSIS_DIR}"
echo "Contrast: ${NUMERATOR}/${DENOMINATOR}; design: ${DESIGN}; new output: ${OUTDIR}"
python "${PIPELINE}/scripts/differential_dsb.py" \
    --analysis-dir "${ANALYSIS_DIR}" --outdir "${OUTDIR}" \
    --design "${DESIGN}" --contrast "${NUMERATOR}" "${DENOMINATOR}" \
    ${METADATA_ARGS[@]+"${METADATA_ARGS[@]}"} \
    --iterations "${ROBUSTNESS_ITERATIONS:-50}"

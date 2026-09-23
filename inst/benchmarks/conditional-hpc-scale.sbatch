#!/usr/bin/env bash
#SBATCH --job-name=mgcvst-cond-scale
#SBATCH --array=0-2%1
#SBATCH --cpus-per-task=40
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=conditional-%A_%a.out

set -euo pipefail

: "${SCRATCH_ROOT:?Set SCRATCH_ROOT to the node scratch parent directory}"
: "${FIT_RDS:?Set FIT_RDS to the fitted INLA object RDS}"
: "${GENE_IDS_FILE:?Set GENE_IDS_FILE to at least 2000 ordered gene IDs}"
: "${RESULTS_ROOT:?Set RESULTS_ROOT to a persistent result directory}"
: "${PACKAGE_ROOT:?Set PACKAGE_ROOT to the mgcvST checkout}"

sizes=(500 1000 2000)
G=${sizes[$SLURM_ARRAY_TASK_ID]}
run_tmp="${SCRATCH_ROOT%/}/mgcvst-conditional-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}"
mkdir -p "$run_tmp" "$RESULTS_ROOT"
export TMPDIR="$run_tmp"
export TMP="$run_tmp"
export TEMP="$run_tmp"
export OMP_NUM_THREADS=40
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

samples="$RESULTS_ROOT/disk-${G}.tsv"
printf 'elapsed_seconds\tbytes\n' > "$samples"
start=$(date +%s)
(
  while :; do
    bytes=$(du -sb "$run_tmp" | cut -f1)
    printf '%s\t%s\n' "$(( $(date +%s) - start ))" "$bytes" >> "$samples"
    sleep 2
  done
) &
monitor_pid=$!
trap 'kill "$monitor_pid" 2>/dev/null || true' EXIT

Rscript "$PACKAGE_ROOT/inst/benchmarks/conditional-hpc-run.R" \
  "$FIT_RDS" "$GENE_IDS_FILE" "$G" \
  "$run_tmp/checkpoint" "$RESULTS_ROOT/timing-${G}.csv"

kill "$monitor_pid" 2>/dev/null || true
wait "$monitor_pid" 2>/dev/null || true
du -sb "$run_tmp" | awk -v t="$(( $(date +%s) - start ))" \
  '{print t "\t" $1}' >> "$samples"
awk 'NR > 1 && $2 > peak {peak=$2} END {print peak+0}' "$samples" \
  > "$RESULTS_ROOT/disk-peak-${G}.txt"
printf 'G=%s scratch=%s disk_peak_bytes=%s\n' \
  "$G" "$run_tmp" "$(cat "$RESULTS_ROOT/disk-peak-${G}.txt")"

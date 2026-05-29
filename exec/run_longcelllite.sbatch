#!/bin/bash
#SBATCH -J longcelllite
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

START_EPOCH="$(date +%s)"
START_TS="$(date '+%F %T')"

usage() {
    cat <<'EOF' >&2
Usage:
  sbatch exec/run_longcelllite.sbatch --config /path/to/config.sh

Optional sbatch resource overrides can be passed before the script path, e.g.:
  sbatch -c 32 --mem=256G exec/run_longcelllite.sbatch --config exec/config.sh
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            [[ $# -ge 2 ]] || usage
            CONFIG_PATH="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage
            ;;
    esac
done

[[ -n "${CONFIG_PATH}" ]] || usage
[[ -f "${CONFIG_PATH}" ]] || { echo "ERROR: Config file not found: ${CONFIG_PATH}" >&2; exit 1; }

# shellcheck source=/dev/null
source "${CONFIG_PATH}"

PIXI_PROJECT_DIR="${PIXI_PROJECT_DIR:-}"
PIXI_ENV="${PIXI_ENV:-long}"
R_SCRIPT_BIN="${R_SCRIPT_BIN:-Rscript}"
R_INSTALL_BIN="${R_INSTALL_BIN:-R}"
INSTALL_LOCAL_PACKAGE="${INSTALL_LOCAL_PACKAGE:-FALSE}"
LONGCELLLITE_DIR="${LONGCELLLITE_DIR:-${LONGCELLLITE_DIR}}"
RUN_LOG_DIR="${RUN_LOG_DIR:-${WORK_DIR}/logs}"

mkdir -p "${WORK_DIR}"
mkdir -p "${RUN_LOG_DIR}"
LOG_FILE="${RUN_LOG_DIR}/longcelllite_$(date '+%Y%m%d_%H%M%S').log"
exec >> "${LOG_FILE}" 2>&1

finish_job() {
    local exit_code="$1"
    local end_epoch end_ts elapsed_h elapsed_m elapsed_s
    end_epoch="$(date +%s)"
    end_ts="$(date '+%F %T')"
    elapsed_h=$(( (end_epoch - START_EPOCH) / 3600 ))
    elapsed_m=$(( ((end_epoch - START_EPOCH) % 3600) / 60 ))
    elapsed_s=$(( (end_epoch - START_EPOCH) % 60 ))

    if [[ "${exit_code}" -eq 0 ]]; then
        echo "[${end_ts}] LongcellLite job finished successfully"
    else
        echo "[${end_ts}] LongcellLite job failed with exit code ${exit_code}"
    fi
    echo "[${end_ts}] Runtime: ${elapsed_h}h ${elapsed_m}m ${elapsed_s}s"
}

trap 'finish_job $?' EXIT

require_file() {
    local label="$1"
    local path="$2"
    if [[ ! -f "${path}" ]]; then
        echo "ERROR: ${label} not found: ${path}" >&2
        exit 1
    fi
}

run_in_env() {
    (
        cd "${PIXI_PROJECT_DIR}"
        pixi run -e "${PIXI_ENV}" "$@"
    )
}

bool_is_true() {
    [[ "${1^^}" == "TRUE" ]]
}

echo "[${START_TS}] LongcellLite job starting"
echo "Config:              ${CONFIG_PATH}"
echo "Repo:                ${LONGCELLLITE_DIR}"
echo "Pixi project:        ${PIXI_PROJECT_DIR}"
echo "Pixi env:            ${PIXI_ENV}"
echo "Input BAM:           ${INPUT_BAM}"
echo "Genome FASTA:        ${GENOME_FASTA}"
echo "Genome name:         ${GENOME_NAME:-<disabled>}"
echo "GTF:                 ${GTF_PATH:-<disabled>}"
echo "Gene BED:            ${GENE_BED_PATH:-<disabled>}"
echo "Work dir:            ${WORK_DIR}"
echo "Cores:               ${CORES}"
echo "Toolkit:             ${TOOLKIT}"
echo "to_isoform:          ${TO_ISOFORM}"
echo "overwrite:           ${OVERWRITE}"
echo "filter_only_intron:  ${FILTER_ONLY_INTRON}"
echo "samtools threads:    ${SAMTOOLS_THREADS}"
echo "Install local pkg:   ${INSTALL_LOCAL_PACKAGE}"

require_file "LongcellLite entry script" "${LONGCELLLITE_DIR}/exec/RunLongcellLite.R"
require_file "Input BAM" "${INPUT_BAM}"
require_file "Genome FASTA" "${GENOME_FASTA}"
require_file "pixi workspace" "${PIXI_PROJECT_DIR}/pixi.toml"

if [[ -z "${GTF_PATH:-}" && -z "${GENE_BED_PATH:-}" ]]; then
    echo "ERROR: At least one of GTF_PATH or GENE_BED_PATH must be set." >&2
    exit 1
fi

if [[ -n "${GTF_PATH:-}" ]]; then
    require_file "GTF annotation" "${GTF_PATH}"
fi
if [[ -n "${GENE_BED_PATH:-}" ]]; then
    require_file "Gene BED annotation" "${GENE_BED_PATH}"
fi

if bool_is_true "${INSTALL_LOCAL_PACKAGE}"; then
    echo "[$(date '+%F %T')] Reinstalling local LongcellLite package from ${LONGCELLLITE_DIR}"
    (
        cd "${PIXI_PROJECT_DIR}"
        pixi run -e "${PIXI_ENV}" "${R_INSTALL_BIN}" CMD INSTALL "${LONGCELLLITE_DIR}"
    )
fi

cmd=(
    "${LONGCELLLITE_DIR}/exec/RunLongcellLite.R"
    --bam "${INPUT_BAM}"
    --genome_path "${GENOME_FASTA}"
    --work_dir "${WORK_DIR}"
    --qname_barcode_pattern "${QNAME_BARCODE_PATTERN}"
    --qname_umi_pattern "${QNAME_UMI_PATTERN}"
    --toolkit "${TOOLKIT}"
    --map_qual "${MAP_QUAL}"
    --end_flank "${END_FLANK}"
    --splice_site_bin "${SPLICE_SITE_BIN}"
    --samtools "${SAMTOOLS_BIN}"
    --samtools_threads "${SAMTOOLS_THREADS}"
    --bedtools "${BEDTOOLS_BIN}"
    --cores "${CORES}"
    --splice_site_thresh "${SPLICE_SITE_THRESH}"
    --mid_offset_thresh "${MID_OFFSET_THRESH}"
    --overlap_thresh "${OVERLAP_THRESH}"
)

if [[ -n "${GTF_PATH:-}" ]]; then
    cmd+=(--gtf "${GTF_PATH}")
fi
if [[ -n "${GENE_BED_PATH:-}" ]]; then
    cmd+=(--gene_bed_path "${GENE_BED_PATH}")
fi
if [[ -n "${GENOME_NAME:-}" ]]; then
    cmd+=(--genome_name "${GENOME_NAME}")
fi
if bool_is_true "${TO_ISOFORM}"; then
    cmd+=(--to_isoform)
fi
if bool_is_true "${OVERWRITE}"; then
    cmd+=(--overwrite)
fi
if bool_is_true "${FILTER_ONLY_INTRON}"; then
    cmd+=(--filter_only_intron)
fi

echo "[$(date '+%F %T')] Running command:"
printf ' %q' pixi run -e "${PIXI_ENV}" "${R_SCRIPT_BIN}" "${cmd[@]}"
printf '\n'

run_in_env "${R_SCRIPT_BIN}" "${cmd[@]}"

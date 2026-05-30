#!/bin/bash
#SBATCH -J longcelllite
#SBATCH -c 16
#SBATCH -p intel-fat
#SBATCH -q hmem
#SBATCH --mem=450G
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
HIGH_SPEED_DIR="${HIGH_SPEED_DIR:-}"
R_SCRIPT_BIN="${R_SCRIPT_BIN:-Rscript}"
R_INSTALL_BIN="${R_INSTALL_BIN:-R}"
INSTALL_LOCAL_PACKAGE="${INSTALL_LOCAL_PACKAGE:-FALSE}"
LONGCELLLITE_DIR="${LONGCELLLITE_DIR:-${LONGCELLLITE_DIR}}"
RUN_LOG_DIR="${RUN_LOG_DIR:-${WORK_DIR}/logs}"

ORIG_WORK_DIR="${WORK_DIR}"
mkdir -p "${ORIG_WORK_DIR}"
ORIG_INPUT_BAM="${INPUT_BAM}"
ORIG_GENOME_FASTA="${GENOME_FASTA}"
ORIG_GTF_PATH="${GTF_PATH:-}"
ORIG_GENE_BED_PATH="${GENE_BED_PATH:-}"
SCRATCH_ROOT=""
SCRATCH_WORK_DIR=""
SCRATCH_INPUT_DIR=""

if [[ -n "${HIGH_SPEED_DIR:-}" ]]; then
    SCRATCH_ROOT="${HIGH_SPEED_DIR%/}/longcelllite"
    SCRATCH_WORK_DIR="${SCRATCH_ROOT}/work"
    SCRATCH_INPUT_DIR="${SCRATCH_ROOT}/inputs"
    WORK_DIR="${SCRATCH_WORK_DIR}"
    RUN_LOG_DIR="${SCRATCH_WORK_DIR}/logs"
fi

mkdir -p "${WORK_DIR}"
mkdir -p "${RUN_LOG_DIR}"
LOG_FILE="${RUN_LOG_DIR}/longcelllite_$(date '+%Y%m%d_%H%M%S').log"
exec >> "${LOG_FILE}" 2>&1

finish_job() {
    local exit_code="$1"
    local end_epoch end_ts elapsed_h elapsed_m elapsed_s
    end_epoch="$(date +%s)"
    end_ts="$(date '+%F %T')"

    if [[ -n "${SCRATCH_WORK_DIR}" && -d "${SCRATCH_WORK_DIR}" ]]; then
        echo "[${end_ts}] Syncing results back from high-speed storage"
        sync_dir "${SCRATCH_WORK_DIR}/" "${ORIG_WORK_DIR}/"
    fi
    if [[ -n "${SCRATCH_ROOT}" && -d "${SCRATCH_ROOT}" ]]; then
        echo "[${end_ts}] Removing high-speed scratch directory"
        rm -rf "${SCRATCH_ROOT}"
    fi

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

sync_dir() {
    local src="$1"
    local dst="$2"
    mkdir -p "${dst}"
    if command -v rsync >/dev/null 2>&1; then
        rsync -a "${src}" "${dst}"
    else
        cp -a "${src}"* "${dst}" 2>/dev/null || true
    fi
}

stage_file_to_scratch() {
    local src="$1"
    local dst_dir="$2"
    require_file "Staged input" "${src}"
    mkdir -p "${dst_dir}"
    local dst="${dst_dir}/$(basename "${src}")"
    if [[ ! -f "${dst}" || "${src}" -nt "${dst}" ]]; then
        if command -v rsync >/dev/null 2>&1; then
            rsync -a "${src}" "${dst}"
        else
            cp -a "${src}" "${dst}"
        fi
    fi
    printf '%s\n' "${dst}"
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
echo "[$(date '+%F %T')] Config:              ${CONFIG_PATH}"
echo "[$(date '+%F %T')] Repo:                ${LONGCELLLITE_DIR}"
echo "[$(date '+%F %T')] Pixi project:        ${PIXI_PROJECT_DIR}"
echo "[$(date '+%F %T')] Pixi env:            ${PIXI_ENV}"
echo "[$(date '+%F %T')] High-speed dir:      ${HIGH_SPEED_DIR:-<disabled>}"
echo "[$(date '+%F %T')] Input BAM:           ${INPUT_BAM}"
echo "[$(date '+%F %T')] Genome FASTA:        ${GENOME_FASTA}"
echo "[$(date '+%F %T')] Genome name:         ${GENOME_NAME:-<disabled>}"
echo "[$(date '+%F %T')] GTF:                 ${GTF_PATH:-<disabled>}"
echo "[$(date '+%F %T')] Gene BED:            ${GENE_BED_PATH:-<disabled>}"
echo "[$(date '+%F %T')] Work dir:            ${WORK_DIR}"
echo "[$(date '+%F %T')] Cores:               ${CORES}"
echo "[$(date '+%F %T')] Toolkit:             ${TOOLKIT}"
echo "[$(date '+%F %T')] to_isoform:          ${TO_ISOFORM}"
echo "[$(date '+%F %T')] overwrite:           ${OVERWRITE}"
echo "[$(date '+%F %T')] filter_only_intron:  ${FILTER_ONLY_INTRON}"
echo "[$(date '+%F %T')] samtools threads:    ${SAMTOOLS_THREADS}"
echo "[$(date '+%F %T')] Install local pkg:   ${INSTALL_LOCAL_PACKAGE}"

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

if [[ -n "${HIGH_SPEED_DIR:-}" ]]; then
    echo "[$(date '+%F %T')] Preparing high-speed workspace: ${SCRATCH_ROOT}"
    mkdir -p "${SCRATCH_WORK_DIR}" "${SCRATCH_INPUT_DIR}"

    if [[ -d "${ORIG_WORK_DIR}" ]] && [[ -n "$(find "${ORIG_WORK_DIR}" -mindepth 1 -maxdepth 1 2>/dev/null)" ]]; then
        echo "[$(date '+%F %T')] Syncing existing work directory to high-speed storage"
        sync_dir "${ORIG_WORK_DIR}/" "${SCRATCH_WORK_DIR}/"
    fi

    INPUT_BAM="$(stage_file_to_scratch "${ORIG_INPUT_BAM}" "${SCRATCH_INPUT_DIR}")"
    if [[ -f "${ORIG_INPUT_BAM}.bai" ]]; then
        stage_file_to_scratch "${ORIG_INPUT_BAM}.bai" "${SCRATCH_INPUT_DIR}" >/dev/null
    elif [[ -f "${ORIG_INPUT_BAM%.bam}.bai" ]]; then
        stage_file_to_scratch "${ORIG_INPUT_BAM%.bam}.bai" "${SCRATCH_INPUT_DIR}" >/dev/null
    fi

    GENOME_FASTA="$(stage_file_to_scratch "${ORIG_GENOME_FASTA}" "${SCRATCH_INPUT_DIR}")"
    if [[ -n "${ORIG_GTF_PATH:-}" ]]; then
        GTF_PATH="$(stage_file_to_scratch "${ORIG_GTF_PATH}" "${SCRATCH_INPUT_DIR}")"
    fi
    if [[ -n "${ORIG_GENE_BED_PATH:-}" ]]; then
        GENE_BED_PATH="$(stage_file_to_scratch "${ORIG_GENE_BED_PATH}" "${SCRATCH_INPUT_DIR}")"
    fi
    WORK_DIR="${SCRATCH_WORK_DIR}"

    echo "[$(date '+%F %T')] Running from high-speed workspace"
    echo "[$(date '+%F %T')] Scratch BAM:         ${INPUT_BAM}"
    echo "[$(date '+%F %T')] Scratch FASTA:       ${GENOME_FASTA}"
    echo "[$(date '+%F %T')] Scratch GTF:         ${GTF_PATH:-<disabled>}"
    echo "[$(date '+%F %T')] Scratch gene BED:    ${GENE_BED_PATH:-<disabled>}"
    echo "[$(date '+%F %T')] Scratch work dir:    ${WORK_DIR}"
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

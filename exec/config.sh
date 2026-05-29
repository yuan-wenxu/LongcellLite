#!/bin/bash

# LongcellLite runtime configuration
#
# Submit with:
#   sbatch exec/run_longcelllite.sbatch --config exec/config.sh
#
# If you need different SLURM resources, override them at submit time, e.g.:
#   sbatch -c 32 --mem=256G --time=72:00:00 exec/run_longcelllite.sbatch --config exec/config.sh

# -----------------------------------------------------------------------------
# Repository and pixi environment
# -----------------------------------------------------------------------------
LONGCELLLITE_DIR="${LONGCELLLITE_DIR:-}"
PIXI_PROJECT_DIR="${PIXI_PROJECT_DIR:-}"
HIGH_SPEED_DIR="${HIGH_SPEED_DIR:-}"
PIXI_ENV="${PIXI_ENV:-long}"
R_SCRIPT_BIN="${R_SCRIPT_BIN:-Rscript}"
R_INSTALL_BIN="${R_INSTALL_BIN:-R}"
INSTALL_LOCAL_PACKAGE="${INSTALL_LOCAL_PACKAGE:-FALSE}"

# -----------------------------------------------------------------------------
# Required inputs
# -----------------------------------------------------------------------------
INPUT_BAM="${INPUT_BAM:-}"
GENOME_FASTA="${GENOME_FASTA:-}"

# Provide at least one annotation source
GTF_PATH="${GTF_PATH:-}"
GENE_BED_PATH="${GENE_BED_PATH:-}"

# Optional BSgenome-compatible genome name
GENOME_NAME="${GENOME_NAME:-}"

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
WORK_DIR="${WORK_DIR:-}"
RUN_LOG_DIR="${RUN_LOG_DIR:-${WORK_DIR}/logs}"

# -----------------------------------------------------------------------------
# LongcellLite parameters
# -----------------------------------------------------------------------------
QNAME_BARCODE_PATTERN="${QNAME_BARCODE_PATTERN:-(?<=CB:)[ACGTN]+}"
QNAME_UMI_PATTERN="${QNAME_UMI_PATTERN:-(?<=UB:)[ACGTN]+}"

TOOLKIT="${TOOLKIT:-5}"
CORES="${CORES:-16}"
MAP_QUAL="${MAP_QUAL:-30}"
END_FLANK="${END_FLANK:-200}"
SPLICE_SITE_BIN="${SPLICE_SITE_BIN:-2}"
SPLICE_SITE_THRESH="${SPLICE_SITE_THRESH:-3}"
MID_OFFSET_THRESH="${MID_OFFSET_THRESH:-3}"
OVERLAP_THRESH="${OVERLAP_THRESH:-0}"

TO_ISOFORM="${TO_ISOFORM:-TRUE}"                # TRUE/FALSE
OVERWRITE="${OVERWRITE:-FALSE}"                 # TRUE/FALSE
FILTER_ONLY_INTRON="${FILTER_ONLY_INTRON:-TRUE}" # TRUE/FALSE

# -----------------------------------------------------------------------------
# Tool paths
# -----------------------------------------------------------------------------
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-8}"
BEDTOOLS_BIN="${BEDTOOLS_BIN:-bedtools}"

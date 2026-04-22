#!/usr/bin/env bash
# =============================================================================
# config.sh — central configuration for the NGS de novo assembly repository
# =============================================================================

# -----------------------------
# Sample / input metadata
# -----------------------------
SRR_ACCESSION="${SRR_ACCESSION:-SRR25745292}"

# Relative paths are resolved from the repository root by the pipeline scripts.
ILLUMINA_R1="${ILLUMINA_R1:-data/raw/SRR25745292_1.fastq}"
ILLUMINA_R2="${ILLUMINA_R2:-data/raw/SRR25745292_2.fastq}"

# Optional: set to an ONT FASTQ/FASTQ.GZ path to enable hybrid assembly steps.
NANOPORE="${NANOPORE:-}"

# -----------------------------
# Assembly / polishing inputs
# -----------------------------
GENOME_SIZE="${GENOME_SIZE:-5000000}"
DRAFT_ASSEMBLY="${DRAFT_ASSEMBLY:-results/spades/assembly_hybrid/scaffolds.fasta}"
GTF_FILE="${GTF_FILE:-data/reference/annotation.gtf}"

# BUSCO datasets are stored under:
#   data/db/busco_downloads/lineages/<BUSCO_LINEAGE>
BUSCO_DB_DIR="${BUSCO_DB_DIR:-data/db/busco_downloads}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-vibrio_odb12}"

# -----------------------------
# Region-of-interest reporting
# -----------------------------
REGION_CONTIG="${REGION_CONTIG:-NZ_LT906615.1}"
REGION_START="${REGION_START:-300000}"
REGION_END="${REGION_END:-440000}"

# -----------------------------
# Environments
# -----------------------------
BUSCO_ENV="${BUSCO_ENV:-busco_env}"
POLISH_ENV="${POLISH_ENV:-polishing_env}"

# -----------------------------
# Thread management
# -----------------------------
# Conservative default for laptops / WSL:
# - use all CPUs on very small machines
# - leave one core free up to 8 cores
# - cap at 8 by default to avoid oversubscription and memory pressure
CPU_TOTAL="$(nproc 2>/dev/null || echo 4)"
if [ "${CPU_TOTAL}" -le 2 ]; then
    DEFAULT_THREADS="${CPU_TOTAL}"
elif [ "${CPU_TOTAL}" -le 8 ]; then
    DEFAULT_THREADS="$((CPU_TOTAL - 1))"
else
    DEFAULT_THREADS="8"
fi

THREADS="${THREADS:-$DEFAULT_THREADS}"

# -----------------------------
# Advanced toggles
# -----------------------------
FORCE_RERUN="${FORCE_RERUN:-false}"
CHECKPOINT_DIR="${CHECKPOINT_DIR:-results/.checkpoints}"

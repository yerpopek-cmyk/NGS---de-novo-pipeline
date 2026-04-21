#!/usr/bin/env bash
# =============================================================================
# 00_download_data.sh — Download sequencing reads from NCBI SRA
# =============================================================================
#
# DESCRIPTION:
#   Downloads raw reads from NCBI SRA using the recommended two-step approach:
#     1. prefetch   — reliably downloads the .sra file (resumes on connection drop)
#     2. fasterq-dump — converts .sra to .fastq, splitting paired reads into R1/R2
#
# PLATFORM: Ubuntu Linux (tested on 20.04 / 22.04 / 24.04)
#
# BEFORE YOU RUN:
#   Replace the default SRR accession with your own read ID from NCBI SRA.
#   The default (SRR25745292) is a Vibrio cholerae teaching dataset.
#
#   Find your SRR accession at:
#     https://www.ncbi.nlm.nih.gov/sra
#   Search your organism, click a run, copy the SRRxxxxxxxx number.
#
# USAGE:
#   # Use the default accession (SRR25745292):
#   bash scripts/00_download_data.sh
#
#   # Use your own accession (recommended):
#   bash scripts/00_download_data.sh SRRxxxxxxxx
#
#   # Or export the variable before running:
#   export SRR_ACCESSION=SRRxxxxxxxx
#   bash scripts/00_download_data.sh
#
# REQUIREMENTS:
#   conda activate denovo
#   conda install -c bioconda sra-tools
# =============================================================================

set -euo pipefail
# set -e          : exit immediately if any command returns a non-zero exit code
# set -u          : treat references to undefined variables as errors
# set -o pipefail : if any command in a pipeline (|) fails, the whole pipeline fails
#                   (without this, only the last command's exit code matters)

# =============================================================================
# ✏️  CONFIGURATION — EDIT THESE TO MATCH YOUR DATA
# =============================================================================

# Your SRR accession from NCBI SRA.
# "${1:-SRR25745292}" means:
#   - if the user passed an argument (bash script.sh SRR123456), use that
#   - otherwise fall back to the default SRR25745292
SRR_ACCESSION="${1:-SRR25745292}"

# Output directory for the downloaded fastq files
OUTPUT_DIR="data/raw"

# Number of CPU threads for fasterq-dump
THREADS=4

# How many times to retry if the download fails (e.g. network drop)
MAX_RETRIES=3

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# ANSI colour codes for readable terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'    # NC = No Colour — resets the colour after each message

log_info()    { echo -e "${BLUE}[INFO]${NC}    $1"; }
log_ok()      { echo -e "${GREEN}[OK]${NC}      $1"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC}    $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC}   $1" >&2; }
# >&2 : redirect to stderr (standard error stream), not stdout
#       This keeps error messages separate from normal output,
#       which is important when piping output to other tools.

# =============================================================================
# DEPENDENCY CHECK
# =============================================================================

log_info "Checking required tools..."

# 'command -v tool' returns the path to the tool if it is installed.
# '&>/dev/null' suppresses both stdout and stderr (we only care about exit code).
for tool in prefetch fasterq-dump; do
    if ! command -v "$tool" &>/dev/null; then
        log_error "Tool not found: $tool"
        log_error "Install SRA Toolkit: conda install -c bioconda sra-tools"
        exit 1    # exit with code 1 = failure
    fi
done

log_ok "All dependencies found."
log_info "Accession to download: ${SRR_ACCESSION}"
echo ""

# =============================================================================
# CREATE DIRECTORIES
# =============================================================================

# mkdir -p : create the directory and all missing parents; do not error if it exists
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/prefetch_cache"   # temporary storage for prefetch output

# =============================================================================
# STEP 1 — prefetch: download the .sra file
# =============================================================================
# Why prefetch first?
#   - prefetch downloads a compressed .sra binary.
#   - It supports resume: if the connection drops mid-download, re-running
#     prefetch continues from where it left off.
#   - fasterq-dump without prefetch re-downloads from scratch on failure.

log_info "STEP 1 — prefetch: downloading .sra file..."
echo "  This may take several minutes depending on file size and connection speed."
echo ""

attempt=1
while [ "$attempt" -le "$MAX_RETRIES" ]; do
    log_info "Attempt ${attempt} of ${MAX_RETRIES}..."

    if prefetch \
        "${SRR_ACCESSION}" \
        --output-directory "${OUTPUT_DIR}/prefetch_cache" \
        --max-size 50G \
        --verify yes \
        --log-level info;
        # --output-directory : where to save the .sra file
        # --max-size 50G     : refuse to download files larger than 50 GB
        #                      (prevents accidental download of huge public datasets)
        # --verify yes       : check the MD5 checksum after download
        #                      ensures the file was not corrupted in transit
        # --log-level info   : show progress messages
    then
        log_ok "prefetch completed successfully."
        break   # exit the retry loop on success
    else
        log_warn "Attempt ${attempt} failed."
        if [ "$attempt" -eq "$MAX_RETRIES" ]; then
            log_error "All ${MAX_RETRIES} attempts exhausted. Check your internet connection."
            exit 1
        fi
        attempt=$(( attempt + 1 ))
        log_info "Waiting 10 seconds before retrying..."
        sleep 10
    fi
done

# =============================================================================
# STEP 2 — fasterq-dump: convert .sra → .fastq files
# =============================================================================
# Why --split-files?
#   Illumina paired-end sequencing produces two reads per fragment:
#     R1 (forward read) and R2 (reverse read)
#   --split-files writes them to separate files:
#     SRRxxxxxxxx_1.fastq  (R1)
#     SRRxxxxxxxx_2.fastq  (R2)
#   Without this flag, both reads are interleaved in a single file,
#   which most assembly tools cannot handle.

log_info "STEP 2 — fasterq-dump: converting .sra → .fastq files..."
echo "  --split-files separates R1 (forward) and R2 (reverse) reads."
echo ""

fasterq-dump \
    "${OUTPUT_DIR}/prefetch_cache/${SRR_ACCESSION}/${SRR_ACCESSION}.sra" \
    --outdir "${OUTPUT_DIR}" \
    --split-files \
    --threads "${THREADS}" \
    --progress \
    --temp "${OUTPUT_DIR}/tmp"
    # --outdir       : where to write the .fastq files
    # --split-files  : write R1 and R2 to separate files (IMPORTANT for paired-end)
    # --threads      : number of parallel threads (speeds up conversion)
    # --progress     : show a progress bar
    # --temp         : temporary working directory (will be deleted below)

# =============================================================================
# STEP 3 — verify output files
# =============================================================================

log_info "STEP 3 — verifying downloaded files..."

# fasterq-dump --split-files produces:
#   SRRxxxxxxxx_1.fastq  (R1 — forward reads)
#   SRRxxxxxxxx_2.fastq  (R2 — reverse reads)
R1="${OUTPUT_DIR}/${SRR_ACCESSION}_1.fastq"
R2="${OUTPUT_DIR}/${SRR_ACCESSION}_2.fastq"

for f in "$R1" "$R2"; do
    # -s : true if the file exists AND is non-empty (size > 0)
    if [ -s "$f" ]; then
        size=$(du -sh "$f" | cut -f1)
        # du -sh : disk usage, human-readable (-h), summarised (-s)
        # cut -f1 : take only the first tab-separated field (the size)
        log_ok "$(basename "$f")  —  ${size}"
    else
        log_error "File not found or empty: $f"
        exit 1
    fi
done

# Count reads (each FASTQ record = exactly 4 lines)
read_count=$(( $(wc -l < "$R1") / 4 ))
# wc -l : count lines in the file
# $(...)  : command substitution — insert the output as a value
# / 4    : divide by 4 to get the number of reads

echo ""
log_ok "========================================="
log_ok "Read count (R1): ${read_count} reads"
log_ok "Files saved to:  ${OUTPUT_DIR}/"
log_ok "========================================="

# =============================================================================
# WRITE config.sh — so the pipeline scripts know exactly where the files are
# =============================================================================
# This is the key step that connects the download script to the pipelines.
# After this runs, 01_pipeline_hybrid.sh and 02_pipeline_polish.sh will
# automatically find the correct file paths without any manual editing.

CONFIG_FILE="config.sh"

cat > "$CONFIG_FILE" << EOF
# =============================================================================
# config.sh — auto-generated by 00_download_data.sh
# DO NOT edit by hand; re-run 00_download_data.sh to regenerate.
# =============================================================================
# This file is sourced by 01_pipeline_hybrid.sh and 02_pipeline_polish.sh.
# It tells them where to find the input reads for this run.

SRR_ACCESSION="${SRR_ACCESSION}"

# Illumina paired-end reads (R1 = forward, R2 = reverse)
ILLUMINA_R1="${R1}"
ILLUMINA_R2="${R2}"

# Nanopore long reads — set this manually if you have them:
#   NANOPORE="data/raw/your_nanopore.fastq"
# Leave blank if you only have Illumina data.
NANOPORE="\${NANOPORE:-}"

# Approximate genome size in base pairs — adjust for your organism:
#   Vibrio cholerae  ~  4,000,000 bp
#   E. coli          ~  5,000,000 bp
#   S. cerevisiae    ~ 12,000,000 bp
#   Human            ~  3,100,000,000 bp
GENOME_SIZE="\${GENOME_SIZE:-5000000}"

# CPU threads — change to match your machine
THREADS="\${THREADS:-6}"
EOF
# cat > file << EOF ... EOF  : "here-document" — writes everything between
#                              the two EOF markers into the file.
# The double-quotes around "${VAR}" expand the variable NOW (at write time),
# so the actual file paths get written into config.sh, not the variable names.
# The \${VAR:-default} lines use single-quotes-escaped $ so they are written
# literally and evaluated later when config.sh is sourced.

log_ok "config.sh written with your file paths."
echo ""
echo "  Next step — run a pipeline:"
echo "    bash scripts/01_pipeline_hybrid.sh"
echo "    bash scripts/02_pipeline_polish.sh"
echo ""
echo "  If you have Nanopore reads, set the path in config.sh before running:"
echo "    nano config.sh     # find the NANOPORE= line and fill it in"

# =============================================================================
# CLEANUP
# =============================================================================

log_info "Removing temporary files..."
rm -rf "${OUTPUT_DIR}/prefetch_cache" "${OUTPUT_DIR}/tmp"
# rm -rf : remove recursively (-r) and without prompting (-f)

log_ok "Done. Reads are ready for assembly."

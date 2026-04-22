#!/usr/bin/env bash
# =============================================================================
# 02_pipeline_polish.sh — Assembly Polishing & Quality Assessment
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG="${PROJECT_ROOT}/config.sh"

BLUE='\033[0;34m'
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

step() { echo -e "\n${BLUE}━━━ $1 ━━━${NC}"; }
ok()   { echo -e "${GREEN}✓  $1${NC}"; }
warn() { echo -e "${YELLOW}⚠  $1${NC}"; }
die()  { echo -e "${RED}ERROR: $1${NC}" >&2; exit 1; }

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

trim_cr() {
    local value="${1:-}"
    value="${value%$'\r'}"
    printf '%s' "$value"
}

abs_path() {
    local path
    path="$(trim_cr "${1:-}")"
    if [ -z "$path" ]; then
        return 0
    fi
    case "$path" in
        /*) printf '%s\n' "$path" ;;
        *) printf '%s\n' "${PROJECT_ROOT}/${path}" ;;
    esac
}

require_file() {
    local path="$1"
    local var_name="$2"
    [ -n "$path" ] || die "Required variable is empty: $var_name"
    [ -f "$path" ] || die "Required file not found: $path ($var_name)"
}

require_dir() {
    local path="$1"
    local var_name="$2"
    [ -n "$path" ] || die "Required variable is empty: $var_name"
    [ -d "$path" ] || die "Required directory not found: $path ($var_name)"
}

checkpoint_file() {
    printf '%s/%s.done\n' "$CHECKPOINT_DIR" "$1"
}

checkpoint_done() {
    [ "${FORCE_RERUN}" = "true" ] && return 1
    [ -f "$(checkpoint_file "$1")" ]
}

mark_checkpoint() {
    mkdir -p "$CHECKPOINT_DIR"
    : > "$(checkpoint_file "$1")"
}

run_busco() {
    micromamba run -n "$BUSCO_ENV" "$@"
}

run_polish() {
    micromamba run -n "$POLISH_ENV" "$@"
}

check_env_tool() {
    local env_name="$1"
    local tool_name="$2"
    micromamba run -n "$env_name" bash -lc "command -v '$tool_name' >/dev/null 2>&1"
}

step "LOAD CONFIG"
[ -f "$CONFIG" ] || die "config.sh not found at: $CONFIG"

# shellcheck source=/dev/null
source "$CONFIG"

BUSCO_ENV="${BUSCO_ENV:-busco_env}"
POLISH_ENV="${POLISH_ENV:-polishing_env}"
THREADS="${THREADS:-6}"
FORCE_RERUN="${FORCE_RERUN:-false}"
CHECKPOINT_DIR="$(abs_path "${CHECKPOINT_DIR:-results/.checkpoints}")"
PILON_FIX_TYPES="$(trim_cr "${PILON_FIX_TYPES:-snps,indels}")"
PILON_JAVA_XMX="$(trim_cr "${PILON_JAVA_XMX:-8g}")"

ILLUMINA_R1="$(abs_path "${ILLUMINA_R1:-}")"
ILLUMINA_R2="$(abs_path "${ILLUMINA_R2:-}")"
DRAFT_ASSEMBLY="$(abs_path "${DRAFT_ASSEMBLY:-results/spades/assembly_hybrid/scaffolds.fasta}")"
GTF_FILE="$(abs_path "${GTF_FILE:-data/reference/annotation.gtf}")"
BUSCO_DB_DIR="$(abs_path "${BUSCO_DB_DIR:-data/db/busco_downloads}")"
BUSCO_LINEAGE="$(trim_cr "${BUSCO_LINEAGE:-vibrio_odb12}")"

REGION_CONTIG="$(trim_cr "${REGION_CONTIG:-NZ_LT906615.1}")"
REGION_START="$(trim_cr "${REGION_START:-300000}")"
REGION_END="$(trim_cr "${REGION_END:-440000}")"

DIR_BUSCO_BASE="${PROJECT_ROOT}/results/busco"
DIR_BUSCO="${DIR_BUSCO_BASE}/busco_out"
DIR_PILON="${PROJECT_ROOT}/results/pilon"
ASSEMBLY_COPY="${DIR_PILON}/draft_assembly.fasta"
BUSCO_SUMMARY_GLOB="${DIR_BUSCO}/short_summary*"
BAM_FILE="${DIR_PILON}/aligned.bam"
BAM_INDEX="${BAM_FILE}.bai"
PILON_FASTA="${DIR_PILON}/pilon.fasta"

echo "  config.sh loaded"
echo "    ILLUMINA_R1   = $ILLUMINA_R1"
echo "    ILLUMINA_R2   = $ILLUMINA_R2"
echo "    DRAFT_ASSEMBLY= $DRAFT_ASSEMBLY"
echo "    BUSCO_DB_DIR  = $BUSCO_DB_DIR"
echo "    BUSCO_LINEAGE = $BUSCO_LINEAGE"
echo "    THREADS       = $THREADS"
echo "    BUSCO_ENV     = $BUSCO_ENV"
echo "    POLISH_ENV    = $POLISH_ENV"
echo "    PILON_FIX     = $PILON_FIX_TYPES"
echo "    PILON_XMX     = $PILON_JAVA_XMX"

step "PRE-FLIGHT CHECKS"
require_command micromamba
require_file "$ILLUMINA_R1" "ILLUMINA_R1"
require_file "$ILLUMINA_R2" "ILLUMINA_R2"
require_file "$DRAFT_ASSEMBLY" "DRAFT_ASSEMBLY"
require_dir "$BUSCO_DB_DIR" "BUSCO_DB_DIR"
require_dir "${BUSCO_DB_DIR}/lineages/${BUSCO_LINEAGE}" "BUSCO_DB_DIR/lineages/${BUSCO_LINEAGE}"

check_env_tool "$BUSCO_ENV" busco || die "Tool 'busco' not found in env '$BUSCO_ENV'"
for tool in bwa samtools bedtools pilon seqkit; do
    check_env_tool "$POLISH_ENV" "$tool" || die "Tool '$tool' not found in env '$POLISH_ENV'"
done

mkdir -p "$DIR_BUSCO_BASE" "$DIR_BUSCO" "$DIR_PILON"
ok "Pre-flight checks passed."

step "STEP 1 — Prepare local working files"
if checkpoint_done "10_prepare_polish_inputs" && [ -f "$ASSEMBLY_COPY" ]; then
    ok "Skipping input preparation (checkpoint found)."
else
    cp -f "$DRAFT_ASSEMBLY" "$ASSEMBLY_COPY"
    mark_checkpoint "10_prepare_polish_inputs"
    ok "Draft assembly copied to $ASSEMBLY_COPY"
fi

step "STEP 2 — BUSCO"
if checkpoint_done "11_busco" && compgen -G "$BUSCO_SUMMARY_GLOB" >/dev/null; then
    ok "Skipping BUSCO (checkpoint found)."
else
    run_busco busco \
        -i "$ASSEMBLY_COPY" \
        --offline \
        --download_path "$BUSCO_DB_DIR" \
        -l "$BUSCO_LINEAGE" \
        -m genome \
        -c "$THREADS" \
        -o busco_out \
        -f \
        --out_path "$DIR_BUSCO_BASE" \
        --opt-out-run-stats
    mark_checkpoint "11_busco"
    ok "BUSCO complete."
fi

SUMMARY="$(find "$DIR_BUSCO" -maxdepth 1 -name 'short_summary*' | head -n 1 || true)"
if [ -n "$SUMMARY" ] && [ -f "$SUMMARY" ]; then
    echo "  BUSCO short summary:"
    grep -E "Complete|Fragmented|Missing|Total" "$SUMMARY" || true
fi

step "STEP 3 — BWA index"
if checkpoint_done "12_bwa_index" && [ -f "${DIR_PILON}/genome_index.bwt" ]; then
    ok "Skipping bwa index (checkpoint found)."
else
    run_polish bwa index -p "${DIR_PILON}/genome_index" "$ASSEMBLY_COPY"
    mark_checkpoint "12_bwa_index"
    ok "BWA index created."
fi

step "STEP 4 — Align Illumina reads and sort BAM"
if checkpoint_done "13_align_sort" && [ -f "$BAM_FILE" ]; then
    ok "Skipping alignment and sort (checkpoint found)."
else
    run_polish bash -lc "
        set -euo pipefail
        bwa mem -t '$THREADS' '${DIR_PILON}/genome_index' '$ILLUMINA_R1' '$ILLUMINA_R2' \
          | samtools sort -@ '$THREADS' -o '$BAM_FILE'
    "
    [ -s "$BAM_FILE" ] || die "Alignment failed: BAM file was not created."
    mark_checkpoint "13_align_sort"
    ok "Read alignment complete."
fi

step "STEP 5 — Index BAM"
if checkpoint_done "14_bam_index" && [ -f "$BAM_INDEX" ]; then
    ok "Skipping BAM indexing (checkpoint found)."
else
    run_polish samtools index "$BAM_FILE"
    mark_checkpoint "14_bam_index"
    ok "BAM index created."
fi

step "STEP 6 — Pilon"
if checkpoint_done "15_pilon" && [ -f "$PILON_FASTA" ]; then
    ok "Skipping Pilon (checkpoint found)."
else
    run_polish bash -lc "
        set -euo pipefail
        export JAVA_TOOL_OPTIONS='-Xmx${PILON_JAVA_XMX}'
        pilon \
            --genome '$ASSEMBLY_COPY' \
            --frags '$BAM_FILE' \
            --output pilon \
            --outdir '$DIR_PILON' \
            --changes \
            --vcf \
            --fix '${PILON_FIX_TYPES}'
    "
    require_file "$PILON_FASTA" "PILON_FASTA"
    mark_checkpoint "15_pilon"
    ok "Pilon polishing complete."
fi

step "STEP 7 — Region annotation"
if [ ! -f "$GTF_FILE" ]; then
    warn "GTF file not found: $GTF_FILE"
    warn "Skipping region-based annotation lookup."
else
    echo "  Gene products in ${REGION_CONTIG}:${REGION_START}-${REGION_END}:"
    echo -e "${REGION_CONTIG}\t${REGION_START}\t${REGION_END}" \
        | run_polish bedtools intersect -a "$GTF_FILE" -b - \
        | grep -oP 'product[^;]*' \
        | sort -u \
        || echo "  (no annotated genes found in this region)"
fi

step "SUMMARY"
echo "  BUSCO report      : $DIR_BUSCO"
echo "  Aligned BAM       : $BAM_FILE"
echo "  Polished assembly : $PILON_FASTA"

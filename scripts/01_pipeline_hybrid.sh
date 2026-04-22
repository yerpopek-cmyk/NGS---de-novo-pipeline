#!/usr/bin/env bash
# =============================================================================
# 01_pipeline_hybrid.sh — Hybrid Genome Assembly & Scaffolding
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

POLISH_ENV="${POLISH_ENV:-polishing_env}"
THREADS="${THREADS:-6}"
GENOME_SIZE="$(trim_cr "${GENOME_SIZE:-5000000}")"
FORCE_RERUN="${FORCE_RERUN:-false}"
CHECKPOINT_DIR="$(abs_path "${CHECKPOINT_DIR:-results/.checkpoints}")"

ILLUMINA_R1="$(abs_path "${ILLUMINA_R1:-}")"
ILLUMINA_R2="$(abs_path "${ILLUMINA_R2:-}")"
NANOPORE="$(abs_path "${NANOPORE:-}")"

DIR_SPADES="${PROJECT_ROOT}/results/spades"
DIR_HYBRID="${DIR_SPADES}/assembly_hybrid"
DIR_ILLUMINA="${DIR_SPADES}/assembly_illumina"
DIR_SIBELIA="${PROJECT_ROOT}/results/sibelia"
DIR_SIBELIA_CMP="${DIR_SIBELIA}/hybrid_vs_illumina"
DIR_SIBELIA_SCF="${DIR_SIBELIA}/hybrid_vs_scaffolded"
DIR_SCAFFOLDING="${PROJECT_ROOT}/results/scaffolding"

FILTERED_HYBRID="${DIR_SIBELIA}/scaffolds.hybrid.fasta"
FILTERED_ILLUMINA="${DIR_SIBELIA}/scaffolds.illumina.fasta"

echo "  config.sh loaded"
echo "    ILLUMINA_R1 = $ILLUMINA_R1"
echo "    ILLUMINA_R2 = $ILLUMINA_R2"
echo "    NANOPORE    = ${NANOPORE:-<not set>}"
echo "    THREADS     = $THREADS"
echo "    POLISH_ENV  = $POLISH_ENV"

step "PRE-FLIGHT CHECKS"
require_command micromamba
require_file "$ILLUMINA_R1" "ILLUMINA_R1"
require_file "$ILLUMINA_R2" "ILLUMINA_R2"

HAVE_NANOPORE=true
if [ -z "${NANOPORE:-}" ]; then
    HAVE_NANOPORE=false
    warn "NANOPORE is not set. Hybrid assembly and long-read scaffolding will be skipped."
elif [ ! -f "$NANOPORE" ]; then
    die "Nanopore file not found: $NANOPORE"
fi

for tool in spades.py seqtk; do
    check_env_tool "$POLISH_ENV" "$tool" || die "Tool '$tool' not found in env '$POLISH_ENV'"
done

if [ "$HAVE_NANOPORE" = true ]; then
    for tool in Sibelia circos longstitch; do
        check_env_tool "$POLISH_ENV" "$tool" || die "Tool '$tool' not found in env '$POLISH_ENV'"
    done
fi

mkdir -p "$DIR_HYBRID" "$DIR_ILLUMINA" "$DIR_SIBELIA" "$DIR_SIBELIA_CMP" "$DIR_SIBELIA_SCF" "$DIR_SCAFFOLDING"
ok "Pre-flight checks passed."

step "STEP 1 — Illumina-only assembly"
if checkpoint_done "01_spades_illumina" && [ -f "${DIR_ILLUMINA}/scaffolds.fasta" ]; then
    ok "Skipping Illumina-only assembly (checkpoint found)."
else
    rm -rf "${DIR_ILLUMINA}/K"* "${DIR_ILLUMINA}/misc" 2>/dev/null || true
    run_polish spades.py \
        -1 "$ILLUMINA_R1" \
        -2 "$ILLUMINA_R2" \
        -o "$DIR_ILLUMINA" \
        -t "$THREADS" \
        --careful
    require_file "${DIR_ILLUMINA}/scaffolds.fasta" "DIR_ILLUMINA/scaffolds.fasta"
    mark_checkpoint "01_spades_illumina"
    ok "Illumina-only assembly complete."
fi

step "STEP 2 — Hybrid assembly"
if [ "$HAVE_NANOPORE" = false ]; then
    warn "Skipping hybrid assembly because NANOPORE is not configured."
else
    if checkpoint_done "02_spades_hybrid" && [ -f "${DIR_HYBRID}/scaffolds.fasta" ]; then
        ok "Skipping hybrid assembly (checkpoint found)."
    else
        rm -rf "${DIR_HYBRID}/K"* "${DIR_HYBRID}/misc" 2>/dev/null || true
        run_polish spades.py \
            -1 "$ILLUMINA_R1" \
            -2 "$ILLUMINA_R2" \
            --nanopore "$NANOPORE" \
            -o "$DIR_HYBRID" \
            -t "$THREADS" \
            --careful
        require_file "${DIR_HYBRID}/scaffolds.fasta" "DIR_HYBRID/scaffolds.fasta"
        mark_checkpoint "02_spades_hybrid"
        ok "Hybrid assembly complete."
    fi
fi

step "STEP 3 — Filter contigs >= 10 kbp"
if checkpoint_done "03_filter_contigs" && [ -f "$FILTERED_ILLUMINA" ] && { [ "$HAVE_NANOPORE" = false ] || [ -f "$FILTERED_HYBRID" ]; }; then
    ok "Skipping contig filtering (checkpoint found)."
else
    run_polish seqtk seq -L 10000 "${DIR_ILLUMINA}/scaffolds.fasta" > "$FILTERED_ILLUMINA"
    if [ "$HAVE_NANOPORE" = true ]; then
        run_polish seqtk seq -L 10000 "${DIR_HYBRID}/scaffolds.fasta" > "$FILTERED_HYBRID"
    fi
    mark_checkpoint "03_filter_contigs"
    ok "Contig filtering complete."
fi

step "STEP 4 — Synteny comparison"
if [ "$HAVE_NANOPORE" = false ]; then
    warn "Skipping Sibelia/Circos comparison because hybrid assembly was not run."
else
    if checkpoint_done "04_sibelia_hybrid_vs_illumina" && [ -f "${DIR_SIBELIA_CMP}/circos/circos.conf" ]; then
        ok "Skipping hybrid vs illumina comparison (checkpoint found)."
    else
        pushd "$DIR_SIBELIA" >/dev/null
        run_polish Sibelia \
            -s fine \
            -o "$DIR_SIBELIA_CMP" \
            "$FILTERED_HYBRID" \
            "$FILTERED_ILLUMINA"
        popd >/dev/null
        mark_checkpoint "04_sibelia_hybrid_vs_illumina"
        ok "Hybrid vs Illumina comparison complete."
    fi

    if checkpoint_done "05_circos_hybrid_vs_illumina" && [ -f "${DIR_SIBELIA_CMP}/circos/circos.svg" -o -f "${DIR_SIBELIA_CMP}/circos/circos.png" ]; then
        ok "Skipping Circos rendering (checkpoint found)."
    else
        run_polish circos --conf "${DIR_SIBELIA_CMP}/circos/circos.conf"
        mark_checkpoint "05_circos_hybrid_vs_illumina"
        ok "Circos rendering complete."
    fi
fi

K_NTLINK=24
W_VALUE=400
SCAFFOLDED="${DIR_SCAFFOLDING}/scaffolds.k${K_NTLINK}.w${W_VALUE}.tigmint-ntLink.longstitch-scaffolds.fa"

step "STEP 5 — LongStitch scaffolding"
if [ "$HAVE_NANOPORE" = false ]; then
    warn "Skipping LongStitch because NANOPORE is not configured."
else
    if checkpoint_done "06_longstitch" && [ -f "$SCAFFOLDED" ]; then
        ok "Skipping LongStitch (checkpoint found)."
    else
        pushd "$DIR_SCAFFOLDING" >/dev/null
        ln -sfn "${DIR_ILLUMINA}/scaffolds.fasta" scaffolds.fa
        ln -sfn "$NANOPORE" nanopore.fq
        run_polish longstitch run \
            draft=scaffolds.fa \
            reads=nanopore \
            G="$GENOME_SIZE" \
            t="$THREADS" \
            gap_fill=True \
            rounds=3 \
            longmap=ont \
            k_ntLink="$K_NTLINK" \
            w="$W_VALUE"
        popd >/dev/null
        require_file "$SCAFFOLDED" "SCAFFOLDED"
        mark_checkpoint "06_longstitch"
        ok "LongStitch scaffolding complete."
    fi
fi

step "SUMMARY"
echo "  Illumina-only assembly : ${DIR_ILLUMINA}/scaffolds.fasta"
if [ "$HAVE_NANOPORE" = true ]; then
    echo "  Hybrid assembly        : ${DIR_HYBRID}/scaffolds.fasta"
    echo "  Hybrid vs Illumina     : ${DIR_SIBELIA_CMP}/circos/"
    echo "  Scaffolded assembly    : ${SCAFFOLDED}"
fi

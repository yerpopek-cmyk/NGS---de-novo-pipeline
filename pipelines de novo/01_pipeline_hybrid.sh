#!/usr/bin/env bash
# =============================================================================
# 01_pipeline_hybrid.sh — Hybrid Genome Assembly & Scaffolding
# =============================================================================
#
# DESCRIPTION:
#   Assembles a genome using two sequencing data types:
#     - Illumina  (short reads: ~150 bp, very accurate)
#     - Nanopore  (long reads: ~10–50 kbp, noisier)
#
#   Then compares the hybrid assembly to an Illumina-only assembly to show
#   students how long reads improve contiguity. Finally, scaffolding with
#   LongStitch joins short contigs into longer, more complete scaffolds.
#
# PIPELINE FLOW:
#   Illumina R1/R2 ──┐
#                    ├──► SPAdes (hybrid)  ──► seqtk filter ──► Sibelia ──► Circos
#   Nanopore reads ──┘
#   Illumina R1/R2 ──────► SPAdes (illumina) ──► seqtk filter ──┘
#                                └──► LongStitch scaffolding ──► Sibelia ──► Circos
#
# PLATFORM: Ubuntu Linux
#
# REQUIREMENTS:
#   conda activate denovo
#   conda install -c bioconda spades seqtk longstitch sibelia circos
#
# USAGE:
#   bash scripts/01_pipeline_hybrid.sh
#
#   Override input files via environment variables:
#   ILLUMINA_R1=my_R1.fastq bash scripts/01_pipeline_hybrid.sh
# =============================================================================

set -euo pipefail
# set -e          : stop on any error
# set -u          : error on undefined variable
# set -o pipefail : fail if any command in a pipe fails (not just the last one)

# =============================================================================
# CONFIGURATION
# =============================================================================
# Paths are loaded automatically from config.sh (written by 00_download_data.sh).
# Run the download script first:
#   bash scripts/00_download_data.sh SRRxxxxxxxx
#
# You can also override any variable before running:
#   THREADS=12 bash scripts/01_pipeline_hybrid.sh

CONFIG="config.sh"
if [ -f "$CONFIG" ]; then
    # source : execute config.sh in the current shell so its variables
    #          become visible here (unlike 'bash config.sh' which runs
    #          in a sub-shell and cannot export back to the parent).
    # shellcheck source=/dev/null
    source "$CONFIG"
    echo "  Loaded settings from config.sh"
    echo "    ILLUMINA_R1 = ${ILLUMINA_R1:-<not set>}"
    echo "    ILLUMINA_R2 = ${ILLUMINA_R2:-<not set>}"
    echo "    NANOPORE    = ${NANOPORE:-<not set — hybrid steps will be skipped>}"
    echo ""
else
    echo ""
    echo "  ERROR: config.sh not found."
    echo "  Run the download script first to generate it:"
    echo "    bash scripts/00_download_data.sh SRRxxxxxxxx"
    echo ""
    echo "  Or set paths manually:"
    echo "    export ILLUMINA_R1=data/raw/SRRxxxxxxxx_1.fastq"
    echo "    export ILLUMINA_R2=data/raw/SRRxxxxxxxx_2.fastq"
    echo "    export NANOPORE=data/raw/your_ont_reads.fastq   # if you have them"
    echo "    bash scripts/01_pipeline_hybrid.sh"
    echo ""
fi

# Fallback defaults in case config.sh did not define a variable
ILLUMINA_R1="${ILLUMINA_R1:-}"
ILLUMINA_R2="${ILLUMINA_R2:-}"
NANOPORE="${NANOPORE:-}"

# Compute resources — override via environment or edit config.sh
THREADS="${THREADS:-6}"
GENOME_SIZE="${GENOME_SIZE:-5000000}"  # bp — 5 Mbp = Vibrio / E. coli default

# --- Output directories ---
DIR_HYBRID="results/spades/assembly_hybrid"
DIR_ILLUMINA="results/spades/assembly_illumina"
DIR_SIBELIA_CMP="results/sibelia/hybrid_vs_illumina"
DIR_SIBELIA_SCF="results/sibelia/hybrid_vs_scaffolded"
DIR_SCAFFOLDING="results/scaffolding"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

BLUE='\033[0;34m'; GREEN='\033[0;32m'; RED='\033[0;31m'; NC='\033[0m'

step()    { echo -e "\n${BLUE}━━━ $1 ━━━${NC}"; }
ok()      { echo -e "${GREEN}✓  $1${NC}"; }
die()     { echo -e "${RED}ERROR: $1${NC}" >&2; exit 1; }

# Check that a required input file exists before we start
require_file() {
    # $1 = file path,  $2 = variable name (for the error message)
    if [ -z "$1" ] || [ ! -f "$1" ]; then
        die "Input file not found: '${1:-<empty>}'\n\
       Fix: run  bash scripts/00_download_data.sh SRRxxxxxxxx\n\
       Or set:  export $2=/path/to/your/file"
    fi
}

# Same check but non-fatal — for optional files like Nanopore reads
require_file_optional() {
    if [ -z "${1:-}" ]; then
        echo -e "  ${YELLOW}SKIP${NC}  $2 is not set — steps requiring it will be skipped."
        return 1   # return 1 = "file missing", caller decides what to do
    elif [ ! -f "$1" ]; then
        die "Variable $2 is set to '$1' but that file does not exist."
    fi
    return 0
}

# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================

step "PRE-FLIGHT CHECKS"

require_file "$ILLUMINA_R1" "ILLUMINA_R1"
require_file "$ILLUMINA_R2" "ILLUMINA_R2"

# Nanopore reads are optional — if absent, hybrid and scaffolding steps are skipped
HAVE_NANOPORE=true
if ! require_file_optional "${NANOPORE:-}" "NANOPORE"; then
    HAVE_NANOPORE=false
    echo "  Illumina-only mode: hybrid assembly and LongStitch will be skipped."
fi

for tool in spades.py seqtk Sibelia circos longstitch; do
    command -v "$tool" &>/dev/null \
        || die "Tool not found: $tool\n       Install with: conda install -c bioconda $tool"
done

mkdir -p "$DIR_HYBRID" "$DIR_ILLUMINA" \
         "$DIR_SIBELIA_CMP" "$DIR_SIBELIA_SCF" \
         "$DIR_SCAFFOLDING" \
         results/sibelia

ok "All inputs and tools are present."

# =============================================================================
# PART 1 — DE NOVO ASSEMBLY
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 1 — Hybrid assembly: Illumina + Nanopore with SPAdes
# ---------------------------------------------------------------------------
if [ "$HAVE_NANOPORE" = true ]; then
    step "STEP 1 — Hybrid assembly (SPAdes: Illumina + Nanopore)"

    echo "  SPAdes builds a De Bruijn Graph from the short, accurate Illumina reads."
    echo "  The Nanopore reads are then used to:"
    echo "    - Bridge contigs separated by long repeats"
    echo "    - Fill gaps that the short reads cannot span"
    echo "  Expected runtime: 15–30 minutes on 6 cores."
    echo ""

    spades.py \
        -1 "$ILLUMINA_R1" \
        -2 "$ILLUMINA_R2" \
        --nanopore "$NANOPORE" \
        -o "$DIR_HYBRID" \
        -t "$THREADS" \
        --careful
        # -1 / -2        : paired-end Illumina reads
        # --nanopore     : long reads — used to resolve repeats in the DBG
        # -o             : output directory
        # -t             : CPU threads
        # --careful      : misassembly correction pass (slower but safer)

    ok "Hybrid assembly done. Output: $DIR_HYBRID/scaffolds.fasta"
else
    step "STEP 1 — SKIPPED (no Nanopore reads)"
    echo "  To run hybrid assembly, add Nanopore reads to config.sh:"
    echo "    NANOPORE=data/raw/your_ont_reads.fastq"
fi

# ---------------------------------------------------------------------------
# STEP 2 — Illumina-only assembly (baseline for comparison)
# ---------------------------------------------------------------------------
step "STEP 2 — Illumina-only assembly (baseline)"

echo "  We assemble the same genome using only short reads."
echo "  Comparing this to the hybrid will clearly show how Nanopore reads"
echo "  reduce fragmentation and resolve repeats."

spades.py \
    -1 "$ILLUMINA_R1" \
    -2 "$ILLUMINA_R2" \
    -o "$DIR_ILLUMINA" \
    -t "$THREADS" \
    --careful
    # Identical to STEP 1 but without --nanopore.
    # The result will typically have a much lower N50 because repeats
    # longer than the Illumina read length cannot be resolved.

ok "Illumina-only assembly done. Output: $DIR_ILLUMINA/scaffolds.fasta"

# =============================================================================
# PART 2 — FILTER SHORT CONTIGS
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 3 — Filter with seqtk (keep only contigs ≥ 10 kbp)
# ---------------------------------------------------------------------------
step "STEP 3 — Filter short contigs with seqtk (min length: 10 kbp)"

echo "  SPAdes outputs thousands of tiny contigs (< 1 kb)."
echo "  These are mostly assembly artefacts, short repeats, or adapter contamination."
echo "  Sibelia's synteny search works best on contiguous, meaningful sequences."
echo "  We keep only contigs ≥ 10,000 bp."

seqtk seq \
    -L 10000 \
    "$DIR_HYBRID/scaffolds.fasta" \
    > results/sibelia/scaffolds.hybrid.fasta
    # seqtk seq  : general-purpose FASTA/FASTQ processor
    # -L 10000   : Minimum Length — discard sequences shorter than 10,000 bp
    # >          : redirect stdout to a file
    #              (seqtk writes the filtered FASTA to stdout by default)

seqtk seq \
    -L 10000 \
    "$DIR_ILLUMINA/scaffolds.fasta" \
    > results/sibelia/scaffolds.illumina.fasta

# Report how many contigs survived the filter
n_hybrid=$(grep -c "^>" results/sibelia/scaffolds.hybrid.fasta   || true)
n_illumina=$(grep -c "^>" results/sibelia/scaffolds.illumina.fasta || true)
# grep -c "^>"  : count lines starting with '>' (FASTA header lines)
# || true       : prevent the script from stopping if grep finds 0 matches

ok "Hybrid assembly:   ${n_hybrid} contigs ≥ 10 kbp"
ok "Illumina-only:     ${n_illumina} contigs ≥ 10 kbp"

# =============================================================================
# PART 3 — ASSEMBLY COMPARISON (Sibelia + Circos)
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 4 — Sibelia: find synteny blocks between the two assemblies
# ---------------------------------------------------------------------------
step "STEP 4 — Sibelia: find synteny blocks (hybrid vs illumina-only)"

echo "  Sibelia identifies collinear blocks: regions present in both assemblies"
echo "  in the same order. This lets us visualise structural differences:"
echo "  rearrangements, missing regions, and fragmentation."

pushd results/sibelia > /dev/null
# pushd : change directory AND remember the previous one (so we can return with popd)
# > /dev/null : suppress the directory name that pushd prints by default

Sibelia \
    -s fine \
    -o "../../$DIR_SIBELIA_CMP" \
    scaffolds.hybrid.fasta \
    scaffolds.illumina.fasta
    # -s fine  : strictness of synteny block detection
    #            fine  = small, precise blocks (best for closely related sequences)
    #            loose = larger blocks with gaps (better for divergent genomes)
    # -o       : output directory
    #            Sibelia writes a ready-made circos.conf inside this directory
    # The two FASTA files are the assemblies to compare.
    # Order matters: the first file will appear as the outer circle in Circos.

popd > /dev/null
# popd : return to the directory we were in before pushd

ok "Sibelia done. Circos config: $DIR_SIBELIA_CMP/circos/circos.conf"

# ---------------------------------------------------------------------------
# STEP 5 — Circos: render the synteny plot
# ---------------------------------------------------------------------------
step "STEP 5 — Circos: render synteny diagram"

echo "  Circos reads the config written by Sibelia and produces an SVG/PNG image."
echo "  The circular plot shows both genomes as arcs; ribbons connect syntenic blocks."

circos \
    --conf "$DIR_SIBELIA_CMP/circos/circos.conf"
    # --conf : path to the Circos configuration file
    #          Sibelia generated this file automatically with correct data paths,
    #          colours, and chromosome sizes.  You can open and edit it to change
    #          colours, font sizes, or other aesthetics.

ok "Circos plot saved in: $DIR_SIBELIA_CMP/circos/"

# =============================================================================
# PART 4 — SCAFFOLDING (LongStitch) — requires Nanopore reads
# =============================================================================

K_NTLINK=24
W_VALUE=400
SCAFFOLDED="scaffolds.k${K_NTLINK}.w${W_VALUE}.tigmint-ntLink.longstitch-scaffolds.fa"

if [ "$HAVE_NANOPORE" = true ]; then

    # ---------------------------------------------------------------------------
    # STEP 6 — LongStitch: scaffold the Illumina-only assembly using Nanopore reads
    # ---------------------------------------------------------------------------
    step "STEP 6 — LongStitch: scaffolding the Illumina-only assembly with Nanopore reads"

    echo "  Problem: the Illumina-only assembly is fragmented — many short contigs"
    echo "           separated by unresolved repeats."
    echo "  Solution: LongStitch uses Nanopore reads as bridges."
    echo "    ntLink selects representative k-mers (minimizers) from reads and contigs."
    echo "    If two contigs share long reads spanning their ends, they are joined."

    pushd "$DIR_SCAFFOLDING" > /dev/null

    DRAFT_NAME="scaffolds.fa"
    READS_PREFIX="nanopore"

    ln -sf "../../$DIR_ILLUMINA/scaffolds.fasta" "$DRAFT_NAME"
    ln -sf "../../$NANOPORE" "${READS_PREFIX}.fq"
    # ln -sf : symbolic link (-s), overwrite if exists (-f)
    # LongStitch expects the reads file as ${READS_PREFIX}.fq — hence the link name

    longstitch run \
        draft="$DRAFT_NAME" \
        reads="$READS_PREFIX" \
        G="$GENOME_SIZE" \
        t="$THREADS" \
        gap_fill=True \
        rounds=3 \
        longmap=ont \
        k_ntLink="$K_NTLINK" \
        w="$W_VALUE"
        # draft=        : input assembly (Illumina-only, fragmented)
        # reads=        : Nanopore reads prefix (LongStitch adds .fq)
        # G=            : genome size estimate — used to assess read coverage
        # t=            : CPU threads
        # gap_fill=True : try to fill N-gaps with long reads after scaffolding
        # rounds=3      : iterations of ntLink (more rounds = more connections)
        # longmap=ont   : Nanopore mode (alternatives: pb, hifi)
        # k_ntLink=24   : k-mer size for minimizer index
        # w=400         : window size for minimizer sampling

    popd > /dev/null
    ok "Scaffolding done. Output: $DIR_SCAFFOLDING/$SCAFFOLDED"

    # ---------------------------------------------------------------------------
    # STEP 7 — Sibelia + Circos: hybrid assembly vs scaffolded Illumina-only
    # ---------------------------------------------------------------------------
    step "STEP 7 — Final comparison: hybrid vs scaffolded Illumina-only assembly"
    echo "  Strategy A (hybrid SPAdes) vs Strategy B (Illumina assembly + LongStitch)."
    echo "  A good scaffolding result should look nearly identical to the hybrid."

    pushd "$DIR_SCAFFOLDING" > /dev/null

    Sibelia \
        -s fine \
        -o "../../$DIR_SIBELIA_SCF" \
        "../../results/sibelia/scaffolds.hybrid.fasta" \
        "$SCAFFOLDED"

    circos --conf "../../$DIR_SIBELIA_SCF/circos/circos.conf"

    popd > /dev/null
    ok "Final Circos plot saved in: $DIR_SIBELIA_SCF/circos/"

else
    step "STEP 6 + 7 — SKIPPED (no Nanopore reads)"
    echo "  Add Nanopore reads to config.sh to enable scaffolding:"
    echo "    NANOPORE=data/raw/your_ont_reads.fastq"
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  PIPELINE 1 COMPLETE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Assemblies:"
echo "    Hybrid  (Illumina + Nanopore) : $DIR_HYBRID/scaffolds.fasta"
echo "    Illumina-only                 : $DIR_ILLUMINA/scaffolds.fasta"
echo "    Scaffolded Illumina-only      : $DIR_SCAFFOLDING/$SCAFFOLDED"
echo ""
echo "  Synteny plots:"
echo "    Hybrid vs Illumina-only       : $DIR_SIBELIA_CMP/circos/"
echo "    Hybrid vs Scaffolded          : $DIR_SIBELIA_SCF/circos/"
echo ""
echo "  Next step: run Pipeline 2 to polish and assess quality"
echo "    bash scripts/02_pipeline_polish.sh"
echo ""

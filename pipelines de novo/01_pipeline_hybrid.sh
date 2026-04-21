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
# ✏️  CONFIGURATION — EDIT THESE TO MATCH YOUR DATA
# =============================================================================

# --- Input files ---
# Illumina paired-end reads.  If the variable is already set in your shell,
# that value is used; otherwise the default path shown here is used.
ILLUMINA_R1="${ILLUMINA_R1:-data/raw/illumina_R1.fastq}"
ILLUMINA_R2="${ILLUMINA_R2:-data/raw/illumina_R2.fastq}"

# Oxford Nanopore long reads
NANOPORE="${NANOPORE:-data/raw/nanopore.fastq}"

# --- Assembly parameters ---
THREADS=6               # CPU threads to use
GENOME_SIZE=5000000     # Approximate genome size in bp (5 Mbp for Vibrio cholerae)
                        # Used by LongStitch to estimate coverage

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
    [ -f "$1" ] || die "Input file not found: $1\n       Set the variable $2 to the correct path."
}

# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================

step "PRE-FLIGHT CHECKS"

require_file "$ILLUMINA_R1" "ILLUMINA_R1"
require_file "$ILLUMINA_R2" "ILLUMINA_R2"
require_file "$NANOPORE"    "NANOPORE"

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
    #                  R1 = forward (5'→3' on the top strand)
    #                  R2 = reverse (5'→3' on the bottom strand)
    #                  SPAdes uses the pair relationship to improve assembly:
    #                  it knows the two reads came from the same DNA fragment.
    # --nanopore     : long reads from Oxford Nanopore Technology
    #                  SPAdes aligns these to the DBG to resolve ambiguous paths
    #                  and connect disconnected subgraphs.
    # -o             : output directory (SPAdes writes many intermediate files here)
    # -t             : number of CPU threads
    # --careful      : enables an extra "misassembly correction" pass —
    #                  SPAdes re-maps reads to the draft assembly and fixes
    #                  any joins that are not supported by the read evidence.
    #                  Slower, but produces fewer chimeric contigs.

ok "Hybrid assembly done. Output: $DIR_HYBRID/scaffolds.fasta"

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
# PART 4 — SCAFFOLDING (LongStitch)
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 6 — LongStitch: scaffold the Illumina-only assembly using Nanopore reads
# ---------------------------------------------------------------------------
step "STEP 6 — LongStitch: scaffolding the Illumina-only assembly with Nanopore reads"

echo "  Problem: the Illumina-only assembly is fragmented — many short contigs"
echo "           separated by unresolved repeats."
echo "  Solution: LongStitch uses Nanopore reads as bridges."
echo "    Internally it runs ntLink, which selects representative k-mers"
echo "    (called minimizers) from both reads and contigs."
echo "    If two contigs share long reads that span their ends, they are joined."

# We work inside the scaffolding output directory to keep LongStitch files tidy.
pushd "$DIR_SCAFFOLDING" > /dev/null

# Symlinks avoid copying large files into the working directory.
# ln -sf : create a symbolic link (-s), overwriting any existing link (-f)
DRAFT_NAME="scaffolds.fa"
READS_PREFIX="nanopore"

ln -sf "../../$DIR_ILLUMINA/scaffolds.fasta" "$DRAFT_NAME"
# Note: LongStitch infers the reads file name as ${reads}.fq, so the link
# must be named exactly ${READS_PREFIX}.fq
ln -sf "../../$NANOPORE" "${READS_PREFIX}.fq"

K_NTLINK=24    # k-mer size for ntLink minimizer index
W_VALUE=400    # window size for minimizer sampling
               # One minimizer is selected from every W consecutive bases.
               # Larger W → fewer minimizers → faster but may miss connections.

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
    # draft=        : the assembly FASTA to scaffold (basename only, with .fa extension)
    # reads=        : prefix of the long-reads file (LongStitch appends .fq automatically)
    # G=            : estimated genome size in bp
    #                 Used to gauge whether coverage is sufficient for scaffolding.
    # t=            : threads
    # gap_fill=True : after scaffolding, try to fill Ns (gaps) using long reads
    # rounds=3      : run 3 iterations of ntLink
    #                 Each round re-aligns reads to the improved assembly,
    #                 potentially linking contigs that were not linkable in round 1.
    # longmap=ont   : read type: ont = Oxford Nanopore Technology
    #                 Alternatives: pb (PacBio CLR), hifi (PacBio HiFi)
    # k_ntLink=24   : k-mer size for the minimizer index
    #                 Smaller k → more connections, but more false positives
    #                 Larger k → fewer connections, but higher specificity
    # w=400         : minimizer window width
    #                 From every 400 bp window, one representative k-mer is chosen.

# The scaffolded assembly is named by LongStitch based on the parameters:
SCAFFOLDED="scaffolds.k${K_NTLINK}.w${W_VALUE}.tigmint-ntLink.longstitch-scaffolds.fa"
# ${var%suffix} : remove 'suffix' from the end of $var
# Result example: scaffolds.k24.w400.tigmint-ntLink.longstitch-scaffolds.fa

popd > /dev/null

ok "Scaffolding done. Output: $DIR_SCAFFOLDING/$SCAFFOLDED"

# =============================================================================
# PART 5 — FINAL COMPARISON: hybrid vs scaffolded
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 7 — Sibelia + Circos: hybrid assembly vs scaffolded Illumina-only
# ---------------------------------------------------------------------------
step "STEP 7 — Final comparison: hybrid vs scaffolded Illumina-only assembly"

echo "  Now we compare two strategies for using long reads:"
echo "    Strategy A: Feed them directly to SPAdes during assembly (hybrid)."
echo "    Strategy B: Assemble short reads first, then scaffold with long reads."
echo "  A good scaffolding result should look nearly identical to the hybrid."

pushd "$DIR_SCAFFOLDING" > /dev/null

Sibelia \
    -s fine \
    -o "../../$DIR_SIBELIA_SCF" \
    "../../results/sibelia/scaffolds.hybrid.fasta" \
    "$SCAFFOLDED"

circos \
    --conf "../../$DIR_SIBELIA_SCF/circos/circos.conf"

popd > /dev/null

ok "Final Circos plot saved in: $DIR_SIBELIA_SCF/circos/"

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

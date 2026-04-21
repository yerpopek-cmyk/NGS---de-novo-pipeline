#!/usr/bin/env bash
# =============================================================================
# 02_pipeline_polish.sh — Assembly Polishing & Quality Assessment
# =============================================================================
#
# DESCRIPTION:
#   Takes a draft genome assembly (e.g. produced by Flye from Nanopore reads)
#   and improves it through four stages:
#     1. BUSCO      — assess gene completeness before polishing
#     2. BWA + samtools — align accurate Illumina reads to the draft
#     3. Pilon      — correct nucleotide-level errors using the alignments
#     4. bedtools   — find gene products in a specified genomic region
#
# WHY POLISH?
#   Long-read assemblers (Flye, Canu) produce contiguous assemblies that span
#   repeats, but Nanopore reads carry ~0.1–1% random errors.
#   That means roughly 1 wrong base per 100–1000 bp — enough to break
#   open reading frames and corrupt protein sequences.
#   Pilon fixes these errors by treating accurate Illumina reads as a reference.
#
# PIPELINE FLOW:
#   Draft assembly ──► BUSCO ──────────────────────────► completeness report
#          │
#          └──► bwa index ──► bwa mem ──► samtools sort/index
#                                                   │
#                                                Pilon ──► polished assembly
#                                                   │
#                                       annotation + bedtools ──► gene products
#
# PLATFORM: Ubuntu Linux
#
# REQUIREMENTS:
#   conda activate denovo
#   conda install -c bioconda busco bwa samtools bedtools pilon
#
# USAGE:
#   bash scripts/02_pipeline_polish.sh
#
#   Override variables at runtime:
#   DRAFT_ASSEMBLY=my_flye/assembly.fasta bash scripts/02_pipeline_polish.sh
# =============================================================================

set -euo pipefail

# =============================================================================
# ✏️  CONFIGURATION — EDIT THESE TO MATCH YOUR DATA
# =============================================================================

# --- Input files ---

# Draft assembly to be polished.
# This is typically the output of Flye (assembly.fasta) or any other assembler.
DRAFT_ASSEMBLY="${DRAFT_ASSEMBLY:-results/flye/assembly.fasta}"

# Illumina paired-end reads used for polishing
ILLUMINA_R1="${ILLUMINA_R1:-data/raw/illumina_R1.fastq}"
ILLUMINA_R2="${ILLUMINA_R2:-data/raw/illumina_R2.fastq}"

# GTF annotation file for the organism.
# Download from NCBI RefSeq:
#   https://www.ncbi.nlm.nih.gov/genome/
# Example for Vibrio cholerae O1 N16961:
#   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/745/\
#        GCF_000006745.1_VibrioCholerae_O1_biovar_eltor_str_N16961/\
#        GCF_000006745.1_VibrioCholerae_O1_biovar_eltor_str_N16961_genomic.gtf.gz
#   gunzip *.gtf.gz
GTF_FILE="${GTF_FILE:-data/reference/annotation.gtf}"

# --- BUSCO settings ---

# Directory containing downloaded BUSCO lineage databases.
# Download a lineage database with:
#   busco --download vibrio_odb12         (Vibrio)
#   busco --download bacteria_odb10       (all bacteria)
#   busco --download enterobacteriales_odb11  (Enterobacteriaceae)
# Then point this variable to the folder that contains the downloaded lineages.
BUSCO_DB_DIR="${BUSCO_DB_DIR:-data/db/busco_downloads}"

# Which lineage to use for BUSCO gene completeness assessment.
# Choose the lineage that best matches your organism:
#   vibrio_odb12            — genus Vibrio
#   bacteria_odb10          — all bacteria (broad, less sensitive)
#   enterobacteriales_odb11 — Enterobacteriaceae (E. coli, Salmonella, ...)
#   fungi_odb10             — fungi
#   vertebrata_odb10        — vertebrates
BUSCO_LINEAGE="${BUSCO_LINEAGE:-vibrio_odb12}"

# --- Gene region search ---
# Contig/chromosome name in the assembly
REGION_CONTIG="${REGION_CONTIG:-NZ_LT906615.1}"
# Coordinates (0-based, half-open interval — standard BED format)
REGION_START="${REGION_START:-300000}"
REGION_END="${REGION_END:-440000}"

# --- Compute resources ---
THREADS=6

# --- Output directories ---
DIR_BUSCO="results/busco/busco_out"
DIR_PILON="results/pilon"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

BLUE='\033[0;34m'; GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; NC='\033[0m'

step()      { echo -e "\n${BLUE}━━━ $1 ━━━${NC}"; }
ok()        { echo -e "${GREEN}✓  $1${NC}"; }
warn()      { echo -e "${YELLOW}⚠  $1${NC}"; }
die()       { echo -e "${RED}ERROR: $1${NC}" >&2; exit 1; }

require_file() { [ -f "$1" ] || die "File not found: $1 (variable: $2)"; }

# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================

step "PRE-FLIGHT CHECKS"

require_file "$DRAFT_ASSEMBLY" "DRAFT_ASSEMBLY"
require_file "$ILLUMINA_R1"    "ILLUMINA_R1"
require_file "$ILLUMINA_R2"    "ILLUMINA_R2"

for tool in busco bwa samtools pilon bedtools; do
    command -v "$tool" &>/dev/null \
        || warn "Tool not found: $tool — the corresponding step will fail."
done

mkdir -p "$DIR_BUSCO" "$DIR_PILON" results/busco

ok "Pre-flight checks passed."

# =============================================================================
# PART 1 — ASSEMBLY COMPLETENESS (BUSCO)
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 1 — Symlink to the BUSCO database (avoid copying gigabytes)
# ---------------------------------------------------------------------------
step "STEP 1 — Preparing BUSCO database"

echo "  We create a symbolic link (shortcut) to the BUSCO lineage database."
echo "  A symlink is a pointer to the original directory — no data is copied."

if [ -d "$BUSCO_DB_DIR" ]; then
    ln -sf "$(realpath "$BUSCO_DB_DIR")" results/busco/busco_downloads
    # ln -sf : create a symbolic link (-s), overwriting any existing one (-f)
    # realpath : converts a relative path to an absolute one
    #            BUSCO needs an absolute path to locate the database correctly
    ok "Symlink created: results/busco/busco_downloads → $BUSCO_DB_DIR"
else
    warn "BUSCO database not found at: $BUSCO_DB_DIR"
    warn "Download it with: busco --download $BUSCO_LINEAGE"
    warn "Then set BUSCO_DB_DIR to the download directory and re-run."
fi

# ---------------------------------------------------------------------------
# STEP 2 — Run BUSCO
# ---------------------------------------------------------------------------
step "STEP 2 — BUSCO: assess assembly gene completeness"

echo "  BUSCO checks how many conserved, single-copy orthologs are present"
echo "  in your assembly.  These genes must be present in all members of the"
echo "  lineage ($BUSCO_LINEAGE) — so a high BUSCO score means a complete assembly."
echo ""
echo "  BUSCO workflow:"
echo "    1. Loads HMM profiles for ~1,000–3,000 genes from the lineage database."
echo "    2. Searches your assembly for each gene using tBlastn / Augustus."
echo "    3. Classifies each gene as Complete, Fragmented, or Missing."

busco \
    -i "$DRAFT_ASSEMBLY" \
    --offline \
    -l "$BUSCO_LINEAGE" \
    -m geno \
    -c "$THREADS" \
    -o "$DIR_BUSCO" \
    -f \
    --opt-out-run-stats
    # -i                  : input file — your assembly in FASTA format
    # --offline           : do NOT download databases from the internet
    #                       use the local symlink we created in STEP 1
    # -l "$BUSCO_LINEAGE" : which lineage database to use
    #                       must match a folder inside busco_downloads/
    # -m geno             : run mode
    #                       geno = genome mode (whole genome FASTA as input)
    #                       tran = transcriptome mode
    #                       prot = protein mode
    # -c "$THREADS"       : number of parallel threads
    # -o "$DIR_BUSCO"     : name of the output directory
    # -f                  : force — overwrite the output directory if it exists
    # --opt-out-run-stats : do not send anonymous usage statistics to BUSCO servers

ok "BUSCO complete. Results: $DIR_BUSCO/"

# Print the short summary immediately so the student can see scores in the terminal
SUMMARY=$(find "$DIR_BUSCO" -name "short_summary*.txt" | head -1)
if [ -f "$SUMMARY" ]; then
    echo ""
    echo "  ── BUSCO short summary ────────────────────────────────"
    grep -E "Complete|Fragmented|Missing|Total" "$SUMMARY" || true
    echo "  ────────────────────────────────────────────────────────"
    echo ""
    echo "  Interpretation:"
    echo "    ≥ 95% Complete  →  excellent assembly"
    echo "    85–95% Complete →  acceptable; check coverage and lineage"
    echo "    < 85% Complete  →  poor assembly or wrong lineage selected"
fi

# =============================================================================
# PART 2 — POLISH WITH PILON
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 3 — bwa index: build a searchable index of the draft assembly
# ---------------------------------------------------------------------------
step "STEP 3 — bwa index: indexing the draft assembly"

echo "  Before aligning reads, BWA needs to build an index of the genome."
echo "  Think of it as a book's index: instead of reading every page to find"
echo "  a word, you jump directly to the right position."
echo ""
echo "  Index files created (all share the prefix 'genome_index'):"
echo "    .bwt  — Burrows-Wheeler Transform of the genome (the core index)"
echo "    .sa   — suffix array for O(1) position lookups"
echo "    .pac  — 2-bit packed representation of the genome"
echo "    .amb  — positions of ambiguous bases (N, R, Y, ...)"
echo "    .ann  — contig names and lengths (used by samtools later)"

pushd "$DIR_PILON" > /dev/null
# pushd : move into the polishing directory; remember where we came from

bwa index \
    -p genome_index \
    "../../$DRAFT_ASSEMBLY"
    # -p genome_index           : prefix for all index files
    #                             every index file will be named genome_index.XXX
    # "../../$DRAFT_ASSEMBLY"   : path to the genome to index
    #                             ../../ is needed because we are inside results/pilon/

ok "Index created: $DIR_PILON/genome_index.*"

# ---------------------------------------------------------------------------
# STEP 4 — bwa mem: align Illumina reads to the draft assembly
# ---------------------------------------------------------------------------
step "STEP 4 — bwa mem: aligning Illumina reads to the draft assembly"

echo "  BWA MEM (Maximal Exact Match) is the standard algorithm for aligning"
echo "  short Illumina reads to a reference genome."
echo ""
echo "  How BWA MEM works:"
echo "    1. Finds exact seed matches (short substrings) between reads and genome."
echo "    2. Extends seeds into full alignments using Smith-Waterman."
echo "    3. Outputs alignments in SAM format."
echo ""
echo "  We pipe directly into 'samtools sort' to avoid writing a large SAM file:"
echo "    SAM  = plain text, human-readable, but 20–100 GB for a typical dataset"
echo "    BAM  = binary compressed SAM; Pilon requires this format"
echo "  Sorting by coordinate is mandatory: Pilon processes the genome"
echo "  position by position and expects reads in order."

bwa mem \
    -t "$THREADS" \
    genome_index \
    "../../$ILLUMINA_R1" \
    "../../$ILLUMINA_R2" \
    | samtools sort \
        -@ "$THREADS" \
        -o aligned.bam
    # bwa mem:
    #   -t "$THREADS"    : number of CPU threads
    #   genome_index     : index prefix (bwa finds genome_index.bwt, .sa, etc.)
    #   R1 and R2        : paired-end reads
    #                      BWA MEM knows they are paired: it expects both reads
    #                      to map to the same contig ~200–1000 bp apart in
    #                      opposite orientations (FR orientation).
    #
    # | (pipe)           : feed stdout of bwa mem directly into samtools sort
    #                      avoids writing an intermediate SAM file (~20 GB)
    #
    # samtools sort:
    #   -@ "$THREADS"    : threads for sorting (note: samtools uses -@, not -t)
    #   -o aligned.bam   : output filename; .bam extension triggers BAM compression

ok "Alignment complete: $DIR_PILON/aligned.bam"

# ---------------------------------------------------------------------------
# STEP 5 — samtools index: create a BAM index for random access
# ---------------------------------------------------------------------------
step "STEP 5 — samtools index: indexing the BAM file"

echo "  Pilon constantly jumps to different positions in the genome:"
echo "    'Give me all reads covering position 500,000'"
echo "    'Now position 750,000'  ..."
echo ""
echo "  Without an index, samtools would have to scan the entire BAM file"
echo "  from the start every time — extremely slow."
echo "  The index (.bai) stores byte offsets: 'chromosome X starts at byte Y',"
echo "  enabling O(1) random access to any genomic region."

samtools index aligned.bam
# Creates: aligned.bam.bai
# This small file (typically a few MB) points into the sorted BAM.

ok "BAM index created: $DIR_PILON/aligned.bam.bai"

# ---------------------------------------------------------------------------
# STEP 6 — Pilon: polish the assembly
# ---------------------------------------------------------------------------
step "STEP 6 — Pilon: correcting nucleotide errors in the draft assembly"

echo "  Pilon uses the aligned Illumina reads to fix errors in the draft:"
echo ""
echo "  For each position in the genome, Pilon:"
echo "    1. Counts how many reads carry each nucleotide (A, T, G, C)."
echo "    2. If the majority disagrees with the draft, it corrects the draft."
echo "    3. Checks for small insertions and deletions (indels) too."
echo ""
echo "  Example at position 1,000:"
echo "    Draft:   ...A..."
echo "    Read 1:  ...A...   Read 2: ...A...   Read 3: ...G...   Read 4: ...A..."
echo "    Vote: A=3, G=1  →  Pilon keeps A (the G is a sequencing error)"
echo ""
echo "  Expected improvement: ~0.1–1% error rate  →  ~0.001% after polishing"

pilon \
    --genome "../../$DRAFT_ASSEMBLY" \
    --frags aligned.bam
    # --genome      : the draft assembly to correct
    # --frags       : BAM file of paired-end (fragment library) alignments
    #                 Pilon uses the fact that reads are paired to detect
    #                 misassembled regions with unusual insert-size distributions.
    #
    # Other useful flags you might add:
    #   --output prefix     : rename output files (default: pilon)
    #   --changes           : write all corrections to pilon.changes
    #   --vcf               : write corrections as a VCF file
    #   --fix all           : fix all error types (default behaviour)
    #   --diploid           : diploid mode (keeps heterozygous variants as IUPAC bases)
    #   --nanopore bam.bam  : also include Nanopore alignments for additional evidence

ok "Polishing complete!"
ok "Polished assembly: $DIR_PILON/pilon.fasta"

popd > /dev/null

# =============================================================================
# PART 3 — FIND GENE PRODUCTS IN A GENOMIC REGION
# =============================================================================

# ---------------------------------------------------------------------------
# STEP 7 — bedtools intersect + grep: extract gene products from a region
# ---------------------------------------------------------------------------
step "STEP 7 — Find gene products in region ${REGION_CONTIG}:${REGION_START}–${REGION_END}"

echo "  A common task in genomics: 'Which genes are in this region?'"
echo ""
echo "  Tools used:"
echo "    echo -e '...'          : create a single BED-format line defining the region"
echo "    bedtools intersect     : find GTF entries that overlap the BED region"
echo "    grep -oP 'product ...' : extract just the gene product names using regex"
echo "    uniq                   : collapse consecutive duplicate lines"
echo ""
echo "  BED format (tab-separated):"
echo "    CHROM  START  END"
echo "    ${REGION_CONTIG}  ${REGION_START}  ${REGION_END}"
echo "    (0-based, half-open: position START is included, END is NOT)"
echo ""

if [ ! -f "$GTF_FILE" ]; then
    warn "GTF file not found: $GTF_FILE"
    warn "Download an annotation for your organism and set GTF_FILE."
    echo ""
    echo "  Example — download Vibrio cholerae O1 N16961 annotation:"
    echo "  wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/745/GCF_000006745.1_VibrioCholerae_O1_biovar_eltor_str_N16961/GCF_000006745.1_VibrioCholerae_O1_biovar_eltor_str_N16961_genomic.gtf.gz'"
    echo "  gunzip *.gtf.gz"
    echo "  mv *.gtf data/reference/annotation.gtf"
else
    echo "  Gene products found in ${REGION_CONTIG}:${REGION_START}–${REGION_END}:"
    echo "  ─────────────────────────────────────────────────────────────────"

    echo -e "${REGION_CONTIG}\t${REGION_START}\t${REGION_END}" \
        | bedtools intersect \
            -a "$GTF_FILE" \
            -b - \
        | grep -oP "product [^;]*" \
        | uniq \
        || echo "  (no annotated genes found in this region)"

    # Full breakdown of the command:
    #
    # echo -e "${REGION_CONTIG}\t${REGION_START}\t${REGION_END}"
    #   -e   : interpret escape sequences (\t = tab character)
    #   \t   : the BED format requires tab-separated columns — NOT spaces
    #   Output: "NZ_LT906615.1<TAB>300000<TAB>440000"
    #
    # | bedtools intersect -a "$GTF_FILE" -b -
    #   bedtools intersect : report features from -a that overlap features from -b
    #   -a "$GTF_FILE"     : the annotation file (GTF/GFF3 — all annotated genes)
    #   -b -               : read the second BED from stdin ('-' = standard input)
    #                        This is the single-line BED from echo above.
    #   Output: every GTF line whose coordinates overlap [300000, 440000)
    #           on contig NZ_LT906615.1
    #
    # | grep -oP "product [^;]*"
    #   -o   : print only the matching part of the line (not the whole line)
    #   -P   : use Perl-compatible regular expressions (extended regex)
    #   Pattern: "product [^;]*"
    #     product   : literal text that appears in GTF attribute fields
    #     [^;]*     : any characters EXCEPT ';', zero or more times
    #                 This stops the match at the first semicolon, which
    #                 separates attributes in GTF.
    #   Example GTF attribute string:
    #     gene_id "VC0001"; product=cholera toxin subunit A; db_xref "UniProt:P01555"
    #   grep extracts:  product=cholera toxin subunit A
    #
    # | uniq
    #   Removes consecutive duplicate lines.
    #   Why needed: each gene has multiple GTF entries (gene, mRNA, CDS, exon, ...)
    #   all carrying the same product attribute.  uniq collapses them to one line.
    #   NOTE: uniq only removes ADJACENT duplicates.  If you need to remove all
    #   duplicates regardless of order, use: sort | uniq

    echo "  ─────────────────────────────────────────────────────────────────"
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  PIPELINE 2 COMPLETE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Outputs:"
echo "    BUSCO report      : $DIR_BUSCO/"
echo "    BWA alignment     : $DIR_PILON/aligned.bam"
echo "    Polished assembly : $DIR_PILON/pilon.fasta"
echo ""
echo "  Useful next steps:"
echo "    # View BUSCO summary:"
echo "    cat \$(find $DIR_BUSCO -name 'short_summary*.txt')"
echo ""
echo "    # Compare contig count before and after polishing:"
echo "    grep -c '^>' $DRAFT_ASSEMBLY"
echo "    grep -c '^>' $DIR_PILON/pilon.fasta"
echo ""
echo "    # Quick assembly statistics (requires seqkit):"
echo "    seqkit stats -a $DIR_PILON/pilon.fasta"
echo ""

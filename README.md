# 🧬 De Novo Genome Assembly — NGS Course

> A hands-on bioinformatics course repository for students learning de novo genome assembly from scratch.  
> Built and tested on **Ubuntu Linux**. All pipelines use command-line tools via conda.

---

## 📚 Table of Contents

- [Why is assembly even a problem?](#-why-is-assembly-even-a-problem)
- [The math behind assembly](#-the-math-behind-assembly)
- [OLC vs DBG algorithms](#-two-algorithms-olc-vs-dbg)
- [Quality metrics](#-quality-metrics-n50-and-busco)
- [Repository structure](#-repository-structure)
- [Quick start](#-quick-start)
- [Pipeline 1: Hybrid assembly](#pipeline-1-hybrid-assembly--scaffolding)
- [Pipeline 2: Polishing & QC](#pipeline-2-long-reads-polishing--quality-assessment)
- [Tools used](#-tools-used)

---

## 🔬 Why Is Assembly Even a Problem?

### The fundamental limitation: read length

Sequencers cannot read an entire chromosome in one go. They produce short fragments called **reads**:

| Platform | Read length | Accuracy |
|----------|-------------|----------|
| Illumina | 100–300 bp | ~99.9% |
| PacBio HiFi | 10–25 kbp | ~99.9% |
| Oxford Nanopore | 10–100+ kbp | ~85–98% |

**Analogy:** Imagine shredding *War and Peace* into a million strips of paper, each with 5–10 words.  
Your task: reconstruct the book **without having the original**.

### Three enemies of every assembler

1. **Sequencing errors** — a read contains a wrong nucleotide.  
   → Any algorithm requiring exact matches will break.

2. **Repeats** — the same sequence appears many times in the genome.  
   → If a read sits entirely inside a repeat, it's **impossible to know which copy it came from**.

3. **Heterozygosity** — diploid organisms have two copies of each chromosome.  
   → The assembler may merge them into an "averaged" sequence, or (ideally) resolve them into two separate haplotypes.

---

## 📐 The Math Behind Assembly

### Coverage formula

```
C = (N × L) / G
```

- **C** — average coverage (how many times each nucleotide was sequenced)
- **N** — number of reads
- **L** — read length
- **G** — genome size

**Example:** 5 million reads × 150 bp / 5 Mbp (Vibrio genome) = **150×**

| Coverage | What happens |
|----------|-------------|
| < 5× | Too many gaps — assembly falls apart |
| 30–60× | Gold standard for Illumina |
| > 100× | Excessive — errors create noise in the graph |

### k-mers — the foundation of DBG assembly

A **k-mer** is a substring of fixed length `k`.

```
Read: ATGCATG   (k = 3)
k-mers: ATG  TGC  GCA  CAT  ATG
```

**Why?** Comparing all reads pairwise is O(N²) — impossible for billions of reads.  
Decomposing reads into k-mers turns the problem into **dictionary lookup** — O(N).

**k-mer spectrum:** a histogram of k-mer frequencies shows two peaks:
- 🔴 **Error peak** (frequency 1–5): random k-mers from sequencing errors
- 🟢 **True genome peak** (frequency ≈ coverage depth): real genomic k-mers

SPAdes uses this to discard noisy k-mers *before* building the assembly graph.

---

## ⚙️ Two Algorithms: OLC vs DBG

### OLC (Overlap–Layout–Consensus) — assembly like a jigsaw puzzle

1. **Overlap**: find which read's tail matches another read's head
2. **Layout**: build a read graph, find the path (Hamiltonian — NP-hard, but the graph is nearly linear in practice)
3. **Consensus**: stack reads, vote on the correct nucleotide at each position

**Used for:** PacBio HiFi, Nanopore (long reads with unique overlaps)  
**Tools:** Canu, Flye, Hifiasm

### DBG (De Bruijn Graph) — the k-mer magic

1. All reads are cut into k-mers
2. Identical k-mers → single node in the graph
3. Assembly = finding an **Eulerian path** (traverses every edge exactly once — solvable in linear time!)

**The key rule:**
- k **shorter** than a repeat → the graph has a bubble → **fragmentation**
- k **longer** than a repeat → the k-mer spanning the repeat boundary is unique → **repeat resolved**

**Used for:** Illumina (short, accurate reads)  
**Tools:** SPAdes, MEGAHIT, Velvet

### When to use which

```
Illumina only (100–300 bp)        →  SPAdes (DBG)
Nanopore / PacBio only            →  Flye / Canu (OLC)
Illumina + Nanopore (hybrid)      →  SPAdes --nanopore  OR  Flye + Pilon
```

---

## 📊 Quality Metrics: N50 and BUSCO

### N50 — contiguity

**Algorithm:**
1. Sort all contigs by length (descending)
2. Accumulate their lengths
3. N50 = the length of the contig at which the running total exceeds **50% of the total assembly size**

**Interpretation:** "50% of the genome is contained in contigs of length ≥ N50"

```
Good bacterial assembly:    N50 > 100 kb
Great assembly:             N50 ≈ chromosome length
T2T (telomere-to-telomere): no N-gaps at all
```

> ⚠️ N50 can be inflated by joining garbage sequences. Always use it together with BUSCO.

### BUSCO — gene completeness

**Core idea:** All organisms in a lineage share a set of **single-copy conserved orthologs** (genes present in every member, in exactly one copy).  
BUSCO checks how many of those your assembly contains.

```
C:98.0%[S:97.5%, D:0.5%], F:1.0%, M:1.0%, n:1200

C = Complete (found in full)
S = Single-copy (one copy — expected)
D = Duplicated (multiple copies — possible assembly artifact)
F = Fragmented (partially found — contig break inside the gene?)
M = Missing (not found — lost or too diverged)
n = total genes in the database
```

**Good bacterial assembly:** C ≥ 95%

---

## 📁 Repository Structure

```
denovo_assembly_course/
├── README.md                    ← you are here
├── scripts/
│   ├── 00_download_data.sh      ← download reads from NCBI SRA
│   ├── 01_pipeline_hybrid.sh    ← Pipeline 1: hybrid assembly + scaffolding
│   └── 02_pipeline_polish.sh    ← Pipeline 2: polishing + QC
├── data/
│   ├── raw/                     ← raw reads (fastq) — added by download script
│   ├── reference/               ← reference genomes / GTF annotation
│   └── db/                      ← BUSCO lineage databases
├── results/
│   ├── spades/                  ← SPAdes assembly outputs
│   ├── flye/                    ← Flye assembly outputs
│   ├── pilon/                   ← Pilon polishing outputs
│   ├── busco/                   ← BUSCO quality reports
│   ├── scaffolding/             ← LongStitch scaffolding outputs
│   └── sibelia/                 ← Sibelia + Circos comparison outputs
├── theory/
│   └── lecture_notes.md         ← full lecture notes on de novo assembly
└── docs/
    └── pipeline2_walkthrough.md ← step-by-step explanation of Pipeline 2
```

---

## 🚀 Quick Start

### 1 — Install dependencies (Ubuntu)

```bash
# Install Miniconda if you don't have it
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create the course environment
conda create -n denovo python=3.10 -y
conda activate denovo

# Install all tools
conda install -c bioconda -c conda-forge -y \
    spades flye bwa samtools busco pilon \
    seqtk sra-tools bedtools longstitch
```

### 2 — Clone this repository

```bash
git clone https://github.com/YOUR_USERNAME/denovo-assembly-course.git
cd denovo-assembly-course
```

### 3 — Download reads

> ✏️ **Replace the SRR accession with your own read ID before running!**

```bash
# Default accession is SRR25745292 (Vibrio cholerae, used as a teaching example)
# To use your own read:
bash scripts/00_download_data.sh SRRxxxxxxxx

# Or just edit the default inside the script:
#   open scripts/00_download_data.sh
#   find the line: SRR_ACCESSION="${1:-SRR25745292}"
#   replace SRR25745292 with your accession
```

### 4 — Run a pipeline

```bash
# Pipeline 1: hybrid assembly + scaffolding
bash scripts/01_pipeline_hybrid.sh

# Pipeline 2: polishing + quality assessment
bash scripts/02_pipeline_polish.sh
```

---

## Pipeline 1: Hybrid Assembly & Scaffolding

**Data flow:**
```
Illumina R1/R2 ──┐
                 ├──► SPAdes (hybrid) ──► seqtk ──► Sibelia ──► Circos
Nanopore reads ──┘         └──► LongStitch ──► Sibelia ──► Circos
```

| Step | Tool | What happens |
|------|------|-------------|
| Hybrid assembly | SPAdes | Builds a De Bruijn Graph from Illumina reads; uses Nanopore reads to resolve repeats and bridge gaps |
| Illumina-only assembly | SPAdes | Same genome assembled without long reads — used as a comparison baseline |
| Filter short contigs | seqtk | Removes contigs shorter than 10 kb (assembly noise) before visualisation |
| Synteny blocks | Sibelia | Finds collinear blocks shared between two assemblies |
| Visualisation | Circos | Generates a circular synteny plot (SVG/PNG) |
| Scaffolding | LongStitch (ntLink) | Joins contigs into longer scaffolds using Nanopore reads as bridges |

**Run:** `bash scripts/01_pipeline_hybrid.sh`

---

## Pipeline 2: Long Reads, Polishing & Quality Assessment

**Data flow:**
```
NCBI SRA ──► prefetch + fasterq-dump ──► raw reads
                                              │
                        Flye draft assembly ──┤◄── BWA index
                                              │
                    Illumina R1/R2 ──► BWA mem ──► samtools sort/index
                                              │
                                           Pilon ──► polished assembly
                                              │
                                           BUSCO ──► completeness report
                                              │
                                    GTF + bedtools ──► gene products in region
```

**Run:** `bash scripts/02_pipeline_polish.sh`

Full command-by-command explanation: [docs/pipeline2_walkthrough.md](docs/pipeline2_walkthrough.md)

---

## 🛠 Tools Used

| Tool | Version | Purpose | Docs |
|------|---------|---------|------|
| SPAdes | ≥ 3.15 | De novo assembly (DBG) | [link](https://ablab.github.io/spades/) |
| Flye | ≥ 2.9 | Long-read assembly (OLC) | [link](https://github.com/fenderglass/Flye) |
| BWA | ≥ 0.7 | Short-read alignment | [link](https://github.com/lh3/bwa) |
| samtools | ≥ 1.15 | BAM/SAM processing | [link](https://www.htslib.org/) |
| BUSCO | ≥ 5.0 | Assembly completeness | [link](https://busco.ezlab.org/) |
| Pilon | ≥ 1.24 | Assembly polishing | [link](https://github.com/broadinstitute/pilon) |
| seqtk | any | FASTA/Q filtering | [link](https://github.com/lh3/seqtk) |
| LongStitch | ≥ 1.0 | Long-read scaffolding | [link](https://github.com/bcgsc/LongStitch) |
| Sibelia | ≥ 3.0 | Synteny block detection | [link](https://github.com/bioinf/Sibelia) |
| Circos | ≥ 0.69 | Circular genome visualisation | [link](http://circos.ca/) |
| bedtools | ≥ 2.30 | Genomic interval operations | [link](https://bedtools.readthedocs.io/) |
| SRA Toolkit | ≥ 3.0 | Download from NCBI SRA | [link](https://github.com/ncbi/sra-tools) |

---

## 📖 Further Reading

- [Lecture notes: de novo assembly theory](theory/lecture_notes.md)
- [Pipeline 2 step-by-step walkthrough](docs/pipeline2_walkthrough.md)
- [SPAdes manual](https://ablab.github.io/spades/)
- [Flye manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
- [BUSCO user guide](https://busco.ezlab.org/busco_userguide.html)

---

## ❓ Questions?

Open an [Issue](../../issues) and include:
1. The exact command you ran
2. The full error message
3. Output of `conda list | grep -E "spades|bwa|busco|pilon"`

> **Tip for beginners:** Read the comments inside the bash scripts before running them.  
> Every flag and every symbol is explained.

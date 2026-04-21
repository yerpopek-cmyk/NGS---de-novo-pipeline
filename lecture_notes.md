# Lecture Notes: De Novo Genome Assembly

> These notes accompany the practical pipelines in this repository.  
> Designed for students taking their first steps in bioinformatics.

---

## 1. Why Is Assembly Even a Problem?

### The fundamental limitation: read length

Sequencers cannot read a whole chromosome at once. They produce **reads** — short fragments:

| Platform | Read length | Error rate |
|----------|-------------|-----------|
| Illumina | 100–300 bp | ~0.1% |
| PacBio HiFi | 10–25 kbp | ~0.1% |
| Oxford Nanopore | 10–100+ kbp | ~1–15% |

**Analogy:** Feed *War and Peace* through a shredder. You get a million strips of paper, each with 5–10 words. Reconstruct the book **without having the original**.

### Three enemies of every assembler

| Enemy | Why it's hard |
|-------|--------------|
| **Sequencing errors** | An algorithm requiring exact matches will break on a wrong nucleotide |
| **Repeats** | A read inside a repeat cannot tell the assembler which copy it came from |
| **Heterozygosity** | Two slightly different haplotypes may be collapsed into a fake "average" |

---

## 2. The Mathematics of Assembly

### Coverage depth

```
C = (N × L) / G

C = coverage depth  (how many times each base was sequenced on average)
N = number of reads
L = read length (bp)
G = genome size (bp)
```

**Example:** 10 million Illumina reads × 150 bp / 5,000,000 bp = **300×**

| Coverage | Outcome |
|----------|---------|
| < 5× | Gaps everywhere — assembly falls apart |
| 30–60× | Gold standard for Illumina short-read assembly |
| > 100× | Diminishing returns; sequencing errors create noise |

### k-mers

A **k-mer** is a substring of fixed length *k*.

```
Read: ATGCATG   k = 3
k-mers: ATG  TGC  GCA  CAT  ATG
```

**Why k-mers?**  
Comparing all *N* reads pairwise is O(N²) — impossible for 10⁹ reads.  
Decomposing reads into k-mers turns the problem into O(N) dictionary lookup.

**k-mer frequency spectrum:**

When you plot how often each unique k-mer appears, you see two peaks:

```
  count
    │   ┌──────── error peak (k-mers from sequencing errors)
    │  ■│
    │  ■│                ┌── genome peak (true k-mers, freq ≈ coverage)
    │  ■│               ■│
    │  ■■              ■■│
    │  ■■■■■■■■■■■■■■■■■■│
    └────────────────────────── frequency
         1  2  3  4  5      30×
```

SPAdes discards k-mers below a frequency threshold before building the graph,  
effectively removing most sequencing errors without touching the reads.

---

## 3. Two Assembly Algorithms

### OLC (Overlap–Layout–Consensus)

| Stage | What happens |
|-------|-------------|
| **Overlap** | Compare each read's end with every other read's start |
| **Layout** | Build a graph of overlapping reads; find a traversal order |
| **Consensus** | Stack reads, vote on each nucleotide |

**Best for:** Long reads (PacBio, Nanopore) — their overlaps span repeats.  
**Tools:** Flye, Canu, Hifiasm.

**Computational cost:** Finding all overlaps naively is O(N²).  
Modern tools use *minimizers* (representative k-mers) to reduce this to near O(N).

### DBG (De Bruijn Graph)

| Stage | What happens |
|-------|-------------|
| **k-merisation** | Cut every read into k-mers |
| **Graph construction** | Each unique k-mer = one node; consecutive k-mers share an edge |
| **Path finding** | Find an **Eulerian path** — visits every edge exactly once |

**The Eulerian path is solvable in linear time** — unlike the Hamiltonian path (NP-hard).  
This is why DBG assemblers are so fast.

**Best for:** Short, accurate Illumina reads.  
**Tools:** SPAdes, MEGAHIT, Velvet.

**The critical k-mer rule:**

```
Repeat length < k  →  the k-mer spanning the repeat boundary is unique
                       → repeat is resolved ✅

Repeat length > k  →  the k-mer sits entirely inside the repeat
                       → assembler cannot tell which path to take
                       → contig is broken ❌
```

This is why choosing the right *k* is important.  
SPAdes automatically tests multiple k values (e.g., 21, 33, 55, 77) and merges the graphs.

### Head-to-head comparison

| Property | OLC | DBG |
|----------|-----|-----|
| Input data | Long reads | Short reads |
| Repeat resolution | Yes (if reads span repeats) | Only if k > repeat length |
| Memory usage | High (stores all reads) | Lower (stores unique k-mers) |
| Speed | Slower (overlap is expensive) | Faster (hash-based) |
| Error tolerance | Handles noisy reads | Sensitive to errors |

---

## 4. Quality Metrics

### N50 — contiguity

**Algorithm:**
1. Sort contigs by length (largest first).
2. Accumulate lengths from the top.
3. N50 = the length of the contig where the cumulative sum first exceeds 50% of the total assembly length.

```
Example:
  Contig lengths: 500, 400, 300, 200, 100  (total = 1500 bp)
  50% threshold = 750 bp
  Running sum:
    500 → 500 (< 750)
    400 → 900 (> 750) ← N50 = 400
```

**Interpretation:** "Half the genome is in contigs of length ≥ N50"

> ⚠️ **Limitation:** N50 measures contiguity, not accuracy or completeness.  
> A bad assembler could produce a few giant chimeric contigs with a great N50.  
> Always pair N50 with BUSCO.

### BUSCO — gene completeness

BUSCO checks how many conserved orthologs (genes present in every member of a lineage) your assembly contains.

```
BUSCO output example:
  C:97.5% [S:97.0%, D:0.5%], F:1.5%, M:1.0%, n:1200

  C = Complete   (found in full length)
  S = Single     (one copy — expected)
  D = Duplicated (> 1 copy — possible assembly duplication artefact)
  F = Fragmented (partial hit — contig may break inside the gene)
  M = Missing    (not found)
  n = total reference genes in the database
```

**What is a good score?**

| Score | Assembly quality |
|-------|----------------|
| C ≥ 98% | Excellent — chromosome-level |
| C ≥ 95% | Good — suitable for most analyses |
| C 85–95% | Acceptable — investigate fragmented genes |
| C < 85% | Poor — check coverage, contamination, lineage choice |

---

## 5. Hybrid Assembly and Polishing

### The problem with long reads

Nanopore assemblies are contiguous (large N50, few gaps) but contain errors:
- Error rate: ~0.1–1% in modern Nanopore chemistry
- That's ~1 wrong base per 100–1000 bp
- Enough to frameshifts open reading frames and corrupt proteins

### Polishing workflow

```
Draft (Flye/Nanopore) ──► rough but complete genome

Illumina reads  ──► BWA MEM ──► aligned.bam ──► Pilon ──► corrected genome
                                                             error < 0.001%
```

**Pilon's logic at each position:**
- Tallies read votes: how many reads support A, T, G, C at this position?
- If the majority disagrees with the draft → correct it.
- Also handles small insertions and deletions.

### The road to T2T (Telomere-to-Telomere)

```
Ultra-long Nanopore (100–200 kbp) ──► Flye ──► chromosome-scale contigs
PacBio HiFi (20 kbp, 99.9%)       ──► Hifiasm ──► accurate assembly
                                            │
                              Hi-C / Pore-C ──► chromosome phasing
                                            │
                                         T2T ──► no Ns, full chromosomes
```

T2T assemblies (like the first complete human genome in 2022) are only possible by combining  
ultra-long reads, high accuracy, and chromosome scaffolding data.

---

## 6. Key Terms Reference

| Term | Definition |
|------|-----------|
| **contig** | Contiguous assembled sequence with no gaps |
| **scaffold** | Contigs joined by estimated gaps (Ns), often using paired-end or long-read evidence |
| **N50** | Contig length at which 50% of the assembly is covered |
| **BUSCO** | Benchmarking Universal Single-Copy Orthologs — gene completeness metric |
| **k-mer** | Substring of fixed length k |
| **coverage (depth)** | Average number of reads covering each base |
| **polishing** | Correcting nucleotide errors in a draft assembly using short reads |
| **scaffolding** | Ordering and orienting contigs using long-range information |
| **synteny** | Conserved gene order between two genomes or assemblies |
| **T2T** | Telomere-to-Telomere — gapless chromosome-scale assembly |
| **OLC** | Overlap–Layout–Consensus assembly strategy (long reads) |
| **DBG** | De Bruijn Graph assembly strategy (short reads) |
| **Eulerian path** | Graph path visiting every edge exactly once — solvable in O(E) |
| **Hamiltonian path** | Graph path visiting every vertex exactly once — NP-hard |

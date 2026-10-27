# HSPlice

**Raw, tblastn-based exon/intron prediction for phylogenetically conserved genes.**

HSPlice turns the tblastn hits of a set of reference proteins into one
predicted protein per gene and genomic contig. It chains the HSPs, cuts them
where tblastn read through an intron or a stop codon, searches for short exons
that tblastn missed, and places every intron on the donor/acceptor pair that
keeps the reading frame and best fits the query.

It is meant for phylogenomics: quickly pulling out orthologs of conserved
genes (PhyloFisher genes, BUSCO-like markers, etc.) from new or poorly annotated
genomes, including organisms with non-standard genetic codes and non-canonical
splice sites.

> **"Raw"** means that only the part of the gene aligned to the query is
> predicted. HSPlice does not search for start or stop codons, UTRs or
> alternative isoforms.

---

## Contents

- [HSPlice](#hsplice)
  - [Contents](#contents)
  - [How it works](#how-it-works)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Quick start](#quick-start)
  - [Input files](#input-files)
    - [Naming the reference proteins](#naming-the-reference-proteins)
    - [Consistency checks](#consistency-checks)
  - [Usage and options](#usage-and-options)
    - [Input and output](#input-and-output)
    - [Gene model](#gene-model)
    - [Run](#run)
  - [Examples](#examples)
  - [Output](#output)
    - [Protein FASTA](#protein-fasta)
    - [GFF3 (`-g`)](#gff3--g)
    - [Run summary](#run-summary)
  - [Non-standard genetic codes](#non-standard-genetic-codes)
  - [Limitations](#limitations)
  - [Exit status](#exit-status)
  - [Citation](#citation)
  - [License](#license)
  - [Authors](#authors)
  - [Acknowledgements](#acknowledgements)

---

## How it works

For every query protein and each of its best hits (contigs):

1. **Chaining** – the HSPs on the strand of the top HSP are combined into the
   best colinear set (dynamic programming over query and genome coordinates,
   respecting the maximum intron length).
2. **Splitting** – HSPs are cut at long query gaps and at stop codons, where
   tblastn may have read through an intron. Redundant HSP tails that overlap
   the neighbouring HSP are discarded.
3. **Missed-exon rescue** – short exons that tblastn could not report are
   searched for between the pieces and at both ends of the model
   (Smith-Waterman in a search space of a few kb, with a Karlin-Altschul
   significance threshold).
4. **Junction placement** – each junction is placed on the donor/acceptor pair
   that preserves the frame, contains no stop codon and gives the best global
   alignment to the query. Junctions are resolved independently, so there is no
   combinatorial search over all introns of a gene.
5. **Assembly and scoring** – the exons are joined, translated and scored
   against the query. A junction that cannot be resolved splits the model; the
   best-scoring part is kept and flagged as truncated.

When several query proteins belong to the same gene (e.g. `Hsap_RPL3`,
`Dmel_RPL3`, `Atha_RPL3`), all their models on a contig are compared, and the
one that fits all queries of that gene best is written.

---

## Requirements

| Software | Version | Notes |
|---|---|---|
| Python | 3.11 – 3.13 (tested) | |
| [Biopython](https://biopython.org) | ≥ 1.80 (tested up to 1.88) | required |
| [tqdm](https://github.com/tqdm/tqdm) | any | optional, for the progress bar |
| [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | 2.x | `makeblastdb` and `tblastn` |

---

## Installation

```
git clone https://github.com/DavidZihala/hsplice.git
cd hsplice
pip install biopython tqdm
chmod +x hsplice.py
./hsplice.py --help
```

---

## Quick start

```bash
# 1. BLAST database of the genome
makeblastdb -in genome.fa -dbtype nucl -out genome

# 2. tblastn with XML output
tblastn -query proteins.fa -db genome -outfmt 5 -evalue 1e-5 -out blast.xml

# 3. gene models
./hsplice.py proteins.fa genome.fa blast.xml Tthe Tthe.fa -t 8 -g Tthe.gff
```

This writes the predicted proteins to `Tthe.fa` and their exon coordinates to
`Tthe.gff`.

---

## Input files

| Argument | Description |
|---|---|
| `PROTEINS` | Reference proteins in FASTA — **the same file used as the tblastn query**. |
| `GENOME` | Genome in FASTA — **the same file the BLAST database was built from**. |
| `BLAST_XML` | tblastn output in XML format (`-outfmt 5`). Several XML results concatenated into one file are accepted. |
| `ORGANISM` | Prefix for the output sequence names (e.g. `Tthe`). |
| `OUTPUT` | Output FASTA file, or `-` for standard output. |

### Naming the reference proteins

Protein names must end with `_<GENE>`; everything after the **last**
underscore is taken as the gene name:

```
>Hsap_RPL3      -> gene RPL3
>Scer_RPL3      -> gene RPL3
>Hsap_EF1A      -> gene EF1A
```

All queries sharing the same `<GENE>` are treated as one gene. Gene names must
therefore not contain underscores (`Hsap_RPL3_A` would be read as gene `A`).

Gaps (`-`, `.`) and terminal `*` in the protein file are removed automatically.

### Consistency checks

HSPlice checks that the inputs belong together: every contig must exist in
`GENOME` with the length reported by BLAST, and the aligned query must match
the sequence in `PROTEINS`. Hits failing these checks are skipped and
reported in the summary; if no hit can be used, the run stops with exit
status 1.

---

## Usage and options

```
hsplice.py [options] PROTEINS GENOME BLAST_XML ORGANISM OUTPUT
```

### Input and output

| Option | Description |
|---|---|
| `-g`, `--gff FILE` | Also write the exon coordinates of the models as GFF3. |

### Gene model

| Option | Default | Description |
|---|---|---|
| `-c`, `--genetic-code CODE` | standard | NCBI translation table (e.g. `6`) or codon reassignments (e.g. `"TAA:Q;TAG:Q"`). |
| `-n`, `--hits INT` | `1` | Number of best hits (contigs) modelled per query protein. |
| `--min-intron NT` | `15` | Shortest intron accepted. |
| `--max-intron NT` | `15000` | Longest intron; also how far terminal exons are searched. |
| `--donors LIST` | `GT` | Donor dinucleotides, comma-separated (e.g. `GT,GC`). All but the first are penalised. |
| `--acceptors LIST` | `AG` | Acceptor dinucleotides, comma-separated. All but the first are penalised. |

### Run

| Option | Default | Description |
|---|---|---|
| `-t`, `--threads INT` | `1` | Worker processes; `0` uses all available cores. |
| `--no-progress` | | Do not draw a progress bar. |
| `-q`, `--quiet` | | Print warnings and errors only. |
| `-v`, `--verbose` | | Also report every model and every skipped hit. |
| `-V`, `--version` | | Show the version and exit. |
| `-h`, `--help` | | Show the help and exit. |

---

## Examples

Standard genetic code, 8 processes, gene structures as GFF3:

```bash
./hsplice.py proteins.fa genome.fa blast.xml Tthe Tthe.fa -t 8 -g Tthe.gff
```

Ciliate code (TAA and TAG encode Q), two best contigs per query — both forms
are equivalent:

```bash
./hsplice.py proteins.fa genome.fa blast.xml Tthe Tthe.fa -c 6 -n 2
./hsplice.py proteins.fa genome.fa blast.xml Tthe Tthe.fa -c "TAA:Q;TAG:Q" -n 2
```

Allow GC donors and long introns (e.g. vertebrates):

```bash
./hsplice.py proteins.fa genome.fa blast.xml Drer Drer.fa \
    --donors GT,GC --max-intron 50000 -t 0
```

Proteins to standard output, messages to a log file:

```bash
./hsplice.py proteins.fa genome.fa blast.xml Tthe - 2> run.log | gzip > Tthe.fa.gz
```

---

## Output

### Protein FASTA

One protein per gene and contig:

```
>Tthe_RPL3@scaffold_12
MSHRKFSAPRHGSLGFLPRKRSSRHRGKVKSFPKDDPSKPVHLTAFLGYKAGMTHIVREV...
```

The header is `>{ORGANISM}_{GENE}@{contig}`. With `-n 1` (default) each
query contributes its single best contig; with `-n N` a gene may be reported
on up to N contigs per query (useful for paralogs or fragmented assemblies).

### GFF3 (`-g`)

One `mRNA` feature per model and one `CDS` feature per exon, in plus-strand
coordinates with correct phases:

```
##gff-version 3
scaffold_12  raw_prediction  mRNA  10234  11873  712  +  .  ID=Tthe_RPL3@scaffold_12;query=Hsap_RPL3
scaffold_12  raw_prediction  CDS   10234  10512  .    +  0  ID=cds.Tthe_RPL3@scaffold_12;Parent=Tthe_RPL3@scaffold_12
scaffold_12  raw_prediction  CDS   10601  11873  .    +  0  ID=cds.Tthe_RPL3@scaffold_12;Parent=Tthe_RPL3@scaffold_12
```

The `mRNA` score is the local alignment score (BLOSUM62, gap 11/1) of the
predicted protein against the query that produced it; `query=` names that
query.

### Run summary

At the end of a run a summary is printed to standard error:

```
Summary
    BLAST queries                               1,520
      with a usable hit                         1,487
    Query/hit pairs modelled                    1,481
    Proteins written                              254  (254 genes on 231 contigs)
      with introns                                198
      with exons missed by tblastn                 37
      truncated at an unresolved junction           6
```

If more than 20 % of the models are truncated, HSPlice warns that the genetic
code or the minimum intron length is probably wrong.

---

## Non-standard genetic codes

For organisms with reassigned stop codons (ciliates, some green algae,
diplomonads, mitochondrial genomes, …) **both** tblastn and HSPlice need the
right code; otherwise the reassigned stops break the HSPs and the models.

```bash
tblastn -query proteins.fa -db genome -outfmt 5 -evalue 1e-5 \
        -db_gencode 6 -out blast.xml
./hsplice.py proteins.fa genome.fa blast.xml Tthe Tthe.fa -c 6
```

`--genetic-code` accepts any [NCBI translation table](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
id, or explicit reassignments such as `"TGA:W"` or `"TAA:Q;TAG:Q"`. The code in
use is printed at the start of every run.

Organisms with very short introns (e.g. some ciliates, ~20 nt) may need a
lower `--min-intron`.

---

## Limitations

- Only the region aligned to the query is predicted; N- and C-terminal parts
  with no similarity to any query are missing (no start/stop codon search).
- Only HSPs on the strand of the top HSP of a hit are used, and each contig
  yields at most one model per gene. Tandem paralogs on the same contig are
  therefore merged or only partially recovered.
- Genes split across contigs are reported as separate partial models.
- Splice sites are placed by sequence similarity and dinucleotide, without a
  splice-site model; junctions in poorly conserved regions may be imprecise.
- Intended for conserved genes; it is not a general-purpose gene finder.

---

## Exit status

| Code | Meaning |
|---|---|
| `0` | success |
| `1` | invalid input file |
| `2` | usage error |
| `130` | interrupted |

---

## Citation

If you use HSPlice in published work, please cite this repository:

```
TBD
```

---

## License

HSPlice is released under the [MIT License](LICENSE).

## Authors
- David Zihala
- Martin Kolisko
- Serafim Nenarokov

## Acknowledgements
This work was supported by the European Union project LERCO (number CZ.10.03.01/00/22_003/0000003).

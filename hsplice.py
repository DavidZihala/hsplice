#!/usr/bin/env python3
"""
Raw (tblastn-only) exon/intron prediction for phylogenetically conserved genes.

usage: raw_gene_prediction.py [options] PROTEINS GENOME BLAST_XML ORGANISM OUTPUT
       raw_gene_prediction.py --help      (all options, with examples)

How a gene model is built for every query protein and each of its best hits:

  1. the HSPs on the strand of the top HSP are chained into the best colinear
     set (dynamic programming),
  2. HSPs are cut at long query gaps and at stop codons, where tblastn may
     have read through an intron,
  3. exons that tblastn missed are searched for between the pieces and at
     both ends (Smith-Waterman in a search space of a few kb),
  4. each junction is placed on the donor/acceptor pair that keeps the frame,
     contains no stop codon and gives the best alignment to the query,
  5. the pieces are joined, translated and scored.

Nothing outside the aligned part of the query is predicted (no search for
start and stop codons), hence "raw".

Query names are expected to end with "_<GENE>" (e.g. Hsap_RPL3).  All queries
sharing <GENE> are treated as one gene; of the models they produce on a
contig, the one that fits these queries best is written as
>{organism}_{GENE}@{contig}.

BLAST must be run with XML output, e.g.

    makeblastdb -in genome.fa -dbtype nucl -out genome
    tblastn -query proteins.fa -db genome -outfmt 5 -evalue 1e-5 \\
            -gapopen 11 -gapextend 2 -out blast.xml

For organisms with a non-standard genetic code give tblastn the matching
-db_gencode and this script --genetic-code, otherwise the reassigned stop
codons break the HSPs and the models.

Requires Biopython >= 1.80; tqdm is optional.  Tested with Python 3.11-3.13
and Biopython 1.80-1.88.
"""

from __future__ import annotations

import argparse
import logging
import math
import multiprocessing as mp
import os
import shutil
import sys
import time
import xml.etree.ElementTree as ET
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
    from tqdm import tqdm as _tqdm
except ImportError:  # the progress bar is a nicety, not a requirement
    _tqdm = None

__version__ = "2.0.0"
PROG = "raw_gene_prediction.py"
log = logging.getLogger("raw_gene_prediction")

# --------------------------------------------------------------------------
# Heuristic constants (the ones that depend on the organism are CLI options)
# --------------------------------------------------------------------------
ANCHOR_AA = 10        # residues on each side of a junction that are re-aligned
SLACK_NT = 15         # how far past an HSP end a splice site may lie
MAX_EXTEND_AA = 30    # max. unaligned query residues an exon is extended by
MAX_INSERT_AA = 100   # max. insertion accepted when joining HSPs without intron
MIN_GAP_RUN = 3       # query-gap run (aa) inside an HSP that may hide an intron
MIN_BLOCK_AA = 5      # smaller pieces of a split HSP are chance matches
NONCANONICAL_PENALTY = 5    # score penalty for every non-first donor/acceptor
MAX_HSPS = 300        # HSPs per hit considered for chaining
MAX_PAIRS = 20000     # donor/acceptor pairs evaluated per junction
MIN_RESCUE_AA = 6     # shortest missed exon that is searched for
RESCUE_FLANK_AA = 25  # HSP ends re-examined when searching for a missed exon
RESCUE_EVALUE = 0.1   # significance required for a missed exon
MAX_RANK_QUERIES = 10  # queries used to choose between models of one gene
GAP_OPEN, GAP_EXTEND = 11, 1            # gap of length k costs 11 + k
KA_LAMBDA, KA_K = 0.267, 0.041          # Karlin-Altschul, BLOSUM62 11/1

_COMPLEMENT = str.maketrans("ACGTRYKMBDHVN", "TGCAYRMKVHDBN")


class InputError(Exception):
    """A problem with the input files that the user has to fix."""


# --------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Config:
    code: dict                       # codon -> amino acid ('*' = stop)
    min_intron: int = 15
    max_intron: int = 15000
    donors: tuple = ("GT",)
    acceptors: tuple = ("AG",)


@dataclass
class Hsp:
    bits: float
    q_start: int      # 0-based, inclusive
    q_end: int        # exclusive
    s_lo: int         # 1-based plus-strand coordinates, s_lo <= s_hi
    s_hi: int
    strand: int       # +1 / -1
    qseq: str         # aligned query (with '-')
    sseq: str         # aligned translated subject (with '-')
    start: int = 0    # 0-based start on the *oriented* contig / window
    end: int = 0      # exclusive

    def is_consistent(self) -> bool:
        """Alignment strings must agree with the reported coordinates."""
        return (len(self.qseq) == len(self.sseq)
                and self.q_end - self.q_start == len(self.qseq) - self.qseq.count("-")
                and self.s_hi - self.s_lo + 1 == 3 * (len(self.sseq) - self.sseq.count("-")))


@dataclass
class Block:
    """An ungapped-ish piece of alignment, i.e. (part of) a putative exon."""
    source: int                   # parental HSP (>= 0) or rescued exon (< 0)
    q_start: int
    s_start: int
    q_end: int = 0
    s_end: int = 0
    cols_q: list = field(default_factory=list)   # aligned columns: query index
    cols_s: list = field(default_factory=list)   # ... and codon start (nt)


@dataclass
class Locus:
    """Self-contained unit of work for one query/hit pair."""
    gene: str
    query_name: str
    query_seq: str
    contig: str
    contig_len: int
    strand: int
    offset: int       # start of `window` on the oriented contig
    window: str       # oriented, upper-case genomic sequence
    chain: list       # HSPs, coordinates relative to `window`


@dataclass
class Prediction:
    gene: str
    query_name: str
    contig: str
    strand: int
    exons: list       # [(start, end)] 1-based inclusive, plus strand, ascending
    protein: str
    score: float      # local alignment score against its own query
    query_seq: str = ""
    partial: bool = False   # truncated at a junction that could not be resolved
    rescued: int = 0        # exons found by the missed-exon search


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------
def build_genetic_code(spec: str | None) -> dict:
    """Standard code, an NCBI table id ("6") or reassignments ("TAA:Q;TAG:Q")."""
    table_id = int(spec) if spec and spec.strip().isdigit() else 1
    try:
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError:
        raise InputError(f"unknown NCBI genetic code table {table_id}") from None
    code = dict(table.forward_table)
    code.update(dict.fromkeys(table.stop_codons, "*"))
    if spec and not spec.strip().isdigit():
        for item in filter(None, (x.strip() for x in spec.split(";"))):
            try:
                codon, aa = (x.strip().upper() for x in item.split(":"))
            except ValueError:
                raise InputError(f'cannot parse --genetic-code item "{item}" (expected '
                                 'e.g. "TAA:Q;TAG:Q" or an NCBI table id)') from None
            codon = codon.replace("U", "T")
            if codon not in code or len(aa) != 1 or aa not in "ACDEFGHIKLMNPQRSTVWYX*":
                raise InputError(f'invalid codon reassignment "{item}"')
            code[codon] = aa
    return code


def read_fasta(path: str) -> dict:
    """{name: sequence}, the name being the first word of the header."""
    with open(path) as handle:
        seqs = {title.split()[0]: seq for title, seq in SimpleFastaParser(handle)
                if title.strip()}
    if not seqs:
        raise InputError(f"no FASTA sequences found in {path}")
    return seqs


def read_proteins(path: str):
    """({header or name: sequence}, [names]); both keys work for the lookup."""
    seqs, names = {}, []
    with open(path) as handle:
        for title, seq in SimpleFastaParser(handle):
            if title.strip():   # BLAST ignores alignment gaps in queries, too
                seq = seq.upper().replace("-", "").replace(".", "").rstrip("*")
                seqs[title] = seqs[title.split()[0]] = seq
                names.append(title.split()[0])
    if not seqs:
        raise InputError(f"no FASTA sequences found in {path}")
    return seqs, names


def parse_blast_xml(path: str, max_hits: int):
    """Stream (query_def, [(hit_id, hit_def, hit_accession, hit_len, [Hsp])]).

    Reads BLAST XML (-outfmt 5), including several results concatenated into
    one file, without loading it into memory.
    """
    def num(node, tag, cast=int):
        return cast(node.findtext(tag))

    def iterations():
        parser, root_seen = None, False
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if parser is None or line.startswith("<?xml"):
                    parser, root_seen = ET.XMLPullParser(("start", "end")), False
                try:
                    parser.feed(line)
                    events = list(parser.read_events())
                except ET.ParseError as err:
                    raise InputError(f"{path} is not valid BLAST XML (-outfmt 5): "
                                     f"{err}") from None
                for event, elem in events:
                    if event == "start":
                        if not root_seen and elem.tag != "BlastOutput":
                            raise InputError(f"{path}: expected BLAST XML (-outfmt 5), "
                                             f"found <{elem.tag}>")
                        root_seen = True
                    elif elem.tag == "BlastOutput_program" and elem.text != "tblastn":
                        raise InputError(f"{path} was produced by {elem.text}; "
                                         "tblastn output is required")
                    elif elem.tag == "Iteration":
                        yield elem
                        elem.clear()

    for elem in iterations():
        query = elem.findtext("Iteration_query-def") or ""
        if not query or query == "No definition line":
            query = elem.findtext("Iteration_query-ID") or ""
        hits = []
        for hit in elem.iterfind("Iteration_hits/Hit"):
            if len(hits) >= max_hits:
                break
            hsps = []
            for h in hit.iterfind("Hit_hsps/Hsp"):
                s_from, s_to = num(h, "Hsp_hit-from"), num(h, "Hsp_hit-to")
                frame = num(h, "Hsp_hit-frame")
                hsps.append(Hsp(
                    bits=num(h, "Hsp_bit-score", float),
                    q_start=num(h, "Hsp_query-from") - 1,
                    q_end=num(h, "Hsp_query-to"),
                    s_lo=min(s_from, s_to), s_hi=max(s_from, s_to),
                    strand=-1 if frame < 0 or s_from > s_to else 1,
                    qseq=h.findtext("Hsp_qseq"), sseq=h.findtext("Hsp_hseq")))
            if hsps:
                hits.append((hit.findtext("Hit_id") or "", hit.findtext("Hit_def") or "",
                             hit.findtext("Hit_accession") or "", num(hit, "Hit_len"), hsps))
        yield query, hits


def resolve_contig(hit_id: str, hit_def: str, hit_acc: str, genome: dict) -> str | None:
    """Find the genome FASTA name for a hit, whatever makeblastdb did to the ids."""
    candidates = [hit_id, hit_def.split()[0] if hit_def.split() else "", hit_acc]
    candidates += [part for part in hit_id.split("|") if part]
    return next((c for c in candidates if c in genome), None)


def query_matches(hsp: Hsp, query: str) -> bool:
    """Is the query in the alignment the sequence we were given? (X = masked)"""
    aligned = hsp.qseq.replace("-", "").upper()
    ref = query[hsp.q_start:hsp.q_end]
    return (len(aligned) == len(ref)
            and sum(x == y or x == "X" for x, y in zip(aligned, ref)) >= 0.9 * len(ref))


def gene_of(query_name: str) -> str:
    return query_name.split()[0].split("_")[-1]


# --------------------------------------------------------------------------
# Sequence helpers
# --------------------------------------------------------------------------
def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(nt: str, code: dict) -> str:
    """Translate complete codons; anything ambiguous becomes X."""
    get = code.get
    return "".join([get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3)])


_CFG: Config | None = None
_GLOBAL = _LOCAL = None      # aligners
_PAIRS: dict = {}            # (aa, aa) -> BLOSUM62 score
_CLEAN: dict = {}            # translation table: unknown letters -> X


def _init_worker(cfg: Config) -> None:
    """Runs once per process (works with fork, spawn and forkserver)."""
    global _CFG, _GLOBAL, _LOCAL, _PAIRS, _CLEAN
    _CFG = cfg
    matrix = substitution_matrices.load("BLOSUM62")
    aligners = []
    for mode in ("global", "local"):
        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.substitution_matrix = matrix
        aligner.open_gap_score = -(GAP_OPEN + GAP_EXTEND)
        aligner.extend_gap_score = -GAP_EXTEND
        aligners.append(aligner)
    _GLOBAL, _LOCAL = aligners
    _PAIRS = {(a, b): matrix[a, b] for a in matrix.alphabet for b in matrix.alphabet}
    _CLEAN = str.maketrans({chr(c): "X" for c in range(65, 91)
                            if chr(c) not in matrix.alphabet})


def _score(a: str, b: str, local: bool = False) -> float:
    if not a or not b:
        n = len(a) + len(b)
        return 0.0 if local or not n else -(GAP_OPEN + GAP_EXTEND * n)
    return (_LOCAL if local else _GLOBAL).score(a, b)


# --------------------------------------------------------------------------
# Step 1: choose a colinear chain of HSPs
# --------------------------------------------------------------------------
def best_chain(hsps: list, contig_len: int, max_intron: int) -> list:
    """Highest-scoring set of HSPs that is colinear in query and subject.

    Only HSPs on the strand of the top HSP are used; HSP coordinates are
    converted to the oriented contig (so that every gene reads left to right).
    """
    strand = hsps[0].strand
    cand = []
    for h in hsps[:MAX_HSPS]:
        if h.strand != strand or not h.is_consistent() or h.s_hi > contig_len:
            continue
        if strand > 0:
            h.start, h.end = h.s_lo - 1, h.s_hi
        else:
            h.start, h.end = contig_len - h.s_hi, contig_len - h.s_lo + 1
        cand.append(h)
    if not cand:
        return []
    cand.sort(key=lambda h: (h.start, h.q_start))
    total = [h.bits for h in cand]
    prev = [-1] * len(cand)
    for j, b in enumerate(cand):
        b_len = b.q_end - b.q_start
        for i in range(j):
            a = cand[i]
            if not (a.q_start < b.q_start and a.q_end < b.q_end
                    and a.start < b.start and a.end < b.end):
                continue
            if b.start - a.end > max_intron:
                continue
            a_len = a.q_end - a.q_start
            q_overlap = max(0, a.q_end - b.q_start)
            if 2 * q_overlap > min(a_len, b_len):
                continue
            if a.end - b.start > 3 * q_overlap + SLACK_NT:
                continue
            penalty = q_overlap * min(a.bits / a_len, b.bits / b_len)
            if total[i] + b.bits - penalty > total[j]:
                total[j] = total[i] + b.bits - penalty
                prev[j] = i
    j = max(range(len(cand)), key=total.__getitem__)
    chain = []
    while j != -1:
        chain.append(cand[j])
        j = prev[j]
    return chain[::-1]


# --------------------------------------------------------------------------
# Step 2: cut the HSPs into blocks
# --------------------------------------------------------------------------
def hsp_blocks(hsp: Hsp, index: int, seq: str) -> list:
    """Split an HSP at long query gaps and at stop codons.

    A long insertion relative to the query may be an in-frame intron that
    tblastn read through; a stop codon cannot be coding.  Whether the two
    sides are joined directly or by an intron is decided in resolve_junction.
    """
    code = _CFG.code
    blocks, cur = [], None
    q, s = hsp.q_start, hsp.start
    qseq, sseq = hsp.qseq, hsp.sseq
    i, n = 0, len(qseq)

    def close():
        nonlocal cur
        if cur is not None and cur.cols_q:
            cur.q_end, cur.s_end = q, s
            blocks.append(cur)
        cur = None

    while i < n:
        if qseq[i] == "-":                       # insertion in the subject
            j = i
            while j < n and qseq[j] == "-":
                j += 1
            end = s + 3 * (j - i)
            if j - i >= MIN_GAP_RUN or any(code.get(seq[p:p + 3]) == "*"
                                           for p in range(s, end, 3)):
                close()
            s, i = end, j
            continue
        aligned = sseq[i] != "-"
        if aligned and code.get(seq[s:s + 3]) == "*":
            close()
        else:
            if cur is None:
                cur = Block(index, q, s)
            if aligned:
                cur.cols_q.append(q)
                cur.cols_s.append(s)
        s += 3 * aligned
        q += 1
        i += 1
    close()
    if len(blocks) > 1:     # tiny pieces (e.g. between a stop codon and an exon) are noise
        blocks = ([b for b in blocks if len(b.cols_q) >= MIN_BLOCK_AA]
                  or [max(blocks, key=lambda b: len(b.cols_q))])
    return blocks


def chain_blocks(chain: list, seq: str) -> list:
    """Blocks of all HSPs in the chain, without redundant HSP tails.

    tblastn often runs past an exon end into the intron (after a gap) and
    aligns residues that really belong to the next exon.  Such a terminal
    block is recognisable because the neighbouring HSP covers the same query
    residues, and it is dropped here.
    """
    def covered(block, others):
        lo, hi = others[0].q_start, others[-1].q_end
        return 2 * (min(block.q_end, hi) - max(block.q_start, lo)) >= block.q_end - block.q_start

    per_hsp = [blocks for blocks in (hsp_blocks(h, i, seq) for i, h in enumerate(chain))
               if blocks]
    for left, right in zip(per_hsp, per_hsp[1:]):
        while len(left) > 1 and covered(left[-1], right):
            left.pop()
        while len(right) > 1 and covered(right[0], left):
            right.pop(0)
    return [block for blocks in per_hsp for block in blocks]


# --------------------------------------------------------------------------
# Step 3: look for exons tblastn missed (too short to seed/extend an HSP)
# --------------------------------------------------------------------------
def columns_score(block: Block, lo: int, hi: int, query: str, seq: str) -> float:
    """Substitution score of the aligned columns lo..hi-1 of a block."""
    code = _CFG.code
    return sum(_PAIRS.get((query[q], code.get(seq[s:s + 3], "X")), -4)
               for q, s in zip(block.cols_q[lo:hi], block.cols_s[lo:hi]))


def find_missed_exon(query, q_lo, q_hi, frames, seq, n_search, uid):
    """Best significant local alignment of query[q_lo:q_hi] to the given frames.

    `frames` holds one (start, end) nt range per reading frame.  Returns
    (block, score, significance threshold) or None.  The search space is tiny
    compared with a genome, so short exons that tblastn cannot report become
    significant here.
    """
    q_seg = query[q_lo:q_hi]
    threshold = math.log(KA_K * len(q_seg) * n_search / RESCUE_EVALUE) / KA_LAMBDA
    best_score, best = threshold, None
    for lo, hi in frames:
        pos = 0
        for orf in translate(seq[lo:hi], _CFG.code).split("*"):
            if len(orf) >= MIN_RESCUE_AA:
                score = _score(q_seg, orf, local=True)
                if score > best_score:
                    best_score, best = score, (lo + 3 * pos, orf)
            pos += len(orf) + 1
    if best is None:
        return None
    origin, orf = best
    block = Block(uid, 0, 0)
    for (q0, q1), (p0, p1) in zip(*_LOCAL.align(q_seg, orf)[0].aligned):
        block.cols_q.extend(range(q_lo + q0, q_lo + q1))
        block.cols_s.extend(range(origin + 3 * p0, origin + 3 * p1, 3))
    block.q_start, block.q_end = block.cols_q[0], block.cols_q[-1] + 1
    block.s_start, block.s_end = block.cols_s[0], block.cols_s[-1] + 3
    return block, best_score, threshold


def rescue_missed_exons(blocks: list, seq: str, query: str) -> list:
    """Add missed exons between the blocks and at both ends of the model.

    tblastn frequently aligns the residues of an exon it missed to whatever
    follows the neighbouring exon, or reads a short exon between two tiny
    introns in the wrong frame.  The search therefore also covers the last
    RESCUE_FLANK_AA residues of the neighbouring blocks (in the other two
    reading frames).  A new exon is accepted when it explains the query
    significantly better than the block ends it replaces.
    """
    max_intron = _CFG.max_intron
    tried, uid, k = set(), -1, 0
    while k <= len(blocks):
        a = blocks[k - 1] if k else None
        b = blocks[k] if k < len(blocks) else None
        k += 1
        if not (a and b) and (b.q_start if b else len(query) - a.q_end) < MIN_RESCUE_AA:
            continue                    # nothing is missing at this end
        # query segment and genomic region to search, flanks included
        if a:
            ia = min(bisect_left(a.cols_q, max(a.q_end - RESCUE_FLANK_AA,
                                               a.cols_q[len(a.cols_q) // 2])),
                     len(a.cols_q) - 1)
            q_lo, s_lo = a.cols_q[ia], a.cols_s[ia]
        else:
            q_lo, s_lo = 0, max(0, b.s_start - max_intron)
        if b:
            ib = bisect_left(b.cols_q, min(b.q_start + RESCUE_FLANK_AA,
                                           b.cols_q[len(b.cols_q) // 2]))
            q_hi, s_hi = b.cols_q[ib], b.cols_s[ib]
        else:
            q_hi, s_hi = len(query), min(len(seq), a.s_end + max_intron)
        gap = (q_lo, q_hi, s_lo, s_hi)
        if gap in tried or q_hi - q_lo < MIN_RESCUE_AA or s_hi - s_lo < 3 * MIN_RESCUE_AA:
            continue
        tried.add(gap)
        frames = []                     # a block's own frame is not searched again
        for frame in range(3):
            lo = a.s_end if a and a.s_start % 3 == frame else s_lo
            hi = b.s_start if b and b.s_start % 3 == frame else s_hi
            frames.append((lo + (frame - lo) % 3, hi))
        found = find_missed_exon(query, q_lo, q_hi, frames, seq, s_hi - s_lo, uid)
        if found is None:
            continue
        new, score, threshold = found
        # columns of the neighbours that the new exon would replace
        cut_a = min(bisect_left(a.cols_q, new.q_start),
                    bisect_right(a.cols_s, new.s_start - 3)) if a else 0
        cut_b = max(bisect_left(b.cols_q, new.q_end),
                    bisect_left(b.cols_s, new.s_end)) if b else 0
        replaced = ((columns_score(a, cut_a, None, query, seq) if a else 0)
                    + (columns_score(b, 0, cut_b, query, seq) if b else 0))
        if score - replaced < threshold:
            continue
        if a:
            del a.cols_q[cut_a:], a.cols_s[cut_a:]
            if a.cols_q:
                a.q_end, a.s_end = a.cols_q[-1] + 1, a.cols_s[-1] + 3
        if b:
            del b.cols_q[:cut_b], b.cols_s[:cut_b]
            if b.cols_q:
                b.q_start, b.s_start = b.cols_q[0], b.cols_s[0]
        blocks.insert(k - 1, new)
        blocks = [blk for blk in blocks if blk.cols_q]
        uid -= 1
        k = 0                           # neighbours changed: scan again
    return blocks


# --------------------------------------------------------------------------
# Step 4: place the junction between two consecutive blocks
# --------------------------------------------------------------------------
FAIL = "fail"
NO_INTRON = "no_intron"


def resolve_junction(a: Block, b: Block, seq: str, query: str, wide: bool = False):
    """Return (donor, acceptor) nt positions, NO_INTRON or FAIL.

    A stretch of ANCHOR_AA residues on each side of the junction is cut out of
    the two blocks; every in-frame, stop-free donor/acceptor combination in
    that region is translated and globally aligned to the corresponding query
    segment; the best-scoring one wins.  Because the frame is enforced
    arithmetically, junctions are independent of each other and no
    combinatorial search over all introns of a gene is needed.
    """
    cfg = _CFG
    margin = ANCHOR_AA * (2 if wide else 1)
    slack = SLACK_NT * (3 if wide else 1)

    # anchors: (query index, codon start) where the re-aligned region begins/ends
    ia = bisect_right(a.cols_q, min(a.q_end, b.q_start) - 1 - margin) - 1
    ia = min(max(ia, len(a.cols_q) // 2), len(a.cols_q) - 1)
    ib = min(bisect_left(b.cols_q, max(a.q_end, b.q_start) + margin), len(b.cols_q) // 2)
    q_a, s_a = a.cols_q[ia], a.cols_s[ia]
    q_b, s_b = b.cols_q[ib], b.cols_s[ib]
    if q_a >= q_b or s_a >= s_b:
        return FAIL
    q_seg = query[q_a:q_b]
    code = cfg.code
    best_score, best = None, FAIL

    # (a) no intron: the two blocks are simply joined
    span = s_b - s_a
    if span % 3 == 0 and (a.source == b.source
                          or span <= 3 * (q_b - q_a + MAX_INSERT_AA)):
        pep = translate(seq[s_a:s_b], code)
        if "*" not in pep:
            best_score, best = _score(q_seg, pep), NO_INTRON

    # (b) every frame-preserving donor ... acceptor pair
    q_gap = max(0, b.q_start - a.q_end)
    reach = slack + 3 * min(q_gap, MAX_EXTEND_AA)
    donors, acceptors = cfg.donors, cfg.acceptors
    d_sites = [d for d in range(s_a, min(a.s_end + reach, s_b) + 1)
               if seq[d:d + 2] in donors]
    a_sites = {0: [], 1: [], 2: []}
    for p in range(max(b.s_start - reach, s_a + cfg.min_intron, 2), s_b + 1):
        if seq[p - 2:p] in acceptors:
            a_sites[p % 3].append(p)
    evaluated = 0
    for d in d_sites:
        left = seq[s_a:d]
        d_pen = NONCANONICAL_PENALTY * (seq[d:d + 2] != donors[0])
        for p in a_sites[(d - s_a + s_b) % 3]:
            if not cfg.min_intron <= p - d <= cfg.max_intron:
                continue
            evaluated += 1
            if evaluated > MAX_PAIRS:
                break
            pep = translate(left + seq[p:s_b], code)
            if "*" in pep:
                continue
            score = (_score(q_seg, pep) - d_pen
                     - NONCANONICAL_PENALTY * (seq[p - 2:p] != acceptors[0]))
            if best_score is None or score > best_score:
                best_score, best = score, (d, p)
    return best


# --------------------------------------------------------------------------
# Step 5: assemble the gene model
# --------------------------------------------------------------------------
def predict_locus(locus: Locus):
    """Return (query name, contig, Prediction | None) for one query/hit pair."""
    cfg, seq = _CFG, locus.window
    label = (locus.query_name.split()[0], locus.contig)
    query = locus.query_seq.translate(_CLEAN)
    blocks = chain_blocks(locus.chain, seq)
    if not blocks:
        return (*label, None)
    blocks = rescue_missed_exons(blocks, seq, query)
    rescued = sum(block.source < 0 for block in blocks)

    # exon lists; a junction that cannot be resolved starts a new partial model
    models, exons, cur = [], [], blocks[0].s_start
    for a, b in zip(blocks, blocks[1:]):
        res = resolve_junction(a, b, seq, query)
        if res == FAIL:
            res = resolve_junction(a, b, seq, query, wide=True)
        if res == FAIL:
            exons.append((cur, a.s_end))
            models.append(exons)
            exons, cur = [], b.s_start
        elif res != NO_INTRON:
            exons.append((cur, res[0]))
            cur = res[1]
    exons.append((cur, blocks[-1].s_end))
    models.append(exons)

    best = None
    for exons in models:
        exons = [(s, e) for s, e in exons if e > s]
        protein = translate("".join(seq[s:e] for s, e in exons), cfg.code)
        if not protein or "*" in protein:
            continue
        score = _score(query, protein, local=True)
        if best is None or score > best[0]:
            best = (score, exons, protein)
    if best is None:
        return (*label, None)

    score, exons, protein = best
    plus = []
    for s, e in exons:            # window -> oriented contig -> plus strand
        s, e = s + locus.offset, e + locus.offset
        plus.append((s + 1, e) if locus.strand > 0
                    else (locus.contig_len - e + 1, locus.contig_len - s))
    return (*label, Prediction(locus.gene, locus.query_name, locus.contig,
                               locus.strand, sorted(plus), protein, score, query,
                               partial=len(models) > 1, rescued=rescued))


def rank_models(task):
    """Index of the protein that fits a set of query proteins best."""
    proteins, queries = task
    totals = [sum(_score(q, p, local=True) for q in queries) for p in proteins]
    return max(range(len(proteins)), key=lambda i: (totals[i], len(proteins[i])))


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def iter_queries(blast_xml, queries, genome, cfg, max_hits, stats):
    """Yield, per BLAST query, the list of self-contained Locus tasks."""
    def skip(reason, query_name, target=""):
        stats[reason] += 1
        log.debug("skipped %s%s: %s", query_name.split()[0],
                  f" on {target}" if target else "", SKIP_REASONS[reason])

    for query_name, hits in parse_blast_xml(blast_xml, max_hits):
        loci = []
        query_seq = queries.get(query_name) or queries.get(query_name.split()[0])
        if hits and not query_seq:
            skip("query_not_in_protein_dataset", query_name)
            hits = []
        for hit_id, hit_def, hit_acc, hit_len, hsps in hits:
            contig = resolve_contig(hit_id, hit_def, hit_acc, genome)
            if contig is None:
                skip("hit_not_in_genome", query_name, hit_def.split()[0] if hit_def.split()
                     else hit_id)
                continue
            if len(genome[contig]) != hit_len:
                skip("contig_length_mismatch", query_name, contig)
                continue
            if not query_matches(hsps[0], query_seq):
                skip("query_mismatch", query_name, contig)
                continue
            chain = best_chain(hsps, hit_len, cfg.max_intron)
            if not chain:
                skip("no_usable_hsp", query_name, contig)
                continue
            strand = chain[0].strand
            lo, hi = chain[0].start, chain[-1].end
            # room to look for terminal exons that tblastn may have missed
            if chain[0].q_start >= MIN_RESCUE_AA:
                lo = max(0, lo - cfg.max_intron)
            if len(query_seq) - chain[-1].q_end >= MIN_RESCUE_AA:
                hi = min(hit_len, hi + cfg.max_intron)
            if strand > 0:
                window = genome[contig][lo:hi].upper()
            else:
                window = reverse_complement(
                    genome[contig][hit_len - hi:hit_len - lo].upper())
            for h in chain:
                h.start -= lo
                h.end -= lo
            loci.append(Locus(gene_of(query_name), query_name, query_seq, contig,
                              hit_len, strand, lo, window, chain))
        stats["queries_with_hits"] += bool(loci)
        yield loci


def predict_query(loci: list) -> list:
    """All models of one BLAST query (one task for the worker pool)."""
    return [predict_locus(locus) for locus in loci]


def count_blast_queries(path: str) -> int:
    """Number of <Iteration> records, counted without parsing the XML."""
    count, tail = 0, b""
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1 << 24)
            if not chunk:
                return count
            data = tail + chunk
            count += data.count(b"<Iteration>")
            tail = data[-10:]


def write_fasta(path: str, organism: str, kept: dict) -> None:
    out = sys.stdout if path == "-" else open(path, "w")
    try:
        for (gene, contig), pred in kept.items():
            out.write(f">{organism}_{gene}@{contig}\n{pred.protein}\n")
    finally:
        if out is not sys.stdout:
            out.close()


def write_gff(path: str, organism: str, kept: dict) -> None:
    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for (gene, contig), p in kept.items():
            ident = f"{organism}_{gene}@{contig}"
            strand = "+" if p.strand > 0 else "-"
            out.write(f"{contig}\traw_prediction\tmRNA\t{p.exons[0][0]}\t{p.exons[-1][1]}"
                      f"\t{p.score:.0f}\t{strand}\t.\tID={ident};"
                      f"query={p.query_name.split()[0]}\n")
            done, rows = 0, []
            for s, e in (p.exons if p.strand > 0 else p.exons[::-1]):
                rows.append((s, e, (3 - done % 3) % 3))
                done += e - s + 1
            for s, e, phase in sorted(rows):
                out.write(f"{contig}\traw_prediction\tCDS\t{s}\t{e}\t.\t{strand}"
                          f"\t{phase}\tID=cds.{ident};Parent={ident}\n")


# --------------------------------------------------------------------------
# Terminal interface
# --------------------------------------------------------------------------
SKIP_REASONS = {
    "no_model": "no stop-free model could be built",
    "no_usable_hsp": "no usable HSP",
    "hit_not_in_genome": "contig is missing from the genome FASTA",
    "contig_length_mismatch": "contig length differs from the genome FASTA",
    "query_mismatch": "aligned query differs from the protein dataset",
    "query_not_in_protein_dataset": "query is missing from the protein dataset",
}
# reasons that point to a wrong input file rather than to biology
INPUT_PROBLEMS = {
    "hit_not_in_genome": "is GENOME the file the BLAST database was built from?",
    "contig_length_mismatch": "is GENOME the file the BLAST database was built from?",
    "query_mismatch": "is PROTEINS the file that was used as the tblastn query?",
    "query_not_in_protein_dataset": "is PROTEINS the file that was used as the "
                                    "tblastn query?",
}

EPILOG = """\
examples:
  # standard genetic code, 8 processes, gene structures as GFF3
  %(prog)s proteins.fa genome.fa blast.xml Tthe Tthe.fa -t 8 -g Tthe.gff

  # ciliate code (TAA and TAG encode Q), the two best contigs per query
  %(prog)s proteins.fa genome.fa blast.xml Tthe Tthe.fa -c 6 -n 2
  %(prog)s proteins.fa genome.fa blast.xml Tthe Tthe.fa -c "TAA:Q;TAG:Q" -n 2

  # proteins to standard output, messages to a log file
  %(prog)s proteins.fa genome.fa blast.xml Tthe - 2> run.log | gzip > Tthe.fa.gz

preparing BLAST_XML:
  makeblastdb -in genome.fa -dbtype nucl -out genome
  tblastn -query proteins.fa -db genome -outfmt 5 -evalue 1e-5 -out blast.xml
  (add -db_gencode N for an organism with a non-standard genetic code)

exit status:
  0  success          1  invalid input file
  2  usage error    130  interrupted
"""


class _HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Keeps the epilog layout and appends defaults that are worth showing."""

    def __init__(self, prog):
        super().__init__(prog, max_help_position=30, width=min(
            100, shutil.get_terminal_size((100, 24)).columns))

    def _get_help_string(self, action):
        text, default = action.help or "", action.default
        if (action.option_strings and action.nargs != 0 and not text.endswith("]")
                and default is not None and default is not argparse.SUPPRESS):
            shown = ",".join(default) if isinstance(default, tuple) else default
            text += f" [{shown}]"
        return text


class _LogHandler(logging.Handler):
    """Writes to stderr without breaking an active progress bar."""

    def emit(self, record):
        try:
            prefix = f"{record.levelname}: " if record.levelno >= logging.WARNING else ""
            stamp = time.strftime("%H:%M:%S", time.localtime(record.created))
            line = f"[{stamp}] {prefix}{record.getMessage()}"
            if _tqdm is not None:
                _tqdm.write(line, file=sys.stderr)
            else:
                print(line, file=sys.stderr, flush=True)
        except Exception:
            self.handleError(record)


def setup_logging(verbosity: int) -> None:
    """verbosity: -1 = warnings and errors only, 0 = normal, 1 = every model."""
    log.handlers[:] = [_LogHandler()]
    log.propagate = False
    log.setLevel({-1: logging.WARNING, 0: logging.INFO}.get(verbosity, logging.DEBUG))


class Progress:
    """A progress bar on a terminal, otherwise a log line every 30 seconds."""

    def __init__(self, total: int, what: str, unit: str, bar: bool):
        self.total, self.unit, self.done, self.last = total, unit, 0, time.time()
        self.bar = None
        if bar and _tqdm is not None and sys.stderr.isatty():
            size = shutil.get_terminal_size((80, 24))    # some ptys report 0 x 0
            self.bar = _tqdm(total=total, desc=what, unit=f" {unit}", file=sys.stderr,
                             ncols=min(100, size.columns), nrows=size.lines, leave=False)

    def update(self) -> None:
        self.done += 1
        if self.bar is not None:
            self.bar.update()
        elif time.time() - self.last >= 30:
            self.last = time.time()
            log.info("  %s of %s %s done (%.0f %%)", f"{self.done:,}", f"{self.total:,}",
                     self.unit, 100 * self.done / max(1, self.total))

    def close(self) -> None:
        if self.bar is not None:
            self.bar.close()


def describe_code(code: dict) -> str:
    standard = build_genetic_code(None)
    changes = [f"{codon}={aa}" for codon, aa in sorted(code.items())
               if standard[codon] != aa]
    stops = sorted(codon for codon, aa in code.items() if aa == "*")
    text = "standard" if not changes else (
        f"standard with {', '.join(changes)}" if len(changes) <= 6
        else f"{len(changes)} codons differ from the standard code")
    return f"{text}; stop codons: {', '.join(stops) or 'none'}"


def parse_args(argv=None):
    def positive(text):
        value = int(text)
        if value < 1:
            raise argparse.ArgumentTypeError("must be a positive integer")
        return value

    def process_count(text):
        value = int(text)
        if value < 0:
            raise argparse.ArgumentTypeError("must be 0 (all cores) or a positive integer")
        if value == 0:
            value = (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity")
                     else os.cpu_count() or 1)
        return value

    def splice_sites(text):
        sites = tuple(s.strip().upper().replace("U", "T") for s in text.split(",") if s.strip())
        if not sites or any(len(s) != 2 or set(s) - set("ACGT") for s in sites):
            raise argparse.ArgumentTypeError("expected dinucleotides such as GT or GT,GC")
        return sites

    def genetic_code(text):
        try:
            return build_genetic_code(text)
        except InputError as err:
            raise argparse.ArgumentTypeError(str(err)) from None

    p = argparse.ArgumentParser(
        prog=PROG, add_help=False, formatter_class=_HelpFormatter, epilog=EPILOG,
        usage="%(prog)s [options] PROTEINS GENOME BLAST_XML ORGANISM OUTPUT",
        description="Raw (tblastn-based) exon/intron prediction of phylogenetically "
                    "conserved genes.\nBuilds one protein per gene and contig from the "
                    "tblastn hits of a set of\nreference proteins; see the top of this "
                    "file for the method.")

    g = p.add_argument_group("input and output")
    g.add_argument("protein_dataset", metavar="PROTEINS",
                   help='reference proteins used as tblastn query (FASTA); names end '
                        'with "_<GENE>", e.g. Hsap_RPL3')
    g.add_argument("genome", metavar="GENOME",
                   help="genome the BLAST database was built from (FASTA)")
    g.add_argument("blast_output", metavar="BLAST_XML",
                   help="tblastn result in XML format (-outfmt 5)")
    g.add_argument("organism_name", metavar="ORGANISM",
                   help="prefix of the output names: >ORGANISM_GENE@contig")
    g.add_argument("output", metavar="OUTPUT",
                   help='predicted proteins (FASTA); "-" for standard output')
    g.add_argument("-g", "--gff", metavar="FILE",
                   help="also write the exon coordinates of the models as GFF3")

    g = p.add_argument_group("gene model")
    g.add_argument("-c", "--genetic-code", metavar="CODE", type=genetic_code,
                   help='NCBI translation table (e.g. 6) or codon reassignments '
                        '(e.g. "TAA:Q;TAG:Q") [standard code]')
    g.add_argument("-n", "--hits", metavar="INT", type=positive, default=1,
                   help="best hits (contigs) to model per query protein")
    g.add_argument("--min-intron", metavar="NT", type=positive, default=15,
                   help="shortest intron accepted")
    g.add_argument("--max-intron", metavar="NT", type=positive, default=15000,
                   help="longest intron; also how far terminal exons are searched")
    g.add_argument("--donors", metavar="LIST", type=splice_sites, default=("GT",),
                   help="donor dinucleotides, comma-separated; all but the first "
                        "are penalised, e.g. GT,GC")
    g.add_argument("--acceptors", metavar="LIST", type=splice_sites, default=("AG",),
                   help="acceptor dinucleotides, comma-separated")

    g = p.add_argument_group("run")
    g.add_argument("-t", "--threads", metavar="INT", type=process_count, default=1,
                   help="worker processes; 0 uses all available cores")
    g.add_argument("--no-progress", dest="progress", action="store_false",
                   help="do not draw a progress bar")
    verbosity = g.add_mutually_exclusive_group()
    verbosity.add_argument("-q", "--quiet", dest="verbosity", action="store_const",
                           const=-1, default=0, help="print warnings and errors only")
    verbosity.add_argument("-v", "--verbose", dest="verbosity", action="store_const",
                           const=1, help="also report every model and every skipped hit")
    g.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}",
                   help="show the version and exit")
    g.add_argument("-h", "--help", action="help", help="show this help and exit")

    args = p.parse_args(argv)
    if args.max_intron < args.min_intron:
        p.error("--max-intron must not be smaller than --min-intron")
    inputs = [args.protein_dataset, args.genome, args.blast_output]
    for path in inputs:
        if not os.path.isfile(path):
            p.error(f"input file not found: {path}")
        if not os.access(path, os.R_OK):
            p.error(f"input file is not readable: {path}")
    for path in filter(None, (args.output if args.output != "-" else None, args.gff)):
        folder = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(folder):
            p.error(f"output directory does not exist: {folder}")
        if any(os.path.exists(path) and os.path.samefile(path, i) for i in inputs):
            p.error(f"output would overwrite an input file: {path}")
    if args.gff and args.gff == args.output:
        p.error("OUTPUT and --gff must be different files")
    if args.genetic_code is None:
        args.genetic_code = build_genetic_code(None)
    return args


def run(args) -> None:
    started = time.time()
    cfg = Config(args.genetic_code, args.min_intron, args.max_intron,
                 args.donors, args.acceptors)
    log.info("%s %s", PROG, __version__)
    log.info("Genetic code  %s", describe_code(cfg.code))
    log.info("Introns       %d-%d nt, donor %s, acceptor %s", cfg.min_intron,
             cfg.max_intron, "/".join(cfg.donors), "/".join(cfg.acceptors))

    queries, names = read_proteins(args.protein_dataset)
    log.info("Proteins      %s sequences, %s genes (%s)", f"{len(names):,}",
             f"{len({gene_of(n) for n in names}):,}", args.protein_dataset)
    genome = read_fasta(args.genome)
    log.info("Genome        %s contigs, %.2f Mb (%s)", f"{len(genome):,}",
             sum(map(len, genome.values())) / 1e6, args.genome)
    n_queries = count_blast_queries(args.blast_output)
    if not n_queries:
        raise InputError(f"{args.blast_output} contains no BLAST queries "
                         "(is it tblastn output in XML format, -outfmt 5?)")
    log.info("BLAST         %s queries (%s)", f"{n_queries:,}", args.blast_output)
    log.info("Building gene models for the best %s per query with %d process%s",
             "hit" if args.hits == 1 else f"{args.hits} hits",
             args.threads, "" if args.threads == 1 else "es")

    skipped, stats = Counter(), Counter()   # filled by iter_queries / by the results
    pool = None
    if args.threads > 1:
        pool = mp.Pool(args.threads, initializer=_init_worker, initargs=(cfg,))
    else:
        _init_worker(cfg)

    def apply(func, tasks):
        return pool.imap(func, tasks, chunksize=8) if pool else map(func, tasks)

    try:
        # 1. one model per query/hit pair
        models = defaultdict(list)            # (gene, contig) -> [Prediction]
        gene_queries = defaultdict(dict)      # gene -> {query sequence: None}
        tasks = iter_queries(args.blast_output, queries, genome, cfg, args.hits, skipped)
        show_bar = args.progress and args.verbosity >= 0
        progress = Progress(n_queries, "Building gene models", "queries", show_bar)
        try:
            for results in apply(predict_query, tasks):
                progress.update()
                for name, contig, pred in results:
                    stats["pairs"] += 1
                    if pred is None:
                        stats["no_model"] += 1
                        log.debug("skipped %s on %s: %s", name, contig,
                                  SKIP_REASONS["no_model"])
                        continue
                    models[pred.gene, pred.contig].append(pred)
                    gene_queries[pred.gene][pred.query_seq] = None
                    log.debug("%s on %s: %d exon%s, %d aa%s", name, contig,
                              len(pred.exons), "" if len(pred.exons) == 1 else "s",
                              len(pred.protein), ", truncated" if pred.partial else "")
        finally:
            progress.close()

        # 2. per gene and contig keep the model that fits all queries best
        kept, undecided = {}, []
        for key, preds in models.items():
            distinct = list({p.protein: p for p in preds}.values())
            kept[key] = distinct[0]
            if len(distinct) > 1:
                undecided.append((key, distinct))
        if undecided:
            log.info("Choosing between alternative models of %s loci",
                     f"{len(undecided):,}")
        tasks = [([p.protein for p in preds],
                  list(gene_queries[key[0]])[:MAX_RANK_QUERIES])
                 for key, preds in undecided]
        progress = Progress(len(tasks), "Choosing models", "loci", show_bar)
        try:
            for (key, preds), best in zip(undecided, apply(rank_models, tasks)):
                progress.update()
                kept[key] = preds[best]
        finally:
            progress.close()
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()

    # ---- report
    stats.update(skipped)
    genes = {gene for gene, _ in kept}
    rows = [("BLAST queries", n_queries, ""),
            ("  with a usable hit", stats["queries_with_hits"], ""),
            ("Query/hit pairs modelled", stats["pairs"] - stats["no_model"], ""),
            ("Proteins written", len(kept),
             f"{len(genes):,} gene{'' if len(genes) == 1 else 's'} on "
             f"{len({c for _, c in kept}):,} contig"
             f"{'' if len({c for _, c in kept}) == 1 else 's'}"),
            ("  with introns", sum(len(p.exons) > 1 for p in kept.values()), ""),
            ("  with exons missed by tblastn", sum(p.rescued > 0 for p in kept.values()), ""),
            ("  truncated at an unresolved junction",
             sum(p.partial for p in kept.values()), "")]
    rows += [(f"Skipped: {text}", stats[key], "")
             for key, text in SKIP_REASONS.items() if stats[key]]
    width = max(len(label) for label, _, _ in rows)
    table = "\n".join(f"    {label:<{width}}  {value:>9,}{'  (' + note + ')' if note else ''}"
                      for label, value, note in rows)
    log.info("Summary\n%s", table)
    for key, hint in INPUT_PROBLEMS.items():
        if stats[key]:
            log.warning("%s hit%s skipped: %s - %s", f"{stats[key]:,}",
                        "" if stats[key] == 1 else "s", SKIP_REASONS[key], hint)
    truncated = sum(p.partial for p in kept.values())
    if kept and truncated > 0.2 * len(kept):
        log.warning("%d %% of the models are truncated; check --genetic-code "
                    "(and tblastn -db_gencode) and --min-intron",
                    100 * truncated / len(kept))
    if not kept:
        if any(stats[key] for key in INPUT_PROBLEMS):
            raise InputError("none of the BLAST hits could be used: the input files do "
                             "not belong together; nothing was written")
        log.warning("no protein could be predicted")
    write_fasta(args.output, args.organism_name, kept)
    if args.gff:
        write_gff(args.gff, args.organism_name, kept)
    outputs = ("standard output" if args.output == "-" else args.output) + (
        f" and {args.gff}" if args.gff else "")
    log.info("Wrote %s protein%s to %s in %s", f"{len(kept):,}",
             "" if len(kept) == 1 else "s", outputs, _elapsed(time.time() - started))


def _elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f} s"
    minutes, seconds = divmod(int(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes:02d} min" if hours else f"{minutes} min {seconds:02d} s"


def main(argv=None) -> int:
    args = parse_args(argv)
    setup_logging(args.verbosity)
    try:
        run(args)
    except InputError as err:
        log.error("%s", err)
        return 1
    except KeyboardInterrupt:
        log.error("interrupted")
        return 130
    except BrokenPipeError:          # e.g. "... - | head"
        sys.stderr.close()
        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
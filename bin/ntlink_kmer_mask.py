#!/usr/bin/env python3
"""
Coordinate-based minimizer masking for ntLink.

ntLink minimizers are ntHash *hash values*, not k-mer sequences, so a meryl
k-mer set cannot be joined to them directly. This module instead masks
minimizers by assembly coordinate: the user supplies a BED file of regions
whose minimizers should be excluded (or exclusively kept), typically derived
from read k-mer copy number (meryl-lookup) or read depth (mosdepth).

Typical use: mask assembly regions whose read-derived k-mer copy number
exceeds the 1-copy/2-copy crossover, i.e. collapsed repeats that look unique
in the assembly but are multi-copy in the genome.
"""

__author__ = "ntLink precision-scaffolding extension"

import sys
import gzip
from collections import defaultdict
import bisect
import numpy as np


def read_bed_regions(bed_filename):
    """
    Read a BED file into {contig: ([starts], [ends])} with merged, sorted
    intervals. BED is 0-based half-open [start, end).

    Returns an empty dict if bed_filename is falsy.
    """
    if not bed_filename:
        return {}

    raw = defaultdict(list)
    with open(bed_filename, 'r') as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                continue
            contig = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end <= start:
                continue
            raw[contig].append((start, end))

    regions = {}
    for contig, intervals in raw.items():
        intervals.sort()
        merged = []
        cur_start, cur_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= cur_end:          # overlapping or abutting
                cur_end = max(cur_end, end)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((cur_start, cur_end))
        regions[contig] = ([iv[0] for iv in merged], [iv[1] for iv in merged])

    return regions


def in_regions(regions, contig, position):
    """Return True if position falls within any interval of contig."""
    if contig not in regions:
        return False
    starts, ends = regions[contig]
    idx = bisect.bisect_right(starts, position) - 1
    return idx >= 0 and position < ends[idx]


def load_kmer_set(path, k):
    """
    Load a canonical k-mer set. Accepts either k-mer sequences (as `meryl print`
    emits, optionally with a trailing count column) or the packed canonical
    integers used elsewhere in this pipeline. Plain or gzipped.

    Canonical packing matches ntlink_complement: 2 bits per base, A=0 C=1 G=2 T=3,
    value = min(forward, reverse complement).
    """
    if not path:
        return None
    opener = gzip.open if path.endswith(".gz") else open
    codes = {"A": 0, "C": 1, "G": 2, "T": 3}
    out = set()
    with opener(path, "rt") as handle:
        for line in handle:
            token = line.split()[0] if line.strip() else ""
            if not token:
                continue
            if token[0] in "ACGTacgt":
                token = token.upper()
                if len(token) != k or any(base not in codes for base in token):
                    continue
                fwd = rev = 0
                for i, base in enumerate(token):
                    value = codes[base]
                    fwd = (fwd << 2) | value
                    rev |= (3 - value) << (2 * i)
                out.add(min(fwd, rev))
            else:
                try:
                    out.add(int(token))
                except ValueError:
                    continue
    return out


def _canonical_at(codes, positions, k):
    """Canonical packed k-mers at the given start positions. Returns (values, valid)."""
    if len(positions) == 0:
        return np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=bool)
    window = positions[:, None] + np.arange(k)
    sub = codes[window]
    valid = (sub != 255).all(axis=1)
    safe = np.where(sub == 255, 0, sub).astype(np.uint64)
    fwd = np.zeros(len(positions), dtype=np.uint64)
    rev = np.zeros(len(positions), dtype=np.uint64)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | safe[:, i]
        rev = rev | ((np.uint64(3) - safe[:, i]) << np.uint64(2 * i))
    return np.minimum(fwd, rev), valid


def filter_mx_info_by_kmers(mx_info, fasta_path, k, exclude_kmers=None,
                            include_kmers=None, verbose=True):
    """
    Filter minimizers by the identity of the k-mer they represent, not by
    coordinate. This is what read-derived copy number should drive: masking a
    BED interval removes every k-mer overlapping it, including the individually
    unique ones, whereas this removes exactly the k-mers named.

    Minimizers are ntHash values, so the k-mer sequence is recovered from the
    assembly at each minimizer's (contig, position).
    """
    if not exclude_kmers and not include_kmers:
        return mx_info

    from collections import defaultdict
    by_contig = defaultdict(list)
    for mx, info in mx_info.items():
        by_contig[info.contig].append((info.position, mx))

    n_before = len(mx_info)
    drop = set()
    n_missing = 0

    for name, seq in _iter_fasta(fasta_path):
        entries = by_contig.get(name)
        if not entries:
            continue
        codes = _encode(seq)
        positions = np.fromiter((p for p, _ in entries), dtype=np.int64, count=len(entries))
        in_range = positions + k <= len(codes)
        canon, valid = _canonical_at(codes, positions[in_range], k)
        usable = np.flatnonzero(in_range)
        for local, (value, ok) in enumerate(zip(canon, valid)):
            _, mx = entries[usable[local]]
            if not ok:
                drop.add(mx)
                continue
            value = int(value)
            if include_kmers is not None and value not in include_kmers:
                drop.add(mx)
            elif exclude_kmers and value in exclude_kmers:
                drop.add(mx)
        n_missing += int((~in_range).sum())

    kept = {mx: info for mx, info in mx_info.items() if mx not in drop}

    if verbose:
        print("Minimizer k-mer masking:", file=sys.stdout)
        print("\tminimizers before masking: {}".format(n_before), file=sys.stdout)
        print("\tdropped by k-mer identity: {}".format(len(drop)), file=sys.stdout)
        if n_missing:
            print("\tpositions past contig end (skipped): {}".format(n_missing),
                  file=sys.stdout)
        print("\tminimizers retained: {} ({:.2f}%)".format(
            len(kept), 100.0 * len(kept) / n_before if n_before else 0.0),
              file=sys.stdout)
    return kept


def _iter_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def _encode(seq):
    lut = np.full(256, 255, dtype=np.uint8)
    for base, code in zip(b"ACGT", [0, 1, 2, 3]):
        lut[base] = code
    return lut[np.frombuffer(seq.upper().encode(), dtype=np.uint8)]
def filter_mx_info_by_regions(mx_info, exclude_bed=None, include_bed=None,
                              verbose=True):
    """
    Filter an mx -> Minimizer dict by assembly coordinate.

    exclude_bed: minimizers inside these regions are dropped (blacklist).
    include_bed: only minimizers inside these regions are kept (whitelist).
    Both may be combined; exclude is applied after include.

    Returns the filtered dict. A Minimizer must expose .contig and .position.
    """
    if not exclude_bed and not include_bed:
        return mx_info

    exclude = read_bed_regions(exclude_bed)
    include = read_bed_regions(include_bed)

    n_before = len(mx_info)
    kept = {}
    n_dropped_include, n_dropped_exclude = 0, 0

    for mx, info in mx_info.items():
        if include and not in_regions(include, info.contig, info.position):
            n_dropped_include += 1
            continue
        if exclude and in_regions(exclude, info.contig, info.position):
            n_dropped_exclude += 1
            continue
        kept[mx] = info

    if verbose:
        print("Minimizer coordinate masking:", file=sys.stdout)
        print("\tminimizers before masking: {}".format(n_before), file=sys.stdout)
        if include:
            print("\tdropped (outside --mx-include-bed): {}".format(n_dropped_include),
                  file=sys.stdout)
        if exclude:
            print("\tdropped (inside --mx-exclude-bed): {}".format(n_dropped_exclude),
                  file=sys.stdout)
        print("\tminimizers retained: {} ({:.2f}%)".format(
            len(kept), 100.0 * len(kept) / n_before if n_before else 0.0),
              file=sys.stdout)
        if n_before and len(kept) / n_before < 0.05:
            print("\tWARNING: <5% of minimizers retained; scaffolding will likely "
                  "produce very few joins. Check that BED coordinates match the "
                  "target assembly.", file=sys.stderr)

    return kept
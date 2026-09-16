#!/usr/bin/env python3
"""
Complement seeding for ntLink.

Minimizer sampling picks k-mers by hash rank, which is blind to uniqueness.
In regions where assembly-unique k-mers are scarce (repeat interiors made
unique only by local divergence), the 2/(w+1) sampling rate frequently picks
none of them, leaving the region unanchored even though usable k-mers exist.

This module finds those anchor-starved windows and emits the unique k-mers
inside them as a supplementary anchor set. Because the set is defined by the
assembly, reads are matched against it by direct lookup rather than by their
own window competition, so a read contributes every complement k-mer it
carries error-free.

Optionally intersect with an allowed k-mer list (read-derived copy number,
e.g. meryl: read 1-copy, or read 2-copy confirmed count-1 in both haplotypes).
"""

__author__ = "ntLink precision-scaffolding extension"

import argparse
import sys
import numpy as np

COMPLEMENT_LUT = None


def load_fasta(path):
    name, chunks = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def encode(seq):
    lut = np.full(256, 255, dtype=np.uint8)
    for base, code in zip(b"ACGT", [0, 1, 2, 3]):
        lut[base] = code
    return lut[np.frombuffer(seq.upper().encode(), dtype=np.uint8)]


def kmer_ints(codes, k):
    """Canonical 2-bit packed k-mers (k<=32) plus validity mask."""
    n = len(codes)
    m = n - k + 1
    if m <= 0:
        empty = np.zeros(0, dtype=np.uint64)
        return empty, empty, np.zeros(0, dtype=bool)
    valid_base = codes != 255
    safe = np.where(valid_base, codes, 0).astype(np.uint64)
    fwd = np.zeros(m, dtype=np.uint64)
    rev = np.zeros(m, dtype=np.uint64)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | safe[i:i + m]
        rev = rev | ((np.uint64(3) - safe[i:i + m]) << np.uint64(2 * i))
    cs = np.concatenate(([0], np.cumsum(~valid_base)))
    valid = (cs[k:] - cs[:-k]) == 0
    return fwd, rev, valid


def splitmix64(x):
    z = x + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


def minimizer_mask(h, w):
    """Positions that are the minimum of at least one window of w."""
    from collections import deque
    n = len(h)
    picked = np.zeros(n, dtype=bool)
    if n == 0:
        return picked
    if n < w:
        picked[int(np.argmin(h))] = True
        return picked
    dq = deque()
    for i in range(n):
        while dq and h[dq[-1]] >= h[i]:
            dq.pop()
        dq.append(i)
        if dq[0] <= i - w:
            dq.popleft()
        if i >= w - 1:
            picked[dq[0]] = True
    return picked


def analyse_contig(seq, k, w):
    """Return canonical k-mers, uniqueness mask, minimizer mask, validity."""
    codes = encode(seq)
    fwd, rev, valid = kmer_ints(codes, k)
    canon = np.minimum(fwd, rev)
    del fwd, rev, codes

    cv = canon[valid]
    if len(cv):
        _, inverse, counts = np.unique(cv, return_inverse=True, return_counts=True)
        uniq_valid = counts[inverse] == 1
    else:
        uniq_valid = np.zeros(0, dtype=bool)
    is_unique = np.zeros(len(canon), dtype=bool)
    is_unique[valid] = uniq_valid

    h = splitmix64(canon)
    h[~valid] = np.iinfo(np.uint64).max
    picked = minimizer_mask(h, w) & valid
    return canon, is_unique, picked, valid


def select_complement(canon, is_unique, anchors, bin_bp, min_anchors, max_per_bin):
    """
    Flag unique k-mers inside bins that hold fewer than `min_anchors`
    minimizer anchors. Returns a boolean mask over k-mer positions.
    """
    n = len(canon)
    take = np.zeros(n, dtype=bool)
    if n == 0:
        return take
    nbin = n // bin_bp
    trim = nbin * bin_bp

    windows = []
    if nbin:
        per_bin = anchors[:trim].reshape(nbin, bin_bp).sum(axis=1)
        windows = [(b * bin_bp, (b + 1) * bin_bp)
                   for b in np.where(per_bin < min_anchors)[0]]

    # The trailing partial window is a contig end, which is exactly where
    # anchors decide whether a join can form, so it is evaluated too. Its
    # threshold is scaled to its length so a short tail is not called starved
    # merely for being short.
    if trim < n:
        tail = n - trim
        tail_anchors = int(anchors[trim:].sum())
        scaled = max(1, int(round(min_anchors * tail / bin_bp)))
        if tail_anchors < scaled:
            windows.append((trim, n))

    for lo, hi in windows:
        idx = np.flatnonzero(is_unique[lo:hi]) + lo
        if len(idx) == 0:
            continue
        if len(idx) > max_per_bin:          # even spacing keeps coverage uniform
            sel = np.linspace(0, len(idx) - 1, max_per_bin).astype(int)
            idx = idx[sel]
        take[idx] = True
    return take


def main():
    ap = argparse.ArgumentParser(
        description="Emit complement anchors for anchor-starved regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("fasta", help="Target assembly FASTA")
    ap.add_argument("-k", type=int, default=32, help="k-mer size (<=32)")
    ap.add_argument("-w", type=int, default=250, help="Minimizer window")
    ap.add_argument("--bin", type=int, default=10000, help="Window size for starvation test")
    ap.add_argument("--min-anchors", type=int, default=2,
                    help="A bin with fewer anchors than this is starved")
    ap.add_argument("--max-per-bin", type=int, default=200,
                    help="Cap on complement k-mers added per starved bin")
    ap.add_argument("--allowed-kmers", default=None,
                    help="Optional file of allowed canonical k-mers (sequences from "
                         "`meryl print` or packed integers). Complement anchors are "
                         "restricted to this set")
    ap.add_argument("-o", "--out", default="complement.tsv",
                    help="Output TSV: contig, position, canonical k-mer, strand")
    ap.add_argument("--stats", default=None, help="Optional npz of per-bin statistics")
    args = ap.parse_args()

    allowed = None
    if args.allowed_kmers:
        import ntlink_kmer_mask
        kmer_set = ntlink_kmer_mask.load_kmer_set(args.allowed_kmers, args.k)
        allowed = np.sort(np.fromiter(kmer_set, dtype=np.uint64, count=len(kmer_set)))
        print(f"Allowed k-mer list: {len(allowed):,} entries", file=sys.stderr)

    out = open(args.out, "w")
    all_stats = []
    tot_anchor = tot_comp = tot_bins = tot_starved = 0

    for name, seq in load_fasta(args.fasta):
        canon, is_unique, picked, valid = analyse_contig(seq, args.k, args.w)
        anchors = picked & is_unique
        take = select_complement(canon, is_unique, anchors,
                                 args.bin, args.min_anchors, args.max_per_bin)

        if allowed is not None and take.any():
            idx = np.flatnonzero(take)
            pos_in = np.searchsorted(allowed, canon[idx])
            pos_in = np.clip(pos_in, 0, len(allowed) - 1)
            keep = allowed[pos_in] == canon[idx]
            take[idx[~keep]] = False

        take &= ~anchors                     # never duplicate an existing anchor
        # Emit the canonical k-mer integer itself as the identifier: it is
        # identical on the assembly and the read side, which is what makes
        # direct lookup work. Strand is '+' when the forward k-mer is the
        # canonical one, matching indexlr's convention.
        codes_f = encode(seq)
        fwd_all, _, _ = kmer_ints(codes_f, args.k)
        for p in np.flatnonzero(take):
            strand = "+" if fwd_all[p] == canon[p] else "-"
            out.write(f"{name}\t{p}\t{int(canon[p])}\t{strand}\n")
        del codes_f, fwd_all

        nbin = len(canon) // args.bin
        if nbin:
            trim = nbin * args.bin
            per_bin = anchors[:trim].reshape(nbin, args.bin).sum(axis=1)
            comp_bin = take[:trim].reshape(nbin, args.bin).sum(axis=1)
            uniq_bin = is_unique[:trim].reshape(nbin, args.bin).sum(axis=1)
            val_bin = valid[:trim].reshape(nbin, args.bin).sum(axis=1)
            all_stats.append((per_bin, comp_bin, uniq_bin, val_bin))
            tot_bins += nbin
            tot_starved += int((per_bin < args.min_anchors).sum())

        tot_anchor += int(anchors.sum())
        tot_comp += int(take.sum())
        print(f"{name}: anchors={anchors.sum():,} complement={take.sum():,}",
              file=sys.stderr)

    out.close()

    if args.stats and all_stats:
        np.savez(args.stats,
                 anchors=np.concatenate([s[0] for s in all_stats]),
                 complement=np.concatenate([s[1] for s in all_stats]),
                 unique=np.concatenate([s[2] for s in all_stats]),
                 valid=np.concatenate([s[3] for s in all_stats]))

    print(f"\nTotal minimizer anchors : {tot_anchor:,}", file=sys.stderr)
    print(f"Total complement anchors: {tot_comp:,} "
          f"(+{100.0*tot_comp/max(tot_anchor,1):.2f}%)", file=sys.stderr)
    print(f"Bins: {tot_bins:,}  starved: {tot_starved:,} "
          f"({100.0*tot_starved/max(tot_bins,1):.2f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
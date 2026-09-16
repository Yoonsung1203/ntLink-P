#!/usr/bin/env python3
"""
Detect telomeric repeat arrays at contig termini and emit a BED mask for ntLink.

Rationale
---------
A contig terminus carrying a telomere array is, by definition, a chromosome end:
nothing should be scaffolded past it. Any join made there is a false positive.
Masking the minimizers in that terminal region removes the anchors that such
joins are built from, at no cost to real scaffolding, because the contig's
*other* end is untouched and still anchors normally.

Mask width
----------
ntLink computes gap = est_distance - a - b, where `a` is the distance from the
anchor to the contig terminus. If every surviving anchor sits at least L bases
inside the terminus, then a >= L, while est_distance <= read length. Choosing
L >= read length therefore forces gap <= -b, which ntLink rejects via
`abs(gap_est) > length_read` (ntlink_pair.py:336) or `filter_pairs_distances`
(ntlink_pair.py:265). Set --mask-width at or above the long-read N50.

Only termini where a telomere motif is actually detected are masked, so contigs
that genuinely need extension are unaffected.
"""

__author__ = "ntLink precision-scaffolding extension"

import argparse
import re
import sys

MOTIF_FWD = "TTAGGG"     # vertebrate telomere repeat
MOTIF_REV = "CCCTAA"     # reverse complement


def read_fasta(path):
    """Yield (name, sequence). Name is the first whitespace-delimited token."""
    name, chunks = None, []
    with open(path) as handle:
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


def motif_density(segment, motif_fwd=MOTIF_FWD, motif_rev=MOTIF_REV):
    """Fraction of bases in `segment` covered by either telomere motif."""
    if not segment:
        return 0.0
    covered = 0
    for motif in (motif_fwd, motif_rev):
        # non-overlapping scan is enough for tandem arrays and is O(n)
        covered += len(re.findall(motif, segment)) * len(motif)
    return min(covered / len(segment), 1.0)


def telomere_extent(seq, from_start, scan_window, density_threshold, max_scan):
    """
    Walk `scan_window`-sized blocks inward from one terminus while the telomere
    motif density stays above threshold. Returns the array length in bp
    (0 if the terminus is not telomeric).
    """
    extent = 0
    limit = min(max_scan, len(seq))
    while extent + scan_window <= limit:
        if from_start:
            block = seq[extent:extent + scan_window]
        else:
            block = seq[len(seq) - extent - scan_window:len(seq) - extent]
        if motif_density(block.upper()) < density_threshold:
            break
        extent += scan_window
    return extent


def main():
    parser = argparse.ArgumentParser(
        description="Emit a BED mask covering telomeric contig termini.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", help="Target assembly FASTA")
    parser.add_argument("-o", "--out", default="-", help="Output BED ('-' for stdout)")
    parser.add_argument("--mask-width", type=int, default=30000,
                        help="Total bases to mask inward from a telomeric terminus; "
                             "set at or above the long-read N50")
    parser.add_argument("--scan-window", type=int, default=1000,
                        help="Block size used to walk the telomere array")
    parser.add_argument("--density", type=float, default=0.4,
                        help="Minimum motif density for a block to count as telomeric")
    parser.add_argument("--max-scan", type=int, default=500000,
                        help="Maximum distance to walk inward while extending the array")
    parser.add_argument("--min-array", type=int, default=1000,
                        help="Minimum telomere array length to call a terminus telomeric")
    parser.add_argument("--report", default=None,
                        help="Optional TSV summary of detected termini")
    args = parser.parse_args()

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    report = open(args.report, "w") if args.report else None
    if report:
        report.write("contig\tlength\tend\tarray_bp\tmask_start\tmask_end\n")

    n_contigs = n_masked_termini = total_masked = 0

    for name, seq in read_fasta(args.fasta):
        n_contigs += 1
        length = len(seq)
        if length == 0:
            continue

        intervals = []
        for from_start in (True, False):
            array_bp = telomere_extent(seq, from_start, args.scan_window,
                                       args.density, args.max_scan)
            if array_bp < args.min_array:
                continue

            # mask the array itself plus enough flank to reach mask_width
            span = min(max(args.mask_width, array_bp), length)
            if from_start:
                intervals.append((0, span, "start", array_bp))
            else:
                intervals.append((length - span, length, "end", array_bp))

        # A short contig telomeric at both ends can yield overlapping intervals;
        # merge so each masked base is emitted exactly once.
        intervals.sort()
        merged = []
        for start, end, label, array_bp in intervals:
            if merged and start <= merged[-1][1]:
                prev = merged[-1]
                merged[-1] = (prev[0], max(prev[1], end),
                              prev[2] + "+" + label, max(prev[3], array_bp))
            else:
                merged.append((start, end, label, array_bp))

        for start, end, label, array_bp in merged:
            out.write(f"{name}\t{start}\t{end}\n")
            n_masked_termini += 1
            total_masked += end - start
            if report:
                report.write(f"{name}\t{length}\t{label}\t{array_bp}\t{start}\t{end}\n")

    if out is not sys.stdout:
        out.close()
    if report:
        report.close()

    print(f"Contigs scanned            : {n_contigs:,}", file=sys.stderr)
    print(f"Telomeric termini detected : {n_masked_termini:,}", file=sys.stderr)
    print(f"Total bases masked         : {total_masked:,}", file=sys.stderr)


if __name__ == "__main__":
    main()

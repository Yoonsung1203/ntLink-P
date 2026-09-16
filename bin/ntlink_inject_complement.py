#!/usr/bin/env python3
"""
Inject complement anchors into ntLink's minimizer TSVs.

ntLink matches a read minimizer to the assembly by string equality of the
`mx` token (ntlink_pair.py:208, 380). The token is opaque to ntLink, so a
complement k-mer can carry its own identifier as long as the SAME identifier
is written on both sides. We use "C<canonical k-mer integer>", which is
identical by construction wherever that k-mer occurs.

Assembly side: complement positions are already known, so entries are merged
into the target TSV.
Read side: reads are scanned and every k-mer belonging to the complement set
is emitted, with no window competition. That is the point of the exercise:
a read contributes every complement k-mer it carries error-free, instead of
only those that happen to win a minimizer window.

Formats (ntlink_pair.py):
  target TSV : contig <TAB> mx:pos:strand mx:pos:strand ...
  read   TSV : name <TAB> length <TAB> mx:pos:strand ...
"""

__author__ = "ntLink precision-scaffolding extension"

import argparse
import gzip
import sys
import numpy as np

import ntlink_complement as ncomp

PREFIX = "C"       # namespaces complement ids away from ntHash tokens


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_sequences(path):
    """Yield (name, seq) from FASTA or FASTQ, plain or gzipped."""
    with open_maybe_gzip(path) as handle:
        first = handle.readline()
        if not first:
            return
        handle.seek(0)
        if first.startswith("@"):                     # FASTQ
            while True:
                name = handle.readline().strip()
                if not name:
                    break
                seq = handle.readline().strip()
                handle.readline()
                handle.readline()
                yield name[1:].split()[0], seq
        else:                                         # FASTA
            name, chunks = None, []
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


def load_complement(path):
    """Return (set of canonical ints, {contig: [(pos, canon, strand), ...]})."""
    kmers = set()
    by_contig = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            contig, pos, canon, strand = fields[0], int(fields[1]), int(fields[2]), fields[3]
            kmers.add(canon)
            by_contig.setdefault(contig, []).append((pos, canon, strand))
    for contig in by_contig:
        by_contig[contig].sort()
    return kmers, by_contig


def augment_target(target_tsv, by_contig, out_path):
    """Merge complement entries into the target minimizer TSV, sorted by position."""
    n_added = 0
    with open(target_tsv) as fin, open(out_path, "w") as fout:
        for record in fin:
            fields = record.rstrip("\n").split("\t")
            if len(fields) < 2:
                fout.write(record)
                continue
            contig = fields[0]
            entries = []
            for token in fields[1].split(" "):
                if not token:
                    continue
                mx, pos, strand = token.split(":")
                entries.append((int(pos), mx, strand))
            for pos, canon, strand in by_contig.get(contig, []):
                entries.append((pos, f"{PREFIX}{canon}", strand))
                n_added += 1
            entries.sort()
            fout.write(contig + "\t" +
                       " ".join(f"{mx}:{pos}:{strand}" for pos, mx, strand in entries) + "\n")
    return n_added


def scan_reads(read_files, kmers, k):
    """{read name: [(pos, canon, strand), ...]} for complement k-mers present."""
    if not kmers:
        return {}
    sorted_kmers = np.sort(np.fromiter(kmers, dtype=np.uint64, count=len(kmers)))
    hits = {}
    for path in read_files:
        for name, seq in read_sequences(path):
            if len(seq) < k:
                continue
            codes = ncomp.encode(seq)
            fwd, rev, valid = ncomp.kmer_ints(codes, k)
            canon = np.minimum(fwd, rev)
            idx = np.searchsorted(sorted_kmers, canon)
            np.clip(idx, 0, len(sorted_kmers) - 1, out=idx)
            match = (sorted_kmers[idx] == canon) & valid
            positions = np.flatnonzero(match)
            if len(positions) == 0:
                continue
            hits[name] = [(int(p), int(canon[p]), "+" if fwd[p] == canon[p] else "-")
                          for p in positions]
    return hits


def augment_reads(read_tsv, hits, out_path):
    """Merge complement hits into the read minimizer TSV, sorted by position."""
    n_added = n_reads = 0
    with open_maybe_gzip(read_tsv) as fin, open(out_path, "w") as fout:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                fout.write(line)
                continue
            name = fields[0]
            entries = []
            for token in fields[2].split(" "):
                if not token:
                    continue
                mx, pos, strand = token.split(":")
                entries.append((int(pos), mx, strand))
            extra = hits.get(name, [])
            existing = {(pos, mx) for pos, mx, _ in entries}
            for pos, canon, strand in extra:
                token_id = f"{PREFIX}{canon}"
                if (pos, token_id) in existing:
                    continue
                entries.append((pos, token_id, strand))
                n_added += 1
            if extra:
                n_reads += 1
            entries.sort()
            fout.write(f"{name}\t{fields[1]}\t" +
                       " ".join(f"{mx}:{pos}:{strand}" for pos, mx, strand in entries) + "\n")
    return n_added, n_reads


def main():
    ap = argparse.ArgumentParser(
        description="Inject complement anchors into ntLink minimizer TSVs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--complement", required=True,
                    help="Complement TSV from ntlink_complement.py")
    ap.add_argument("--target-tsv", required=True, help="Target minimizer TSV from indexlr")
    ap.add_argument("--read-tsv", required=True, help="Read minimizer TSV from indexlr")
    ap.add_argument("--reads", required=True, nargs="+",
                    help="Read FASTA/FASTQ files (may be gzipped)")
    ap.add_argument("-k", type=int, default=32, help="k-mer size, must match indexlr")
    ap.add_argument("--out-target", required=True, help="Augmented target TSV")
    ap.add_argument("--out-reads", required=True, help="Augmented read TSV")
    args = ap.parse_args()

    kmers, by_contig = load_complement(args.complement)
    print(f"Complement k-mers loaded : {len(kmers):,}", file=sys.stderr)

    n_target = augment_target(args.target_tsv, by_contig, args.out_target)
    print(f"Target entries added     : {n_target:,}", file=sys.stderr)

    hits = scan_reads(args.reads, kmers, args.k)
    n_read, n_reads = augment_reads(args.read_tsv, hits, args.out_reads)
    print(f"Read entries added       : {n_read:,} across {n_reads:,} reads", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate synthetic test data for the ntLink precision pipeline.

Five scenarios, each isolating one behaviour:

  S1 clean        ctg1a/ctg1b  unique flanks, spanning reads
                  -> should join under every configuration

  S2 starved      ctg2a ends in a repeat flank with sparse unique k-mers
                  -> strict settings fail without complement, succeed with it

  S3 telomere     ctg3 ends in a TTAGGG array; chimeric reads link it to ctg4
                  -> false join unless the telomeric terminus is masked

  S4 collapsed    a 5 kb segment exists twice in the true genome but once in
                  the assembly (end of ctg5); reads from the second locus carry
                  ctg6 sequence -> false join unless that region is blacklisted

  S5 whitelist    ctg7a/ctg7b behave like S1 but sit outside the whitelist
                  -> join disappears when --mx-include-bed restricts anchors

Outputs
  target.fa          assembly (one haplotype)
  reads.fa.gz        HiFi-like reads, 0.2 % substitution error
  telomere.bed       terminus mask for S3 (what a user would curate)
  blacklist.bed      collapsed-repeat mask for S4
  whitelist.bed      anchor whitelist for S5
  allowed_kmers.txt  canonical k-mers permitted as complement anchors
  truth.tsv          expected outcome per scenario and configuration
"""

import gzip
import os
import sys
import numpy as np

K = 32
READ_LEN = 20000
READS_PER_JUNCTION = 40
ERROR_RATE = 0.002
MASK_WIDTH = 25000            # >= read length, see gap sanity argument
TELO_UNIT = "TTAGGG"

rng = np.random.default_rng(2024)
BASES = np.array(list("ACGT"))


def rand(n):
    return "".join(rng.choice(BASES, n))


def repeat_flank(n, unit_len=150, snv_rate=0.0015):
    """Tandem repeat with sparse substitutions: few unique k-mers, fewer minimizers."""
    unit = rand(unit_len)
    out = []
    for _ in range(n // unit_len):
        arr = np.array(list(unit))
        idx = np.flatnonzero(rng.random(unit_len) < snv_rate)
        if len(idx):
            arr[idx] = rng.choice(BASES, len(idx))
        out.append("".join(arr))
    return "".join(out)


def mutate(seq):
    arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    idx = np.flatnonzero(rng.random(len(arr)) < ERROR_RATE)
    if len(idx):
        arr[idx] = np.frombuffer(b"ACGT", dtype=np.uint8)[rng.integers(0, 4, len(idx))]
    return arr.tobytes().decode()


def span_reads(left, right, gap, n=READS_PER_JUNCTION, jitter=500):
    """Reads centred on the junction between two sequences."""
    joined = left + gap + right
    centre = len(left) + len(gap) // 2
    out = []
    for _ in range(n):
        start = centre - READ_LEN // 2 + int(rng.integers(-jitter, jitter))
        start = max(0, min(start, len(joined) - READ_LEN))
        out.append(mutate(joined[start:start + READ_LEN]))
    return out


contigs = {}
reads = []

# ---------------------------------------------------------------- S1 clean
c1a = rand(40000)
c1b = rand(40000)
gap1 = rand(300)
contigs["ctg1a"] = c1a
contigs["ctg1b"] = c1b
reads += span_reads(c1a, c1b, gap1)

# -------------------------------------------------------------- S2 starved
c2a = rand(30000) + repeat_flank(12000)      # starved right flank
c2b = rand(40000)
gap2 = rand(300)
contigs["ctg2a"] = c2a
contigs["ctg2b"] = c2b
reads += span_reads(c2a, c2b, gap2)

# ------------------------------------------------------------- S3 telomere
sub_telo = rand(35000)                        # subtelomeric unique sequence
telo = TELO_UNIT * 700                        # 4.2 kb array at the right end
c3 = sub_telo + telo
c4 = rand(40000)
contigs["ctg3"] = c3
contigs["ctg4"] = c4
# Chimeric reads: end of ctg3 spliced onto the start of ctg4. Nothing should be
# scaffolded past a telomere, so any ctg3->ctg4 edge is a false positive.
for _ in range(READS_PER_JUNCTION):
    left = c3[-12000:]
    right = c4[:8000]
    reads.append(mutate(left + right))

# ------------------------------------------------------------ S4 collapsed
collapsed = rand(5000)                        # present twice in the true genome
c5 = rand(35000) + collapsed                  # assembly keeps only one copy
c6 = rand(40000)
contigs["ctg5"] = c5
contigs["ctg6"] = c6
# Reads from the *other* locus: same collapsed sequence, different neighbour.
for _ in range(READS_PER_JUNCTION):
    left = collapsed
    right = c6[:15000]
    reads.append(mutate(left + right))

# ------------------------------------------------------------ S5 whitelist
c7a = rand(40000)
c7b = rand(40000)
gap7 = rand(300)
contigs["ctg7a"] = c7a
contigs["ctg7b"] = c7b
reads += span_reads(c7a, c7b, gap7)

# --------------------------------------------------------------- write data
with open("target.fa", "w") as fh:
    for name, seq in contigs.items():
        fh.write(f">{name}\n{seq}\n")

with gzip.open("reads.fa.gz", "wt") as fh:
    for i, seq in enumerate(reads):
        fh.write(f">read{i}\n{seq}\n")

# telomere mask: the terminus of ctg3 that carries the array
with open("telomere.bed", "w") as fh:
    start = max(0, len(c3) - MASK_WIDTH)
    fh.write(f"ctg3\t{start}\t{len(c3)}\n")

# collapsed-repeat mask: what read depth or meryl copy number would flag
with open("blacklist.bed", "w") as fh:
    start = len(c5) - len(collapsed)
    fh.write(f"ctg5\t{start}\t{len(c5)}\n")

# whitelist: everything except the S5 pair
with open("whitelist.bed", "w") as fh:
    for name, seq in contigs.items():
        if name.startswith("ctg7"):
            continue
        fh.write(f"{name}\t0\t{len(seq)}\n")

# combined mask a user would normally curate and pass as --exclude-bed
with open("exclude.bed", "w") as fh:
    fh.write(open("telomere.bed").read())
    fh.write(open("blacklist.bed").read())

# ------------------------------------------- allowed k-mers for complement
NTLINK_BIN = os.environ.get("NTLINK_BIN", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")))
sys.path.insert(0, NTLINK_BIN)
import ntlink_complement as nc

blocked = set()
for name, region in (("ctg3", (len(c3) - MASK_WIDTH, len(c3))),
                     ("ctg5", (len(c5) - len(collapsed), len(c5)))):
    seq = contigs[name]
    codes = nc.encode(seq)
    fwd, rev, valid = nc.kmer_ints(codes, K)
    canon = np.minimum(fwd, rev)
    lo, hi = region
    for p in range(max(0, lo), min(hi, len(canon))):
        if valid[p]:
            blocked.add(int(canon[p]))

allowed = []
for name, seq in contigs.items():
    codes = nc.encode(seq)
    fwd, rev, valid = nc.kmer_ints(codes, K)
    canon = np.minimum(fwd, rev)
    cv = canon[valid]
    _, inverse, counts = np.unique(cv, return_inverse=True, return_counts=True)
    keep = cv[counts[inverse] == 1]
    allowed.append(keep)
allowed = np.unique(np.concatenate(allowed))
allowed = np.array([v for v in allowed if int(v) not in blocked], dtype=np.uint64)
np.savetxt("allowed_kmers.txt", allowed, fmt="%d")

# k-mer level blacklist for the BASE minimizer seeding: exactly the k-mers of
# the collapsed segment and the telomere array, written as sequences the way
# `meryl print` emits them. Unlike the BED masks these name individual k-mers,
# so unique k-mers that merely sit nearby are untouched.
def kmer_seqs(seq, lo, hi):
    out = []
    for p in range(max(0, lo), min(hi, len(seq) - K + 1)):
        sub = seq[p:p + K].upper()
        if set(sub) <= set("ACGT"):
            out.append(sub)
    return out

with open("blacklist_kmers.txt", "w") as fh:
    for kmer in kmer_seqs(c5, len(c5) - len(collapsed), len(c5)):
        fh.write(kmer + "\t20\n")          # second column mimics meryl counts
    for kmer in kmer_seqs(c3, len(c3) - MASK_WIDTH, len(c3)):
        fh.write(kmer + "\t20\n")

# k-mer whitelist: every unique k-mer outside the S5 pair
with open("whitelist_kmers.txt", "w") as fh:
    for name, seq in contigs.items():
        if name.startswith("ctg7"):
            continue
        codes = nc.encode(seq)
        fwd, rev, valid = nc.kmer_ints(codes, K)
        canon_all = np.minimum(fwd, rev)
        cv = canon_all[valid]
        _, inv, cnt = np.unique(cv, return_inverse=True, return_counts=True)
        for value in np.unique(cv[cnt[inv] == 1]):
            if int(value) not in blocked:      # also drop telomere / collapsed k-mers
                fh.write(f"{int(value)}\n")

# ------------------------------------------------------------------- truth
with open("truth.tsv", "w") as fh:
    fh.write("scenario\tedge\ttruth\tno_masks\twith_masks\twith_masks_complement\n")
    fh.write("S1 clean\tctg1a-ctg1b\ttrue join\tjoin\tjoin\tjoin\n")
    fh.write("S2 starved\tctg2a-ctg2b\ttrue join\tmiss (strict)\tmiss (strict)\tjoin\n")
    fh.write("S3 telomere\tctg3-ctg4\tFALSE join\tjoin (bad)\tblocked\tblocked\n")
    fh.write("S4 collapsed\tctg5-ctg6\tFALSE join\tjoin (bad)\tblocked\tblocked\n")
    fh.write("S5 whitelist\tctg7a-ctg7b\ttrue join\tjoin\tjoin\tdropped if whitelist used\n")

total = sum(len(s) for s in contigs.values())
print(f"contigs        : {len(contigs)}  total {total:,} bp")
for name, seq in contigs.items():
    print(f"  {name:8} {len(seq):>8,} bp")
print(f"reads          : {len(reads)} x {READ_LEN:,} bp, error {ERROR_RATE}")
print(f"telomere mask  : ctg3 last {MASK_WIDTH:,} bp")
print(f"blacklist      : ctg5 last {len(collapsed):,} bp (collapsed repeat)")
print(f"whitelist      : everything except ctg7a/ctg7b")
print(f"allowed k-mers : {len(allowed):,}")

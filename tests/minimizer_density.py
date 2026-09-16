#!/usr/bin/env python3
"""
CHM13 chr21 analysis: how much unique-k-mer anchoring capacity does
minimizer sampling (k=32, w=250) lose, as a function of local repeat content?

Definitions used (matching ntLink semantics):
  - k-mer uniqueness: canonical 32-mer occurring exactly once in the sequence
    analysed. ntLink's read_minimizers() drops any minimizer seen >=2 times in
    the target, so "usable anchor" == minimizer whose k-mer count == 1.
  - minimizer: per window of w consecutive k-mers, the one with minimum hash.

Caveat: uniqueness is computed on chr21 alone, not the whole genome, so it
OVERSTATES uniqueness. Both the sampled and the all-k-mer track share this
bias, so their ratio (the quantity of interest) is largely unaffected.
"""
import sys
import numpy as np

K = 32
W = 250
WINDOW_BP = 10000   # analysis bin size

# ---------------------------------------------------------------- load
def load_fasta(path):
    chunks = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                chunks.append(line.strip())
    return "".join(chunks).upper()


def encode(seq):
    """2-bit encode; non-ACGT -> 255 sentinel."""
    lut = np.full(256, 255, dtype=np.uint8)
    for base, code in zip(b"ACGT", [0, 1, 2, 3]):
        lut[base] = code
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    return lut[arr]


def kmer_ints(codes, k=K):
    """
    Rolling 2-bit packed k-mers (k=32 fits exactly in uint64).
    Returns (fwd, rev, valid) arrays of length n-k+1.
    """
    n = len(codes)
    m = n - k + 1
    valid_base = codes != 255
    safe = np.where(valid_base, codes, 0).astype(np.uint64)

    fwd = np.zeros(m, dtype=np.uint64)
    rev = np.zeros(m, dtype=np.uint64)
    # build by shifting in each of the k positions (vectorised over all windows)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | safe[i:i + m]
        # reverse complement: complement base, placed in reverse order
        comp = (np.uint64(3) - safe[i:i + m])
        rev = rev | (comp << np.uint64(2 * i))

    # validity: window contains no sentinel
    cs = np.concatenate(([0], np.cumsum(~valid_base)))
    valid = (cs[k:] - cs[:-k]) == 0
    return fwd, rev, valid


def splitmix64(x):
    """Vectorised finalizer; stands in for ntHash (uniformity is what matters)."""
    x = x.copy()
    x += np.uint64(0x9E3779B97F4A7C15)
    z = x
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


def sliding_argmin_mask(h, w):
    """
    Boolean mask: position i is a minimizer if h[i] is the minimum of at least
    one window of w consecutive values. Monotonic-deque sweep, O(n).
    """
    n = len(h)
    picked = np.zeros(n, dtype=bool)
    if n < w:
        if n:
            picked[int(np.argmin(h))] = True
        return picked
    from collections import deque
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


def main(fasta):
    print("Loading sequence ...", flush=True)
    seq = load_fasta(fasta)
    n = len(seq)
    print(f"  length = {n:,} bp")

    codes = encode(seq)
    del seq
    print("Encoding k-mers ...", flush=True)
    fwd, rev, valid = kmer_ints(codes)
    del codes
    canon = np.minimum(fwd, rev)
    del fwd, rev
    print(f"  k-mer positions = {len(canon):,} (valid: {valid.sum():,})")

    print("Counting k-mer multiplicity (whole chr21) ...", flush=True)
    cv = canon[valid]
    uniq_vals, inverse, counts = np.unique(cv, return_inverse=True, return_counts=True)
    is_unique_valid = counts[inverse] == 1
    del uniq_vals, inverse, counts, cv
    is_unique = np.zeros(len(canon), dtype=bool)
    is_unique[valid] = is_unique_valid
    del is_unique_valid
    print(f"  unique (count==1) k-mers: {is_unique.sum():,} "
          f"({100.0*is_unique.sum()/valid.sum():.2f}% of valid)")

    print("Hashing + minimizer selection (w=%d) ..." % W, flush=True)
    h = splitmix64(canon)
    del canon
    h[~valid] = np.iinfo(np.uint64).max      # never pick invalid windows
    picked = sliding_argmin_mask(h, W)
    del h
    picked &= valid
    print(f"  minimizers selected: {picked.sum():,} "
          f"(density {picked.sum()/valid.sum():.4f}, theory {2/(W+1):.4f})")

    seeded = picked & is_unique      # ntLink-usable anchors
    print(f"  of which unique (ntLink anchors): {seeded.sum():,}")
    print(f"  minimizers dropped by dedup: {picked.sum()-seeded.sum():,} "
          f"({100.0*(picked.sum()-seeded.sum())/picked.sum():.2f}%)")

    # ------------------------------------------------ per-window statistics
    print("Binning ...", flush=True)
    m = len(picked)
    nbin = m // WINDOW_BP
    trim = nbin * WINDOW_BP
    seeded_bin = seeded[:trim].reshape(nbin, WINDOW_BP).sum(axis=1)
    unique_bin = is_unique[:trim].reshape(nbin, WINDOW_BP).sum(axis=1)
    valid_bin = valid[:trim].reshape(nbin, WINDOW_BP).sum(axis=1)

    np.save("seeded_bin.npy", seeded_bin)
    np.save("unique_bin.npy", unique_bin)
    np.save("valid_bin.npy", valid_bin)
    print(f"  bins = {nbin:,} of {WINDOW_BP} bp")

    ok = valid_bin > WINDOW_BP * 0.9
    uf = np.zeros(nbin)
    uf[ok] = unique_bin[ok] / valid_bin[ok]     # unique-kmer fraction (repeat proxy)

    print("\n" + "=" * 78)
    print("Minimizer anchor density vs local unique-k-mer fraction")
    print("(unique fraction low = repeat-rich)")
    print("=" * 78)
    print(f"{'unique-kmer frac':>18} | {'bins':>7} | {'anchors/10kb':>12} | "
          f"{'unique kmers/10kb':>18} | {'captured':>9} | {'bins<2 anchors':>14}")
    print("-" * 78)

    edges = [0.0, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 1.001]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = ok & (uf >= lo) & (uf < hi)
        cnt = int(sel.sum())
        if cnt == 0:
            continue
        a = seeded_bin[sel].mean()
        u = unique_bin[sel].mean()
        frac_low = 100.0 * (seeded_bin[sel] < 2).mean()
        cap = (a / u * 100) if u > 0 else float("nan")
        rows.append((lo, hi, cnt, a, u, cap, frac_low))
        print(f"{lo:>8.2f}-{hi:<8.2f} | {cnt:>7,} | {a:>12.1f} | "
              f"{u:>18,.0f} | {cap:>8.2f}% | {frac_low:>13.1f}%")

    print("-" * 78)
    tot_ok = int(ok.sum())
    print(f"{'ALL':>18} | {tot_ok:>7,} | {seeded_bin[ok].mean():>12.1f} | "
          f"{unique_bin[ok].mean():>18,.0f} | "
          f"{100.0*seeded_bin[ok].sum()/unique_bin[ok].sum():>8.2f}% | "
          f"{100.0*(seeded_bin[ok]<2).mean():>13.1f}%")

    # ---------------- the decision-relevant quantity -----------------
    print("\n" + "=" * 78)
    print("COMPLEMENT POTENTIAL: bins with <2 anchors but unique k-mers available")
    print("=" * 78)
    starved = ok & (seeded_bin < 2)
    recoverable = starved & (unique_bin >= 2)
    print(f"  bins analysed                         : {tot_ok:,}")
    print(f"  bins with <2 minimizer anchors        : {int(starved.sum()):,} "
          f"({100.0*starved.sum()/tot_ok:.2f}%)")
    print(f"  ...of which have >=2 unique k-mers    : {int(recoverable.sum()):,} "
          f"({100.0*recoverable.sum()/tot_ok:.2f}% of all bins)")
    if starved.sum():
        print(f"  ...i.e. {100.0*recoverable.sum()/starved.sum():.1f}% of starved bins "
              f"are recoverable by complementing")
        print(f"  truly anchor-less (no unique k-mer)   : "
              f"{int((starved & (unique_bin < 2)).sum()):,}")
    if recoverable.sum():
        print(f"  median unique k-mers in recoverable bins: "
              f"{np.median(unique_bin[recoverable]):,.0f}")

    return rows


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "chr21.fa")

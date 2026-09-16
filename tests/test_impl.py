#!/usr/bin/env python3
"""Validate encoding, canonical k-mers, and sliding-window minimizer selection."""
import sys
import numpy as np
import os
NTLINK_ROOT = os.environ.get("NTLINK_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
NTLINK_BIN = os.path.join(NTLINK_ROOT, "bin")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import minimizer_density as md

fails = []
def check(name, got, exp):
    ok = got == exp
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: got={got} expected={exp}")
    if not ok:
        fails.append(name)

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# --- canonical k-mer correctness on a short sequence, k=4 (brute force) ---
md_K = 4
seq = "ACGTTTGACAGTNACGTACGT"
codes = md.encode(seq)
fwd, rev, valid = md.kmer_ints(codes, k=md_K)

def pack(s):
    v = 0
    for ch in s:
        v = (v << 2) | "ACGT".index(ch)
    return v

exp_canon, exp_valid = [], []
for i in range(len(seq) - md_K + 1):
    sub = seq[i:i+md_K]
    if set(sub) <= set("ACGT"):
        exp_canon.append(min(pack(sub), pack(revcomp(sub))))
        exp_valid.append(True)
    else:
        exp_canon.append(None)
        exp_valid.append(False)

got_canon = np.minimum(fwd, rev)
check("validity mask matches brute force", list(valid), exp_valid)
check("canonical k-mers match brute force",
      [int(got_canon[i]) if exp_valid[i] else None for i in range(len(exp_valid))],
      exp_canon)

# --- canonical invariance: revcomp of sequence gives same canonical multiset ---
s2 = "ACGTTTGACAGTACGTACGTAGGCTA"
c1 = np.minimum(*md.kmer_ints(md.encode(s2), k=md_K)[:2])
c2 = np.minimum(*md.kmer_ints(md.encode(revcomp(s2)), k=md_K)[:2])
check("canonical set invariant under revcomp",
      sorted(c1.tolist()), sorted(c2.tolist()))

# --- sliding window minimizer: brute force comparison ---
rng = np.random.default_rng(0)
h = rng.integers(0, 10**6, size=500).astype(np.uint64)
for w in (1, 2, 5, 37, 250):
    got = md.sliding_argmin_mask(h, w)
    exp = np.zeros(len(h), dtype=bool)
    for start in range(0, len(h) - w + 1):
        exp[start + int(np.argmin(h[start:start+w]))] = True
    check(f"sliding minimizer mask w={w}", got.tolist(), exp.tolist())

# --- density sanity: random sequence should approach 2/(w+1) ---
rng2 = np.random.default_rng(1)
rand_seq = "".join(rng2.choice(list("ACGT"), size=300000))
f, r, v = md.kmer_ints(md.encode(rand_seq), k=32)
hh = md.splitmix64(np.minimum(f, r))
for w in (100, 250):
    picked = md.sliding_argmin_mask(hh, w)
    dens = picked.sum() / len(picked)
    theory = 2 / (w + 1)
    ok = abs(dens - theory) / theory < 0.10
    print(f"[{'PASS' if ok else 'FAIL'}] density w={w}: {dens:.5f} vs theory {theory:.5f}")
    if not ok:
        fails.append(f"density w={w}")

print()
if fails:
    print("FAILED:", fails); sys.exit(1)
print("All implementation tests passed.")

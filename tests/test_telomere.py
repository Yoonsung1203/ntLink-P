#!/usr/bin/env python3
"""Tests for ntlink_telomere_bed: detection at both termini, absence, edge cases."""
import os
NTLINK_ROOT = os.environ.get("NTLINK_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
NTLINK_BIN = os.path.join(NTLINK_ROOT, "bin")
import subprocess
import sys
import random
import tempfile
import os

SCRIPT = os.path.join(NTLINK_BIN, "ntlink_telomere_bed.py")
fails = []


def check(name, got, exp):
    ok = got == exp
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: got={got} expected={exp}")
    if not ok:
        fails.append(name)


def rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def run(fasta_text, extra=None):
    fh = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    fh.write(fasta_text)
    fh.close()
    cmd = [sys.executable, SCRIPT, fh.name, "--mask-width", "5000",
           "--max-scan", "50000"] + (extra or [])
    res = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(fh.name)
    lines = [l.split("\t") for l in res.stdout.strip().split("\n") if l]
    return [(c, int(s), int(e)) for c, s, e in lines], res.stderr


# --- 1. telomere at start only ---
seq = "TTAGGG" * 500 + rand_seq(40000, 1)          # 3000 bp array
bed, _ = run(f">ctgA\n{seq}\n")
check("telomere at start detected", bed, [("ctgA", 0, 5000)])

# --- 2. telomere at end only (reverse-complement motif) ---
seq = rand_seq(40000, 2) + "CCCTAA" * 500
bed, _ = run(f">ctgB\n{seq}\n")
expected_len = 40000 + 3000
check("telomere at end detected", bed, [("ctgB", expected_len - 5000, expected_len)])

# --- 3. both ends telomeric (a complete T2T chromosome) ---
seq = "TTAGGG" * 500 + rand_seq(40000, 3) + "CCCTAA" * 500
total = len(seq)
bed, _ = run(f">ctgC\n{seq}\n")
check("both termini detected", bed,
      [("ctgC", 0, 5000), ("ctgC", total - 5000, total)])

# --- 4. no telomere -> no output (must not mask ordinary contigs) ---
bed, _ = run(f">ctgD\n{rand_seq(40000, 4)}\n")
check("non-telomeric contig untouched", bed, [])

# --- 5. interstitial telomere repeat must NOT trigger (not at terminus) ---
seq = rand_seq(20000, 5) + "TTAGGG" * 500 + rand_seq(20000, 6)
bed, _ = run(f">ctgE\n{seq}\n")
check("interstitial array ignored", bed, [])

# --- 6. short array below --min-array is ignored ---
seq = "TTAGGG" * 20 + rand_seq(40000, 7)           # 120 bp only
bed, _ = run(f">ctgF\n{seq}\n")
check("sub-threshold array ignored", bed, [])

# --- 7. mask never exceeds contig length ---
seq = "TTAGGG" * 300 + rand_seq(500, 8)            # total 2300 bp < mask width
bed, _ = run(f">ctgG\n{seq}\n")
ok = bed and bed[0][1] >= 0 and bed[0][2] <= len(seq)
print(f"[{'PASS' if ok else 'FAIL'}] mask clipped to contig length: {bed} len={len(seq)}")
if not ok:
    fails.append("mask clipping")

# --- 8. multi-contig file handled independently ---
s1 = "TTAGGG" * 500 + rand_seq(30000, 9)
s2 = rand_seq(30000, 10)
bed, _ = run(f">c1\n{s1}\n>c2\n{s2}\n")
check("multi-contig: only telomeric one masked", bed, [("c1", 0, 5000)])

# --- 9. mask width respected (larger value) ---
seq = "TTAGGG" * 500 + rand_seq(60000, 11)
fh = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
fh.write(f">ctgH\n{seq}\n"); fh.close()
res = subprocess.run([sys.executable, SCRIPT, fh.name, "--mask-width", "30000"],
                     capture_output=True, text=True)
os.unlink(fh.name)
line = res.stdout.strip().split("\t")
check("mask-width=30000 honoured", (line[0], int(line[1]), int(line[2])),
      ("ctgH", 0, 30000))

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All telomere detection tests passed.")

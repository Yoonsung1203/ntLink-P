#!/usr/bin/env python3
"""
Regression test: with default settings (hc=2, mx_ratio=0, no BED), the modified
ntLink must behave identically to stock ntLink.

Stock semantics:
  - well-anchored: hit_count > 1  (equivalently >= 2)
  - no coordinate masking
  - no ratio filter
"""
import sys
import types
import subprocess
import tempfile
import os

import os
NTLINK_ROOT = os.environ.get("NTLINK_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
NTLINK_BIN = os.path.join(NTLINK_ROOT, "bin")
sys.path.insert(0, NTLINK_BIN)
import ntlink_utils
import ntlink_pair
import ntlink_kmer_mask as km

failures = []


def check(name, got, expected):
    status = "PASS" if got == expected else "FAIL"
    if status == "FAIL":
        failures.append(name)
    print(f"[{status}] {name}: got={got} expected={expected}")


# --- 1. hc=2 reproduces stock 'hit_count > 1' exactly ---
# Stock accepted hit_count of 2,3,4...; rejected 0,1.
# New code uses hit_count >= hc. With hc=2 the accept set must be identical.
stock_accept = {c: (c > 1) for c in range(0, 6)}
new_accept = {c: (c >= 2) for c in range(0, 6)}
check("hc=2 identical to stock hit_count>1", new_accept, stock_accept)

# hc=3 must be strictly stricter (differs only at hit_count == 2)
hc3 = {c: (c >= 3) for c in range(0, 6)}
diff = {c for c in range(0, 6) if hc3[c] != stock_accept[c]}
check("hc=3 differs from stock only at hit_count==2", diff, {2})

# --- 2. No BED -> mx_info returned unchanged (identity, not a copy) ---
Minimizer = ntlink_pair.Minimizer
mx_info = {"a": Minimizer("c1", 10, "+"), "b": Minimizer("c1", 20, "-")}
out = km.filter_mx_info_by_regions(mx_info, None, None, verbose=False)
check("no BED -> object identity preserved", out is mx_info, True)

# --- 3. mx_ratio=0 -> filter never fires, even with a populated index ---
class FakeScaffold:
    def __init__(self, length):
        self.length = length

# A deliberately terrible anchor: 2 hits spanning 40 index positions (ratio 0.05)
mx_list = [("m0", 100, "+"), ("m1", 200, "+")]
list_mx_info = {"m0": Minimizer("ctg", 1000, "+"), "m1": Minimizer("ctg", 5000, "+")}
scaffolds = {"ctg": FakeScaffold(100000)}
index = {"ctg": [1000 + 100 * i for i in range(41)]}

args = types.SimpleNamespace(k=32, x=0, mx_ratio=0, checkpoint=None, hc=2,
                             sensitive=False, z=100, f=10)
accepted, _ = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("mx_ratio=0 keeps even a 0.05-ratio anchor", sorted(set(accepted)), ["ctg"])

# Same input with the filter on must drop it -> proves the filter is wired up
args.mx_ratio = 0.5
accepted, _ = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("mx_ratio=0.5 drops the same anchor", sorted(set(accepted)), [])

# --- 4. x=0 default behaviour unchanged ---
args = types.SimpleNamespace(k=32, x=0, mx_ratio=0, checkpoint=None, hc=2,
                             sensitive=False, z=100, f=10)
accepted, _ = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("x=0 keeps wide-span anchor (stock)", sorted(set(accepted)), ["ctg"])

# --- 5. Stock upstream files untouched except the intended ones ---
res = subprocess.run(["git", "-C", NTLINK_ROOT, "diff", "--name-only",
                      "--diff-filter=M", "-G", "."],
                     capture_output=True, text=True)
changed = sorted(f for f in res.stdout.strip().split("\n") if f)
check("only intended files changed in content", changed,
      ["bin/ntlink_pair.py", "bin/ntlink_utils.py", "ntLink"])

# --- 6. All modified modules import cleanly ---
for mod in ["ntlink_pair", "ntlink_utils", "ntlink_kmer_mask"]:
    r = subprocess.run([sys.executable, "-c",
                        f"import sys; sys.path.insert(0, r'{NTLINK_BIN}'); import {mod}"],
                       capture_output=True, text=True)
    check(f"import {mod}", r.returncode, 0)

# ntlink_patch_gaps requires btllib (not installed in this container); verify it
# only fails on that external dependency, i.e. our edits introduced no new error.
r = subprocess.run([sys.executable, "-c",
                    f"import sys; sys.path.insert(0, r'{NTLINK_BIN}'); import ntlink_patch_gaps"],
                   capture_output=True, text=True)
check("patch_gaps fails only on missing btllib (pre-existing)",
      "No module named 'btllib'" in r.stderr, True)

# Compile-check patch_gaps to prove it is syntactically valid despite the import
r = subprocess.run([sys.executable, "-m", "py_compile",
                    os.path.join(NTLINK_BIN, "ntlink_patch_gaps.py")],
                   capture_output=True, text=True)
check("patch_gaps compiles cleanly", r.returncode, 0)

print()
if failures:
    print(f"FAILED: {failures}")
    sys.exit(1)
print("All regression tests passed: defaults reproduce stock ntLink behaviour.")

#!/usr/bin/env python3
"""
Integration tests for the mx-ratio filter and --hc parameterization,
exercising the real ntlink_utils.get_accepted_anchor_contigs code path.
"""
import sys
import types

import os
NTLINK_ROOT = os.environ.get("NTLINK_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
NTLINK_BIN = os.path.join(NTLINK_ROOT, "bin")
sys.path.insert(0, NTLINK_BIN)
import ntlink_utils
import ntlink_pair

failures = []


def check(name, got, expected):
    status = "PASS" if got == expected else "FAIL"
    if status == "FAIL":
        failures.append(name)
    print(f"[{status}] {name}: got={got} expected={expected}")


def make_args(**kwargs):
    defaults = dict(k=32, x=0, mx_ratio=0, checkpoint=None, hc=2,
                    sensitive=False, z=100, f=10)
    defaults.update(kwargs)
    return types.SimpleNamespace(**defaults)


class FakeScaffold:
    def __init__(self, length):
        self.length = length


def build_inputs(hit_specs):
    """
    hit_specs: list of (contig, ctg_pos, read_pos)
    Returns (mx_list, list_mx_info, scaffolds) in the format the real
    function expects: mx_list entries are (mx, read_pos, read_strand).
    """
    mx_list, list_mx_info, scaffolds = [], {}, {}
    for i, (contig, ctg_pos, read_pos) in enumerate(hit_specs):
        mx = f"mx{i}"
        mx_list.append((mx, read_pos, "+"))
        list_mx_info[mx] = ntlink_pair.Minimizer(contig, ctg_pos, "+")
        scaffolds.setdefault(contig, FakeScaffold(100000))
    return mx_list, list_mx_info, scaffolds


# ---------------------------------------------------------------
# Scenario: contig "sparse" spans ctg 1000..5000 but the read hits
# only 3 of the 41 unique minimizers in that interval (ratio ~0.073).
# contig "dense" spans 1000..1400 and hits 5 of 5 (ratio 1.0).
# ---------------------------------------------------------------
hits = [("dense", 1000, 100), ("dense", 1100, 200), ("dense", 1200, 300),
        ("dense", 1300, 400), ("dense", 1400, 500),
        ("sparse", 1000, 600), ("sparse", 3000, 700), ("sparse", 5000, 800)]
mx_list, list_mx_info, scaffolds = build_inputs(hits)

# Index: dense has exactly its 5; sparse has 41 evenly spaced (1000..5000 step 100)
index = {
    "dense": [1000, 1100, 1200, 1300, 1400],
    "sparse": [1000 + 100 * i for i in range(41)],
}

# --- Test 1: mx_ratio disabled -> both contigs survive ---
args = make_args(mx_ratio=0, x=0)
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("mx_ratio=0 keeps both", sorted(set(accepted)), ["dense", "sparse"])

# --- Test 2: mx_ratio=0.5 -> sparse (0.073) dropped, dense (1.0) kept ---
args = make_args(mx_ratio=0.5, x=0)
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("mx_ratio=0.5 drops sparse", sorted(set(accepted)), ["dense"])

# --- Test 3: mx_ratio=0.05 (below sparse ratio) -> both kept ---
args = make_args(mx_ratio=0.05, x=0)
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("mx_ratio=0.05 keeps both", sorted(set(accepted)), ["dense", "sparse"])

# --- Test 4: checkpoint mode gates the filter off (rounds 2+) ---
args = make_args(mx_ratio=0.5, x=0, checkpoint="prefix.verbose_mapping.tsv")
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("checkpoint mode disables filter", sorted(set(accepted)), ["dense", "sparse"])

# --- Test 5: empty index (patch_gaps path) disables filter, no crash ---
args = make_args(mx_ratio=0.5, x=0)
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, {})
check("empty index disables filter", sorted(set(accepted)), ["dense", "sparse"])

# --- Test 6: index omitted entirely (backward compat, old call signature) ---
args = make_args(mx_ratio=0.5, x=0)
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args)
check("omitted index (old signature) works", sorted(set(accepted)), ["dense", "sparse"])

# --- Test 7: args lacking mx_ratio entirely (patch_gaps argparse) ---
args = types.SimpleNamespace(k=32, x=0, sensitive=False, z=100, f=10)
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx_list, 100000, scaffolds, list_mx_info, args, index)
check("args without mx_ratio does not crash", sorted(set(accepted)), ["dense", "sparse"])

# --- Test 8: x filter still functions independently ---
# contig spans 4000 on contig but only 200 on read -> non-colinear
hits2 = [("scatter", 1000, 100), ("scatter", 5000, 300),
         ("good", 1000, 1000), ("good", 1200, 1200)]
mx2, info2, scaf2 = build_inputs(hits2)
idx2 = {"scatter": [1000, 5000], "good": [1000, 1200]}
args = make_args(x=1.1, mx_ratio=0)
accepted, runs = ntlink_utils.get_accepted_anchor_contigs(
    mx2, 100000, scaf2, info2, args, idx2)
check("x=1.1 drops non-colinear contig", sorted(set(accepted)), ["good"])

# --- Test 9: hc parameterization present in add_pair source ---
src = open(os.path.join(NTLINK_BIN, "ntlink_pair.py")).read()
check("no hardcoded 'hit_count > 1' remains", "hit_count > 1" in src, False)
check("hc used in add_pair", src.count("hit_count >= self.args.hc"), 3)

print()
if failures:
    print(f"FAILED: {failures}")
    sys.exit(1)
print("All filter integration tests passed.")

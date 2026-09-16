#!/usr/bin/env python3
"""Unit tests for ntlink_kmer_mask: BED parsing, interval lookup, masking."""
import sys
import os
import tempfile
from collections import namedtuple

import os
NTLINK_ROOT = os.environ.get("NTLINK_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
NTLINK_BIN = os.path.join(NTLINK_ROOT, "bin")
sys.path.insert(0, NTLINK_BIN)
import ntlink_kmer_mask as km

Minimizer = namedtuple("Minimizer", ["contig", "position", "strand"])

failures = []


def check(name, got, expected):
    status = "PASS" if got == expected else "FAIL"
    if status == "FAIL":
        failures.append(name)
    print(f"[{status}] {name}: got={got} expected={expected}")


def write_bed(content):
    fh = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
    fh.write(content)
    fh.close()
    return fh.name


# --- Test 1: BED parsing with merging of overlapping/abutting intervals ---
bed = write_bed(
    "# comment line\n"
    "track name=test\n"
    "ctg1\t100\t200\n"
    "ctg1\t150\t250\n"      # overlaps previous -> merge to 100-250
    "ctg1\t250\t300\n"      # abuts -> merge to 100-300
    "ctg1\t500\t600\n"      # separate
    "ctg2\t0\t50\n"
    "ctg3\t10\t10\n"        # zero-length, skipped
    "ctg3\t30\t20\n"        # end<start, skipped
    "malformed_line\n"
)
regions = km.read_bed_regions(bed)
check("merge overlapping+abutting", regions["ctg1"], ([100, 500], [300, 600]))
check("second contig parsed", regions["ctg2"], ([0], [50]))
check("degenerate intervals skipped", "ctg3" in regions, False)

# --- Test 2: interval membership, half-open [start, end) ---
check("pos inside first interval", km.in_regions(regions, "ctg1", 150), True)
check("pos at start boundary (inclusive)", km.in_regions(regions, "ctg1", 100), True)
check("pos at end boundary (exclusive)", km.in_regions(regions, "ctg1", 300), False)
check("pos just before end", km.in_regions(regions, "ctg1", 299), True)
check("pos in gap between intervals", km.in_regions(regions, "ctg1", 400), False)
check("pos before all intervals", km.in_regions(regions, "ctg1", 50), False)
check("unknown contig", km.in_regions(regions, "ctgX", 150), False)
check("pos 0 in zero-start interval", km.in_regions(regions, "ctg2", 0), True)

# --- Test 3: exclude (blacklist) semantics ---
mx_info = {
    "h1": Minimizer("ctg1", 150, "+"),   # inside excluded -> dropped
    "h2": Minimizer("ctg1", 400, "+"),   # in gap -> kept
    "h3": Minimizer("ctg1", 550, "-"),   # inside excluded -> dropped
    "h4": Minimizer("ctg2", 10, "+"),    # inside excluded -> dropped
    "h5": Minimizer("ctgX", 999, "+"),   # contig not in BED -> kept
}
out = km.filter_mx_info_by_regions(mx_info, exclude_bed=bed, verbose=False)
check("exclude keeps correct set", sorted(out.keys()), ["h2", "h5"])

# --- Test 4: include (whitelist) semantics ---
out = km.filter_mx_info_by_regions(mx_info, include_bed=bed, verbose=False)
check("include keeps correct set", sorted(out.keys()), ["h1", "h3", "h4"])

# --- Test 5: no BED given -> identity (same object, no copy) ---
out = km.filter_mx_info_by_regions(mx_info, verbose=False)
check("no bed returns input unchanged", out is mx_info, True)

# --- Test 6: empty BED file -> nothing excluded ---
empty_bed = write_bed("")
out = km.filter_mx_info_by_regions(mx_info, exclude_bed=empty_bed, verbose=False)
check("empty bed drops nothing", len(out), 5)

# --- Test 7: whitespace-separated (non-tab) BED still parses ---
ws_bed = write_bed("ctg1 100 200\n")
check("whitespace-delimited bed", km.read_bed_regions(ws_bed)["ctg1"], ([100], [200]))

# --- Test 8: combined include + exclude ---
inc = write_bed("ctg1\t0\t1000\n")
exc = write_bed("ctg1\t100\t300\n")
out = km.filter_mx_info_by_regions(mx_info, exclude_bed=exc, include_bed=inc, verbose=False)
check("include then exclude", sorted(out.keys()), ["h2", "h3"])

for f in [bed, empty_bed, ws_bed, inc, exc]:
    os.unlink(f)

print()
if failures:
    print(f"FAILED: {failures}")
    sys.exit(1)
print("All ntlink_kmer_mask tests passed.")

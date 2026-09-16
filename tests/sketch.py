#!/usr/bin/env python3
"""Stand-in for `indexlr --long --pos --strand`: emit ntLink-format minimizer TSVs."""
import sys, numpy as np
import os
sys.path.insert(0, os.environ.get("NTLINK_BIN", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))))
import ntlink_complement as nc
import ntlink_inject_complement as inj

fasta, k, w, mode, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
with open(out,"w") as fh:
    for name, seq in inj.read_sequences(fasta):
        codes = nc.encode(seq)
        fwd, rev, valid = nc.kmer_ints(codes, k)
        canon = np.minimum(fwd, rev)
        h = nc.splitmix64(canon); h[~valid] = np.iinfo(np.uint64).max
        picked = nc.minimizer_mask(h, w) & valid
        toks = [f"{int(h[p])}:{p}:{'+' if fwd[p]==canon[p] else '-'}"
                for p in np.flatnonzero(picked)]
        if mode == "target":
            fh.write(f"{name}\t{' '.join(toks)}\n")
        else:
            fh.write(f"{name}\t{len(seq)}\t{' '.join(toks)}\n")

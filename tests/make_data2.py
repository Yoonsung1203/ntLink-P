#!/usr/bin/env python3
"""Redesign: flank must be genuinely starved right up to contig A's terminus."""
import sys, numpy as np
import os
sys.path.insert(0, os.environ.get("NTLINK_BIN", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))))
import ntlink_complement as nc
rng=np.random.default_rng(11); B=np.array(list("ACGT"))
def rand(n): return "".join(rng.choice(B,n))

def starved_flank(n, unit=150, snv=0.0015):
    """Tandem repeat, very sparse SNVs -> some unique k-mers, few minimizers."""
    u=rand(unit); out=[]
    for _ in range(n//unit):
        a=np.array(list(u)); idx=np.flatnonzero(rng.random(unit)<snv)
        if len(idx): a[idx]=rng.choice(B,len(idx))
        out.append("".join(a))
    return "".join(out)

coreA=rand(30000); flankA=starved_flank(12000)
GAP=rand(300); coreB=rand(30000)
A=coreA+flankA; Bc=coreB; genome=A+GAP+Bc
open("target2.fa","w").write(f">ctgA\n{A}\n>ctgB\n{Bc}\n")

c,u,p,v=nc.analyse_contig(A,32,250); anc=p&u
pos=np.flatnonzero(anc)
print(f"ctgA len={len(A)} anchors={len(pos)} last_anchor={pos.max()} overhang_a={len(A)-pos.max()}")
print(f"flank unique kmers={u[30000:].sum()}  flank anchors={(pos>=30000).sum()}")

READ,N,ERR=20000,40,0.002
with open("reads2.fa","w") as f:
    for i in range(N):
        centre=len(A)+len(GAP)//2
        start=max(0,min(centre-READ//2+int(rng.integers(-500,500)),len(genome)-READ))
        s=np.frombuffer(genome[start:start+READ].encode(),dtype=np.uint8).copy()
        idx=np.flatnonzero(rng.random(len(s))<ERR)
        if len(idx): s[idx]=np.frombuffer(b"ACGT",dtype=np.uint8)[rng.integers(0,4,len(idx))]
        f.write(f">read{i}\n{s.tobytes().decode()}\n")
print(f"reads: {N} x {READ}")

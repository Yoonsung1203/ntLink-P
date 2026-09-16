#!/usr/bin/env python3
"""Correctness tests for complement selection."""
import sys, numpy as np, random
import os
NTLINK_ROOT = os.environ.get("NTLINK_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
NTLINK_BIN = os.path.join(NTLINK_ROOT, "bin")
sys.path.insert(0,NTLINK_BIN)
import ntlink_complement as nc

fails=[]
def check(name, got, exp):
    ok = got==exp
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: got={got} expected={exp}")
    if not ok: fails.append(name)

# --- canonical kmer vs brute force (k=4) ---
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def pack(s):
    v=0
    for ch in s: v=(v<<2)|"ACGT".index(ch)
    return v
seq="ACGTTTGACAGTNACGTACGT"
f,r,v = nc.kmer_ints(nc.encode(seq),4)
canon=np.minimum(f,r)
exp=[min(pack(seq[i:i+4]),pack(rc(seq[i:i+4]))) if set(seq[i:i+4])<=set("ACGT") else None
     for i in range(len(seq)-3)]
check("canonical kmers", [int(canon[i]) if v[i] else None for i in range(len(v))], exp)

# --- minimizer mask vs brute force ---
rng=np.random.default_rng(0)
h=rng.integers(0,10**6,size=400).astype(np.uint64)
for w in (3,10,50):
    got=nc.minimizer_mask(h,w)
    e=np.zeros(len(h),dtype=bool)
    for s in range(len(h)-w+1): e[s+int(np.argmin(h[s:s+w]))]=True
    check(f"minimizer mask w={w}", got.tolist(), e.tolist())

# --- select_complement: starved bins only ---
n=30000
canon=np.arange(n,dtype=np.uint64)
is_unique=np.zeros(n,dtype=bool)
anchors=np.zeros(n,dtype=bool)
# bin0 (0-10k): 5 anchors -> not starved; has 100 unique
anchors[[100,200,300,400,500]]=True
is_unique[1000:1100]=True
# bin1 (10k-20k): 0 anchors -> starved; 50 unique
is_unique[15000:15050]=True
# bin2 (20k-30k): 1 anchor -> starved (min=2); 0 unique
anchors[25000]=True
take=nc.select_complement(canon,is_unique,anchors,bin_bp=10000,min_anchors=2,max_per_bin=200)
check("complement only in starved bins (bin0 excluded)", int(take[:10000].sum()), 0)
check("complement taken in starved bin1", int(take[10000:20000].sum()), 50)
check("no unique kmers in bin2 -> nothing", int(take[20000:].sum()), 0)

# --- max_per_bin cap + even spacing ---
is_unique2=np.zeros(n,dtype=bool); is_unique2[10000:19000]=True
anchors2=np.zeros(n,dtype=bool)
take2=nc.select_complement(canon,is_unique2,anchors2,10000,2,max_per_bin=100)
check("max_per_bin cap honoured", int(take2[10000:20000].sum()), 100)

# --- trailing partial window (contig end) must be evaluated ---
n=25000
canon3=np.arange(n,dtype=np.uint64)
uniq3=np.zeros(n,dtype=bool); uniq3[20000:21000]=True      # unique kmers in the tail
anc3=np.zeros(n,dtype=bool); anc3[[100,200,300,5000,15000,15001]]=True
t3=nc.select_complement(canon3,uniq3,anc3,10000,2,200)
check("tail window evaluated (was skipped before fix)", int(t3.sum()), 200)
check("tail complement lands in the tail", bool(t3[20000:].all() or t3[20000:21000].sum()==200), True)

# short tail must not be called starved merely for being short
anc4=np.zeros(n,dtype=bool); anc4[[100,200,15000,15001,20500]]=True
t4=nc.select_complement(canon3,uniq3,anc4,10000,2,200)
check("short tail with scaled threshold not starved", int(t4.sum()), 0)

# a full-length contig (no remainder) behaves as before
m=20000
canon5=np.arange(m,dtype=np.uint64)
uniq5=np.zeros(m,dtype=bool); uniq5[12000:12500]=True
anc5=np.zeros(m,dtype=bool); anc5[[10,20,30]]=True          # bin1 (10k-20k) starved
t5=nc.select_complement(canon5,uniq5,anc5,10000,2,200)
check("exact multiple of bin size unchanged", int(t5.sum()), 200)

# --- density sanity on random sequence ---
_rng=random.Random(7)
rs="".join(_rng.choice("ACGT") for _ in range(200000))
c,u,p,val = nc.analyse_contig(rs,32,250)
dens=p.sum()/val.sum(); theory=2/251
ok=abs(dens-theory)/theory<0.10
print(f"[{'PASS' if ok else 'FAIL'}] minimizer density {dens:.5f} vs theory {theory:.5f}")
if not ok: fails.append("density")
ok2 = u.sum()/val.sum() > 0.99
print(f"[{'PASS' if ok2 else 'FAIL'}] random seq nearly all-unique: {u.sum()/val.sum():.4f}")
if not ok2: fails.append("uniqueness")

print()
if fails: print("FAILED:",fails); sys.exit(1)
print("All complement tests passed.")

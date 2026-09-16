#!/usr/bin/env python3
"""k-mer identity filtering: loader formats, and BED-vs-kmer difference."""
import sys, os, tempfile, gzip
from collections import namedtuple
import os
NTLINK_ROOT = os.environ.get("NTLINK_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
NTLINK_BIN = os.path.join(NTLINK_ROOT, "bin")
sys.path.insert(0,NTLINK_BIN)
import ntlink_kmer_mask as km

fails=[]
def check(n,g,e):
    ok=g==e; print(f"[{'PASS' if ok else 'FAIL'}] {n}: got={g} expected={e}")
    if not ok: fails.append(n)

Minimizer=namedtuple("Minimizer",["contig","position","strand"])
K=8
def canon(s):
    c={"A":0,"C":1,"G":2,"T":3}; f=r=0
    for i,b in enumerate(s):
        v=c[b]; f=(f<<2)|v; r|=(3-v)<<(2*i)
    return min(f,r)

def wr(txt,suffix=".txt",gz=False):
    fh=tempfile.NamedTemporaryFile("wb" if gz else "w",suffix=suffix,delete=False)
    if gz: fh.write(gzip.compress(txt.encode()))
    else: fh.write(txt)
    fh.close(); return fh.name

# --- loader: sequences, integers, meryl two-column, gzip, wrong length ---
f=wr("ACGTACGT\nTTTTAAAA\n")
check("loader: sequences", km.load_kmer_set(f,K)=={canon("ACGTACGT"),canon("TTTTAAAA")}, True)
f2=wr("ACGTACGT\t17\nTTTTAAAA\t3\n")
check("loader: meryl 2-column", km.load_kmer_set(f2,K)=={canon("ACGTACGT"),canon("TTTTAAAA")}, True)
f3=wr(f"{canon('ACGTACGT')}\n")
check("loader: packed integers", km.load_kmer_set(f3,K)=={canon("ACGTACGT")}, True)
f4=wr("ACGTACGT\n", gz=True, suffix=".txt.gz")
check("loader: gzipped", km.load_kmer_set(f4,K)=={canon("ACGTACGT")}, True)
f5=wr("ACGT\nACGTACGTAA\n")
check("loader: wrong-length skipped", km.load_kmer_set(f5,K)==set(), True)
check("loader: no path -> None", km.load_kmer_set(None,K), None)
# canonical: revcomp maps to same value
check("canonical revcomp equal", canon("ACGTACGT")==canon("ACGTACGT"[::-1].translate(str.maketrans("ACGT","TGCA"))), True)

# --- build a tiny assembly and minimizer dict ---
seq="AAAACCTCGAGGTATTACGTACCTCTCCGGAG"          # 32 bp, three distinct canonical k-mers
fa=wr(f">ctg1\n{seq}\n",suffix=".fa")
# three "minimizers" at positions 0, 8, 16
mx_info={"h0":Minimizer("ctg1",0,"+"),
         "h8":Minimizer("ctg1",8,"+"),
         "h16":Minimizer("ctg1",16,"+")}
k0,k8,k16 = seq[0:8], seq[8:16], seq[16:24]

# --- exclude by k-mer identity ---
ex=wr(f"{k8}\n")
out=km.filter_mx_info_by_kmers(mx_info,fa,K,exclude_kmers=km.load_kmer_set(ex,K),verbose=False)
check("exclude drops only the named k-mer", sorted(out), ["h0","h16"])

# --- include (whitelist) ---
inc=wr(f"{k0}\n{k16}\n")
out=km.filter_mx_info_by_kmers(mx_info,fa,K,include_kmers=km.load_kmer_set(inc,K),verbose=False)
check("include keeps only whitelisted", sorted(out), ["h0","h16"])

# --- no filter -> identity ---
out=km.filter_mx_info_by_kmers(mx_info,fa,K,verbose=False)
check("no filter returns input object", out is mx_info, True)

# --- KEY: BED removes everything in an interval, k-mer filter does not ---
bed=wr("ctg1\t6\t18\n",suffix=".bed")
out_bed=km.filter_mx_info_by_regions(mx_info,exclude_bed=bed,verbose=False)
out_kmer=km.filter_mx_info_by_kmers(mx_info,fa,K,exclude_kmers=km.load_kmer_set(wr(f"{k8}\n"),K),verbose=False)
check("BED interval drops both pos 8 and 16", sorted(out_bed), ["h0"])
check("k-mer filter drops only pos 8", sorted(out_kmer), ["h0","h16"])
print("   -> 좌표 마스킹은 구간 내 고유 k-mer까지 제거, k-mer 필터는 지정된 것만 제거")

for f in [f,f2,f3,f4,f5,fa,ex,inc,bed]:
    try: os.unlink(f)
    except: pass
print()
if fails: print("FAILED:",fails); sys.exit(1)
print("All k-mer filter tests passed.")

# ntLink-P

**Precision-oriented anchor control for ntLink.**

ntLink-P is a modified build of [ntLink](https://github.com/bcgsc/ntLink) v1.3.11
aimed at haplotype-resolved long-read assemblies, where a false join is more
costly than a missed one. It does not change how ntLink scaffolds. It changes
which minimizers are allowed to act as anchors, and adds anchors back where
minimizer sampling leaves a region unanchored.

> **This code was written with Claude (Anthropic).** Every change in this fork,
> the helper tools, the test suites, and this document were produced in a
> working session with Claude. Treat it accordingly: read the diff, run the
> tests, and validate on your own data before trusting any result. Upstream
> ntLink is unaffiliated with this work.

---

## Why

Three properties of stock ntLink motivated the changes.

**Minimizer selection is blind to uniqueness.** A minimizer wins its window by
hash rank, not by whether its k-mer is unique in the assembly. In regions where
unique k-mers are scarce, such as repeat interiors made unique only by local
variants, the sampling rate of `2/(w+1)` frequently picks none of them, and the
region is left unanchored even though usable k-mers exist. Measured on CHM13
chr21 at k=32, w=250: the fraction of unique k-mers captured is a flat 0.70 to
0.82 % regardless of local repeat content, and 7.9 % of 10 kb windows hold
fewer than two anchors.

**Anchor uniqueness is judged from the assembly alone.** A repeat that is
collapsed in the assembly looks unique, so reads from every genomic copy anchor
to the single assembled copy. Because the resulting false joins are supported
consistently by many reads, the multi-read agreement that normally protects
ntLink reinforces them instead.

**Some anchors should never be used.** A contig terminus carrying a telomere
array is a chromosome end; anything scaffolded past it is a false positive.

---

## What is added

| Area | Option | Effect |
|---|---|---|
| Anchor quality | `--hc` | Minimum minimizer hits for a read to count as well-anchored. Parameterises a value hard-coded as `> 1` upstream; the default of 2 reproduces stock behaviour. |
| Anchor quality | `--mx-ratio` | Fraction of a contig's unique minimizers inside the anchored interval that the read must hit. Round 1 only. Off by default; see caveats. |
| Masking (interval) | `--mx-exclude-bed`, `--mx-include-bed` | Drop or keep minimizers by assembly coordinate. Intended for things that genuinely are intervals, such as telomeric termini. |
| Masking (k-mer) | `--mx-exclude-kmers`, `--mx-include-kmers` | Drop or keep minimizers by k-mer identity. Accepts `meryl print` output or packed canonical integers, plain or gzipped. Unlike an interval mask this removes exactly the k-mers named, so individually unique k-mers inside a flagged region survive. |
| Complement seeding | `ntlink_complement.py`, `ntlink_inject_complement.py` | Find windows that minimizer sampling left starved, emit the unique k-mers inside them as supplementary anchors, and inject them into both the target and the read sketches. |
| Helper | `ntlink_telomere_bed.py` | Detect telomeric termini and emit a draft BED mask for curation. |
| Runner | `run_ntlink_precision.sh` | End-to-end pipeline with a proper `--outdir`. |

All defaults reproduce stock ntLink. The regression suite asserts this.

### How complement seeding works

Minimizer sampling and direct lookup are different operations, and the fork uses
both. Assembly-side, windows with fewer than `--comp-min-anchors` anchors are
flagged and the unique k-mers inside them are emitted (capped and evenly spaced
by `--comp-max-per-bin`). Read-side, reads are scanned for those exact k-mers
with **no window competition**, so a read contributes every complement k-mer it
carries error-free rather than only those that happen to win a window.

ntLink matches a read minimizer to the assembly by string equality of the `mx`
token, and the token is opaque to ntLink, so complement k-mers carry their own
identifier `C<canonical k-mer integer>`, identical by construction wherever that
k-mer occurs.

On a synthetic contig whose terminal flank is anchor-starved, under strict
settings (`-n 2 -a 3 --hc 3 -x 1.1`) the baseline produced no edge at all while
complement seeding produced a join with a gap estimate of 332 bp against a true
gap of 300 bp.

---

## Install

Requires a working ntLink installation (`indexlr` from btllib, `abyss-scaffold`,
`MergeContigs`) plus Python 3 with `numpy` and `python-igraph`.

```bash
git clone https://github.com/bcgsc/ntLink.git
cd ntLink
git checkout v1.3.11

# apply the modifications to the three upstream files
git apply /path/to/changes.patch

# add the new tools
cp /path/to/bin/ntlink_kmer_mask.py \
   /path/to/bin/ntlink_complement.py \
   /path/to/bin/ntlink_inject_complement.py \
   /path/to/bin/ntlink_telomere_bed.py bin/
cp /path/to/run_ntlink_precision.sh .
chmod +x bin/*.py run_ntlink_precision.sh
```

`bin/ntlink_kmer_mask.py` is imported by `ntlink_pair.py`; without it the
patched `ntlink_pair.py` will not start.

---

## Usage

### 1. Build the k-mer lists with meryl

Count with **the same k** that ntLink will use.

```bash
meryl count k=32 memory=64 threads=16 reads.fq.gz output reads.meryl
meryl count k=32 memory=64 threads=16 hap1.fa   output hap1.meryl
meryl count k=32 memory=64 threads=16 hap2.fa   output hap2.meryl

meryl histogram reads.meryl > reads.hist    # read the peaks before going on
```

With haploid coverage `c`, heterozygous k-mers peak at `c` and homozygous ones at
`2c`. The bands below assume `c = 15`; use your own histogram.

```bash
# read-level copy-number bands
meryl greater-than  7 less-than 23 reads.meryl output reads.1copy.meryl
meryl greater-than 22 less-than 45 reads.meryl output reads.2copy.meryl

# exactly one copy in each haplotype assembly
meryl equal-to 1 hap1.meryl output hap1.uniq.meryl
meryl equal-to 1 hap2.meryl output hap2.uniq.meryl

# haplotype-specific: read 1-copy and unique in hap1
meryl intersect reads.1copy.meryl hap1.uniq.meryl output allowed.1copy.meryl

# homozygous, confirmed: read 2-copy and unique in BOTH haplotypes
meryl intersect hap1.uniq.meryl hap2.uniq.meryl  output hom.uniq.meryl
meryl intersect reads.2copy.meryl hom.uniq.meryl output allowed.2copy.meryl

meryl union allowed.1copy.meryl allowed.2copy.meryl output allowed.meryl
meryl print allowed.meryl | gzip > allowed_kmers.txt.gz

# collapsed repeats: far above the homozygous peak
meryl greater-than 45 reads.meryl output collapsed.meryl
meryl print collapsed.meryl | gzip > blacklist_kmers.txt.gz
```

Requiring one copy in *each* haplotype is what separates a true homozygous
k-mer from a segmental duplication carried by a single haplotype.

### 2. Curate the interval mask

Telomeric termini are an interval, not a k-mer set. Generate a draft and check
it by eye before use.

```bash
ntlink_telomere_bed.py hap1.fa --mask-width 30000 \
  --report hap1_telo.tsv -o hap1_telo.bed
```

Set `--mask-width` at or above the long-read N50. If every surviving anchor sits
at least one read length inside the terminus, the gap arithmetic forces any join
past the telomere to a negative gap that ntLink rejects.

### 3. Run

```bash
run_ntlink_precision.sh \
  --target hap1.fa --reads reads.fq.gz \
  --outdir hap1_scaffold \
  -k 32 -w 250 -z 2000 -n 2 -a 2 -x 1.1 --hc 2 \
  --exclude-bed hap1_telo.bed \
  --exclude-kmers blacklist_kmers.txt.gz \
  --allowed-kmers allowed_kmers.txt.gz \
  --comp-bin 2000 --comp-min-anchors 3 --comp-max-per-bin 100 \
  --repeat-filter
```

Run each haplotype separately. Merging `hap1.fa` and `hap2.fa` into one FASTA
makes every homozygous k-mer appear twice, so ntLink's duplicate filter removes
it and only heterozygous k-mers survive: roughly 3 % of anchors at k=32.

---

## Parameters

Values that differ from upstream defaults are marked.

| Option | Default | Notes |
|---|---|---|
| `-k` | 32 | Must match the meryl k. HiFi tolerates 32 to 40; noisy ONT needs less. |
| `-w` | 250 | Minimizer window. Density is about `2/(w+1)`. Independent of `k`. |
| `-z` | 2000 *(was 1000)* | Minimum contig length considered. |
| `-n` | 2 *(was 1)* | Minimum edge weight, counting every supporting read. Also the start of the abyss-scaffold sweep, so it applies twice. |
| `-a` | 2 *(was 1)* | Minimum well-anchored reads per edge, counting only reads with at least `--hc` hits on **both** contigs. Always `a <= n`, so at equal values `-a` is the stricter filter. |
| `-x` | 1.1 *(was 0)* | Colinearity: drop a contig whose span on the contig exceeds `x * read_span + k`. **Never use a value below 1** — it imposes a hard cap of `k/(1-x)` on read span and discards most legitimate anchors. |
| `--hc` | 2 | 2 equals stock behaviour. |
| `--mx-ratio` | 0 (off) | See caveats. |
| `--comp-bin` | 10000 | Window for the starvation test. 2000 localises it to flank scale. |
| `--comp-min-anchors` | 2 | At `w=250` a 2 kb window holds about 16 anchors, so 3 targets genuinely sparse windows. |
| `--comp-max-per-bin` | 200 | With `--comp-bin 2000`, a cap of 100 gives roughly one complement anchor per 20 bp. |
| `--repeat-filter` | off | Stock ntLink option. Removes minimizers repeated **within a single read**, per read, not globally. Useful against sequencing errors that duplicate a k-mer inside one read; complement anchors are subject to it as well. |
| `-g` | 20 | abyss-scaffold `--min-gap`: the floor on the N run written between contigs. Affects output length only, never which contigs join. |

A suggested starting point for HiFi: `-k 32 -w 250 -z 2000 -n 2 -a 2 -x 1.1
--hc 2 --repeat-filter`, with `--mx-ratio 0`.

---

## Caveats

**Overlap trimming is not run.** The runner drives pairing directly and hands the
graph to abyss-scaffold, so `ntlink_overlap_sequences.py` never executes. That
step is skipped deliberately: when the trims it computes for the two ends of a
contig are incompatible it drops the contig from the path entirely, which is
sequence loss. The cost is that a join with a negative gap estimate keeps the
overlapping sequence on both contigs, separated by a short N run. Curate those
junctions afterwards, for example by aligning adjacent contig ends with
`minimap2 -c`, and do not report duplication metrics before doing so.

**`--mx-ratio` is a diagnostic, not a filter.** Its value tracks roughly
`identity^k`, so it reports read accuracy more than mapping quality: on HiFi all
genuine anchors cluster near 0.9, and on ONT they are low and variable. It also
inherits an interval-inflation bias, since a single outlier hit widens the
denominator. Log it if useful; leave the threshold at 0.

**Complement candidates are chosen before masking.** `ntlink_complement.py` does
not see the BED or k-mer masks, so masked anchors are selected, scanned for, and
only then discarded in `read_minimizers`. The result is correct but wasteful.
Excluding masked k-mers from `--allowed-kmers` avoids the wasted work.

**Verification status.** Sketching, complement selection, injection, pairing, and
the graph-level differences were all exercised end to end. The stitch and
`MergeContigs` stages were not: the development container had no
`abyss-scaffold` or `MergeContigs`, and the stand-ins used do not emit the
gap-token path format that `ntlink_stitch_paths.py` parses. Check those stages on
your first real run.

---

## Tests

```bash
python3 tests/test_regression.py    # defaults reproduce stock ntLink
python3 tests/test_kmer_mask.py     # BED parsing, interval membership
python3 tests/test_kmer_filter.py   # k-mer identity masking, loader formats
python3 tests/test_complement.py    # starved-window selection, tail windows
python3 tests/test_telomere.py      # telomere detection
python3 tests/test_filters.py       # mx-ratio and hc in the real code path
python3 tests/test_impl.py          # canonical k-mers, minimizer density
```

`testdata/` carries a synthetic assembly with five scenarios that separate the
features: a clean join, an anchor-starved flank, a telomeric false join, a
collapsed-repeat false join, and a whitelist case. `testdata/run_scenarios.sh`
runs seven configurations and prints which edges survive each one;
`testdata/truth.tsv` is the expected outcome.

---

## Licence and attribution

ntLink is distributed under the GNU General Public License v3, and this fork
inherits it. Upstream ntLink is developed at the BC Cancer Genome Sciences
Centre; please cite the original work when using it
([bcgsc.ca/resources/software/ntlink](https://www.bcgsc.ca/resources/software/ntlink)):

> Coombe L, Li JX, Lo T, Wong J, Nikolic V, Warren RL and Birol I. LongStitch:
> high-quality genome assembly correction and scaffolding using long reads.
> *BMC Bioinformatics* 22, 534 (2021). doi:10.1186/s12859-021-04451-7
>
> Coombe L, Warren RL, Wong J, Nikolic V and Birol I. ntLink: A toolkit for de
> novo genome assembly scaffolding and mapping using long reads. *Current
> Protocols* 3, e733 (2023). doi:10.1002/cpz1.733

The modifications in this fork carry no such claim. They were written with
Claude (Anthropic), validated against synthetic data and CHM13 chr21, and have
not been benchmarked on a production assembly.
#!/usr/bin/env bash
#
# ntLink precision scaffolding: full pipeline
#
#   sketch -> [telomere mask] -> [complement selection] -> [injection]
#          -> pairing -> abyss-scaffold sweep -> stitch -> gap fill
#
# Why pairing is run directly instead of purely through the ntLink Makefile:
# the Makefile pipes the read sketch straight into ntlink_pair.py (ntLink:246-248)
# and never writes it to disk, so there is no file to inject complement anchors
# into. We therefore materialise both sketches, inject, run pairing ourselves,
# and hand the resulting .dot graph back to the Makefile for the remaining steps.
#
# Requirements on PATH: indexlr (btllib), abyss-scaffold, ntLink's bin/
set -euo pipefail

# ------------------------------------------------------------------ defaults
TARGET=""              # assembly to scaffold (one haplotype at a time)
READS=""               # long reads, gzipped fastq/fasta
PREFIX=""
OUTDIR="."             # every generated file lands here
K=32
W=250
T=8
Z=2000                 # min contig size considered for scaffolding
N=2                    # min edge weight / start of abyss-scaffold sweep
MAX_N=10
A=2                    # min well-anchored reads per edge
X=1.1                  # colinearity tolerance
HC=2                   # min minimizer hits for a read to be well-anchored
MX_RATIO=0             # interval coverage filter (0 = off)
G=20                   # min gap size for abyss-scaffold
EXCLUDE_BED=""         # telomere / collapsed-repeat mask (curated by the user)
INCLUDE_BED=""         # optional whitelist: keep anchors only inside these regions
EXCLUDE_KMERS=""       # k-mer blacklist for base seeding (meryl copy number)
INCLUDE_KMERS=""       # k-mer whitelist for base seeding
ALLOWED_KMERS=""       # meryl-curated k-mer whitelist for complement
REPEAT_FILTER=0        # 1 = drop minimizers repeated within a single read
COMPLEMENT=1           # 1 = run complement seeding, 0 = skip
COMP_BIN=10000
COMP_MIN_ANCHORS=2
COMP_MAX_PER_BIN=200
NTLINK_BIN=""          # defaults to the directory holding ntlink_pair.py

usage() {
  cat <<EOF
Usage: $0 --target hap1.fa --reads reads.fq.gz [options]

Required:
  --target FILE          assembly FASTA (run each haplotype separately)
  --reads  FILE          long reads, gzipped

Common options:
  --outdir DIR           write every output here  [$OUTDIR]
  --prefix STR           output basename         [<target basename>.k<K>.w<W>.z<Z>]
  -k INT                 k-mer size              [$K]
  -w INT                 minimizer window        [$W]
  -t INT                 threads                 [$T]
  -z INT                 min contig size         [$Z]
  -n INT                 min edge weight         [$N]
  --max-n INT            end of abyss sweep      [$MAX_N]
  -a INT                 min well-anchored reads [$A]
  -x FLOAT               colinearity tolerance   [$X]
  --hc INT               min hits per anchor     [$HC]
  --mx-ratio FLOAT       interval coverage filter[$MX_RATIO]
  --exclude-bed FILE     curated mask (telomeres, collapsed repeats)
  --include-bed FILE     whitelist; anchors kept only inside these regions
  --exclude-kmers FILE   k-mer blacklist for base minimizer seeding
  --include-kmers FILE   k-mer whitelist for base minimizer seeding
  --allowed-kmers FILE   meryl k-mer whitelist for complement
  --repeat-filter        drop minimizers that repeat within one read (stock ntLink
                         option; off by default)
  --no-complement        skip complement seeding
  --comp-bin INT         starvation window       [$COMP_BIN]
  --comp-min-anchors INT starved if fewer than   [$COMP_MIN_ANCHORS]
  --comp-max-per-bin INT cap per starved window  [$COMP_MAX_PER_BIN]
  --ntlink-bin DIR       ntLink bin directory
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target) TARGET="$2"; shift 2;;
    --reads) READS="$2"; shift 2;;
    --prefix) PREFIX="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -k) K="$2"; shift 2;;
    -w) W="$2"; shift 2;;
    -t) T="$2"; shift 2;;
    -z) Z="$2"; shift 2;;
    -n) N="$2"; shift 2;;
    --max-n) MAX_N="$2"; shift 2;;
    -a) A="$2"; shift 2;;
    -x) X="$2"; shift 2;;
    --hc) HC="$2"; shift 2;;
    --mx-ratio) MX_RATIO="$2"; shift 2;;
    -g) G="$2"; shift 2;;
    --exclude-bed) EXCLUDE_BED="$2"; shift 2;;
    --include-bed) INCLUDE_BED="$2"; shift 2;;
    --exclude-kmers) EXCLUDE_KMERS="$2"; shift 2;;
    --include-kmers) INCLUDE_KMERS="$2"; shift 2;;
    --allowed-kmers) ALLOWED_KMERS="$2"; shift 2;;
    --repeat-filter) REPEAT_FILTER=1; shift;;
    --no-complement) COMPLEMENT=0; shift;;
    --comp-bin) COMP_BIN="$2"; shift 2;;
    --comp-min-anchors) COMP_MIN_ANCHORS="$2"; shift 2;;
    --comp-max-per-bin) COMP_MAX_PER_BIN="$2"; shift 2;;
    --ntlink-bin) NTLINK_BIN="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown option: $1" >&2; usage;;
  esac
done

[[ -n "$TARGET" && -n "$READS" ]] || usage
[[ -f "$TARGET" ]] || { echo "Target not found: $TARGET" >&2; exit 1; }
[[ -f "$READS"  ]] || { echo "Reads not found: $READS" >&2; exit 1; }

if [[ -z "$NTLINK_BIN" ]]; then
  NTLINK_BIN="$(dirname "$(command -v ntlink_pair.py)")" ||
    { echo "ntlink_pair.py not on PATH; use --ntlink-bin" >&2; exit 1; }
fi
[[ -z "$PREFIX" ]] && PREFIX="$(basename "$TARGET").k${K}.w${W}.z${Z}"

# Inputs may live anywhere; only OUTDIR is written to. A prefix containing a
# slash is treated as a path the caller chose and is left alone.
mkdir -p "$OUTDIR"
[[ -w "$OUTDIR" ]] || { echo "Output directory not writable: $OUTDIR" >&2; exit 1; }
case "$PREFIX" in
  */*) ;;                       # caller gave an explicit path
  *) PREFIX="${OUTDIR%/}/${PREFIX}";;
esac
mkdir -p "$(dirname "$PREFIX")"

for tool in indexlr abyss-scaffold; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 1; }
done

# The target sketch is written next to the other outputs, never beside the
# input assembly, so a shared or read-only assembly directory still works.
TARGET_TSV="${OUTDIR%/}/$(basename "$TARGET").k${K}.w${W}.tsv"
READ_TSV="${PREFIX}.reads.k${K}.w${W}.tsv"

echo "### ntLink precision scaffolding"
echo "  target=$TARGET"
echo "  reads =$READS"
echo "  outdir=$OUTDIR  prefix=$PREFIX"
echo "  k=$K w=$W z=$Z n=$N a=$A x=$X hc=$HC mx_ratio=$MX_RATIO"
echo "  complement=$COMPLEMENT repeat_filter=$REPEAT_FILTER"

# ------------------------------------------------------- 1. sketch both sides
echo
echo "### [1/6] Sketching minimizers"
if [[ ! -s "$TARGET_TSV" ]]; then
  indexlr --long --pos --strand -k "$K" -w "$W" -t "$T" "$TARGET" > "$TARGET_TSV"
fi
if [[ ! -s "$READ_TSV" ]]; then
  gzip -cd "$READS" | indexlr --long --pos --strand --len -k "$K" -w "$W" -t "$T" - > "$READ_TSV"
fi
echo "  target sketch: $TARGET_TSV"
echo "  read sketch  : $READ_TSV"

# --------------------------------------------- 2. complement selection/inject
PAIR_TARGET_TSV="$TARGET_TSV"
PAIR_READ_TSV="$READ_TSV"

if [[ "$COMPLEMENT" -eq 1 ]]; then
  echo
  echo "### [2/6] Complement seeding"
  COMP_TSV="${PREFIX}.complement.tsv"
  ALLOWED_ARG=()
  [[ -n "$ALLOWED_KMERS" ]] && ALLOWED_ARG=(--allowed-kmers "$ALLOWED_KMERS")

  "$NTLINK_BIN/ntlink_complement.py" "$TARGET" \
    -k "$K" -w "$W" --bin "$COMP_BIN" \
    --min-anchors "$COMP_MIN_ANCHORS" --max-per-bin "$COMP_MAX_PER_BIN" \
    "${ALLOWED_ARG[@]}" -o "$COMP_TSV" --stats "${PREFIX}.complement_stats.npz"

  if [[ -s "$COMP_TSV" ]]; then
    "$NTLINK_BIN/ntlink_inject_complement.py" \
      --complement "$COMP_TSV" \
      --target-tsv "$TARGET_TSV" --read-tsv "$READ_TSV" \
      --reads "$READS" -k "$K" \
      --out-target "${PREFIX}.target.aug.tsv" --out-reads "${PREFIX}.reads.aug.tsv"
    PAIR_TARGET_TSV="${PREFIX}.target.aug.tsv"
    PAIR_READ_TSV="${PREFIX}.reads.aug.tsv"
    echo "  injected: $PAIR_TARGET_TSV / $PAIR_READ_TSV"
  else
    echo "  WARNING: no complement anchors selected; continuing with plain sketches."
    echo "  If this is unexpected, --comp-bin ($COMP_BIN) may exceed your contig sizes,"
    echo "  or no window fell below --comp-min-anchors ($COMP_MIN_ANCHORS)."
  fi
else
  echo
  echo "### [2/6] Complement seeding skipped"
fi

# ----------------------------------------------------------- 3. build pairing
echo
echo "### [3/6] Pairing and scaffold graph"
PAIR_OPTS=(--verbose)
[[ -n "$EXCLUDE_BED" ]] && PAIR_OPTS+=(--mx-exclude-bed "$EXCLUDE_BED")
[[ -n "$INCLUDE_BED" ]] && PAIR_OPTS+=(--mx-include-bed "$INCLUDE_BED")
[[ -n "$EXCLUDE_KMERS" ]] && PAIR_OPTS+=(--mx-exclude-kmers "$EXCLUDE_KMERS")
[[ -n "$INCLUDE_KMERS" ]] && PAIR_OPTS+=(--mx-include-kmers "$INCLUDE_KMERS")
[[ "$REPEAT_FILTER" -eq 1 ]] && PAIR_OPTS+=(--repeat-filter)

"$NTLINK_BIN/ntlink_pair.py" \
  -p "$PREFIX" -n "$N" -m "$PAIR_TARGET_TSV" -s "$TARGET" \
  -k "$K" -a "$A" -z "$Z" -f 10 -x "$X" \
  --hc "$HC" --mx-ratio "$MX_RATIO" \
  "${PAIR_OPTS[@]}" "$PAIR_READ_TSV"

DOT="${PREFIX}.n${N}.scaffold.dot"
[[ -s "$DOT" ]] || { echo "Pairing produced no graph: $DOT" >&2; exit 1; }
EDGES=$(grep -c ' -> ' "$DOT" || true)
echo "  graph: $DOT  edges=$EDGES"
if [[ "$EDGES" -eq 0 ]]; then
  echo
  echo "No edges passed the filters, so there is nothing to scaffold."
  echo "Loosen -n ($N), -a ($A), --hc ($HC) or -x ($X), or enable complement seeding."
  exit 0
fi

# ------------------------------------------------ 4. abyss-scaffold weight sweep
echo
echo "### [4/6] abyss-scaffold sweep (n=$N..$MAX_N)"
PATHS=()
for nn in $(seq "$N" "$MAX_N"); do
  out="${PREFIX}.n${nn}.abyss-scaffold.path"
  # Mirrors ntLink:255 — note -k2 (not the minimizer k) and the .sterr log,
  # which ntlink_stitch_paths.py parses to pick the best n (ntlink_stitch_paths.py:376).
  abyss-scaffold -k2 -n "$nn" -s"$Z" --min-gap "$G" - "$DOT" 1> "$out" 2> "${out}.sterr"
  PATHS+=("$out")
done

# ------------------------------------------------------------- 5. stitch paths
echo
echo "### [5/6] Stitching paths"
STITCH="${PREFIX}.stitch.path"
"$NTLINK_BIN/ntlink_stitch_paths.py" \
  --min_n "$N" --max_n "$MAX_N" -p "$PREFIX" -g "$DOT" \
  --conservative -o "$STITCH" "${PATHS[@]}"

# ------------------------------------------------------------ 6. build scaffolds
echo
echo "### [6/6] Writing scaffolds"
if command -v MergeContigs >/dev/null; then
  abyss-scaffold --version >/dev/null 2>&1 || true
  MergeContigs -k"$K" -o "${PREFIX}.scaffolds.fa" "$TARGET" "$DOT" "$STITCH"
  echo "  scaffolds: ${PREFIX}.scaffolds.fa"
else
  echo "  MergeContigs not found; final assembly step skipped."
  echo "  Path file ready for merging: $STITCH"
fi

echo
echo "### Done"
echo "  graph     : $DOT"
echo "  paths     : $STITCH"
[[ "$COMPLEMENT" -eq 1 && -s "${PREFIX}.complement.tsv" ]] && \
  echo "  complement: ${PREFIX}.complement.tsv ($(wc -l < "${PREFIX}.complement.tsv") anchors)"
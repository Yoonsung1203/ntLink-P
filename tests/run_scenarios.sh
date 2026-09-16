#!/usr/bin/env bash
# Exercise every masking mode and tabulate which edges survive.
#
#   A  no masking, no complement          -> false joins present
#   B  BED masks (telomere + collapsed)   -> false joins removed
#   C  k-mer masks instead of BED         -> same removal, finer granularity
#   D  BED masks, strict (a3 hc3)         -> S2 lost
#   E  BED masks, strict + complement     -> S2 recovered
#   F  BED masks + BED whitelist          -> S5 dropped
#   G  BED masks + k-mer whitelist        -> S5 dropped
set -uo pipefail
ROOT="${NTLINK_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RUN="$ROOT/run_ntlink_precision.sh"
BIN="$ROOT/bin"
BASE="-k 32 -w 250 -z 2000 --max-n 3 -x 1.1 --ntlink-bin $BIN"

run() {
  local name="$1"; shift
  rm -rf "out_$name"; mkdir -p "out_$name"
  $RUN --target target.fa --reads reads.fa.gz --outdir "out_$name" \
       --prefix "$name" $BASE "$@" > "out_$name/log.txt" 2>&1
  grep -oE '"ctg[0-9a-b]+\+" -> "ctg[0-9a-b]+\+"' "out_$name/${name}.n2.scaffold.dot" 2>/dev/null \
    | sed 's/"//g; s/+//g; s/ -> /-/' | sort -u | tr '\n' ' '
  echo
}

printf "%-42s %s\n" "A no masking" "$(run A -n 2 -a 2 --hc 2 --no-complement)"
printf "%-42s %s\n" "B BED masks" "$(run B -n 2 -a 2 --hc 2 --no-complement --exclude-bed exclude.bed)"
printf "%-42s %s\n" "C k-mer masks" "$(run C -n 2 -a 2 --hc 2 --no-complement --exclude-kmers blacklist_kmers.txt)"
printf "%-42s %s\n" "D BED masks, strict" "$(run D -n 2 -a 3 --hc 3 --no-complement --exclude-bed exclude.bed)"
printf "%-42s %s\n" "E BED masks, strict + complement" "$(run E -n 2 -a 3 --hc 3 --exclude-bed exclude.bed --allowed-kmers allowed_kmers.txt --comp-bin 2000 --comp-min-anchors 3 --comp-max-per-bin 100)"
printf "%-42s %s\n" "F BED masks + BED whitelist" "$(run F -n 2 -a 2 --hc 2 --no-complement --exclude-bed exclude.bed --include-bed whitelist.bed)"
printf "%-42s %s\n" "G BED masks + k-mer whitelist" "$(run G -n 2 -a 2 --hc 2 --no-complement --include-kmers whitelist_kmers.txt)"

#!/usr/bin/env bash
#
# wgs_target_prefilter.sh: reduce an aligned WGS BAM/CRAM to the reads an
# exome-mode `immuno` run uses, emitted as complete FASTQ read pairs.
#
# Region filtering an aligned BAM with `samtools view -L` keeps a record only if
# it overlaps a region, so any pair with one mate on target and one mate off
# target loses a mate. Picard SamToFastq (sequenceAlignAndTag / sequenceToFastq
# BAM mode) then fails, and bwa-mem2 -p mis-pairs the interleaved stream. This
# script avoids that by:
#
#   1. fetching every primary record that overlaps the regions, plus its mate
#      wherever the mate maps (`samtools view -P`, absent from samtools 1.16.1, present in 1.21);
#   2. grouping records by read name (`samtools collate`);
#   3. writing R1/R2 FASTQs per read group (`samtools fastq`), which routes any
#      record without its mate to a separate singleton file. Singletons are
#      counted, reported, and not passed on.
#
# Regions kept (union, then merged):
#   - every -L interval file (Picard interval_list or BED), padded by -p bases.
#     Pass the same bait_intervals and target_intervals the immuno run uses.
#   - the MHC region and every HLA / chr6 alt contig in the header. The immuno
#     HLA typers (OptiType, HLA-HD, PHLAT) extract chr6:29836259-33148325, the
#     HLA alt contigs and all unmapped reads from the CRAM, so these must
#     survive the filter.
#   - unplaced unmapped pairs (the BAM's `*` bin), for the same reason.
#     Disable with --no-unplaced.
#
# The immuno WDL realigns the FASTQs to its own reference, so the input BAM may
# come from any GRCh38 alignment whose contig names match the interval files.
#
# Requires: samtools with `view --fetch-pairs` (tested on 1.21) on PATH (or set SAMTOOLS), awk, sort.
# gs:// inputs work if samtools is built with GCS support and GCS_OAUTH_TOKEN
# is set, e.g. `export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)`.

set -euo pipefail

SAMTOOLS=${SAMTOOLS:-samtools}
PADDING=250
THREADS=4
HLA_REGION="chr6:29836259-33148325"
HLA_CONTIG_REGEX='HLA|HSCHR6_MHC|^chr6_.*_alt$'
INCLUDE_UNPLACED=1
MAX_SINGLETON_FRACTION=0.001
INPUT=""
INDEX=""
REFERENCE=""
OUT_PREFIX=""
TMPDIR_ROOT="${TMPDIR:-/tmp}"
INTERVAL_FILES=()

usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT.bam|cram -L INTERVALS [-L INTERVALS ...] -o OUT_PREFIX [options]

Required:
  -i FILE      Coordinate-sorted, indexed WGS BAM or CRAM (local path or gs://).
  -L FILE      Interval file: Picard interval_list (1-based) or BED (0-based).
               Repeat for bait_intervals and target_intervals.
  -o PREFIX    Output prefix; may include a directory, which must exist.

Options:
  -x FILE      Index file, if not at INPUT.bai / INPUT.crai (needed for some gs:// inputs).
  -T FILE      Reference FASTA (required for CRAM).
  -p INT       Padding added to each interval, in bases (default ${PADDING}).
  -@ INT       Threads (default ${THREADS}).
  --hla-region REGION       MHC region to keep (default ${HLA_REGION}).
  --hla-contig-regex REGEX  Contigs kept whole (default '${HLA_CONTIG_REGEX}').
  --no-unplaced             Do not keep unplaced unmapped pairs.
  --max-singleton-fraction F
               Exit non-zero if singleton reads exceed this fraction of paired
               reads written (default ${MAX_SINGLETON_FRACTION}).
  --tmpdir DIR Scratch directory (default \$TMPDIR or /tmp).
  -h           Show this help.

Outputs:
  PREFIX.<RG>.R1.fastq.gz, PREFIX.<RG>.R2.fastq.gz   one pair per read group
  PREFIX.regions.bed                                  merged regions used
  PREFIX.sequences.yaml                               immuno SequenceData entries
  PREFIX.prefilter_qc.tsv                             per-read-group counts
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT=$2; shift 2 ;;
    -x) INDEX=$2; shift 2 ;;
    -T) REFERENCE=$2; shift 2 ;;
    -L) INTERVAL_FILES+=("$2"); shift 2 ;;
    -o) OUT_PREFIX=$2; shift 2 ;;
    -p) PADDING=$2; shift 2 ;;
    -@) THREADS=$2; shift 2 ;;
    --hla-region) HLA_REGION=$2; shift 2 ;;
    --hla-contig-regex) HLA_CONTIG_REGEX=$2; shift 2 ;;
    --no-unplaced) INCLUDE_UNPLACED=0; shift ;;
    --max-singleton-fraction) MAX_SINGLETON_FRACTION=$2; shift 2 ;;
    --tmpdir) TMPDIR_ROOT=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%H:%M:%S')] $*" >&2; }

[[ -n "$INPUT" ]] || die "-i is required"
[[ -n "$OUT_PREFIX" ]] || die "-o is required"
[[ ${#INTERVAL_FILES[@]} -gt 0 ]] || die "at least one -L interval file is required"
[[ "$PADDING" =~ ^[0-9]+$ ]] || die "-p must be a non-negative integer"
for f in "${INTERVAL_FILES[@]}"; do [[ -r "$f" ]] || die "cannot read interval file $f"; done

# -P (--fetch-pairs) is absent from samtools 1.16.1; tested on 1.21.
"$SAMTOOLS" view --help 2>&1 | grep -q -- '--fetch-pairs' \
  || die "$SAMTOOLS lacks 'view --fetch-pairs' (tested on samtools 1.21)"

if [[ "$INPUT" == *.cram ]]; then
  [[ -n "$REFERENCE" ]] || die "CRAM input needs -T reference.fa"
fi

READ_ARGS=()
[[ -n "$REFERENCE" ]] && READ_ARGS+=(-T "$REFERENCE")
# samtools view takes an explicit index with -X, after which the input and
# index are the first two positional arguments.
INPUT_ARGS=("$INPUT")
if [[ -n "$INDEX" ]]; then
  READ_ARGS+=(-X)
  INPUT_ARGS=("$INPUT" "$INDEX")
fi

WORK=$(mktemp -d "${TMPDIR_ROOT%/}/wgs_prefilter.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

log "Reading header"
"$SAMTOOLS" view -H ${READ_ARGS[@]+"${READ_ARGS[@]}"} "${INPUT_ARGS[@]}" > "$WORK/header.sam"
grep -q '^@HD.*SO:coordinate' "$WORK/header.sam" \
  || die "input header does not declare SO:coordinate; sort and index it first"
awk -F'\t' '$1=="@SQ"{n="";l="";for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4);if($i~/^LN:/)l=substr($i,4)}print n"\t"l}' \
  "$WORK/header.sam" > "$WORK/contigs.tsv"
[[ -s "$WORK/contigs.tsv" ]] || die "no @SQ lines in input header"

# ---- Build the region BED --------------------------------------------------
log "Building regions (padding ${PADDING} bp)"
: > "$WORK/raw.bed"
for f in "${INTERVAL_FILES[@]}"; do
  if grep -q -m1 '^@SQ' "$f"; then
    # Picard interval_list: 1-based closed -> BED 0-based half-open.
    awk -F'\t' '!/^@/ && NF>=3 {print $1"\t"($2-1)"\t"$3}' "$f" >> "$WORK/raw.bed"
  else
    awk -F'\t' '!/^(#|track|browser)/ && NF>=3 {print $1"\t"$2"\t"$3}' "$f" >> "$WORK/raw.bed"
  fi
done

hla_chrom=${HLA_REGION%%:*}
hla_range=${HLA_REGION#*:}
hla_start=${hla_range%-*}
hla_end=${hla_range#*-}
printf '%s\t%s\t%s\n' "$hla_chrom" "$((hla_start - 1))" "$hla_end" >> "$WORK/raw.bed"

# Pad, clamp to contig length, drop contigs absent from the header, add whole
# HLA/alt contigs, then sort in header order and merge.
awk -F'\t' -v OFS='\t' -v pad="$PADDING" -v re="$HLA_CONTIG_REGEX" -v dropped="$WORK/dropped.txt" '
  FNR==NR { len[$1]=$2; ord[$1]=FNR; if ($1 ~ re) print FNR, $1, 0, $2; next }
  {
    if (!($1 in len)) { missing[$1]++; next }
    s = $2 - pad; if (s < 0) s = 0
    e = $3 + pad; if (e > len[$1]) e = len[$1]
    if (e > s) print ord[$1], $1, s, e
  }
  END { for (c in missing) print c"\t"missing[c] > dropped }
' "$WORK/contigs.tsv" "$WORK/raw.bed" \
  | sort -k1,1n -k3,3n \
  | awk -F'\t' -v OFS='\t' '
      $2!=c || $3>e { if (c!="") print c, s, e; c=$2; s=$3; e=$4; next }
      $4>e { e=$4 }
      END { if (c!="") print c, s, e }
    ' > "$OUT_PREFIX.regions.bed"

if [[ -s "$WORK/dropped.txt" ]]; then
  log "WARNING: intervals on contigs absent from the input header were skipped:"
  sort "$WORK/dropped.txt" | head -20 >&2
fi
region_bp=$(awk '{s+=$3-$2} END{print s+0}' "$OUT_PREFIX.regions.bed")
genome_bp=$(awk '{s+=$2} END{print s+0}' "$WORK/contigs.tsv")
log "Regions: $(wc -l < "$OUT_PREFIX.regions.bed" | tr -d ' ') intervals, ${region_bp} bp ($(awk -v a="$region_bp" -v b="$genome_bp" 'BEGIN{printf "%.2f", 100*a/b}')% of the reference)"

# ---- Extract complete pairs ------------------------------------------------
# -F 0x900 drops secondary and supplementary records; samtools fastq would drop
# them anyway, and their mates are fetched through the primary records.
log "Fetching region reads and their mates"
"$SAMTOOLS" view -@ "$THREADS" -u -h -M -P -F 0x900 -L "$OUT_PREFIX.regions.bed" \
  ${READ_ARGS[@]+"${READ_ARGS[@]}"} -o "$WORK/regions.bam" "${INPUT_ARGS[@]}"
parts=("$WORK/regions.bam")
if [[ $INCLUDE_UNPLACED -eq 1 ]]; then
  log "Fetching unplaced unmapped reads"
  "$SAMTOOLS" view -@ "$THREADS" -u -h -F 0x900 ${READ_ARGS[@]+"${READ_ARGS[@]}"} \
    -o "$WORK/unplaced.bam" "${INPUT_ARGS[@]}" '*'
  parts+=("$WORK/unplaced.bam")
fi

log "Grouping by read name"
"$SAMTOOLS" cat "${parts[@]}" \
  | "$SAMTOOLS" collate -@ "$THREADS" -u -O - "$WORK/collate" \
  | "$SAMTOOLS" view -@ "$THREADS" -b -1 -o "$WORK/collated.bam" -

# Every record is fetched once even when the regions overlap or a mate falls
# inside another region; a duplicate record would surface as a third read with
# the same name. Check rather than assume.
# (Flag bits 0x40/0x80 are extracted arithmetically; BSD awk has no and().)
dup_names=$("$SAMTOOLS" view "$WORK/collated.bam" \
  | awk -F'\t' '{r=(int($2/64)%2)+2*(int($2/128)%2); k=$1":"r} k==prev{d++} {prev=k} END{print d+0}')
[[ "$dup_names" -eq 0 ]] || die "${dup_names} records were emitted twice; refusing to write FASTQs"

# ---- FASTQs per read group -------------------------------------------------
RG_LINES=()
while IFS= read -r line; do RG_LINES+=("$line"); done < <(grep '^@RG' "$WORK/header.sam" || true)
printf 'read_group\tpairs\tsingletons\tother\n' > "$OUT_PREFIX.prefilter_qc.tsv"
: > "$OUT_PREFIX.sequences.yaml"
total_pairs=0
total_single=0

write_fastqs() {  # $1 = RG id ("" for no read groups), $2 = file tag
  local rg=$1 tag=$2 select=()
  [[ -n "$rg" ]] && select=(-r "$rg")
  local r1="$OUT_PREFIX.$tag.R1.fastq.gz" r2="$OUT_PREFIX.$tag.R2.fastq.gz"
  # -O writes original qualities (OQ tag) when present, i.e. pre-BQSR values.
  "$SAMTOOLS" view -@ "$THREADS" -u ${select[@]+"${select[@]}"} "$WORK/collated.bam" \
    | "$SAMTOOLS" fastq -@ "$THREADS" -O -n -c 6 \
        -1 "$r1" -2 "$r2" -0 "$WORK/$tag.other.fq.gz" -s "$WORK/$tag.single.fq.gz" - \
        2> "$WORK/$tag.fastq.log"
  local pairs single other
  pairs=$(( $(gzip -dc "$r1" | wc -l) / 4 ))
  local pairs2=$(( $(gzip -dc "$r2" | wc -l) / 4 ))
  [[ $pairs -eq $pairs2 ]] || die "read group $tag: R1 has $pairs reads, R2 has $pairs2"
  single=$(( $(gzip -dc "$WORK/$tag.single.fq.gz" | wc -l) / 4 ))
  other=$(( $(gzip -dc "$WORK/$tag.other.fq.gz" | wc -l) / 4 ))
  printf '%s\t%s\t%s\t%s\n' "$tag" "$pairs" "$single" "$other" >> "$OUT_PREFIX.prefilter_qc.tsv"
  total_pairs=$((total_pairs + pairs))
  total_single=$((total_single + single))
  log "  $tag: $pairs pairs, $single singletons, $other unflagged"
}

yaml_entry() {  # $1 = file tag, $2 = @RG line or ""
  local tag=$1 rgline=$2 rg_escaped
  {
    echo "- sequence:"
    echo "    fastq1: \"$OUT_PREFIX.$tag.R1.fastq.gz\""
    echo "    fastq2: \"$OUT_PREFIX.$tag.R2.fastq.gz\""
    if [[ -n "$rgline" ]]; then
      rg_escaped=$(printf '%s' "$rgline" | sed 's/\t/\\t/g; s/"/\\"/g')
      echo "  readgroup: \"$rg_escaped\""
    else
      echo "  readgroup: \"@RG\\tID:1\\tSM:REPLACE_ME\\tLB:REPLACE_ME\\tPL:Illumina\""
    fi
  } >> "$OUT_PREFIX.sequences.yaml"
}

log "Writing FASTQs"
if [[ ${#RG_LINES[@]} -eq 0 ]]; then
  write_fastqs "" "all"
  yaml_entry "all" ""
else
  for line in "${RG_LINES[@]}"; do
    id=$(printf '%s' "$line" | tr '\t' '\n' | awk -F: '$1=="ID"{sub(/^ID:/,""); print; exit}')
    tag=$(printf '%s' "$id" | tr -c 'A-Za-z0-9._-' '_')
    write_fastqs "$id" "$tag"
    yaml_entry "$tag" "$line"
  done
fi

# Reads without a read-group tag are not selected by -r; count them so they
# cannot vanish silently.
if [[ ${#RG_LINES[@]} -gt 0 ]]; then
  no_rg=$("$SAMTOOLS" view -c "$WORK/collated.bam" \
    | awk -v want=$((2 * total_pairs + total_single)) '{print $1 - want}')
  [[ "$no_rg" -le 0 ]] || log "WARNING: ${no_rg} records carry no read group listed in the header and were not written"
fi

frac=$(awk -v s="$total_single" -v p="$total_pairs" 'BEGIN{ if (p==0) print 1; else printf "%.6f", s/(2*p) }')
log "Total: ${total_pairs} pairs, ${total_single} singletons (fraction ${frac})"
[[ $total_pairs -gt 0 ]] || die "no read pairs written"
awk -v f="$frac" -v m="$MAX_SINGLETON_FRACTION" 'BEGIN{exit !(f<=m)}' \
  || die "singleton fraction ${frac} exceeds ${MAX_SINGLETON_FRACTION}; the input may be truncated or not paired-end"
log "Done. Inputs for immuno: $OUT_PREFIX.sequences.yaml"

version 1.0

# Builds the region BED that restrictToRegions keeps when exome-design
# workflows run on whole-genome data: the union of the given interval lists,
# each interval padded, plus the MHC region and every HLA / chr6 alt contig in
# the reference. The HLA typers (OptiType, HLA-HD, PHLAT) read those contigs
# and chr6:29836259-33148325 from the alignments, so they must be kept.
task wgsRegionsBed {
  input {
    Array[File] interval_lists  # Picard interval_list (1-based) or BED (0-based)
    File reference_dict
    Int padding = 250
    String hla_region = "chr6:29836259-33148325"
    String hla_contig_regex = "HLA|HSCHR6_MHC|^chr6_.*_alt$"
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "mgibio/samtools-cwl:1.16.1"
    disks: "local-disk 10 HDD"
  }

  command <<<
    set -euo pipefail

    awk -F'\t' '$1=="@SQ"{n="";l="";for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4);if($i~/^LN:/)l=substr($i,4)}print n"\t"l}' \
      "~{reference_dict}" > contigs.tsv

    : > raw.bed
    for f in ~{sep=" " interval_lists}; do
      before=$(wc -l < raw.bed)
      if grep -q -m1 '^@SQ' "$f"; then
        awk -F'\t' '!/^@/ && NF>=3 {print $1"\t"($2-1)"\t"$3}' "$f" >> raw.bed
      else
        awk -F'\t' '!/^(#|track|browser)/ && NF>=3 {print $1"\t"$2"\t"$3}' "$f" >> raw.bed
      fi
      # A gzipped or space-delimited file parses to nothing; fail, not skip.
      [[ $(wc -l < raw.bed) -gt $before ]] || { echo "no tab-delimited intervals in $f" >&2; exit 1; }
    done

    HLA_REGION="~{hla_region}"
    [[ "$HLA_REGION" =~ ^[^:]+:[0-9]+-[0-9]+$ ]] || { echo "hla_region must be CHROM:START-END" >&2; exit 1; }
    hla_chrom=${HLA_REGION%%:*}
    hla_range=${HLA_REGION#*:}
    printf '%s\t%s\t%s\n' "$hla_chrom" "$(( ${hla_range%-*} - 1 ))" "${hla_range#*-}" >> raw.bed

    : > missing.txt
    awk -F'\t' -v OFS='\t' -v pad=~{padding} -v re='~{hla_contig_regex}' '
      FNR==NR { len[$1]=$2; ord[$1]=FNR; if ($1 ~ re) print FNR, $1, 0, $2; next }
      {
        if (!($1 in len)) { missing[$1]++; next }
        s = $2 - pad; if (s < 0) s = 0
        e = $3 + pad; if (e > len[$1]) e = len[$1]
        if (e > s) print ord[$1], $1, s, e
      }
      END { for (c in missing) print c"\t"missing[c] > "missing.txt" }
    ' contigs.tsv raw.bed \
      | sort -k1,1n -k3,3n \
      | awk -F'\t' -v OFS='\t' '
          $2!=c || $3>e { if (c!="") print c, s, e; c=$2; s=$3; e=$4; next }
          $4>e { e=$4 }
          END { if (c!="") print c, s, e }
        ' > regions.bed

    # A contig-name mismatch ("1" vs "chr1") would drop the whole exome while
    # the MHC/HLA lines keep regions.bed non-empty, so fail instead.
    if [[ -s missing.txt ]]; then
      echo "intervals on contigs absent from reference_dict (contig, count):" >&2
      sort missing.txt >&2
      exit 1
    fi
    [[ -s regions.bed ]] || { echo "no regions produced" >&2; exit 1; }
  >>>

  output {
    File regions_bed = "regions.bed"
  }
}

workflow wf {
  input {
    Array[File] interval_lists
    File reference_dict
    Int padding = 250
  }
  call wgsRegionsBed {
    input:
    interval_lists=interval_lists,
    reference_dict=reference_dict,
    padding=padding
  }
}

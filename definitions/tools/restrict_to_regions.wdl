version 1.0

# Keeps every read template (all records sharing a read name: both mates,
# supplementary and secondary alignments) that has a primary or supplementary
# record overlapping regions_bed, plus unmapped pairs if keep_unmapped_pairs.
#
# Filtering record by record (`samtools view -L`) splits pairs whose mates map
# on and off the regions, and the orphaned records break Picard SamToFastq,
# pair-aware duplicate marking and mate-aware callers. Selecting by read name
# keeps templates whole.
#
# The input must be grouped by read name, as sequenceAlignAndTag writes it
# (bwa-mem2 | samblaster). `samtools view -N` preserves that order, so
# mergeBams and markDuplicatesAndSort see the same layout as without filtering.
task restrictToRegions {
  input {
    File bam
    File regions_bed
    # The HLA typers read every unmapped read (`samtools view -f 4`), and an
    # exome run's BAM contains the unmapped pairs from the captured library.
    Boolean keep_unmapped_pairs = true
  }

  Int cores = 4
  Int space_needed_gb = 10 + ceil(1.5 * size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    cpu: cores
    docker: "mgibio/samtools-cwl:1.16.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outname = "restricted.bam"
  command <<<
    set -euo pipefail

    # Pass 1: names of templates with a primary or supplementary record in the
    # regions (a capture library contains such chimeric fragments too).
    # Without an index, -L tests each record against the BED while streaming.
    /usr/bin/samtools view -@ ~{cores} -F 0x100 -L "~{regions_bed}" "~{bam}" | cut -f 1 > names.txt
    if ~{keep_unmapped_pairs}; then
      /usr/bin/samtools view -@ ~{cores} -f 12 -F 0x900 "~{bam}" | cut -f 1 >> names.txt
    fi
    sort -u -S 8G names.txt -o names.txt
    [[ -s names.txt ]] || { echo "no reads overlap the regions" >&2; exit 1; }

    # Pass 2: every record of those templates, in input order.
    /usr/bin/samtools view -@ ~{cores} -b -N names.txt -o "~{outname}" "~{bam}"

    in_n=$(/usr/bin/samtools view -c -@ ~{cores} -F 0x900 "~{bam}")
    out_n=$(/usr/bin/samtools view -c -@ ~{cores} -F 0x900 "~{outname}")
    echo "restrictToRegions: kept ${out_n} of ${in_n} primary records ($(wc -l < names.txt) templates)" >&2
  >>>

  output {
    File restricted_bam = outname
  }
}

workflow wf {
  input {
    File bam
    File regions_bed
    Boolean keep_unmapped_pairs = true
  }
  call restrictToRegions {
    input:
    bam=bam,
    regions_bed=regions_bed,
    keep_unmapped_pairs=keep_unmapped_pairs
  }
}

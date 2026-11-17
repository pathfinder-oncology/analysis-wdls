version 1.0

# Keeps every read template (all records sharing a read name: both mates,
# supplementary and secondary alignments) that has a primary or supplementary
# record overlapping regions_bed, plus every unmapped read and its mate.
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
  }

  Int cores = 4
  Int space_needed_gb = 10 + ceil(1.5 * size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    cpu: cores
    docker: "mgibio/samtools-cwl:1.16.1"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  String outname = "restricted.bam"
  command <<<
    set -euo pipefail
    # Byte-order collation: deterministic, and faster than a UTF-8 locale.
    export LC_ALL=C

    # Pass 1: names of templates with a primary or supplementary record in the
    # regions (a capture library contains such chimeric fragments too).
    # Without an index, -L tests each record against the BED while streaming.
    /opt/samtools/bin/samtools view -@ ~{cores} -F 0x100 -L "~{regions_bed}" "~{bam}" | cut -f 1 > names.txt
    # Every unmapped read and its mate, wherever the mate maps: the HLA typers
    # read all unmapped reads (`samtools view -f 4`), and divergent HLA reads
    # end up here.
    /opt/samtools/bin/samtools view -@ ~{cores} -f 4 -F 0x900 "~{bam}" | cut -f 1 >> names.txt
    sort -u -S 8G names.txt -o names.txt
    [[ -s names.txt ]] || { echo "no reads overlap the regions" >&2; exit 1; }
    echo "restrictToRegions: keeping $(wc -l < names.txt) templates" >&2

    # Pass 2: every record of those templates, in input order.
    /opt/samtools/bin/samtools view -@ ~{cores} -b -N names.txt -o "~{outname}" "~{bam}"
  >>>

  output {
    File restricted_bam = outname
  }
}

workflow wf {
  input {
    File bam
    File regions_bed
  }
  call restrictToRegions {
    input:
    bam=bam,
    regions_bed=regions_bed
  }
}

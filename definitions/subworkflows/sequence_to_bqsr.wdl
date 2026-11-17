version 1.0

import "../types.wdl"  # !UnusedImport
import "../tools/sequence_align_and_tag.wdl" as saat
import "../tools/restrict_to_regions.wdl" as rtr
import "../tools/merge_bams.wdl" as mb
import "../tools/mark_duplicates_and_sort.wdl" as mdas
import "../tools/index_bam.wdl" as ib
import "../tools/bqsr.wdl" as b

workflow sequenceToBqsr {
  input {
    Array[SequenceData] unaligned

    Array[File] bqsr_known_sites
    Array[File] bqsr_known_sites_tbi

    TrimmingOptions? trimming

    File reference
    File reference_fai
    File reference_dict
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_0123

    String final_name = "final"

    # WGS mode: keep only read templates with an alignment (primary or supplementary) in these
    # regions (see tools/restrict_to_regions.wdl). Unset = no filtering.
    File? restrict_to_bed
  }

  scatter(seq_data in unaligned) {
    call saat.sequenceAlignAndTag {
      input:
      unaligned=seq_data,
      trimming=trimming,
      reference=reference,
      reference_alt=reference_alt,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_0123=reference_0123
    }

    if (defined(restrict_to_bed)) {
      call rtr.restrictToRegions {
        input:
        bam=sequenceAlignAndTag.aligned_bam,
        regions_bed=select_first([restrict_to_bed])
      }
    }
    File aligned_bam = select_first([restrictToRegions.restricted_bam, sequenceAlignAndTag.aligned_bam])
  }

  call mb.mergeBams {
    input:
    bams=aligned_bam,
    name=final_name
  }

  call mdas.markDuplicatesAndSort {
    input: bam=mergeBams.merged_bam
  }

  call b.doBqsr {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=markDuplicatesAndSort.sorted_bam,
    bam_bai=markDuplicatesAndSort.sorted_bam_bai,
    known_sites=bqsr_known_sites,
    known_sites_tbi=bqsr_known_sites_tbi,
    output_name=final_name
  }

  output {
    File final_bam =  doBqsr.bqsr_bam
    File final_bai =  doBqsr.bqsr_bai
    File final_bam_bai =  doBqsr.bqsr_bam_bai
    File mark_duplicates_metrics_file = markDuplicatesAndSort.metrics_file
  }
}

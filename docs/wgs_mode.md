# Running `immuno` on whole-genome DNA

`immuno` is designed for exome DNA. Its variant callers, CNVkit and QC are
restricted to `bait_intervals` / `target_intervals`. On whole-genome input,
however, every whole-BAM step (merge, duplicate marking, BQSR, CRAM
conversion, HLA extraction and metrics) processes the full genome, which is
17 times the padded bait and target territory (185 Mb of 3.2 Gb).

## Why region filtering breaks the workflow

Restricting an aligned BAM with `samtools view -L regions.bed` keeps each
record only if that record overlaps a region. When one mate of a pair maps
inside the regions and the other maps outside, the pair is split. On a 30x
genome, filtered to padded xGen targets, 7,474 of 107,810 templates (6.9%)
lost a mate. Picard `SamToFastq`, which `sequenceAlignAndTag` and
`sequenceToFastq` use for BAM input, then fails:

```
SAMFormatException: SAM validation error: ERROR::MATE_NOT_FOUND:Found 7474 unpaired mates
```

The fix is to select whole read templates (all records sharing a read name)
instead of individual records.

## Option 1: `wgs_mode` in the workflow

Set `immuno.wgs_mode: true`. The workflow then:

1. builds a region BED (`tools/wgs_regions_bed.wdl`): the union of
   `bait_intervals` and `target_intervals`, each padded by
   `wgs_region_padding` (default 250 bp), plus the MHC
   (chr6:29836259-33148325) and every HLA / chr6 alt contig in the reference;
2. aligns each read group as usual, then keeps every template that has a
   primary or supplementary alignment in those regions, plus fully unmapped pairs
   (`tools/restrict_to_regions.wdl`). The aligner output is grouped by read
   name, so the step reads it twice by streaming: once to collect read names
   with `-L`, then once to keep those names with `-N`. It needs no sort or
   index, and it preserves the input order;
3. continues unchanged: merge, duplicate marking, BQSR, calling, HLA typing.

The step applies to the tumor and normal alignments in `somaticExome` and to
the normal alignment in `germlineExome`. FASTQ and BAM inputs are handled
the same way. Alignment still runs on the full genome; every step after it
sees exome-sized data.

### Regions that must be kept

- **Padded baits and targets.** `target_interval_padding` (default 100) sets
  the calling window. Any read that overlaps the window overlaps the padded
  region, so the default 250 bp keeps every read a caller can see.
- **MHC, HLA contigs and unmapped reads.** OptiType, HLA-HD and PHLAT extract
  `chr6:29836259-33148325`, the HLA alt contigs and all unmapped reads
  (`samtools view -f 4`) from the CRAM. HLA reads that are too divergent to
  align are among the unmapped reads.

### What differs from an exome run

- BQSR is trained on the reads that survive the filter, the same as for an
  exome, instead of on the whole genome.
- The unaligned FDA metrics still describe the complete input FASTQs/BAMs. The
  aligned metrics describe the filtered alignments.
- Coverage around the targets is flat, without the bait-shaped falloff of a
  capture library. CNVkit off-target bins receive reads only from the padding.
- A template whose only in-region record is a secondary alignment is dropped.
  `wgs_mode` keeps templates with a supplementary alignment in the regions.
  The prefilter script (option 2) does not: their primary record lies outside
  the regions, `--fetch-pairs` returns the mate but not the primary, and the
  script reports the leftover read as a singleton.

## Option 2: `scripts/wgs_target_prefilter.sh` before the workflow

If the genome is already aligned (a vendor BAM/CRAM), run the prefilter before
`immuno` to skip full-genome alignment:

```bash
scripts/wgs_target_prefilter.sh \
  -i tumor_wgs.bam \
  -L IDT_xGen_Lockdown_Exome_v1.baits.plusEnsemblCDS.2.interval_list \
  -L IDT_xGen_Lockdown_Exome_v1.targets.plusEnsemblCDS.2.interval_list \
  -o out/tumor -@ 8
```

For each read group the script writes complete R1/R2 FASTQs, plus
`out/tumor.sequences.yaml`, which holds ready-to-paste `tumor_sequence`
entries with the original `@RG` lines. It uses `samtools view
--fetch-pairs`, so it needs a coordinate-sorted, indexed input and a samtools
newer than 1.16.1 (tested on 1.21). It exits non-zero if R1 and R2 disagree,
a record appears twice, or singletons exceed `--max-singleton-fraction`.
Because mate fetching makes random reads, run it against a local copy of the
BAM. Run against `gs://` URLs, it was impractically slow.

The FASTQs are realigned by `immuno`, so the source alignment can use any
GRCh38 reference whose contig names match the interval lists.

## Validation

Both approaches were checked on a 30x WGS normal slice: a chr17 window, part
of the MHC and 20,000 unplaced unmapped reads, 275,413 complete pairs in all.
With the xGen bait and target lists:

| Check | Result |
|---|---|
| Templates kept (script and `restrictToRegions`) | 127,810, identical sets |
| Templates with a primary record in the regions that were missing | 0 |
| Kept templates explained by neither the regions nor the unmapped bin | 0 |
| Incomplete pairs in the output | 0 |
| R1/R2 read names, line by line (script) | identical |
| Read-name grouping and order (`restrictToRegions`) | preserved |
| Picard `SamToFastq` (strict) on the `restrictToRegions` output | succeeds |
| Picard `SamToFastq` on the naive `samtools view -L` output | fails, 7,474 unpaired mates |

The padded region set covers 185 Mb, 5.75% of the reference.

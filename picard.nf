nextflow.enable.dsl=2
params.picardMarkDuplicates = "java -jar ~/working/packages/picard.jar MarkDuplicates"

// NOTE: by default the workflow uses `samtools` to remove duplicate reads
// change lines 12-13 of `ATACseq_worklow.nf` (import REMOVEDUPLICATEREADS) if
// you want to use the picard version below

process REMOVEDUPLICATEREADS {
    tag { sample }
    input:
     tuple val(sample), path(bam_dupes_ch)
     tuple val(sample), path(bai_dupes_ch)

    output:
      tuple val(sample), path("${bam_dupes_ch.baseName}.noDuplicates.bam"), emit: bam_nodupes_ch
      tuple val(sample), path("${bai_dupes_ch.baseName}.noDuplicates.bai"), emit: bai_nodupes_ch

    script:
    """
    ${params.picardMarkDuplicates} QUIET=true INPUT=${bam_dupes_ch} OUTPUT=${bam_dupes_ch.baseName}.noDuplicates.bam METRICS_FILE=${bam_dupes_ch}.duplicates.metrics REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
    """
}
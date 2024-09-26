nextflow.enable.dsl=2

process SORTBAM {
  tag { sample }
  input:
    tuple val(sample), path(unsorted_bam_ch)
    tuple val(sample), path(unsorted_bai_ch)

  output:
    tuple val(sample), path("${unsorted_bam_ch.baseName}.bam"), emit: sorted_bam_ch
    tuple val(sample), path("${unsorted_bai_ch.baseName}.bai"), emit: sorted_bai_ch

  script:
  """
  #!/usr/bin/env bash
  samtools sort -O bam ${unsorted_bam_ch} -o ${unsorted_bam_ch.baseName}.tmp.bam
  samtools index ${unsorted_bam_ch.baseName}.tmp.bam -o ${unsorted_bam_ch.baseName}.bai
  mv ${unsorted_bam_ch.baseName}.tmp.bam ${unsorted_bam_ch.baseName}
  """    
}

process REMOVEMITOREADS{
  tag { sample }
  input:
    tuple val(sample), path(bam_mito_ch)
    tuple val(sample), path(bai_mito_ch)

  output:
    tuple val(sample), path("${bam_mito_ch.baseName}.noChrM.bam"), emit:bam_noChrM_ch
    tuple val(sample), path("${bai_mito_ch.baseName}.noChrM.bai"), emit:bai_noChrM_ch

  script:
  """
  #!/usr/bin/env bash
  samtools view -h ${bam_mito_ch} | grep -v chrM | samtools sort -O bam -o ${bam_mito_ch.baseName}.noChrM.bam -T .
  samtools index ${bam_mito_ch.baseName}.noChrM.bam -o ${bam_mito_ch.baseName}.noChrM.bai
  """
}

process ADDREADGROUPS{
  tag { sample }
  input:
    tuple val(sample), path(bam_noRG_ch)
    tuple val(sample), path(bai_noRG_ch)

  output:
    tuple val(sample), path("${bam_noRG_ch.baseName}.RGadded.bam"), emit:bam_RGadded_ch
    tuple val(sample), path("${bai_noRG_ch.baseName}.RGadded.bai"), emit:bai_RGadded_ch

  script:
  """
  #!/usr/bin/env bash
  samtools addreplacerg -r "@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina\\tLB:Library.fa" -o ${bam_noRG_ch.baseName}.RGadded.bam ${bam_noRG_ch.baseName}.bam
  samtools index ${bam_noRG_ch.baseName}.RGadded.bam -o ${bai_noRG_ch.baseName}.RGadded.bai
  """
}


process REMOVEDUPLICATEREADS {
  tag { sample }
  input:
    tuple val(sample), path(bam_dupes_ch)
    tuple val(sample), path(bai_dupes_ch)

  output:
    tuple val(sample), path("${bam_dupes_ch.baseName}.noDuplicates.bam"), emit:bam_nodupes_ch
    tuple val(sample), path("${bai_dupes_ch.baseName}.noDuplicates.bai"), emit:bai_nodupes_ch

  script:
  """
  #!/usr/bin/env bash
  samtools view -h -b -f 2 -F 1548 -q 30 ${bam_dupes_ch} | samtools sort -O bam -o ${bam_dupes_ch.baseName}.noDuplicates.bam -T .
  samtools index ${bam_dupes_ch.baseName}.noDuplicates.bam -o ${bam_dupes_ch.baseName}.noDuplicates.bai
  """    
}
nextflow.enable.dsl=2

process ALIGNTOGENOME {
  tag { sample }
  input:
    tuple val(sample), path(raw_reads)
  
  output:
    tuple val(sample), path("*.bam"), emit: unsorted_bam_ch
    tuple val(sample), path("*.bai"), emit: unsorted_bai_ch

  script:
  def (read1, read2) = raw_reads
  """
  #!/usr/bin/env bash
  ${params.alignmentProgram} ${params.alignmentParams} ${params.referenceGenome} -1 ${read1} -2 ${read2} | samtools view -bS - > ${sample}.bam 

  samtools sort ${sample}.bam -o ${sample}.tmp.bam
  samtools index ${sample}.tmp.bam -o ${sample}.bai
  mv ${sample}.tmp.bam ${sample}.bam
  """
}
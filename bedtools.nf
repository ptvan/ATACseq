nextflow.enable.dsl=2

process FILTERBLACKLISTREGIONS {
  tag { sample }
  input:
    tuple val(sample), path(input_bam_ch)
    tuple val(sample), path(input_bai_ch)
    
    path(bed_file)

  output:      
    tuple val(sample), path("${input_bam_ch.baseName}.blacklisted.bam"), emit: blacklisted_bam_ch
    tuple val(sample), path("${input_bai_ch.baseName}.blacklisted.bai"), emit: blacklisted_bai_ch

  script:
  """
  bedtools intersect -nonamecheck -v -abam ${input_bam_ch} -b ${bed_file} | samtools sort -o ${input_bam_ch.baseName}.blacklisted.bam
  samtools index ${input_bam_ch.baseName}.blacklisted.bam -o ${input_bam_ch.baseName}.blacklisted.bai
  """    
}

process BAMTOBED {
  tag { sample }
  input:
    tuple val(sample), path(input_bam_ch)
    
  output:
    path("${bed_ch.baseName}.bed"), emit: bed_ch

  script:
  """
  bedtools bamtobed -i ${input_bam_ch} > ${bed_ch.baseName}.bed
  """    
}

process BAMTOBEDPE {
  tag { sample }
  input:
    tuple val(sample), path(input_bam_ch)
    tuple val(sample), path(input_bai_ch)
    
  output:
    tuple val(sample), path("${sample}.bedPE"), emit: bedpe_ch

  script:
  """
  bedtools bamtobed -i ${input_bam_ch} -bedpe > ${sample}.bedPE
  """    
}
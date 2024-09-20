nextflow.enable.dsl=2

process RANDSAMPLE {
  publishDir "${params.output}", mode:"copy", overwrite: true
  tag { sample }
  input:
    tuple val(sample), path(bam_ch)

  output:
    tuple val(sample), path("*.bed"), emit: randsampleoutput_ch
  
  script:
  """
  macs3 randsample \
  -i ${bam_ch} \
  -f BAMPE \
  -p 100 \
  -o ${sample}.bed

  """
}

process CALLPEAKS {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(bedPE_ch)

    output:
      tuple val(sample), path("${sample}_macs/*"), emit: callpeaks_ch

    script:
    """
    macs3 callpeak \
    -f BEDPE \
    --nomodel \
    --shift -37 \
    --extsize 73 \
    -g ${params.effectiveGenomeSize} \
    -B --broad \
    --keep-dup all \
    --cutoff-analysis -n ${sample} \
    -t ${bedPE_ch} \
    --outdir ${sample}_macs 2> ${sample}_macs3.log
    """
}

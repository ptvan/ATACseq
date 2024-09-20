nextflow.enable.dsl=2

process ALIGNMENTSIEVE {
    publishDir "${params.output}", mode:"copy", overwrite: true
    tag { sample }
    input:
      tuple val(sample), path(bam_ch)
      tuple val(sample), path(bai_ch)

    output:
      tuple val(sample), path("${bam_ch.baseName}.shifted.bam"), emit: shiftedbam_ch
      tuple val(sample), path("${bam_ch.baseName}.shifted.bai"), emit: shiftedbai_ch

    script:
    """
    alignmentSieve \
    --verbose \
    --ATACshift \
    --blackListFileName ${params.genomeBlacklist} \
    --bam ${bam_ch.baseName}.bam \
    -o ${bam_ch.baseName}.shifted.bam
    
    samtools sort ${bam_ch.baseName}.shifted.bam -o ${bam_ch.baseName}.shifted.tmp.bam
    samtools index ${bam_ch.baseName}.shifted.tmp.bam -o ${bam_ch.baseName}.shifted.bai
    mv ${bam_ch.baseName}.shifted.tmp.bam ${bam_ch.baseName}.shifted.bam 
    """
}

process BAMCOVERAGE {
    tag { sample }
    input:
      tuple val(sample), path(bam_ch)
      tuple val(sample), path(bai_ch)

    output:
      tuple val(sample), path("*.bw"), emit: coveragebigwig_ch

    script:
    """
    bamCoverage \
    --numberOfProcessors ${params.nCPUs} \
    --binSize 10 \
    --normalizeUsing BPM \
    --effectiveGenomeSize ${params.effectiveGenomeSize} \
    --bam ${bam_ch} \
    -o ${sample}_coverage_BPM.bw
    """
}

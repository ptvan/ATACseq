#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// params.suffix = "*_R{1,2}.trimmed.fastq.gz"
params.alignmentProgram = "bowtie2"
params.alignmentParams = "--local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x "
params.referenceGenome = "$HOME/working/Databases/GCRh38_ATACseq"
params.genomeBlacklist = "$HOME/working/raw_data/hg38.blacklist.bed.gz"
params.effectiveGenomeSize = 2862010578
params.nCPUs = 8
params.help = false

include { ALIGNTOGENOME } from './bowtie2.nf'
include { SORTBAM; REMOVEMITOREADS; ADDREADGROUPS; REMOVEDUPLICATEREADS } from './samtools.nf'
// include { REMOVEDUPLICATEREADS } from './picard.nf'
include { FILTERBLACKLISTREGIONS; BAMTOBED; BAMTOBEDPE } from './bedtools.nf'
include { ALIGNMENTSIEVE; BAMCOVERAGE } from './deeptools.nf'
include { CALLPEAKS } from './MACS.nf'

if (params.help || params.csv == null || params.output == null){
    helpMessage()
    exit 1
}

def helpMessage() {
    log.info"""
    Usage:
    nextflow run ATACseq_workflow.nf --csv <> --output <>
    
    Required Arguments:
      --csv          Path to samplesheet CSV
      --output       Path to store output files 
    
    """.stripIndent()
}



workflow {
    raw_reads = Channel.fromPath( params.csv )
        .splitCsv( header: true)
        .map { row -> tuple( row.sampleID, file(row.read1), file(row.read2) ) } 
        
    aligned_reads = ALIGNTOGENOME(raw_reads)
    noChrM_reads = REMOVEMITOREADS(aligned_reads)
    readgroup_reads = ADDREADGROUPS(noChrM_reads)
    noduplicates_reads = REMOVEDUPLICATEREADS(readgroup_reads)
    blacklisted_reads = FILTERBLACKLISTREGIONS(noduplicates_reads, "${params.genomeBlacklist}")
    bedPE = BAMTOBEDPE(blacklisted_reads)

    coordshifted_reads = ALIGNMENTSIEVE(blacklisted_reads)
    bamcoverage = BAMCOVERAGE(blacklisted_reads)
    
    peaks = CALLPEAKS(bedPE)
}


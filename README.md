# ATACSeq Nextflow workflow

This workflow takes FASTQs, performs alignment, QC, and peak-calling for ATACSeq analysis. Additional annotation
and analyses are found in [ATACSeq_postpeakcalling_workflow.R](https://github.com/ptvan/workflows/blob/master/R/ATACSeq_postpeakcalling_workflow.R)

## Prerequisites

- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://www.htslib.org)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [deeptools](https://deeptools.readthedocs.io/en/latest/)
- [MACS](https://github.com/macs3-project/MACS) currently `macs3` at the time of this workflow's completion
- [picard](https://broadinstitute.github.io/picard/) _optional_ if you prefer to use it instead of `samtools` for removing duplicate reads, see `picard.nf`

## Installation & setup

1. Install all tools above and verify that they are callable at the command prompt
2. Install bowtie2 index appropriate for your species after creating them from scratch or downloading from [here](https://benlangmead.github.io/aws-indexes/bowtie)  
3. Similarly, install blacklist files appropriate for your species from [here](https://github.com/Boyle-Lab/Blacklist)
4. Change the `params.xxx` section of `ATACseq_workflow.nf` to change tool parameters (or pass them in when calling)
5. Create `FASTQ_directory` where your input FASTQs reside and `output_directory` where you want outputs to go

## Running the workflow

`nextflow run ATACseq_workflow.nf -resume --input <path_to_FASTQ_directory> --output <output_directory>`
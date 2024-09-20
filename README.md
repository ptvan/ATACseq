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

## Preparation

1. Install bowtie2 index appropriate for your species after creating them from scratch or downloading from [here](https://benlangmead.github.io/aws-indexes/bowtie)  
2. Similarly, install blacklist files appropriate for your species from [here](https://github.com/Boyle-Lab/Blacklist) in a `Databases/` directory
3. Create `FASTQ_directory/` where your input FASTQs reside and `output_directory/` where you want outputs to go

After this you have two options:

## OPTION A : Run with tools installed on your system:

4a. Install all tools above locally and verify that they are callable from the command prompt 

5a. Change the `params.xxx` section of `ATACseq_workflow.nf` to change tool parameters (or pass them in when calling `nextflow run` in step 6a below)

6a. Run the workflow: 
```
nextflow run https://github.com/ptvan/ATACseq -r main --input <local_input_dir_containing_FASTQs> --output ~/<local_output_dir>
```

## OPTION B: Run with dockerized tools:

4b. Pull the docker containing all the tools from: [nulzilla/atacseq-dkr](https://hub.docker.com/repository/docker/nulzilla/atacseq-dkr/general): 
```
docker pull nulzilla/atacseq-dkr
```

5a. Start up the docker with your `FASTQ_directory/` and `Databases/` mounted : 

```
docker run -it -v ./:/pipeline/nextflow:ro -v <FASTQ_directory>:/root/working/raw_data -v <your_Databases_dir>:/root/working/Databases:ro <$DOCKER_IMAGE_ID> /bin/bash
```

6b. Run the workflow: 
```
nextflow run /pipeline/main.nf --input ~/working/raw_data/input/ --output ~/working/raw_data/output/
```
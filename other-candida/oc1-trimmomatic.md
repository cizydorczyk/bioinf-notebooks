# oc1 - Trimming raw reads

**Date:** 2025-06-22
**System:** Local
**Project:** other-candida

## Objective
Trim sequencing reads with Trimmomatic (v0.39). This also included running FastQC and multiqc on 
the trimmed reads.
 

## Inputs
- raw fastq files: `/media/user/DATA/other-candida/raw-fastq/`


## Outputs
- trimmed fastq files: `/media/user/DATA/other-candida/results/trimmed-fastq/`
- MultiQC report: `/media/user/DATA/other-candida/results/multiqc/multiqc_report.html`

## Commands/Pipeline
Trimmomatic was run through nextflow (although I am moving away from Nextflow). The nextflow 
script `/media/user/DATA/other-candida/run-qc.nf` was used (see below for its code). A few initial 
samples were tested with FastQC to determine the adapters present (NexteraPE-PE.fa).


## Notes
Initially, fastp was used to trim reads (run through nextflow). However, I found that when trying 
to use the trimmed reads for sourmash (to determine species), several files were truncated or had 
unexpected or missing EOF stuff going on...I do not know what happened, but indeed even `seqkit` 
`stats` was failing on some samples, suggesting real issues with the trimmed reads. Rather than 
waste time trying to figure out why fastp is messing up my reads, I just re-trimmed them with good 
ol' Trimmomatic.


## Scripts
- `run-qc.nf`:
```
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline parameters
 */

params.trimmomatic_adapters = "${projectDir}/trimmomatic-adapters/NexteraPE-PE.fa"
params.raw_fastq = "${projectDir}/raw-fastq-files/illumina-fastq/*{_1,_2}.fastq.gz"


samples_ch = Channel.fromFilePairs(params.raw_fastq)

process trimmomatic {
    
    publishDir "${projectDir}/results/trimmed-fastq/", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_t_1.fastq.gz"), path("${sample}_t_2.fastq.gz"), emit: trimmomatic_raw_ch


    script:
    """
    trimmomatic PE \\
        -threads ${task.cpus} \\
        ${reads[0]} ${reads[1]} \\
        ${sample}_t_1.fastq.gz ${sample}_u_1.fastq.gz \\
        ${sample}_t_2.fastq.gz ${sample}_u_2.fastq.gz \\
        CROP:150 ILLUMINACLIP:${params.trimmomatic_adapters}:2:30:10:8:TRUE SLIDINGWINDOW:4:5 MINLEN:31
    """

}

process fastqc {
    
    input:
    tuple val(sample), path(reads)

    output:
    path "${sample}_raw_fastqc", emit: trimmed_fastqc_ch

    script:
    """
    mkdir ${sample}_raw_fastqc
    fastqc -t ${task.cpus} -o ${sample}_raw_fastqc ${reads}
    """
}

process multiqc {

    publishDir "${projectDir}/results/multiqc/", mode: 'copy'

    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}


workflow {

    trimmomatic(samples_ch)
    trimmomatic.out.map { sample, reads1, reads2 -> [sample, [reads1, reads2]] }.set { trimmomatic_ch }

    fastqc(trimmomatic_ch)
    
    multiqc(fastqc.out.collect())

}
```


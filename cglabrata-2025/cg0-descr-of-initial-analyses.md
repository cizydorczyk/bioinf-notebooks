# pa00 - Blank Template

**Date:** 2025-06-22
**System:** Local or Cluster
**Project:** cglabrata-2025

## Objective
This is a description of the initial analyses that have been performed for C. glabrata.


## FastQC:

Input: a sample of raw fastq files
Output: not saved

FastQC (v0.12.1) was run on a few samples (not recorded) to determine read lengths and which 
adapters required trimming.

Reads were 150bp (requiring cropping of 151st base)
Adapters were Nextera transposase sequence


## Trimmomatic:

Input: all raw fastq files
Output: trimmed fastq files in ../c1-c2-trimmed-fastq/

Trimmomatic (v0.39) was run on all raw fastq files. Parameters were:

CROP:150
ILLUMINACLIP:/path/to/NexteraPE-PE.fa:2:30:10:8:true
SLIDINGWINDOW:4:5
MINLEN:31

e.g.:
trimmomatic PE \
    -threads 4 \
    /home/conrad/c1-raw-fastq/35699_R1.fastq.gz \
    /home/conrad/c1-raw-fastq/35699_R2.fastq.gz \
    /home/conrad/c1-trimmed-fastq/35699_t_R1.fastq.gz \
    /home/conrad/c1-trimmed-fastq/35699_u_R1.fastq.gz \
    /home/conrad/c1-trimmed-fastq/35699_t_R2.fastq.gz \
    /home/conrad/c1-trimmed-fastq/35699_u_R2.fastq.gz \
    CROP:150 \
    ILLUMINACLIP:/home/conrad/software/NexteraPE-PE.fa:2:30:10:8:true \
    SLIDINGWINDOW:4:5 \
    MINLEN:31


## SKESA:

Input: trimmed fastq ../c1-c2-trimmed-fastq/
Output: ../skesa-assemblies/

SKESA was run to quickly assemble genomes. SKESA was run with default settings. It was run on ARC 
and assemblies transferred to ../skesa-assemblies/.


## Kraken2:

Input: trimmed reads in ../c1-c2-trimmed-fastq/
Output: ../kraken-reports/

Kraken2 was run to confirm species as C. glabrata and check for contamination. 
Kraken2 was run using the PlusPFP-16 pre-built Kraken2 database downloaded from 
https://benlangmead.github.io/aws-indexes/k2.

All settings were otherwise default (e.g., kmer size).


## BUSCO

Input: skesa assemblies in ../skesa-assemblies/
Output: ../busco/

BUSCO was run to test assembly completeness, since SKESA is not a commonly used assembler for fungi. 
It seems to have done a decent job.

e.g.:
busco \
    -i ../skesa-assemblies/35687-skesa-assembly.fasta \
    -m genome \
    -l saccharomycetaceae_odb12 \
    --metaeuk \
    --out busco-output-dirs/35687 \
    --cpu 8


## QUAST

Input: ../skesa-assemblies/
Output: ../quast-skesa/

QUAST (v5.2.0) was run to collect assembly metrics.


## MultiQC

Input: kraken-reports/, quast-skesa/, busco/
Output: ../mlutiqc_data/

MultiQC (v1.27.1) was run to create a single report for all kraken, quast, and busco data.


## MLST

Input: ../c1-c2-trimmed-fastq/
Output: ../mlst/

Three MLST tools were run: PyMLST (v2.1.6, release 2.1.0), stringMLST (v0.6.3), and ARIBA (2.14.6). 
All outputs are in the parent 'mlst/' folder. All tools downloaded the Candida glabrata MLST scheme 
from PubMLST.

All tools were run with default settings.

Consensus was taken as 2/3 tools for an ST call.

e.g.:
### PyMLST:
claMLST search2 \
    --output pymlst-fastq-outputs/35606-pymlst-fastq.txt \
    --paired cglabrata \
    ../c1-c2-trimmed-fastq/35606_t_R1.fastq.gz \
    ../c1-c2-trimmed-fastq/35606_t_R2.fastq.gz

### stringMLST:
stringMLST.py --predict \
    -1 /home/user/candida/candida-glabrata/c1-c2-trimmed-fastq/35446_t_R1.fastq.gz \
    -2 /home/user/candida/candida-glabrata/c1-c2-trimmed-fastq/35446_t_R2.fastq.gz \
    -P cg_mlst/cglabrata \
    -p \
    -o /home/user/candida/candida-glabrata/stringmlst/results/35446-stringmlst-results.txt

### ARIBA:
ariba run \
    cg_mlst/ref_db/ \
    ../c1-c2-trimmed-fastq/35687_t_R1.fastq.gz \
    ../c1-c2-trimmed-fastq/35687_t_R2.fastq.gz \
    results/35687/


## PopPUNK

PopPUNK was run using SKESA assemblies to get an idea of diversity/etc. without running SNP calling 
on all samples against a single reference. It was run using the 'dbscan' model and the 'refine' model, 
as the partitioning of distances by just the dbscan model looked a bit off.

e.g.:
(poppunk-env) user@user-ThinkStation-P3-Tower:~$ poppunk --create-db --output test-db --r-files poppunk-input.txt --threads 32 --min-k 15 --max-k 31 --sketch-size 100000 --plot-fit 3

(poppunk-env) user@user-ThinkStation-P3-Tower:~$ poppunk --fit-model dbscan --ref-db test-db --threads 32

(poppunk-env) user@user-ThinkStation-P3-Tower:~$ poppunk --fit-model refine --ref-db test-db --model-dir test-db

(poppunk-env) user@user-ThinkStation-P3-Tower:~$ poppunk_visualise --ref-db test-db --output test-viz --microreact

# oc3 - Sourmash to get species

**Date:** 2025-07-18
**System:** Local
**Project:** other-candida

## Objective
Run sourmash to get species for each sample, because I do not have metadata at the moment. Also, 
good to confirm...Kraken2 seems to fail so trying sourmash.
 
Sourmash runs only on R1 or R2 files (could not figure out how to make run on both). So, here I am 
running it only on R1 files.


## Inputs
- Trimmed R1 fastq: `/media/user/DATA/other-candida/results/trimmed-fastq/`


## Outputs
- Sourmash output folder: `/media/user/DATA/other-candida/oc3-sourmash/`
- Sourmash classifications: `/media/user/DATA/other-candida/oc3-sourmash/tax/`
- Sourmash species summary: `/media/user/DATA/other-candida/oc3-sourmash/SPECIES-summary.txt`


## Commands/Pipeline
Commands were run using parallel. Below are some example commands:

### sourmash sketch
```
sourmash sketch dna -p scaled=10000,k=31 -o 
/media/user/DATA/other-candida/sourmash/sketches/25879.sig /media/user/DATA/other-candida/results/trimmed-fastq/25879_t_1.fastq.gz
```

### sourmash gather
```
sourmash gather /media/user/DATA/other-candida/sourmash/sketches/25879.sig 
/media/user/DATA/other-candida/sourmash/genbank-2022.03-fungi-k31.zip -o /media/user/DATA/other-candida/sourmash/gathers/25879-gather.csv
```

### sourmash tax
```
sourmash tax genome --gather-csv gathers/25879-gather.csv --taxonomy-csv 
genbank-2022.03-fungi.lineages.csv -o tax/25879-tax
```


## Scripts
Parallel commands for the above three commands (sketch, gather, tax) are found in:
- `/media/user/DATA/other-candida/oc3-sourmash/run-sourmash-sketch.txt`
- `/media/user/DATA/other-candida/oc3-sourmash/run_sourmash_gather.txt`
- `/media/user/DATA/other-candida/oc3-sourmash/run_sourmash_tax.txt`

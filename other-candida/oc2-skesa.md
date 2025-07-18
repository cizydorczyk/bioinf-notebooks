# oc2 - SKESA assemblies (for QC)

**Date:** 2025-07-18
**System:** Local
**Project:** other-candida

## Objective
Assemble all 486 (485 + 1 duplicate) genomes with SKESA (v2.5.1). Also run QUAST (v5.2.0).


## Inputs
- Trimmed fastq: `/media/user/DATA/other-candida/results/trimmed-fastq/`


## Outputs
- SKESA assemblies: `/media/user/DATA/other-candida/oc2-skesa-assemblies/`
- QUAST data: `/media/user/DATA/other-candida/oc2-quast/`

## Commands/Pipeline

### SKESA
SKESA was run using parallel; all parallel commands are found in 
`/media/user/DATA/other-candida/run-skesa.txt`.

An example command is:
```
skesa --reads /media/user/DATA/other-candida/results/trimmed-fastq/25879_t_1.fastq.gz,/media/user/DATA/other-candida/results/trimmed-fastq/25879_t_2.fastq.gz --contigs_out skesa-assemblies/25879-skesa-assembly.fasta --cores 4 --memory 16
```


### QUAST
```
(quast-5.2.0) user@user-ThinkStation-P3-Tower:/media/user/DATA/other-candida$ quast -m 1000 -t 32 --no-plots --no-icarus --no-html -o oc2-quast/ oc2-skesa-assemblies/*.fasta
```


## Notes
37515-S45 assembled poorly and appears to have little data; its partner sample 37515-S245 will be 
the one I continue with. I think this one is just some sort of mislabeling...

**This leaves 485 samples**

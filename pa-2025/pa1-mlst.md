# pa0 - Initial Entry

**Date:** 2025-06-22
**System:** Local & Cluster
**Project:** pa-2025

## Objective

Rerun of MLST for PA data from trimmed reads using `ARIBA` and `stringMLST`. `pyMLST` data will be re-used from 
a prior run in January 2025 (not documented I don't think...).

## Inputs

- Trimmed reads: `/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/`

### stringMLST

- stringMLST database is in `/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/pa_mlst/`
- `run_stringmlst.txt` contains individual stringMLST commands for all isolates
- `run_stringmlst.slurm` is a job script that runs the commands in the `.txt` file

### ARIBA (run locally and uploaded to arc)

- ARIBA database is in `/work/parkins_lab/project/conrad/pa-2025/mlst/ap1-ariba/pa_mlst/`
- `run_pa_ariba_mlst.txt` contains individual commands to run ARIBA for all PA genomes

### pyMLST

- *Not sure where the pyMLST db is but it would have been the PA db...*

## Outputs

### stringMLST

- General output directory: `/work/parkins_lab/projec/conrad/pa-2025/mlst/pa1-stringmlst/`
- Raw stringMLST output directory: 
`/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/stringmlst-output/`
- Summary file: `/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/stringmlst-summary.txt`

### ARIBA

- General output directory: `/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/`
- Raw ARIBA output directory: `/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/ariba-output/`
- Summary file: `/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/ariba-mlst-summary.txt`

### pyMLST

- General output directory: `/work/parkins_lab/project/conrad/pa-2025/mlst/pymlst-jan2025/`
- Summary file: `/work/parkins_lab/project/conrad/pa-2025/mlst/pymlst-jan2025/pymlst-summary.txt`

## Commands/Pipeline

### stringMLST

- Database generated using:
```
stringMLST.py --getMLST -P /work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/pa_mlst/ --species 
"Pseudomonas aeruginosa"
```

- Example run command:
```
stringMLST.py --predict --prefix 
/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/pa_mlst/Pseudomonas_aeruginosa -1 
/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/A001A074-P0137-10-05-2004_1.fastq.gz -2 
/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/A001A074-P0137-10-05-2004_2.fastq.gz -o 
/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/stringmlst-output/A001A074-P0137-10-05-2004-stringmlst-output.txt
```

- A summary file was generated from the outputs like so:
```
for i in $(cat /work/parkins_lab/project/conrad/pa-2025/good-pa-list-clean.txt); do echo -e "$(cat 
/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/stringmlst-output/${i}-stringmlst-output.txt | sed 
-n '2p')" >> /work/parkins_lab/project/conrad/pa-2025/mlst/pa1-stringmlst/stringmlst-summary.txt; done
```

### ARIBA

- ARIBA db was generated like so:
```
ariba pubmlstget "Pseudomonas aeruginosa" /work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/pa_mlst
```

- Example of an ARIBA command for one isolate:
```
ariba run --threads 2 /work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/pa_mlst/ref_db 
/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/A001A074-P0137-10-05-2004_1.fastq.gz 
/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/A001A074-P0137-10-05-2004_2.fastq.gz 
/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/ariba-output/A001A074-P0137-10-05-2004/
```

- A summary file was generated using a command similar to the following:
```
for i in $(cat /work/parkins_lab/project/conrad/pa-2025/good-pa-list.clean.txt); do echo -e "${i}\t$(cat 
/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/ariba-output/${i}/mlst_report.tsv | sed -n '2p')" >> 
/work/parkins_lab/project/conrad/pa-2025/mlst/pa1-ariba/ariba-mlst-summary.txt
```

### pyMLST

- *pyMLST commands were not recorded...see website for general commands; I assume standard commands were used.*

## Notes

## Related Analyses

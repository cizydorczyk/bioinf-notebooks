# pa00 - Blank Template

**Date:** 2025-06-22
**System:** Cluster
**Project:** pa-2025

## Objective
Assemble PA genomes with Unicycler v0.5.1.
 

## Inputs
- Trimmed reads: `/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/`


## Outputs
- Assemblies: `/work/parkins_lab/project/conrad/pa-2025/de-novo-assemblies/pa2-unicycler/`
- Assembly copies: 
`/work/parkins_lab/project/conrad/pa-2025/de-novo-assemblies/pa2-unicycler-assembly-copies/`
- Job log files: `/work/parkins_lab/project/conrad/pa-2025/de-novo-assemblies/pa2-logs/`


## Commands/Pipeline
Unicycler was run with default settings. See script below for details.

Assemblies were first copied with `cp` to include the isolate name in the assembly filename, and then these 
copies were moved to the assembly copies folder. This was done to facilitate using these for whatever other 
analyses.


## Notes



## Scripts

### Unicycler job array script:
```
#!/bin/bash
#SBATCH --job-name=unicycler_array
#SBATCH --partition=cpu2017-bf05  # FILL IN YOUR PARTITION HERE
#SBATCH --array=1-1427            # Adjust range based on number of samples
#SBATCH --cpus-per-task=28       # Adjust based on your needs
#SBATCH --mem=0               # Adjust memory as needed
#SBATCH --time=5:00:00         # Adjust time limit
#SBATCH --output=logs/unicycler_%A_%a.out
#SBATCH --error=logs/unicycler_%A_%a.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Load modules (adjust based on your system)

source activate unicycler-0.5.1

# Define input and output directories
INPUT_DIR="/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/"
OUTPUT_DIR="/work/parkins_lab/project/conrad/pa-2025/de-novo-assemblies/pa2-unicycler/"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Get sample list
SAMPLE_LIST="/work/parkins_lab/project/conrad/pa-2025/good-pa-list-clean.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

# Define input files
SHORT_R1="${INPUT_DIR}/${SAMPLE}_1.fastq.gz"
SHORT_R2="${INPUT_DIR}/${SAMPLE}_2.fastq.gz"

# Create sample-specific output directory
SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p ${SAMPLE_OUTPUT}

# Run Unicycler
unicycler -1 ${SHORT_R1} -2 ${SHORT_R2} \
          -o ${SAMPLE_OUTPUT} \
          --threads ${SLURM_CPUS_PER_TASK}
```


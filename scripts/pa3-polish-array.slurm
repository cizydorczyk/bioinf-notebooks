#!/bin/bash
#SBATCH --job-name=polish_array
#SBATCH --output=logs/polish_%A_%a.out
#SBATCH --error=logs/polish_%A_%a.err
#SBATCH --array=0-46
#SBATCH --cpus-per-task=14
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --partition=cpu2017-bf05  # <-- Set this to your actual partition

set -euo pipefail

# Load modules if necessary
# module load bwa-mem2 polypolish

source activate polypolish-env

# Read sample name
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /work/parkins_lab/project/conrad/pa-2025/snp-calling/reference-genomes/multi-patient-references.txt)

BWA=/home/conrad.izydorczyk/bwa-mem2-2.2.1_x64-linux/bwa-mem2-2.2.1_x64-linux/bwa-mem2
ASSEMBLY=/work/parkins_lab/project/conrad/pa-2025/de-novo-assemblies/pa2-unicycler-assembly-copies/${SAMPLE}-assembly.fasta
READ1=/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/${SAMPLE}_1.fastq.gz
READ2=/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/${SAMPLE}_2.fastq.gz
OUTDIR=/work/parkins_lab/project/conrad/pa-2025/snp-calling/reference-genomes/polished-assemblies/

mkdir -p "$OUTDIR"

# Index
$BWA index "$ASSEMBLY"

# Align
$BWA mem -t 14 -a "$ASSEMBLY" "$READ1" > "$OUTDIR/${SAMPLE}-R1.sam"
$BWA mem -t 14 -a "$ASSEMBLY" "$READ2" > "$OUTDIR/${SAMPLE}-R2.sam"

# Filter
polypolish filter \
  --in1 "$OUTDIR/${SAMPLE}-R1.sam" \
  --in2 "$OUTDIR/${SAMPLE}-R2.sam" \
  --out1 "$OUTDIR/${SAMPLE}-FILT-R1.sam" \
  --out2 "$OUTDIR/${SAMPLE}-FILT-R2.sam"

# Polish
polypolish polish "$ASSEMBLY" \
  "$OUTDIR/${SAMPLE}-FILT-R1.sam" "$OUTDIR/${SAMPLE}-FILT-R2.sam" \
  > "$OUTDIR/${SAMPLE}-polished-assembly.fasta"

# Cleanup
rm "$OUTDIR/${SAMPLE}-R1.sam" "$OUTDIR/${SAMPLE}-R2.sam" "$OUTDIR/${SAMPLE}-FILT-R1.sam" "$OUTDIR/${SAMPLE}-FILT-R2.sam"

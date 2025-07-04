#!/bin/bash

# Function to print usage/help/docstring
usage() {
    cat << EOF

NAME
    process-reference.sh - Assemble reference genome with SPAdes and run MitoFinder

SYNOPSIS
    process_fastq.sh -i <forward_reads> -j <reverse_reads> -o <output_dir> -p <prefix>

DESCRIPTION
    This script prepares a reference for SNP calling using mycosnp.
    It runs SPAdes to assemble the genome, removes short contigs, runs 
    MitoFinder, and filters out mitochondrial contigs from the assembly.

OPTIONS
    -i <file>       Forward reads FASTQ file
    -j <file>       Reverse reads FASTQ file
    -o <dir>        Output directory within which subdirectories will be
    		    created as required. Will be created if does not exist.
    -p <string>     Prefix for output files
    -t <num>        Number of threads to use. Required.
    -m <num>        Max memory for SPAdes in Gb. Required.
    -r <file>       Reference mitochondrioa seq for mitofinder in gb format.
    -h              Show this help message

EXAMPLE
    ./process_fastq.sh -i R1.fastq.gz -j R2.fastq.gz -o results/ -p sample

EOF
    exit 1
}

# Initialize variables
freads=""
rreads=""
output=""
prefix=""
threads=4
mem=32
MITOFINDER="/media/user/DATA/candida-glabrata/cg1-snp-calling/mitofinder_v1.4.2.sif"
SCRIPTDIR=$(dirname "$0")

# Parse options
while getopts ":i:j:o:p:s:t:m:r:h" opt; do
  case ${opt} in
    i)
      freads=$OPTARG
      ;;
    j)
      rreads=$OPTARG
      ;;
    o)
      output=$OPTARG
      ;;
    p)
      prefix=$OPTARG
      ;;
    t)
      threads=$OPTARG
      ;;
    m)
      mem=$OPTARG
      ;;
    r)
      ref=$OPTARG
      ;;
    h)
      usage
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Shift past processed options
shift $((OPTIND -1))

# Show help if no arguments were passed
if [[ "$OPTIND" -eq 1 ]]; then
    usage
fi

# Validate text arguments
if [[ -z "$freads" || -z "$rreads" || -z "$output" || -z "$prefix" || -z "$ref" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

# Validate numerical arguments
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: Threads must be a non-negative integer." >&2
    usage
fi

if ! [[ "$mem" =~ ^[0-9]+$ ]]; then
    echo "Error: Max memory must be a non-negative integer." >&2
    usage
fi

# Input files exist
[[ ! -f "$freads" ]] && echo "Error: Forward reads file not found: $freads" >&2 && usage
[[ ! -f "$rreads" ]] && echo "Error: Reverse reads file not found: $rreads" >&2 && usage

# Now proceed with your logic...
echo -e "\nForward reads: $freads"
echo "Reverse reads: $rreads"
echo "Output dir: $output"
echo "Prefix: $prefix"
echo "Using ${threads} thread(s)"
echo "Max SPAdes memory: ${mem}GB"

# Run SPAdes assembly:
SPADESOUT="${output}${prefix}-spades/"
SPADESASSEMBLY="${SPADESOUT}/${prefix}-filt-assembly.fasta"

spades.py --isolate -1 ${freads} -2 ${rreads} -t ${threads} -m ${mem} -o ${SPADESOUT}

# Filter out contigs <500 bp
bioawk -c fastx '{ if(length($seq) > 500) { print ">"$name; print $seq }}' ${SPADESOUT}/scaffolds.fasta > ${SPADESASSEMBLY}


# Run Polypolish:
POLYPOLISHDIR="${output}${prefix}-polypolish/"
FSAM="${output}${prefix}-polypolish/f_aln.sam"
RSAM="${output}${prefix}-polypolish/r_aln.sam"
FILTFSAM="${output}${prefix}-polypolish/filt_f_aln.sam"
FILTRSAM="${output}${prefix}-polypolish/filt_r_aln.sam"
POLISHEDFASTA="${output}${prefix}-polypolish/${prefix}-polished.fasta"

mkdir -p ${POLYPOLISHDIR}

bwa index ${SPADESASSEMBLY}
bwa mem -t ${threads} -a ${SPADESASSEMBLY} ${freads} > ${FSAM}
bwa mem -t ${threads} -a ${SPADESASSEMBLY} ${rreads} > ${RSAM}
polypolish filter --in1 ${FSAM} --in2 ${RSAM} --out1 ${FILTFSAM} --out2 ${FILTRSAM}
polypolish polish ${SPADESASSEMBLY} ${FILTFSAM} ${FILTRSAM} > ${POLISHEDFASTA}
rm ${FSAM} ${RSAM} ${FILTFSAM} ${FILTRSAM}

# Run MitoFinder

singularity run ${MITOFINDER} -j ${prefix} -a ${POLISHEDFASTA} -r ${ref} -o 3 -p ${threads} --new-genes

echo "Removing mitochondrial contigs..."
for i in ${prefix}/${prefix}_MitoFinder_mitfi_Final_Results/${prefix}_mtDNA_contig_[0-9].fasta; do grep '^>' ${i} | cut -c2- >> ${output}${prefix}-mito-contigs.txt; done

bioawk -c fastx -v ids_to_remove="${output}${prefix}-mito-contigs.txt" 'BEGIN{while((getline id < ids_to_remove) > 0) remove[id]=1} !(remove[$name]){print ">"$name; print$seq}' ${POLISHEDFASTA} > ${output}${prefix}-final-ref.fasta

rm -rf ${prefix}/ ${prefix}_MitoFinder.log




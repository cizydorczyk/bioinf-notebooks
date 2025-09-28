# pa4 - SNP Calling

**Date:** 2025-09-28
**System:** Cluster
**Project:** pa-2025

## Objective
Call SNPs per ST on both multi- and single-patient STs.

SNP calling steps include:
- snippy
- snippy-core
- snippy-clean
- iq-tree
- cfml
- maskrc-svg
- snp-dists

The `bacterial-genomics-pipeline-2.sh` contains these steps.
 

## Inputs
- Trimmed reads: `/work/parkins_lab/project/conrad/pa-2025/fastq-files/trimmed-fastq/`
- Reference genomes: 
`/work/parkins_lab/project/conrad/pa-2025/snp-calling/reference-genomes/final-filtered-assemblies/`
- List of multi-patient reference genomes: 
`/work/parkins_lab/project/conrad/pa-2025/snp-calling/reference-genomes/multi-patient-references.txt`
- List of single-patient reference genomes: 
`/work/parkins_lab/project/conrad/pa-2025/snp-calling/reference-genomes/single-patient-refrences.txt`


## Outputs
- Output folder per ST: `/work/parkins_lab/project/conrad/snp-calling/${ST}-snp-calling/`


## Commands/Pipeline
The `bacterial-genomics-pipeline-2.sh` bash script was used to run SNP calling.


## Notes


## Scripts
The `bacterial-genomics-pipeline-2.sh` script was used to run all SNP calling. It is copied below.

```
#!/bin/bash

# Bacterial Genomics Pipeline - Single ST Version
# Bash implementation

#set -euo pipefail

BOLD_CYAN='\e[1;36m'
BOLD_GREEN='\e[1;32m'
BOLD_RED='\e[1;31m'
NC='\e[0m' # Reset color

# Default parameters
SAMPLESHEET=""
OUTDIR="results"
ST_NAME=""
REFERENCE=""
SNIPPY_MINFRAC=0.8
SNIPPY_CLEANUP="yes"
SNIPPY_UNMAPPED="no"
IQTREE_MODEL="GTR+I+R"
IQTREE_BB=10000
REMOVE_REFERENCE=false
SNIPPY=false
IQTREE=false
CFML=false
CPUS=$(nproc)

# Help message
show_help() {
    cat << EOF
Bacterial Genomics Pipeline - Single ST Version (Bash)

Usage:
    $0 --samplesheet samples.csv --st_name ST1 --reference /path/to/reference.fasta --outdir results

Required parameters:
    --samplesheet       CSV file with columns: sample_id,reads_1,reads_2
    --st_name          Name of the sequence type (e.g., ST1, ST131, etc.)
    --reference        Path to reference genome FASTA file
	
Pipeline control options:
	--snippy			Run snippy (+ snippy-core and snippy-clean) (default: false)
	--iqtree			Run IQ-Tree (default: false)
	--cfml				Run CFML (default: false)
    
Optional parameters:
    --outdir            Output directory (default: results)
    --snippy_minfrac    Snippy minimum fraction (default: 0.8)
    --snippy_cleanup    Clean up intermediate files (default: yes)
    --snippy_unmapped   Keep unmapped reads (default: no)
    --iqtree_model      IQ-TREE model (default: GTR+I+R)
    --iqtree_bb         Bootstrap replicates (default: 10000)
    --remove_reference  Remove reference sequence from alignments (default: false)
    --cpus              Number of CPUs to use (default: all available)

Samplesheet format:
    sample_id,reads_1,reads_2
    sample001,/path/to/sample001_1.fastq.gz,/path/to/sample001_2.fastq.gz
    sample002,/path/to/sample002_1.fastq.gz,/path/to/sample002_2.fastq.gz
    sample003,/path/to/sample003_1.fastq.gz,/path/to/sample003_2.fastq.gz
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --samplesheet)
            SAMPLESHEET="$2"
            shift 2
            ;;
        --st_name)
            ST_NAME="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --snippy_minfrac)
            SNIPPY_MINFRAC="$2"
            shift 2
            ;;
        --snippy_cleanup)
            SNIPPY_CLEANUP="$2"
            shift 2
            ;;
        --snippy_unmapped)
            SNIPPY_UNMAPPED="$2"
            shift 2
            ;;
        --iqtree_model)
            IQTREE_MODEL="$2"
            shift 2
            ;;
        --iqtree_bb)
            IQTREE_BB="$2"
            shift 2
            ;;
        --remove_reference)
            REMOVE_REFERENCE=true
            shift 1
            ;;
        --snippy)
            SNIPPY=true
            shift 1
            ;;
        --iqtree)
            IQTREE=true
            shift 1
            ;;
        --cfml)
            CFML=true
            shift 1
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check required parameters
if [[ -z "$SAMPLESHEET" ]]; then
    echo "Error: Please provide a samplesheet with --samplesheet"
    exit 1
fi

if [[ -z "$ST_NAME" ]]; then
    echo "Error: Please provide a sequence type name with --st_name"
    exit 1
fi

if [[ -z "$REFERENCE" ]]; then
    echo "Error: Please provide a reference genome with --reference"
    exit 1
fi

# Check if at least one pipeline step is selected
if [[ "$SNIPPY" == false && "$IQTREE" == false && "$CFML" == false ]]; then
    echo "Error: Please select at least one pipeline step (--snippy, --iqtree, or --cfml)"
    exit 1
fi

# Check if required files exist
if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "Error: Samplesheet file not found: $SAMPLESHEET"
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: Reference file not found: $REFERENCE"
    exit 1
fi

# Create output directories
mkdir -p "$OUTDIR"/{raw-snippy-output,snippy-core,iqtree,cfml,maskrc-svg}

printf "${BOLD_CYAN}Starting bacterial genomics pipeline for $ST_NAME${NC}\n"
echo "Using $CPUS CPUs"

# Show which steps will be run
printf "Pipeline steps selected:\n"
if [[ "$SNIPPY" == true ]]; then
    printf "  - ${BOLD_CYAN}SNIPPY${NC} (+ snippy-core + snippy-clean)\n"
fi
if [[ "$IQTREE" == true ]]; then
    printf "  - ${BOLD_CYAN}IQ-TREE${NC}\n"
fi
if [[ "$CFML" == true ]]; then
    printf "  - ${BOLD_CYAN}ClonalFrameML${NC} (+ maskrc-svg + snp-dists)\n"
fi
echo ""

# Create log file:
LOGFILE="${OUTDIR}/${ST_NAME}.log"
printf "Logging commands to: ${LOGFILE}\n\n"

# Set up snippy parameters
CLEANUP_FLAG=""
UNMAPPED_FLAG=""
if [[ "$SNIPPY_CLEANUP" == "yes" ]]; then
    CLEANUP_FLAG="--cleanup"
fi
if [[ "$SNIPPY_UNMAPPED" == "yes" ]]; then
    UNMAPPED_FLAG="--unmapped"
fi

# Array to store sample directories for snippy-core
SAMPLE_DIRS=()

# Count total samples first (only if snippy is enabled)
if [[ "$SNIPPY" == true ]]; then
    TOTAL_SAMPLES=0
    echo "Counting valid samples in samplesheet..."

    {
        read -r header  # Skip header
        while IFS=',' read -r sample_id reads_1 reads_2; do
            [[ -z "$sample_id" ]] && continue

            # Clear the previous output (4 lines)
            if [[ $TOTAL_SAMPLES -gt 0 ]]; then
                for _ in {1..4}; do
                    printf "\033[1A"  # Move cursor up
                    printf "\033[2K"  # Clear entire line
                done
            fi

            echo "Checking sample: $sample_id"
            echo "  Reads 1: $reads_1"
            echo "  Reads 2: $reads_2"

            if [[ -f "$reads_1" && -f "$reads_2" ]]; then
                TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
                echo "  Valid sample found (total: $TOTAL_SAMPLES)"
            else
                echo "  Sample files not found"
            fi

            sleep 0.05  # Optional: slow down so you can see it
        done
        echo ""  # Final newline after loop
    } < "$SAMPLESHEET"

    if [[ $TOTAL_SAMPLES -eq 0 ]]; then
        echo "Error: No valid samples found in samplesheet"
        exit 1
    fi

    echo "Found $TOTAL_SAMPLES valid samples"
fi

# SNIPPY SECTION
if [[ "$SNIPPY" == true ]]; then
    # Process each sample with snippy (parallelized)
    printf "\n${BOLD_CYAN}SNIPPY:${NC}\t\t\tRunning snippy on individual samples...\n"

    # Calculate number of parallel jobs (assuming CPUS is a multiple of 4)
    PARALLEL_JOBS=$((CPUS / 4))

    # Create a function to process a single sample
    process_sample() {
        local sample_id="$1"
        local reads_1="$2"
        local reads_2="$3"
        local reference="$4"
        local outdir="$5"
        local minfrac="$6"
        local cleanup_flag="$7"
        local unmapped_flag="$8"
        local logfile="$9"
        
        # Create sample output directory
        sample_outdir="$outdir/raw-snippy-output/$sample_id"
        
        # Run snippy (suppress output)
        snippy --R1 "$reads_1" --R2 "$reads_2" --ref "$reference" \
               --outdir "$sample_outdir" --cpus 4 \
               --minfrac "$minfrac" $cleanup_flag $unmapped_flag &>> "${logfile}"
        
        # Output the sample directory for collection
        echo "$sample_outdir"
    }

    # Export the function and variables for parallel
    export -f process_sample
    export REFERENCE OUTDIR SNIPPY_MINFRAC CLEANUP_FLAG UNMAPPED_FLAG LOGFILE

    # Create a temporary file to store sample directories
    TEMP_SAMPLE_DIRS=$(mktemp)

    # Process samples in parallel with progress tracking
    {
        read -r header  # Skip header line
        while IFS=',' read -r sample_id reads_1 reads_2; do
            # Skip empty lines
            if [[ -z "$sample_id" ]]; then
                continue
            fi
            
            # Output the sample info for parallel processing
            echo "$sample_id,$reads_1,$reads_2"
        done
    } < "$SAMPLESHEET" | \
    parallel -j "$PARALLEL_JOBS" --colsep ',' --bar \
        'process_sample {1} {2} {3} "$REFERENCE" "$OUTDIR" "$SNIPPY_MINFRAC" "$CLEANUP_FLAG" "$UNMAPPED_FLAG" "$LOGFILE"' \
        > "$TEMP_SAMPLE_DIRS"

    # Read the sample directories back into the array
    while IFS= read -r sample_dir; do
        SAMPLE_DIRS+=("$sample_dir")
    done < "$TEMP_SAMPLE_DIRS"

    # Clean up temporary file
    rm "$TEMP_SAMPLE_DIRS"

    # Check if we have any valid samples
    if [[ ${#SAMPLE_DIRS[@]} -eq 0 ]]; then
        echo "Error: No valid samples found or processed"
        exit 1
    fi

    printf "${BOLD_CYAN}SNIPPY:${NC}\t\t\tProcessed ${#SAMPLE_DIRS[@]} samples\n"

    # Run snippy-core
    printf "${BOLD_CYAN}SNIPPY-CORE:${NC}\t\tProcessing core genome alignment...\n"
    cd "$OUTDIR/snippy-core"
    snippy-core --prefix "$ST_NAME" --ref "$REFERENCE" "${SAMPLE_DIRS[@]}" &>> "${LOGFILE}"
    cd - > /dev/null
    printf "\r${BOLD_CYAN}SNIPPY-CORE:${NC}\t\tCompleted core genome alignment\n"

    # Check if snippy-core output exists
    if [[ ! -f "$OUTDIR/snippy-core/${ST_NAME}.full.aln" ]]; then
        echo "Error: snippy-core failed to produce alignment"
        exit 1
    fi

    # Clean alignments
    printf "${BOLD_CYAN}SNIPPY-CLEAN:${NC}\t\tCleaning alignments...\n"
    snippy-clean_full_aln "$OUTDIR/snippy-core/${ST_NAME}.full.aln" > "$OUTDIR/snippy-core/${ST_NAME}.clean.full.aln"
    printf "\r${BOLD_CYAN}SNIPPY-CLEAN:${NC}\t\tCompleted cleaning alignments\n"

    # Choose alignment for downstream analysis
    ALIGNMENT_FILE="$OUTDIR/snippy-core/${ST_NAME}.clean.full.aln"

    # Remove reference if requested
    if [[ "$REMOVE_REFERENCE" == true ]]; then
        printf "${BOLD_CYAN}REMOVE-REF:${NC}\t\tRemoving reference sequence from alignment...\n"
        python3 << EOF
with open("$ALIGNMENT_FILE", 'r') as infile, open("$OUTDIR/snippy-core/${ST_NAME}.no_ref.clean.full.aln", 'w') as outfile:
    skip_seq = False
    for line in infile:
        if line.startswith('>Reference'):
            skip_seq = True
            continue
        elif line.startswith('>') and skip_seq:
            skip_seq = False
        
        if not skip_seq:
            outfile.write(line)
EOF
        ALIGNMENT_FILE="$OUTDIR/snippy-core/${ST_NAME}.no_ref.clean.full.aln"
        printf "\r${BOLD_CYAN}REMOVE-REF:${NC}\t\tCompleted removing reference sequence\n"
    fi
fi

# IQ-TREE SECTION
if [[ "$IQTREE" == true ]]; then
    # Check if we have an alignment file (either from snippy or user should provide one)
    if [[ -z "$ALIGNMENT_FILE" ]]; then
        # Look for existing alignment file
        if [[ -f "$OUTDIR/snippy-core/${ST_NAME}.clean.full.aln" ]]; then
            ALIGNMENT_FILE="$OUTDIR/snippy-core/${ST_NAME}.clean.full.aln"
        elif [[ -f "$OUTDIR/snippy-core/${ST_NAME}.no_ref.clean.full.aln" ]]; then
            ALIGNMENT_FILE="$OUTDIR/snippy-core/${ST_NAME}.no_ref.clean.full.aln"
        else
            echo "Error: No alignment file found for IQ-TREE. Run --snippy first or provide an alignment file."
            exit 1
        fi
    fi

    # Run IQ-TREE
    printf "${BOLD_CYAN}IQTREE:${NC}\t\t\tRunning...\n"
    cd "$OUTDIR/iqtree"

    # Count sequences in alignment
    NUM_SEQS=$(grep -c "^>" "$ALIGNMENT_FILE")

    # Set bootstrap parameter based on number of sequences
    if [[ $NUM_SEQS -ge 4 ]]; then
        BB_PARAM="-bb $IQTREE_BB"
    else
        BB_PARAM=""
    fi

    # Run IQ-TREE (suppress output)
    iqtree -s "$ALIGNMENT_FILE" -m "$IQTREE_MODEL" $BB_PARAM -nt "$CPUS" -pre "$ST_NAME" &>> "${LOGFILE}"

    cd - > /dev/null
    printf "\r${BOLD_CYAN}IQTREE:${NC}\t\t\tCompleted\n"
fi

# CFML SECTION
if [[ "$CFML" == true ]]; then
    # Check if we have an alignment file
    if [[ -z "$ALIGNMENT_FILE" ]]; then
        # Look for existing alignment file
        if [[ -f "$OUTDIR/snippy-core/${ST_NAME}.clean.full.aln" ]]; then
            ALIGNMENT_FILE="$OUTDIR/snippy-core/${ST_NAME}.clean.full.aln"
        elif [[ -f "$OUTDIR/snippy-core/${ST_NAME}.no_ref.clean.full.aln" ]]; then
            ALIGNMENT_FILE="$OUTDIR/snippy-core/${ST_NAME}.no_ref.clean.full.aln"
        else
            echo "Error: No alignment file found for CFML. Run --snippy first or provide an alignment file."
            exit 1
        fi
    fi

    # Run ClonalFrameML
    printf "${BOLD_CYAN}CFML:${NC}\t\t\tRunning...\n"
    cd "$OUTDIR/cfml"

    # Find the best tree file (prefer .contree, fallback to .treefile)
    TREE_FILE=""
    if [[ -f "$OUTDIR/iqtree/${ST_NAME}.treefile" ]]; then
        TREE_FILE="$OUTDIR/iqtree/${ST_NAME}.treefile"
    fi

    if [[ -z "$TREE_FILE" ]]; then
        echo -e "\nError: No tree file found from IQ-TREE. Run --iqtree first."
        exit 1
    fi

    # Run ClonalFrameML (suppress output)
    ClonalFrameML "$TREE_FILE" "$ALIGNMENT_FILE" "$ST_NAME" -em true -show_progress false &>> "${LOGFILE}"

    cd - > /dev/null
    printf "\r${BOLD_CYAN}CFML:${NC}\t\t\tCompleted\n"

    # Run maskrc-svg
    printf "${BOLD_CYAN}MASKRC-SVG:${NC}\t\tMasking recombinant regions...\n"
    cd "$OUTDIR/maskrc-svg"

    # Copy ClonalFrameML output files to current directory
    cp "$OUTDIR/cfml/${ST_NAME}".* . 2>/dev/null || true

    # Run maskrc-svg (suppress output)
    maskrc-svg.py --symbol N --aln "$ALIGNMENT_FILE" \
                  --out "${ST_NAME}.rc_masked.clean.full.aln" \
                  --regions "${ST_NAME}.rc_masked.regions.txt" \
                  "$ST_NAME" &>> "${LOGFILE}"

    cd - > /dev/null
    printf "\r${BOLD_CYAN}MASKRC-SVG:${NC}\t\tCompleted masking recombinant regions\n"

    # Run snp-dists
    printf "${BOLD_CYAN}SNP-DISTS:${NC}\t\tCalculating SNP distance matrix...\n"
    cd "$OUTDIR/maskrc-svg"

    snp-dists "${ST_NAME}.rc_masked.clean.full.aln" > "${ST_NAME}.rc_masked.clean.full.aln.snpdists_matrix.txt" 2>/dev/null

    cd - > /dev/null
    printf "\r${BOLD_CYAN}SNP-DISTS:${NC}\t\tCompleted SNP distance matrix\n\n"
fi

printf "${BOLD_GREEN}Pipeline completed successfully!${NC}\n"
printf "${BOLD_GREEN}Results are in: $OUTDIR${NC}\n"

# Show relevant output files based on what was run
if [[ "$SNIPPY" == true ]]; then
    printf "${BOLD_GREEN}Core genome alignment: $OUTDIR/snippy-core/${ST_NAME}.clean.full.aln${NC}\n"
fi
if [[ "$IQTREE" == true ]]; then
    printf "${BOLD_GREEN}Phylogenetic tree: $OUTDIR/iqtree/${ST_NAME}.treefile${NC}\n"
fi
if [[ "$CFML" == true ]]; then
    printf "${BOLD_GREEN}SNP distance matrix: $OUTDIR/maskrc-svg/${ST_NAME}.rc_masked.clean.full.aln.snpdists_matrix.txt${NC}\n"
fi
```

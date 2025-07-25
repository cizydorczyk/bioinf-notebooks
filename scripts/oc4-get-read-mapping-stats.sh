#!/usr/bin/env bash
set -euo pipefail

# Usage check
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <sorted.bam> <sample_name> <output_file.tsv>"
    exit 1
fi

BAM="$1"
SAMPLE="$2"
OUTFILE="$3"
DEPTH_FILE=$(mktemp)

COL='\033[38;5;222m'
RESET='\033[0m'

printf "\n${COL}Beginning analysis for sample: ${SAMPLE}${RESET}\n\n"
printf "${COL}\tInput BAM:\t${BAM}${RESET}\n"
printf "${COL}\tOutput file:\t${OUTFILE}${RESET}\n\n"

# Generate per-base depth file
printf "${COL}\tCalculating genome sequencing depth...${RESET}\n"
~/software/bedtools2/bin/bedtools genomecov -d -ibam "$BAM" > "$DEPTH_FILE"

# Compute summary stats per contig (min, Q1, median, mean, Q3, max)
printf "${COL}\tCalculating genome depth statistics...${RESET}\n"
STATS_FILE=$(mktemp)
datamash -g 1 min 3 q1 3 median 3 mean 3 q3 3 max 3 < "$DEPTH_FILE" > "$STATS_FILE"

# Compute % ≥1x, % ≥10x, and number of zero-coverage positions per contig
printf "${COL}\tCalculating genome coverage statistics...${RESET}\n"
COV_FILE=$(mktemp)
awk '{
    total[$1]++
    if ($3 == 0) zero[$1]++
    if ($3 >= 1) cov1[$1]++
    if ($3 >= 10) cov10[$1]++
}
END {
    for (c in total) {
        pct1 = 100 * cov1[c] / total[c]
        pct10 = 100 * cov10[c] / total[c]
        z = (zero[c] ? zero[c] : 0)
        printf "%s\t%.4f\t%.4f\t%d\n", c, pct1, pct10, z
    }
}' "$DEPTH_FILE" > "$COV_FILE"

printf "${COL}\tCreating final report...${RESET}\n"

# Join stats and coverage
JOINED_FILE=$(mktemp)
join -t $'\t' <(sort "$STATS_FILE") <(sort "$COV_FILE") > "$JOINED_FILE"

# Add sample name column
FINAL_FILE=$(mktemp)
awk -v sample="$SAMPLE" 'BEGIN { OFS="\t" } { print sample, $0 }' "$JOINED_FILE" > "$FINAL_FILE"

# Add header if output file doesn't exist
if [[ ! -f "$OUTFILE" ]]; then
    echo -e "Sample\tContig\tMin\tQ1\tMedian\tMean\tQ3\tMax\tPct_≥1x\tPct_≥10x\tZero_depth_positions" > "$OUTFILE"
fi

# Append results to output file
cat "$FINAL_FILE" >> "$OUTFILE"

# Clean up
printf "${COL}\tCleaning up...${RESET}\n"
rm -f "$DEPTH_FILE" "$STATS_FILE" "$COV_FILE" "$JOINED_FILE" "$FINAL_FILE"

printf "${COL}Finished analysis for sample: ${SAMPLE}${RESET}\n"


# cg1 - Post-SNP Calling Processing

**Date:** 2025-06-22
**System:** Local
**Project:** cgla 2025

## Objective
This entry records post snp calling analyses to interpret SNP calling results.

 
## Inputs

### Build new trees with IQ-Tree
- SNP alignments in: 
`/media/user/DATA/candida-glabrata/cg1-snp-calling/mycosnp-output/${ST}-mycosnp/combined/vcftofasta/vcf-to-fasta.fasta`

### Process SNP distances
- SNP distances in: 
`/media/user/DATA/candida-glabrata/cg1-snp-calling/mycosnp-output/${ST}-mycosnp/combined/snpdists/combined.tsv`
- Phylogeines: `/media/user/DATA/candida-glabrata/cg1-snp-calling/iqtree-reruns/${ST}.treefile`

### Plot phylogenies with metadata
- Phylogenies: `/media/user/DATA/candida-glabrata/cg1-snp-calling/iqtree-reruns/${ST}.treefile`
- Public genomes metadata: `/media/user/DATA/candida-glabrata-public-genomes/public-genomes-metadata.txt`
- Private genomes metatdata: `/media/user/DATA/candida-glabrata/CIHR_Cgla_extracted.xlsx`
- Public genomes MLST data: `/media/user/DATA/candida-glabrata-public-genomes/public-genomes-mlst-summary.txt`
- Public genomes ST lists: `/media/user/DATA/candida-glabrata-public-genomes/public-st-lists`
- Private genomes ST lists: `/media/user/DATA/candida-glabrata/st-lists`


## Outputs

### IQ-Tree
- Phylogenies in: `/media/user/DATA/candida-glabrata/cg1-snp-calling/iqtree-reruns/`

### Process SNP distances
- Long format SNP distances files in: `../cg1-snp-calling/long-formatted-snp-distances/`
    - This includes two sets of files:
        - ALL snp distances
        - SNP distances <=100
- HTML Rmd report: `../cg1-snp-calling/parse-snp-distance-R/cg1-parse-snp-distances.nb.html`
    - This report provides for each ST:
        - A heatmap of SNP distances
        - SNP distance distribution
        - Table of SNP distances <= 100
        - Quick phylogeny

### Plot phylogenies with metadata
- Annotated phylogenies in: `../cg1-snp-calling/iqtree-reruns/${ST}_tree.pdf`


## Commands/Pipeline

### IQ-Tree
IQ-Tree was re-run on SNP alignments to build new trees, as those output by mycosnp have a lot of 
polytomies. Some of these can be resolved by running IQ-Tree properly (i.e., with the +ASC 
option). Here, IQ-Tree was run with the substitution model set to `-m MFP+ASC+R`.

For most STs, this was done in a loop:
```
(i18-snp-pipeline-env) user@user-ThinkStation-P3-Tower:/media/user/DATA/candida-glabrata/cg1-snp-calling/iqtree-reruns$ for i in $(cat ../snp-calling-st-list.txt); do iqtree -s ../mycosnp-output/${i}-mycosnp/combined/vcf-to-fasta/vcf-to-fasta.fasta -m MFP+ASC+R -nt 16 -bb 10000 --prefix  ${i}; done
```

For a few STs, this initial run failed, with IQ-Tree complaining there were invariant sites in 
the SNP alignment. I am not sure why this would be the case, except for maybe a case where all 
isolates have a base except a few have 'N'? Unsure. Anyway, IQ-Tree output a variant sites 
alignment that I then used to run it again, this time successfully. A sample command follows:
```
(i18-snp-pipeline-env) user@user-ThinkStation-P3-Tower:/media/user/DATA/candida-glabrata/cg1-snp-calling/iqtree-reruns$ iqtree -s ../mycosnp-output/st47-mycosnp/combined/vcf-to-fasta/vcf-to-fasta.fasta -m MFP+ASC+R -nt 16 --prefix st47
```

### Process SNP distances
SNP distances were processed using the `../cg1-snp-calling/parse-snp-distances-R/cg1-parse-snp-distances.Rmd` script.

### Plot phylogenies with metadata
Phylogenies were plotted using the 
`../cg1-snp-calling/parse-snp-distance-R/cg1-plot-phylogenies.Rmd` R script.

## Notes



## Scripts

**cg1-parse-snp-distance.Rmd**: `/media/user/DATA/candida-glabrata/cg1-snp-calling/parse-snp-distances-R/cg1-parse-snp-distances.Rmd`

**cg1-plot-phylogenies.Rmd**: `/media/user/DATA/candida-glabrata/cg1-snp-calling/parse-snp-distances-R/cg1-plot-phylogenies.Rmd`

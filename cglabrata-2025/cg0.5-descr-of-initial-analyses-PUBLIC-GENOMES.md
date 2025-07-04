# pa00 - Blank Template

**Date:** 2025-06-22
**System:** Local
**Project:** cgla 2025

## Objective
Process public C. glabrata genomes to prepare for SNP calling.
 
The main project directory for this entry is the entire folder 
`/media/user/DATA/candida-glabrata-public-genomes/`.

## Inputs
- Public genomes metadata: `/media/user/DATA/candida-glabrata-public-genomes/SraRunTable.csv`
- Public genomes accessions list: 
`/media/user/DATA/candida-glabrata-public-genomes/SRR_Acc_List.csv`

## Commands/Pipeline
This was inspired by this paper: https://www.nature.com/articles/s41564-023-01547-z#Sec9

NCBI SRA run selector was used to find and download an accessions list of 1439 Candida glabrata (Nakaseomyces glabratus) fastsq read sets.
(I don't remember the exact search and filter options...I am not getting 1439 when I try it today but am getting 1425...(June 3, 2025)...not good but not the end of the world, I suppose)...
- Today's filters (and assumed previous filters as well...)
    - DATASTORE TYPE fastq
    - Library Layout PAIRED
    - Library Source GENOMIC
    - Organism  Nakaseomyces glabratus (none of the other options)
This was downloaded May 5, 2025. Either way, I have the accessions and run selector table downloaded from May 5, so can at least point to it and say these were the runs/genomes I recovered that day.

I then used the `nf-core/fetchngs` tool to download these genomes. **I successfully downloaded** 
**1435 genomes fastq sets.**
- 4 failed to download and fall out here.

I then trimmed all downloaded reads using fastp (`run-fastp.nf` file). Output is in 
`/media/user/DATA/candida-glabrata-public-genomes/results/trimmed-fastq/`.

I then ran ARIBA mlst (`run-ariba.nf` file). Output is in 
`/media/user/DATA/candida-glabrata-public-genomes/results/ARIBA-MLST/`.
- 5 genomes failed ARIBA due to misformatted fastq files. So these 5 genomes fall out here.

**This leaves me with 1430 Candida glabrata genomes that had MLST attempted by ARIBA.**
**Not all have MLST assigned...but most do.**

Nothing further has been done with these reads as of this entry.

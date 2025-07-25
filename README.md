# Bioinf Notebooks

This will be a new way of keeping track of analyses.

It avoids having to log in all the time to Benchling and copy and paste code there. Rather, I can 
simply create a folder for each project here and store `.md` files in it.

## Scripts directory
This directory contains useful scripts (mostly bash) that can be re-used in the future, and is 
simply a good place to keep track of scripts I have used.

These scripts are generally pipeline scripts, where tools are chained together; nothing fancy, but 
useful nonetheless, and it keeps track of what steps were performed, how, and their order.

### oc4-get-read-mapping-stats.sh
This script takes as input a sorted/indexed bam file, isolate/sample name, and an output file (can 
exist already; will be appended to).

It first runs `bedtools genomecov` (edit script to add correct path to bedtools) to get sequencing 
depths against the reference.

Next, it uses `datamash` (must be in path) and `awk` to calculate depth and coverage statistics, 
respectively.

It then puts these into a single output file, labeled with isolate/sample so that multiple samples 
can share a single output file for downstream analysis.

Notably, the script removes all intermediate files, so no space will be taken up by large bedtools 
genomecov intermediates!


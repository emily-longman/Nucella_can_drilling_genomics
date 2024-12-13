# Nucella_can_drilling_genomics

## Project Goal

Assemble and annotate a genome for Nucella canaliculata. Identify candidate loci associated with drilling phenotypes in Nucella canaliculata by utilizing a population of dogwhelks that contains a mix of drilling phenotypes.

## File Structure

The files in this project are organized in the following structure:
 - data/
     - raw/
     - processed/
 - src/
 - results/
     - stats/
     - figures/
     - tables/
 - docs/
 - scratch/


## Genome Assembly 

21) Finalize assembly - run Quast and BUSCO 

22) Mask repeats

23) Annotation (Braker)
Run braker in two ways: "ab initio" using just the genome and supplying braker with a bam file from a cDNA Nanopore data. 
Download the gtf files, and move them to the resulst (rename files since braker output is just "braker.gtf") 

## Genomics of Drilling 

1) 


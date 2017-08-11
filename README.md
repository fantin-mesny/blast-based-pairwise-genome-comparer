# blast-based-pairwise-genome-comparer
A blast-based programme to compare two similar genomes, and retrieve functional data (paralog families with network representations, id-to-GO/Pfam/COG associations, reciprocal best matches list...)

## Installation
Download the .py file and make it executable using a "chmod +x .py" command.

## Simple run


## Input files
This programme has been developed to analyse genomes annotated by PROKKA. Therefore, it accepts two types of input:
- either two PROKKA folders (one per genome) including .faa, .ffn, .fna and .gff files.
- or separated predicted transcripts and proteins multi-fasta files, genomes sequences and .gff-type annotation files.


## Data it retrieves
- all-against-all blastp outputs
- all-against-all tblastx outputs
-

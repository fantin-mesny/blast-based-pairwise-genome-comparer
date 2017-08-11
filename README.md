# blast-based-pairwise-genome-comparer
A blast-based programme to compare two similar genomes, and retrieve functional data (paralog families with network representations, id-to-GO/Pfam/COG associations, reciprocal best matches list...)

## Installation
Download the .py file and make it executable using a "chmod +x .py" command.

## Input files
This programme has been developed to analyse genomes annotated by PROKKA. Therefore, it accepts two types of input:
- either two PROKKA folders (one per genome) including .faa, .ffn, .fna and .gff files.
- or separated predicted transcripts and proteins multi-fasta files, genomes sequences and .gff-type annotation files.

## Simple run
    $ app.py -folders "data/PROKKA1,data/PROKKA2"

## Required arguments
- **-folders** Two comma-separated folders including transcripts (.ffn), proteins (.faa) and genome (.fna) fasta files
OR
- **-genomes** Two comma-separated genomes (fasta-files)
- **-prot** Two comma-separated predicted proteins multi-fasta files
- **-rna** Two comma-separated predicted transcripts multi-fasta files

## Optional arguments
- **-o** Output folder, default: creates a folder in the programme's directory 
- **-minLrap** minLrap significance threshold for blast hits, default=0.7
- **-pident** Identity percentage significance threshold for blast hits, default=70
- **-graph** Ouput type for the plot showing genomes similarity. 'img' (default) or 'show' to display using Python, default='img'
- **-cog** Perform COG rpsblasts ? 'yes'(default) or 'no'
- **-pfam** Perform Pfam rpsblasts and get protein-GO associations ? 'yes'(default) or 'no'
- **-paralog_pid** Identity percentage significance threshold to consider proteins as homologous, default=60


## Data it retrieves
- all-against-all blastp original outputs + "hitsTables" giving supplementary informations (reciprocity, minLrap, maxLrap)
- all-against-all tblastx outputs + "hitsTables" giving supplementary informations (reciprocity, minLrap, maxLrap)
- "CDStables", showing for each gene the result of the blastp and tblastx (cds start and end, %id, minLrap and reciprocity)
- a genomes similarity graph (.png) showing the %id for each cds along the genome
- a reciprocal best-matches list, based on chosen %id and minLrap thresholds
- a "paralog" folder containing a list of numbered paralog families and a network weighted graph for each one (.png) 
- rpsblast on Pfam database output
- rpsblast on COG database output
- a list of Pfam numbers per protein id for each organism
- a list of GO numbers per protein id for each organism
- a list of COG numbers per protein id for each organism



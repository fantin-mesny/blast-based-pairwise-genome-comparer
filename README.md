# blast-based-pairwise-genome-comparer
A blast-based programme to compare two similar genomes, and retrieve functional data (paralog families with network representations, id-to-GO/Pfam/COG associations, reciprocal best matches list...)
It has been developed to study two ecotypes of a single bacterial species. Therefore, it is adapted to the comparison of two highly similar genomes.

## Installation
Download the .py file and make it executable using a "chmod +x .py" command.
Requires the following softwares to be installed:
   - python (http://python.org)
   - blast (https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/)
   - argparse (https://pypi.python.org/pypi/argparse)
   - matplotlib (http://matplotlib.org)
   - networkx (http://networkx.github.io)

## Input files
This programme has been developed to analyse genomes annotated by PROKKA. Therefore, it accepts two types of input:
- either two PROKKA folders (one per genome) including .faa, .ffn, .fna and .gff files.
- or separated predicted transcripts and proteins multi-fasta files, genomes sequences and .gff-type annotation files.

## How to run
### Required arguments
- **-folders** Two comma-separated folders including transcripts (.ffn), proteins (.faa) and genome (.fna) fasta files    
    
OR

- **-genomes** Two comma-separated genomes (fasta-files)
- **-prot** Two comma-separated predicted proteins multi-fasta files
- **-rna** Two comma-separated predicted transcripts multi-fasta files
- **-gff** Two comma-separated .gff annotation files

NB: please make sure to follow the same order organism1,organism2 in each of these arguments

    $ app.py -folders "data/PROKKA1,data/PROKKA2
    
    OR
    
    $ app.py -genomes "data/genome1.fna,data/genome2.fna" -prot "data/p1.faa,data/p2.faa" 
        -rna "data/t1.ffn,data/t2.ffn" -gff "data/1.gff,data/2.gff"

### Optional arguments
- **-o** Output folder, default: creates a folder in the programme's directory 
- **-minLrap** minLrap significance threshold for blast hits, default=0.7
- **-pident** Identity percentage significance threshold for blast hits, default=70
- **-graph** Ouput type for the plot showing genomes similarity. 'img' (default) or 'show' to display using Python, default='img'
- **-cog** Perform COG rpsblasts ? 'yes' or 'no' (default)
- **-pfam** Perform Pfam rpsblasts and get protein-GO associations ? 'yes'(default) or 'no' (default)
- **-paralog_pid** Identity percentage significance threshold to consider proteins as homologous, default=60


## Data it retrieves
- all-against-all blastp original outputs + "hitsTables" giving supplementary informations (reciprocity, minLrap, maxLrap)
- all-against-all tblastx outputs + "hitsTables" giving supplementary informations (reciprocity, minLrap, maxLrap)
- "CDStables", showing for each gene the result of the blastp and tblastx (cds start and end, %id, minLrap and reciprocity)
- a genomes similarity graph (.png) showing the %id for each cds along the genome
- a reciprocal best-matches list, based on chosen %id and minLrap thresholds
- a "paralog" folder containing a list of numbered paralog families and a network weighted graph for each one (.png) 

If you wish to get functional data using COG, Pfam and GO annotation systems, uncomment the last lines of the script. You would then obtain:

- rpsblast on Pfam database output
- rpsblast on COG database output
- a list of Pfam numbers per protein id for each organism
- a list of GO numbers per protein id for each organism
- a list of COG numbers per protein id for each organism



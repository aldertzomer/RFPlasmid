# RFPlasmid
Predicting plasmid contigs from assemblies using single copy marker genes, plasmid genes, kmers


Introduction: 
Antimicrobial resistant (AMR) genes in bacteria are often carried on plasmids. Since
these plasmids can spread the AMR genes between bacteria, it is important to know if the genes
are located on highly transferable plasmids or in the more stable chromosomes. Whole genome
sequence (WGS) analysis makes it easy to determine if a strain contains a resistance gene,
however, it is not easy to determine if the gene is located on the chromosome or on a plasmid as
genome sequence assembly generally results in 50-300 DNA fragments (contigs). With our newly
developed prediction tool, we analyze the composition of these contigs to predict their likely
source, plasmid or chromosomal. This information can be used to determine if a resistant gene is
chromosomally located or on a plasmid. The tool is optimized for 19 different bacterial species,
including Campylobacter, E. coli, and Salmonella, and can also be used for metagenomic
assemblies.

Methods: The tool identifies the number of chromosomal marker genes, plasmid replication genes
and plasmid typing genes using CheckM and DIAMOND Blast, and determines pentamer
frequencies and contig sizes per contig. A prediction model was trained using Random Forest on
an extensive set of plasmids and chromosomes from 19 different bacterial species and validated
on separate test sets of known chromosomal and plasmid contigs of the different bacteria.
Results: Prediction of plasmid contigs was nearly perfect when calculated based on number of
correctly predicted bases, with up to 99% specificity and 99% sensitivity. Prediction of small
contigs remains difficult, since these contigs consists primarily of repeated sequences present in
both plasmid and chromosome, e.g. transposases.

Conclusion: The newly developed tool is able to determine if contigs are chromosomal or plasmid
with a very high specificity and sensitivity (up to 99%) and can be very useful to analyze WGS
data of bacterial genomes and their antimicrobial resistance genes.

Requirements

CheckM

Python 3 with pandas ( https://pandas.pydata.org/)

RandomForest package in R

DIAMOND

Optional: Jelly

Usage: 

python3 plaspred_total_nov2018v2.py [-h] --species SPECIES --input INPUT
                                   [--training] [--specieslist] [--jelly]
                                   [--out OUT] [--debug] [--threads THREADS]

--jelly required a functional jelly install. Greatly speeds up the analysis

read specieslist.txt for species specific models.

Plasmid databases can be downloaded from: http://klif.uu.nl/download/plasmid_db/


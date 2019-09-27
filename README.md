# RFPlasmid
Predicting plasmid contigs from assemblies using single copy marker genes, plasmid genes, kmers

Linda van der Graaf-van Bloois, Jaap Wagenaar, Aldert Zomer

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

A rudimentary web-interface to test single fasta files is available here: http://klif.uu.nl/rfplasmid/


Requirements

Python 3 with pandas ( https://pandas.pydata.org/)
```
$ sudo pip3 install pandas
```
CheckM ( https://ecogenomics.github.io/CheckM/ )
Unfortunately CheckM is still Python 2 while RFplasmid requires Python 3.
```
$ sudo pip2 install numpy
$ sudo pip2 install checkm-genome
```

RandomForest package in R ( https://cran.r-project.org/web/packages/randomForest/index.html )
```
$ sudo R
> install.packages("randomForest")
```

DIAMOND ( https://github.com/bbuchfink/diamond )
```
$ wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
$ tar xzf diamond-linux64.tar.gz
$ sudo cp diamond /usr/local/bin/diamond
```

Recommended: Jellyfish ( http://www.genome.umd.edu/jellyfish.html )
```
$ wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-linux
$ sudo cp jellyfish-linux /usr/local/bin/jellyfish
$ sudo chmod +x /usr/local/bin/jellyfish
```

Getting the software:
```
$ git clone https://github.com/aldertzomer/RFPlasmid.git
$ cd RFPlasmid
$ bash getdb.sh # downloads and formats the plasmid DBs
```

Usage: 
```
$ python3 rfplasmid.py [-h] --species SPECIES --input INPUT
                                   [--training] [--specieslist] [--jelly]
                                   [--out OUT] [--debug] [--threads THREADS]
# example
python3 rfplasmid.py --species Campylobacter --input example --jelly --threads 8 --out output
# compare output with the folder example_out
```

--jelly requires a functional jellyfish install. Greatly speeds up the analysis. Strongly recommended as our kmer profiling method in Python is slow

Read specieslist.txt for species specific models. We have a general Enterobacteriaceae model instead of a species model. All others are species except for the "Bacteria" model which can be used for unknown or metagenomics samples.

Plasmid databases can be downloaded from: http://klif.uu.nl/download/plasmid_db/

Data used for training can be downloaded here: http://klif.uu.nl/download/plasmid_db/trainingsets2/

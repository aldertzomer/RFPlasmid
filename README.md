[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# RFPlasmid
Predicting plasmid contigs from assemblies

## Webinterface

A web-interface to test single fasta files is available here: http://klif.uu.nl/rfplasmid/

## Abstract

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

## Getting the software:

Using Conda (thanks to https://github.com/rpetit3 )
```
conda create -n rfplasmid -c conda-forge -c bioconda rfplasmid
```

Using Pip. Installs most requirements except DIAMOND and JellyFish and R (see below). You need to download additional databases for CheckM if you have never installed it. 

```
$  pip3 install rfplasmid
$  export PATH=$PATH:~/.local/bin # pip installs in ~/.local/bin and it should be in your path but some distros don't have this set (even though they should).
$  rfplasmid --initialize #We makes use of a bash helper script to locate the rfplasmid.py file and to download the plasmid files
$  rfplasmid
```

### Required if you have never installed CheckM before: 

CheckM relies on a number of precalculated data files which can be downloaded from https://data.ace.uq.edu.au/public/CheckM_databases/. Decompress the file to an appropriate folder and run the following to inform CheckM of where the files have been placed. The example below uses wget to download the an archive of file and installs them in your homedir.
```
$  cd ~
$  mkdir checkm_data
$  cd checkm_data
$  wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
$  tar xzvf checkm_data_2015_01_16.tar.gz
$  checkm data setRoot ~/checkm_data
```

### Dependencies you need to install 
RandomForest package in R ( https://cran.r-project.org/web/packages/randomForest/index.html ) (likely already installed). 
```
$  R
> install.packages("randomForest")
```

DIAMOND ( https://github.com/bbuchfink/diamond ) 
```
$ wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
$ tar xzf diamond-linux64.tar.gz
$ cp diamond ~/bin/diamond
```

Strongly recommended: Jellyfish ( http://www.genome.umd.edu/jellyfish.html )
```
$ wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-linux
$ cp jellyfish-linux ~/bin/jellyfish
$ chmod +x ~/bin/jellyfish
```

### This is for advanced users. 

You can get the source and using git and run from the folder you downloaded it to. You will need to install the requirements by hand as well
```
$ git clone https://github.com/aldertzomer/RFPlasmid.git
$ cd RFPlasmid
$ bash getdb.sh # downloads and formats the plasmid DBs
$ python3 rfplasmid.py
```

Installing requirements. Assumes you have ~/bin/ in your PATH. Depending on your setup you may need to follow the systemwide version (see further below)

Python 3 with pandas ( https://pandas.pydata.org/)
```
$  pip3 install pandas
```
CheckM ( https://ecogenomics.github.io/CheckM/ ). According to the github page of CheckM:
```
$  pip3 install numpy
$  pip3 install scipy
$  pip3 install pysam
$  pip3 install checkm-genome
$  cd ~
$  mkdir checkm_data
$  cd checkm_data
$  wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
$  tar xzvf checkm_data_2015_01_16.tar.gz
$  checkm data setRoot ~/checkm_data
```

RandomForest package in R ( https://cran.r-project.org/web/packages/randomForest/index.html )
```
$  R
> install.packages("randomForest")
```

DIAMOND ( https://github.com/bbuchfink/diamond )
```
$ wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
$ tar xzf diamond-linux64.tar.gz
$ cp diamond ~/bin/diamond
```

Strongly recommended: Jellyfish ( http://www.genome.umd.edu/jellyfish.html )
```
$ wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-linux
$ cp jellyfish-linux ~/bin/jellyfish
$ chmod +x ~/bin/jellyfish
```

## Usage: 
```
$ python3 rfplasmid.py [-h] --species SPECIES --input INPUT
                                   [--training] [--specieslist] [--jelly]
                                   [--out OUT] [--debug] [--threads THREADS]
# example
python3 rfplasmid.py --species Campylobacter --input example --jelly --threads 8 --out output
# compare output with the folder example_out
```

A folder containing .fasta file is required as input.

--jelly requires a functional jellyfish install. Greatly speeds up the analysis. Strongly recommended as our kmer profiling method in Python is slow

Read specieslist.txt for species specific models. We have a general Enterobacteriaceae model instead of a species model. All others are species except for the "Bacteria" model which can be used for unknown or metagenomics samples.

## Systemwide install
Only if you are a system administrator and you know what you are doing. 
```
$ sudo pip3 install rfplasmid
$ sudo rfplasmid --initialize # Should install the databases as well provided all requirements are met
$ rfplasmid
```

Requirements if you want to install the requirements systemwide. Let your system administrator do this. 

Python 3 with pandas ( https://pandas.pydata.org/)
```
$ sudo pip3 install pandas
```
CheckM ( https://ecogenomics.github.io/CheckM/ ). According to the github page of CheckM:
```
$ sudo pip3 install numpy
$ sudo pip3 install scipy
$ sudo pip3 install pysam
$ sudo pip3 install checkm-genome #install the CheckM databases as well. See the CheckM page
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

Strongly recommended: Jellyfish ( http://www.genome.umd.edu/jellyfish.html )
```
$ wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-linux
$ sudo cp jellyfish-linux /usr/local/bin/jellyfish
$ sudo chmod +x /usr/local/bin/jellyfish
```

## Output Files

| File | Description |
| --------- | ----------- |
| prediction.csv | This is the primary output as comma delimited file. c = chromosome, p=plasmid. |
| prediction_full.csv | This is the primary output as comma delimited file including all input data. |
| classification.RData | The R object containing input and output. |
| outputdataframe.csv | Output from the various scripts used as input by the randomForest classifier. |
| other files | Temporary input files. Can be safely removed. |

## Output explained

The file  prediction.csv contains the contig number (column 1), the prediction wether it's chromosomal or plasmid (column 2), the votes for chromosome or plasmids (columns 3 and 4, and the original contig ID (column 5). 

| Contig | Prediction | Votes chromosomal | Votes plasmid | ContigID |
| --------- | ----------- | ----------- | ----------- | ----------- |
| Kp1_1 | p | 0.100 | 0.900 | Kp1_ctg1 |


The file  prediction_full.csv contains the data as before, but also includes the input data.

| Contig     | prediction | Votes chromosomal | Votes plasmid | contigID | Genome | contig length | % SCM | % plasmid genes | %ID plasmidfinder | Number of kmers | SCM genes | plasmid genes | kmer fractions | etc |
|------------|------------|-------------------|---------------|----------|--------|---------------|-------|-----------------|-------------------|-----------------|-----------|---------------|----------------|-----|
| Kp1_1 | p          | 0.101             | 0.899         | Kp1_ctg1 | Kp1    | 50000         | 0.12  | 0.51            | 69.5              | 138845          | 12        | 32            | 0.01477        |     |
|            |            |                   |               |          |        |               |       |                 |                   |                 |           |               |                |     |

All headers are explained below

| Header              | Header contents   | Explanation                                                                                                                |
|---------------------|-------------------|----------------------------------------------------------------------------------------------------------------------------|
|                     | Contig            | Contig number                                                                                                              |
| prediction          | prediction        | Plasmid or chromosome                                                                                                      |
| votes chromosomal   | Votes chromosomal | Random forest votes for the chromosome class (0-1)                                                                         |
| votes plasmid       | Votes plasmid     | Random forest votes for the plasmid class (0-1)                                                                            |
| contigID            | contigID          | Original contig ID                                                                                                         |
| genome              | Genome            | Name of genome                                                                                                             |
| contig_length       | contig length     | length of contig in bases                                                                                                  |
| SCM_genes           | % SCM             | Percentage of genes that is classified as chromosomal marker gene based   on CheckM                                        |
| plasmid_genes_genes | % plasmid genes   | Percentage of genes that is classified as chromosomal marker gene based   on alignment against a plasmid proteins database |
| plasmidcge_id       | %ID plasmidfinder | Highest identity percentage of hits against the PlasmidFinder database                                                     |
| kmer_number         | Number of kmers   | Number of kmers in the contig                                                                                              |
| SCM                 | SCM genes         | Number of genes that is classified as chromosomal marker gene based on   CheckM                                            |
| plasmid_genes       | plasmid genes     | Number of genes that is classified as chromosomal marker gene based on   alignment against a plasmid proteins database     |
| kmer1               | kmer fractions    | Fractions of the 1024 4-mers from AAAA to TTTTT                                                                            |
| etc                 |                   |                                                                                                                            |

## Training your own model

RFPlasmid comes with an option to train your own model if you provide it with contigs of chromosomes and plasmids. The following steps will generate the model for you and make it available for your RFplasmid installation. We will use the genus Lactococcus as an example.

The following steps need to be taken.
1. Get and install RFPlasmid from Github
2. Prepare a folder containing your plasmid contigs and your chromosome contigs. Draft genome assemblies are highly recommended, don't use completed genomes, as the genomes you want to analyze also won't be complete. Make sure that the plasmid contig file(s) start with "p" and the file(s) with the chromosomal contigs start with "c". In the example they are called plasmids.fasta and chromosomes.fasta, but as long as the files start with p or c your are good. They may also be separate files per genome/plasmid. You need about 20 plasmids and 20 chromosomes for a decent classification.
3. Add the name of your genus/family/etc to specieslist.txt. Use the name exactly as it is listed in CheckM for the appropriate taxon level. It is not recommended to use a species level model, else you will have to deal with spaces in names. Use only genus level or higher. 
4. Run RFplasmid in training mode. If you have very large numbers of contigs (>5000) it is recommended to adjust the sampsize option in the training.R script (change it to e.g. 1000, at most to 2/3rds of the smallest number of contigs per class).
5. Move the resulting training.rfo R object in the output folder to the RFplasmid folder and name it appropriately 
6. Test your new plasmid prediction model

Below an example of the commands. Replace Lactococcus with the genus of interest and ofcourse change "/mnt/data/files" to the folder containing the input
```
$ git clone https://github.com/aldertzomer/RFPlasmid.git
$ cd RFPlasmid
$ bash getdb.sh
$ mkdir Lactococcus
$ cat /mnt/data/files/p*.fasta > Lactococcus/plasmids.fasta
$ cat /mnt/data/files/c*.fasta > Lactococcus/chromosomes.fasta
$ checkm taxon_list |grep Lactococcus
$ echo "Lactococcus genus" >> specieslist.txt # this is important. Don't skip this step.
$ python3 rfplasmid.py --training --jelly --input Lactococcus --threads 16 --out Lactococcus.out --species Lactococcus
$ cp Lactococcus.out/training.rfo Lactococcus.rfo
```

It is recommended to load the object in R and explore it to check how well it performs using the OOB output and the confusion matrix. Generally you want to have the number of chromosome contigs that are predicted as being plasmid contigs as low as possible. It is also recommended to check the sizes of the contigs that are incorrectly predicted, as contigs <1 kb are more difficult to predict and at some point no improvements are possible. Modifying the ratio of the contigs sampled in the sampsize option in the training.R script may push the values to more chromosomal or plasmid contigs correctly predicted. In all our models the ratio was kept 1:1. 

```
$ R
> library(randomForest)
> load("Lactococcus.rfo")
> rf
```


Test your model using 

```
python3 rfplasmid.py --jelly --input Lactococcus --threads 16 --out Lactococcus.out2 --species Lactococcus
```

## Training data
Plasmid databases can be downloaded from: http://klif.uu.nl/download/plasmid_db/
Data used for training can be downloaded here: http://klif.uu.nl/download/plasmid_db/trainingsets2/

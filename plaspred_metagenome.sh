#!/bin/bash

SCRIPT=$(readlink -f "$0")
scriptlocation=$(dirname "$SCRIPT")

echo "Usage = plasmidpredictor <species> <folder>"
echo "Use metagenome as species for species-agnostic recognition"

# $1 = species
# $2 = folder
# $3 = training, this is for development purposes

if [ -z "$1" ]
  then
    echo "Error: no species supplied. Exiting"
    exit 1
fi


if [ -z "$2" ] ; then
    echo "Error: no folder with .fasta files supplied. Exiting"
    exit 1
fi

# is er een model beschikbaar van deze species
# TODO: de tooldir aangeven bij species list

available=`cat $scriptlocation/specieslist.txt |grep -w $1`

# also check if $3 contains the word training. If yes, the developmental option for running the training script is activated
if [ -z "$available" ] && [ "$3" != "training" ] ; then
    echo "Error: This species has no model yet. Exiting"
    exit 1
fi

tmpdir="outputplasmidpredictor"_`date +%s%N`
mkdir $tmpdir

# stop de species en de taxonomic level in de environment oproepbaar voor de volgende scripts. Environment variables zijn ook inleesbaar door python en R mocht dat nodig zijn
species=`echo $available |cut -f 1 -d " "`
level=`echo $available |cut -f 2 -d " "`
export PLASMIDPREDICTORSPECIES=$species
export PLASMIDPREDICTORLEVEL=$level

echo 'cleanup_contigs' 
ls $2/*.fasta | while read genome ; do basename=`basename $genome` ; cat $genome |awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' $genome |cut -f 1 |cat -n |while read hit ; do echo $basename $hit |cut -f 1,2,3 -d ' ' |sed 's/.fasta//g' ; done ; done > $tmpdir/contigID.out
ls $2/*.fasta | while read genome ; do basename=`basename $genome` ; cat $genome |awk '/^>/{print ">" ++i; next}{print}' > $tmpdir/$basename ; done 
echo 'cleanup done'

cd $tmpdir

echo 'Start making checkm files'
if [ "$1" != "metagenome" ] ; then
checkm taxonomy_wf $PLASMIDPREDICTORLEVEL $PLASMIDPREDICTORSPECIES . checkm_output -x fasta --nt -t 16 > checkm_output.tsv
fi

if [ "$1" = "metagenome" ] ; then
echo "You have selected the species agnostic model"
checkm taxonomy_wf domain Bacteria . checkm_output -x fasta --nt -t 16 > checkm_output.tsv
fi

echo 'checkm ; done'
         
echo 'Start blast plasmid genes'
find * |grep genes.faa |while read faa; do diamond blastp -d $scriptlocation/plasmiddb_total -q $faa -e 1E-30 -k 1 |awk '{if ($3 > 80)print $0}' |cut -f 1 -d '_' |sort |uniq -c |awk '{t=$1; $1=$2; $2=t; print}' > $faa.count_plasmid_hit.out ; done
echo 'Blast and count hits per contig; done'

echo 'Start blast cge genes'
find * |grep fasta$ | while read fasta ; do diamond blastx -d $scriptlocation/plasmiddb_cge -q $fasta -e 1E-30 | while read hit ; do echo $fasta $hit |cut -f 1,2,4 -d ' ' |sed 's/.fasta//g' ; done | awk 'unique[$1FS$2]<$3{unique[$1FS$2]=$3; next}END{for (i in unique) print i,unique[i]}' ; done > blast_plasmid_cge_outputpercontig.out
echo 'Blast plasmid genes ; done'

echo 'Start count genes per contig'
find * |grep genes.faa |while read faa ; do cat $faa |grep '>' |awk -F '\\_' '{print $1}' |sed 's/>//g' |sed '/#/d' |sort |uniq -c |awk '{t=$1; $1=$2; $2=t; print}' > $faa.genespercontig.csv ; done
echo 'Counting genes per contig ; done'

   
echo 'Start count SCMs per contig'
find * |grep hmmer.analyze.txt |grep -v result |while read hmm ; do cat $hmm |cut -f 1 -d ' ' |sort|uniq |cut -f 1 -d _ |sed '/#/d' |sort |uniq -c |awk '{t=$1; $1=$2; $2=t; print}' > $hmm.results.txt ; done
echo 'Counting SCMs per contig ; done'

echo 'start python script'
python3 $scriptlocation/plaspred_total.py $tmpdir $scriptlocation

echo 'Start R'
if [ "$3" = "training" ] ; then
    echo "You have selected training mode. This will attempt to generate a model for your organism, which will be saved as training.rfo. Rename to $PLASMIDPREDICTORSPECIES.rfo , copy it to $scriptlocation and adjust specieslist.txt to include species"
    R --vanilla < $scriptlocation/training.R
fi

if [ -z "$3" ] ; then
    echo "Running $PLASMIDPREDICTORSPECIES prediction model"
    R --vanilla --args $PLASMIDPREDICTORSPECIES $scriptlocation < $scriptlocation/classification.R 
fi

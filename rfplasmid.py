#!/usr/bin/env python3

import os, sys, glob, argparse, shutil, multiprocessing
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from itertools import groupby
from datetime import datetime
from contextlib import suppress

scriptlocation = os.path.dirname(os.path.realpath(__file__))

def cpu_threads(max_threads):
    if multiprocessing.cpu_count() > max_threads:
        return max_threads
    else:
        return multiprocessing.cpu_count()

# Import arguments
parser = argparse.ArgumentParser()
parser.add_argument("--species", help="define species (required)", required=('--input') in sys.argv and len(sys.argv) != 1)
parser.add_argument("--input", help="directory with input fasta files (required)", required=('--species') in sys.argv and len(sys.argv) != 1)
parser.add_argument("--specieslist", help="list of available species models", action="store_true", default=False)
parser.add_argument("--jelly", help="run jellyfish as kmer-count (faster)", action="store_true", default=False)
parser.add_argument("--out", help="specify output directory", default=False)
parser.add_argument("--debug", help="no cleanup of intermediate files", action="store_true", default=False)
parser.add_argument("--training", help="trainings mode Random Forest", action="store_true", default=False)
parser.add_argument("--threads", help="specify number of threads to be used, default is max available threads up to 16 threads", default=cpu_threads(16), type=int)
parser.add_argument("--version", help="print version number", action="store_true")

args = parser.parse_args()
species_import = args.species
input_directory = args.input

# version
if args.version:
    print('RFPlasmid version 0.0.15')
    sys.exit()

species_file = os.path.join(scriptlocation, "specieslist.txt")
df_species = pd.read_csv(species_file, header=None, sep=' ', names=['species', 'level'])
level_import = next(iter(df_species.loc[df_species['species'] == species_import, 'level']), 'no match')
species_list = df_species['species'].tolist()

if args.specieslist:
    print('Available species models: \n{}'.format(df_species.species.to_csv(index=False, header=False)))
    sys.exit()

if len(sys.argv) <= 1:
    print("Error; no arguments. Required to specificy --input and --species")
    parser.print_help(sys.stderr)
    sys.exit()

# check if input is directory
if not os.path.isdir(input_directory):
    print('%r is not a directory. Please use directory as input.' % input_directory)
    sys.exit()

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit()

# check if species exist in list
if not args.species in species_list:
    print('Error: species classification model is not available')
    print('Available species models: \n{}'.format(df_species.species.to_csv(index=False, header=False)))
    sys.exit()

print('Start RFPlasmid, version 0.0.15')

# make output folder
if args.out:
    print(args.out)
    new_dir = args.out
if not (args.out):
    now = datetime.now().strftime('%Y%m%d_%H%M%S')
    new_dir = 'output_plasmidpredictor' + '_' + now

os.mkdir(new_dir)

# copy original contigIDs
print('copy original contig names')
for fasta_file in glob.glob(os.path.join(input_directory, '*.fasta')):
    basename = os.path.splitext(os.path.basename(fasta_file))[0]
    with open(fasta_file, 'r') as f_fasta:
        i = 0
        records = SeqIO.parse(f_fasta, 'fasta')
        for record in records:
            i = i + 1
            contig_id = basename + '\t' + str(i) + '\t' + record.description
            with open(new_dir+'/contigID.out', 'a') as ID_file:
                ID_file.write(contig_id + '\n')

# cleanup fasta files
print('cleanup contig names')	
for fasta_file in glob.glob(os.path.join(input_directory, '*.fasta')):
    basename = os.path.splitext(os.path.basename(fasta_file))[0]
    new_file = os.path.join(new_dir, basename+'.fasta')
    i = 0
    with open(fasta_file, 'r') as file_fasta, open(new_file, 'w') as new:
        for i, seq_record in enumerate(SeqIO.parse(file_fasta, "fasta")):
            seq_record.id = str(i + 1)
            seq_record.description = ""
            SeqIO.write(seq_record, new, "fasta")

# change working directory to new directory
os.chdir(new_dir)

# Checkm
print('start Checkm')
os.system('checkm taxonomy_wf {} {} . checkm_output -x fasta --nt -t 16 > checkm_output.tsv'.format(level_import, species_import))

# check if Checkm files are made
if not os.path.isfile("./checkm_output/{file}.ms".format(file=species_import)):
    print("Error: Checkm is not working properly")
    sys.exit()

# check if contigs contain ORFs
for fasta_file in glob.glob('*.fasta'):
    basename = os.path.splitext(os.path.basename(fasta_file))[0]
    if not os.path.isfile("./checkm_output/bins/{contig}/hmmer.analyze.txt".format(contig=basename)):
        print('Error: File {} is too short for prediction'.format(basename))
        sys.exit()
print('Checkm done')

# Blast plasmiddb genes
print('start blast plasmiddb')
plasmiddb_total_file = os.path.join(scriptlocation, "plasmiddb_total")

for root, dirs, files in os.walk(os.getcwd()):
    filepath = os.path.join(root, 'genes.faa')
    basename = os.path.basename(os.path.dirname(filepath))
    if os.path.isfile(filepath):
        output_file = os.path.join(root, 'blast_plasmiddb_output.txt')
        os.system('diamond blastp -d {} -q {} -e 1E-30 -k 1 > {}'.format(plasmiddb_total_file, filepath, output_file))
    else:
        continue
	
for root, dirs, files in os.walk(os.getcwd()):
    filepath = os.path.join(root, 'blast_plasmiddb_output.txt')
    if os.path.isfile(filepath):
        with open(filepath, 'r') as f1:
            with suppress(ValueError):
                df1 = pd.read_csv(f1, header=None, delim_whitespace=True, usecols=[0,2], names=['contig_gene', 'hit'])
                df1['hit'] = df1['hit'].round(0).astype(int)
                df1 = df1.loc[df1['hit'] > 80]
                df1['contig_gene'] = df1['contig_gene'].drop_duplicates()
                df1['contig'], df1['gene'] = df1['contig_gene'].str.split('_', 1).str
                df2 = df1.groupby(['contig']).size().reset_index(name='counts')
                output_file = os.path.join(root, 'plasmiddbgenes_per_contig.txt')
                with open(output_file, 'w') as out_file:
                    df2.to_csv(out_file, sep=' ', index=False, header=False)
    else:
        continue
print('blast plasmid genes done')

# Start blast cge genes
print('start blast cge')
cgedb_file = os.path.join(scriptlocation, "plasmiddb_cge")

for fasta_file in glob.glob('*.fasta'):
    basename = os.path.splitext(os.path.basename(fasta_file))[0]
    output_file = basename + '_blast_plasmidcge_output.txt'
    os.system('diamond blastx -d {} -q {} -e 1E-30 -k 1 > {}'.format(cgedb_file, fasta_file, output_file))

for cge_file in glob.glob('*_blast_plasmidcge_output.txt'):
    with open(cge_file, 'r') as f1:
        with suppress(ValueError):
            genome_name = str(cge_file).split('_blast')[:-1][0]
            df1 = pd.read_csv(f1, header=None, delim_whitespace=True, usecols=[0,2], names=['contig', 'hit'])
            df1['hit'] = df1['hit'].astype(float)
            df1 = df1.loc[df1.groupby('contig')['hit'].idxmax()]
            df1['genome'] = genome_name
            df1.to_csv('blast_plasmid_cge_outputpercontig.out', index=False, header=False, sep=' ', mode='a')

print('blast cge genes done')

# Start count genes per contig
for root, dirs, files in os.walk(os.getcwd()):
    filepath = os.path.join(root, 'genes.faa')
    basename = os.path.basename(os.path.dirname(filepath))
    if os.path.isfile(filepath):
        with open(filepath, 'r') as faa_file:
            id = []
            for gene in faa_file:
                if gene.startswith(">"):
                    gene = gene.split(' ')[0]
                    id.append(gene[1:])
            df = pd.DataFrame({'id': id})
            df[['contig', 'gene']] = df['id'].str.split('_', expand=True)
            df2 = df['contig'].value_counts()
            output_file = os.path.join(root, 'genespercontig.txt')
            with open(output_file, 'w') as out_file:
                df2.to_csv(out_file, sep=' ', index=True, header=False)
print('count genes per contig done')

# Start count SCMs per contig
for root, dirs, files in os.walk(os.getcwd()):
    filepath = os.path.join(root, 'hmmer.analyze.txt')
    basename = os.path.basename(os.path.dirname(filepath))
    if os.path.isfile(filepath):
        with open(filepath, 'r') as scm_file:
            with suppress(ValueError):
                df = pd.read_csv(scm_file, header=None, delim_whitespace=True, usecols=[0], names=['contig_gene'])
                df = df[~df['contig_gene'].str.contains("#")]
                df['contig_gene'] = df['contig_gene'].drop_duplicates()
                df['contig'], df['gene'] = df['contig_gene'].str.split('_', 1).str
                df2 = df.groupby(['contig']).size().reset_index(name='counts')
                output_file = os.path.join(root, 'SCMpercontig.txt')
                with open(output_file, 'w') as out_file:
                    df2.to_csv(out_file, sep=' ', index=False, header=False)
print('count SCMs per contig done')

print('Start kmer count')
if args.jelly:
    print('You have selected kmer-counting with jellyfish')
    kmer_dir = 'fasta_split'
    os.mkdir(kmer_dir)
    for fasta_file in glob.glob('*.fasta'):
        basename = os.path.splitext(os.path.basename(fasta_file))[0]
        for record in SeqIO.parse(fasta_file,"fasta"):
            id = record.id
            id_file = basename + '_' + id + '.fasta'
            with open(os.path.join(kmer_dir, id_file), 'w') as out:
                out.write('>'+record.id+'\n'+str(record.seq + 'NN' + record.seq.reverse_complement() + '\n'))
    os.chdir(kmer_dir)
    for split_fasta in glob.glob('*.fasta'):
        basename = os.path.splitext(os.path.basename(split_fasta))[0]
        os.system('cat {} | jellyfish count /dev/fd/0 -m 5 -c 20 -s 1000000 -t 16 -o /dev/fd/1 |jellyfish dump -c -t /dev/fd/0 |while read kmer ; do echo {} $kmer ; done >> ../kmers.out'.format(split_fasta, basename))
    os.chdir("../")
    print('Jellyfish kmer-count done')
if not (args.jelly):
    print('Standard kmer-counting method')
    kmer_dir = 'fasta_split'
    os.mkdir(kmer_dir)

    def count_kmers(read, k):
        counts = defaultdict(list)
        num_kmers = len(read) - k + 1
        for i in range(num_kmers):
            kmer = read[i:i+k]
            if kmer not in counts:
                counts[kmer] = 0
            counts[kmer] += 1
        return(counts)

    for fasta_file in glob.glob('*.fasta'):
        basename = os.path.splitext(os.path.basename(fasta_file))[0]
        with open(fasta_file) as f_fasta:
            for k, g in groupby(f_fasta, lambda x: x.startswith('>')):
                if k:
                    sequence = next(g).strip('>\n')
                else:
                    d1 = list(''.join(line.strip() for line in g))
                    d1 = [item.upper() for item in d1]
                    d2 = ''.join(d1)
                    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'N', 'M': 'N', 'R': 'N', 'W': 'N', 'S': 'N', 'Y': 'N', 'K': 'N', 'V': 'N', 'H': 'N', 'D': 'N', 'B': 'N', 'X': 'N'}
                    reverse_complement = ''.join(complement.get(base, base) for base in reversed(d2))
                    d3 = (d2+'NN'+reverse_complement)
                    counting = count_kmers(d3, 5)
                    for item in counting:
                        output = basename, sequence, item, counting[item]
                        with open('kmer.out', 'a') as text_file:
                            text_file.write("%s %s %s %s\n" % output)
    print('Kmer count standard method ready')

print('Count contig length')
for fasta_file in glob.glob('*.fasta'):
    basename = os.path.splitext(os.path.basename(fasta_file))[0]
    with open(fasta_file) as f_fasta:
        for k, g in groupby(f_fasta, lambda x: x.startswith('>')):
            if k:
                sequence = next(g).strip('>\n')
            else:
                d1 = list(''.join(line.strip() for line in g))
                d2 = ''.join(d1)
                data = sum(1 for ch in d2)
                output = basename, sequence, data
                with open('count_contiglength.out', 'a') as out_file:
                    out_file.write("%s %s %s\n" % output)
print('count contig length ; done')

# next step: Work on files in different folders and merge to separate dfs
print('Start merging contig files')

for root, dirs, files in os.walk(os.getcwd()):
    filepath = os.path.join(root, 'genespercontig.txt')
    if os.path.isfile(filepath):
        with open(filepath, 'r') as f1:
            df1 = pd.read_csv(f1, header=None, delim_whitespace=True, names=['contig', 'genes'])
            df1['genome'] = os.path.basename(os.path.dirname(filepath))
    else:
        continue

    filepath = os.path.join(root, 'SCMpercontig.txt')
    if os.path.isfile(filepath):
       	with open(filepath, 'r') as f2:
            df2 = pd.read_csv(f2, header=None, delim_whitespace=True, names=['contig', 'SCM'])
            df2['genome'] = os.path.basename(os.path.dirname(filepath))
    if not os.path.isfile(filepath):
        df2 = df1
        df2 = df2.rename(columns={'genes': 'SCM'})
        df2['SCM'] = np.NaN

    filepath = os.path.join(root, 'plasmiddbgenes_per_contig.txt')
    if os.path.isfile(filepath):
        with open(filepath, 'r') as f3:
            df3 = pd.read_csv(f3, header=None, delim_whitespace=True, names=['contig', 'plasmid_genes'])
            df3['genome'] = os.path.basename(os.path.dirname(filepath))
    if not os.path.isfile(filepath):
        df3 = df1
        df3 = df3.rename(columns={'genes': 'plasmid_genes'})
        df3['SCM'] = np.NaN

    # merge dataframes for each folder
    df_merge1 = pd.merge(df1, df2, left_on=['genome', 'contig'], right_on=['genome', 'contig'], how='left')
    df_merge2 = pd.merge(df_merge1, df3, on=['genome', 'contig'], how='left')
    df_merge2.to_csv(os.path.join(root,'outputgenesdf.csv'))
print('Create df_files ; done')

# next step: concatenate csv files
dfList = []

for root, dirs, files in os.walk(os.getcwd()):
    filepath = os.path.join(root, 'outputgenesdf.csv')
    if os.path.isfile(filepath):
        with open(filepath, 'r') as fnew:
            frame = pd.read_csv(fnew, converters={'genome': lambda x: str(x)})
            dfList.append(frame)

df1a = pd.concat(dfList)
df1a['contig'] = df1a['contig'].astype(object)

#import contigIDnumbers
df1b = pd.read_csv('contigID.out', header=None, sep='\t', names=['genome', 'contig', 'contigID'])
df1b['contig'] = df1b['contig'].astype(object)
df1 = pd.merge(df1a, df1b, how='right', left_on=['genome','contig'], right_on=['genome','contig'])
df1['contig'] = df1['contig'].astype(int)
df1['source'] = df1.genome.astype(str).str.cat(df1.contig.astype(str), sep='_')
df1 = df1.set_index('source')

print('Start to add blast output to dataframe')

# next step; make dataframe

# import kmers
kmer_file = os.path.join(scriptlocation, "kmer.txt")
df2 = pd.read_csv(kmer_file, delim_whitespace=True, header=None)
# use row as columnheaders
df2.columns = df2.iloc[0]
df2 = df2.loc[:-1]

# samenvoegen dfs
df3 = pd.concat([df1, df2], axis=1)
df3['source'] = df3.genome.astype(str).str.cat(df3.contig.astype(str), sep='_')
df3 = df3.set_index('source')

# import kmers
if args.jelly:
    df4 = pd.read_csv('kmers.out', sep=' ', low_memory=False, header=None, names=['source', 'gene', 'bitscore'])
    df4['genome'], df4['contig'] = df4['source'].str.rsplit('_', 1).str
if not (args.jelly):
    df4 = pd.read_csv('kmer.out', sep=' ', low_memory=False, header=None, names=['genome', 'contig', 'gene', 'bitscore'])
    df4['source'] = df4.genome.astype(str).str.cat(df4.contig.astype(str), sep='_')

# combine df3 and df4
df5 = df4.pivot_table(index='source', columns='gene', values='bitscore').reindex(index=df3.index, columns=df3.columns).fillna(df3)
df5['genome'] = df5['genome'].astype(object)
df5['contig'] = df5['contig'].astype(object)

# import blastcge_outputpercontig
df5aa = pd.read_csv('blast_plasmid_cge_outputpercontig.out', header=None, delim_whitespace=True, names=["contig", "plasmidcge_id", "genome"])
df5aa['source'] = df5aa.genome.astype(str).str.cat(df5aa.contig.astype(str), sep='_')
df5aa = df5aa.set_index('source')
df5aa['plasmidcge_id'] = df5aa['plasmidcge_id'].astype(object)
df5aa['contig'] = df5aa['contig'].astype(object)
df5aa['genome'] = df5aa['genome'].astype(object)
df5a = pd.merge(df5, df5aa, how='outer')
df5a['source'] = df5a.genome.astype(str).str.cat(df5a.contig.astype(str), sep='_')
df5a['contig'] = df5a['contig'].astype(object)

# import contig length
df5bb = pd.read_csv('count_contiglength.out', header=None, delim_whitespace=True, names=["genome", "contig", "contig_length"])
df5bb['source'] = df5bb.genome.astype(str).str.cat(df5bb.contig.astype(str), sep='_')
df5bb = df5bb.set_index('source')
df5bb['genome'] = df5bb['genome'].astype(object)
df5bb['contig'] = df5bb['contig'].astype(object)
df5bb['contig_length'] = df5bb['contig_length'].astype(object)
df5b = pd.merge(df5a, df5bb, how='outer')
df5b.fillna(int(0), inplace=True)

# add column with kmer_numbers
df5b['kmer_number'] = df5b.contig_length - 4
df5b['kmer_number'] = df5b.kmer_number * 2

# add column with genes/SCM en round to 2 decimals
df5b['SCM_genes'] = df5b['SCM']/df5b['genes']
df5b['SCM_genes'] = df5b['SCM_genes'].fillna(0)

# add column with genes/plasmid_genes en round to 2 decimals
df5b['plasmid_genes_genes'] = df5b['plasmid_genes']/df5b['genes']
df5b['plasmid_genes_genes'] = df5b['plasmid_genes_genes'].fillna(0)

# ID contig column to first position
cols = list(df5b.columns.values)
cols.pop(cols.index('source'))
cols.pop(cols.index('contigID'))
cols.pop(cols.index('contig_length'))
cols.pop(cols.index('kmer_number'))
cols.pop(cols.index('plasmid_genes_genes'))
cols.pop(cols.index('SCM_genes'))
cols.pop(cols.index('contig'))
cols.pop(cols.index('genome'))
cols.pop(cols.index('plasmidcge_id'))
df5b = df5b[['source', 'contigID', 'genome', 'contig', 'contig_length', 'SCM_genes', 'plasmid_genes_genes', 'plasmidcge_id', 'kmer_number']+cols]
df5b.drop('Unnamed: 0', axis=1, inplace=True)

df5da = df5b.loc[:, 'source':'plasmidcge_id']
df5db = df5b.loc[:, 'kmer_number':'TTTTT'].astype(int)

df6 = pd.concat([df5da, df5db], axis=1)

# kmer fraction
df6[['kmer{}'.format(i) for i in range(1, 1+(2002-979)+1)]] = df6.loc[:, 'AAAAA':'TTTTT'].div(df6['kmer_number'], axis=0)
df6.loc[:, 'kmer1':'kmer1024'] = df6.loc[:, 'kmer1':'kmer1024'].applymap('{:.2e}'.format)

# delete some columns
df6.drop(df6.columns.to_series()["AAAAA":"TTTTT"], axis=1, inplace=True)
df6.drop(['contig','genes'], axis=1, inplace=True)

df6.to_csv('outputdataframe.csv', index=False)

print('CSV file from dataframe ready in file outputdataframe.csv')

if args.debug:
    print('debug mode, no cleaning up')
if not (args.debug):
    print('Cleaning-up intermediate files')
    shutil.rmtree(kmer_dir)
    shutil.rmtree('checkm_output')
    for file in glob.glob('*.txt'):
        os.remove(file)
    for file in glob.glob('*.out'):
        os.remove(file)
    for file in glob.glob('*.fasta'):
        os.remove(file)
    os.remove('checkm_output.tsv')
	
print('Start R')

if args.training:
    print('You have selected training mode')
    trainings_location = os.path.join(scriptlocation, "training.R")
    os.system('R --vanilla < {}'.format(trainings_location))
    print('training done')
if not (args.training):
    print('prediction mode')
    classification_location = os.path.join(scriptlocation, "classification.R")
    os.system('R --vanilla --args {} {} < {}'.format(species_import, scriptlocation, classification_location))
    print('Prediction done')

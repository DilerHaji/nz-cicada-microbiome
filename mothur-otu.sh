#!/bin/bash
#SBATCH --job-name=grafmag
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mail-user=diler.haji@uconn.edu
#SBATCH -c 32

cd .

module load mothur


mothur "#summary.seqs(fasta=*.fasta); 
screen.seqs(fasta=current, summary=current, group=*.group maxambig=0, maxlength=275);
unique.seqs(fasta=current);
count.seqs(name=current, group=current);
align.seqs(fasta=current, reference=silva.v4.v132.align);
summary.seqs(fasta=current);
screen.seqs(fasta=current, summary=current, summary=current, count=current, start = 13862, end = 23444, maxhomop=8);
summary.seqs(fasta=current);
filter.seqs(fasta=current, vertical=T, trump=.);
unique.seqs(fasta=current, count=current);
summary.seqs(fasta=current);
pre.cluster(fasta=current, count=current, diffs=1);
summary.seqs(fasta=current);
chimera.vsearch(fasta=current, count=current, dereplicate=t);
remove.seqs(fasta=current, accnos=current);
summary.seqs(fasta=current);
classify.seqs(fasta=current, count=current, reference=silva.v4.v132.align, taxonomy=silva.v4.v132.tax, cutoff=80);
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria);
dist.seqs(fasta=current, cutoff = 0.03);
cluster(column=current, count=current);
make.shared(list=current, count=current, label=0.03);
classify.otu(list=current, count=current, taxonomy=current, label=0.03);
get.oturep(column=current, name=current, fasta=current, list=current);"


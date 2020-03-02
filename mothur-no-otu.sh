#!/bin/bash
#SBATCH --job-name=$NAME
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mail-user=diler.haji@uconn.edu
#SBATCH -c 32

cd .

module load mothur

mothur "#summary.seqs(fasta=$NAME.fasta); 
screen.seqs(fasta=current, maxambig=0, maxlength=275);
unique.seqs(fasta=current);
count.seqs(name=current);
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
system(mkdir tempp);
system(mkdir tempp2);
system(mv *pick.pick.fasta tempp/);
system(mv *pick.count_table tempp/);
system(mv *pick.taxonomy tempp/);
system(mv graf* tempp2/);
system(mv tempp/* .);
system(mv tempp2/ files/);
system(rm -r tempp/);"



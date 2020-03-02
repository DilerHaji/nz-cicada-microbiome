#!/bin/bash
#SBATCH --job-name=grafmag
#SBATCH --mail-type=END
#SBATCH --mail-user=diler.haji@uconn.edu
#SBATCH -c 48
#SBATCH --mem=128G

cd .

module load qiime/1.9.1

mkdir paired_reads/
multiple_join_paired_ends.py –i . –o paired_reads/
mkdir nonjoin/
find paired_reads/ -name “fastqjoin.un*” –print –exec mv {}   nonjoin/ \;
multiple_split_libraries_fastq.py -i paired_reads/ -o paired_reads/ --demultiplexing_method sampleid_by_file --include_input_dir_path -p parameters.txt
sed 's/_S.*_L001.*join.fastq//g' paired_reads/seqs.fna >  paired_reads/paired_seqs.fna
rm paired_reads/seqs.fna
cp paired_reads/paired_seqs.fna paired_reads/seqs.fna
pick_open_reference_otus.py -i fastq/paired_reads/seqs.fna -o fastq/otus/ -m uclust -r silva.v4.v132.align -a 

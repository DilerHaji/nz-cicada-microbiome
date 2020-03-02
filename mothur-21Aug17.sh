Mothurstability21Aug17 

	#Have 2 terminal windows open. One for mothur, one for bash


pwd #/Users/dilerhaji
# Going into stability-31942914 folder to change sample names so they make more sense
../ #up a level
../../ #up two levels

# Changing folder names (Kendra's code)
# i%%-* means go only up to the dash
mkdir fastq
for i in */*.gz; do mv $i "fastq""/"${i%%-*}"."`basename $i`; done 

# -lh will give you details about the files 
ls -lh fastq/

# Printing to screen and you can use arrow keys to move up and down
ls -lh fastq | less




##################### Now getting into mothur #####################

mothur

# Need to use system to use things like cd or pwd
system(cd fastq/)

# Even though mothur makes this file for you, still check it so that the names are correct 
# MOTHUR WILL NOT LET YOU HAVE DASHES IN THE NAMES 
make.file(inputdir=fastq, type = gz)

# Renamed and moved stability.files to working directory (stability-#####)
# mv fastq/stability.files stability.files 

# Use forward and reverse reads to make concatenated sequence. if the bases agree, keep it. If disagree, look at quality score. If quality score difference greater than 6, keep the base with hgher score. If less than 6, given an "n"
# Phred score - logorithmic scale, what we talk about is the exponent. Score of 20 = 1 base in 100 is an error. 
# WARNING: Went back and changed file names (both actual files in fastq folder and file names in stability.files) so there aren't any dashes. replaced with period 
# WARNING: Had to take out the 3 deep water samples from stability.files. Now rerunning
system(cd fastq/)
make.contigs(file=stability.files)

# Could be a mismatch in unzip, maybe unzipping all the files would work
# trying to unzip then process 
# unzip fastq.gz to fastq
# Replace out gz in stability.files 
# Redownloaded files. replaced dashes in names with periods in stability.file but not in the file names of the sequences associated with the sample 
# gzip -d fastq/*001.fastq.gz
# make.contigs(file=stability.files)

# rename 's/what you're looking for/ what you're replacing it with/'
# s means start at the beginning
# Not actually using this code because file names already stability, but just in case 
# rename 's/stability/stability' fastq/*

# man will give you the about page 
# man rename

# stability.trim.contigs.fasta 
# you don't have to run summary.seqs, but Kendra says it is useful
# gives you a sumary of your sequences 
# homopolymers are number of repeated single bases - 454 had a problem with that 
# Should be that most sequences fall into 250 bp
# 8 homopolymers is the most that someone has identified in this amplicon 
summary.seqs(fasta = stability.trim.contigs.fasta)

# Now screening sequences 
# stability.contigs.groups is an important file = 2 column file where 1 column is sample name and other is sequence name, keeps track of which sequences were in which folder
# DO NOT REMOVE FASTA FILE NAMES AND REPLACE WITH SAMPLE NAMES- the groups files keeps track of which sequences match to which sample 
# We are throwing out bad sequences 
# We add the summary file to make the process faster 
# screen.seqs(maxambig, maxlength)
# Remove ambiguous bases, too short/long seqs
# We remove any seq with any N, length trim could have phylogenetic bias
# In the code, we are allowing 0 ambiguous bases and the max length of a sequence is 275 (suggested by Pat)
# 275 because there could be insertion in stems within 16S that are real 
# There is a minlength option (doesn't really apply to Illumina data) 
# OUTPUT: good file of fasta and good file of groups 
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, summary=stability.trim.contigs.summary, maxambig=0, maxlength=275)

# threw out ~15% of sequences, but that's fine
summary.seqs(fasta = stability.trim.contigs.good.fasta)

# We are going to align sequences to database, but aligment is computationally expensive 
# Reducing number of sequences
# Compares each sequence to others and picks 1 sequence of sequences that are EXACTLY the same  
unique.seqs(fasta=stability.trim.contigs.good.fasta)

# Gives you a warning that you need to use names file 
summary.seqs(fasta = stability.trim.contigs.good.unique.fasta)

# Adding names 
# names file is which sequences are exactly the same as others 
# Dropped computational load by a 1/3 
# REMINDER: this is not rarefaction, just taking duplicates out 
summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, name=stability.trim.contigs.good.names)

# Making a table of sequence counts 
# reducing computational time by replacing names and groups files which have all the sequence names
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table)

# Aligning sequences to SILVA alignment
# Need to download the SILVA database (we are using the SEED database in this case, but there is a full database too)
# Mac can unzip tgz fine, but if you need to do it in terminal, do following: tar -xvzf Silva.seed_v128.tgz

# Trims out everything before were we tell it to start and after where we tell it to end
# These numbers are for the v4 primers 
# Keepdots means - we don't want to keep anything that started in the middle of v4
pcr.seqs(fasta=silva.seed_v128.align, 
, keepdots = F)

# Don't want to keep the name that mothur gives
# Adding that this is alignment for v4 to make things easier in the future 
mv silva.seed_v128.pcr.align silva.v128.v4.align

# summary.seqs(fasta=silva.seed_v128.align)
summary.seqs(fasta=silva.v128.v4.align)

# For MiSeq, flip paramater is not set to true because it's known that the directionality of the sequences is correct. This is directional library prep
# If the barcodes were ligated on than we need to flip = true 
align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.v128.v4.align, processors = 16)
summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, processors = 2)

# Aligning and screening gets rid of spurious sequences 
# Anything that starts before 1968 and after 11550 will be thrown away 
screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, start=1968, end=11550, maxhomop=8)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table)

# Getting rid of all columns that are all gaps (bad alignment, eukaryotes, etc.)
# You don't really need vertical=T because that's the default 
# use trump in filter.seqs to remove colums that contain at least one NA (more for phylogenetic inference)
filter.seqs(fasta=stability.trim.contigs.good.unique.good.align,vertical=T)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta,count=stability.trim.contigs.good.good.count_table)

# pre clustering 
# decreasing the number of sequencing. Any sequences less than X number of difference 
# We are expanding definition of unique 
# If genus level difference (5% difference)
# If species levele (1% difference)
# Precluster less than half of the divergence your are interested in 
# You can pick the most abundant sequence within the each pre clustered unique
# If you have a smaller dataset, than you don't need to do this
# Diffs means the the number of different bases (2 if 250 bases long for about 1% difference) 
pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.fasta,count=stability.trim.contigs.good.good.count_table,diffs=2)

summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.precluster.fasta,count=stability.trim.contigs.good.unique.good.filter.precluster.count_table)

# Chimera vsearch
# All seqs from sample, makes database, takes every seq and compares to see if it matches different seqs at front and end using a sliding window approach
# THIS DID NOT WORK SO WE ARE MOVING ON WITHOUT REMOVING CHIMERAS> FOR SOME REASON THE VSEARCH IS NOT BEING DETECTED IN THE PATH 
## FIXED THE PROBLEM: Put the uchime, vsearch, and mothur executable files into the working directory that you are currently 
chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.precluster.count_table)

# We will remove the sequences 
# accnos function - accnos file contains the names of the sequences that it flags as chimera
# Dereplicate - chimera checking is only done within each sample. Sometimes sequence is flagged in one in one sample and not others. dereplicate take them both out
# Accnos files are created when you remove seqs
remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.precluster.count_table, accnos=stability.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.accnos)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.precluster.pick.count_table)


###### Kendra does not take out singletons ! You are already taking sample of the community 




#############N NOTE THAT FROM HERE ON OUT, USING THE DATA from pre cluster without filtering for chimeras to see if it's different from everyone else in the class
############# So no "pick" in the file names. In the future, all file names from here on out will have "pick" in the file names 

# We will classify sequences 
# If you want to use archeae just remember that the v4 primers are designed for bacteria and won't amplify all archeae the evenly 
# RDP method for classification - uses bayesian stats.
# Drop cutoff if there is risk of bad classification (60 for cicadas)
classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.precluster.pick.count_table, reference=silva.v128.v4.align, taxonomy=silva.seed_v128.tax, cutoff=60)

# Now we will cluster! 
# Decisions to be made: reference or no reference 
# No reference alignment will be VERY computationally expensive and create files that are terabytes big 
# Heuristic clustering = Kendra is a big fan. 
# Greedy clustering = pairwise clustering with a threshold: choosing type sequence is really important though, which is why Kendra doesn't like this. SO THIS ACTUALLY MAKES THE SEQUENCES WITHIN AN OTU ACTUALLY MORE DIFFERENT THAN THE SIMILARITY YOU PICK. Sometimes people sort by aduncance and the most abundant sequence is the start. The sequence that gets chosen by the type is REALLY IMPORTANT. 
# Pat's solution = classify seqs first than cluster within specified taxa (i.e. proteobacteria). This is kind of a 2-step heuristic. Doesn't work when you have a huge diverse community of proteobacteria. 
# Pat's new solution = starts with all the seqs in 1 pool and pulls out the most different one and that becomes a seed 

get.current()

# Removing based on the taxonomy we choose
# Taxon = everything we want to take out
# REMOVE THE LINEAGES THAT YOU DON'T WANT !!!! 
# remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.precluster.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.precluster.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)	

summary.seqs(fasta=current, count=current)


# Creating distance matrix of all sequences 
dist.seqs(fasta=current, cutoff=0.03)

# Clustering 
# Matthews correlation coefficient is used 
# Average distance can't be more than 3% or nearest neighbor can't be more than 3% or furthest neighbor can't be more than 3% 
# Average neighbor tend to be the one that alot of people use 
# Furthest neighbor make many OTUs
# Nearest neighbor makes big OTUs
# OptiClust is an alternative to the above (dont really know how it works, see OptiClust paper by Pat) - supposed to use heuristic clustering but with the same speed as greedy clustering
# Qimme uses usearch and probably average neighbor 
# DADA2 uses hamming distances instead of MCC (matthews correlation coefficient). 
# Rob Knights group does not merge forward and reverse reads (discard reverse read because more error prone) but Kendra says that this makes you lose a lot of base error information 
# Extremem mock communities - INTERESTING - how well do different clustering methods do when some OTUs dominate dataset (i.e. hodgkinia)

cluster(column=stability.trim.contigs.good.unique.good.filter.precluster.pick.dist, count= stability.trim.contigs.good.unique.good.filter.precluster.pick.count_table)

# Make OTU table 
# This is your shared OTU file
# 
make.shared(list=current)

# Classifying all of the OTU
# Looks at taxonomic identification of every seq inside of the OTU and creates consencus taxonomy in the OTU. 
# If all of the seqs within an OTU are not identified to the same taxa, it will go up a taxa level
classify.otu(list=current, taxonomy=current, threshold=80)


get.oturep(fasta=current, count=current, list=current, method=abundance)
count.groups(shared=current)


# Tell you how many seqs for each sample 
# You can rarefy down - people like 10,000. THIS IS RANDOM SUBSAMPLING IN MOTHUR
# count.groups(shared=stability.trim.contigs.good.unique.good.filter.precluster.opti_mcc.shared)

# Shannon diverisity takes into account abundance? 
# species obsserved 
# Good's coverage : what percentage of the community you have good coverage for. 1 - (number of individuals in species / total number of individuals)
# Plot shannon and shannon's evenness or inverse simpson etc. to see how each differs 
# nseqs
# invsimpson (Kendra's favorite)
summary.single(shared=current, calc=nseqs-sobs-coverage-shannon-shannoneven-invsimpson, subsample=4000)

# Jaccard is a presence absence dissimilarity measure (does not incorporate abundance). Upweights the importance of rares 
# BrayCurtis DOES take into account abundance. Tries to balance the importance of rares
# Thetayc (2011) - like braycurtis but really! punishes abundant and low abundant. REALLY downweights the importance of rares
# Jaccard, braycurtis, thetayc are a continuum of how important rares are 
# 
dist.shared(shared=current, calc=braycurtis-jest-thetayc, subsample=4000)

# Sub sample shared. 
sub.sample(taxonomy=stability.trim.contigs.good.unique.good.filter.precluster.pick.seed_v128.wang.taxonomy, count=stability.trim.contigs.good.unique.good.filter.precluster.pick.count_table, list=stability.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.list, size=4000, persample=true, label=0.03)

summary.tax(taxonomy=current, count=current)

# PARTIAL MATELL 

# Always check the make file before sending off to server 
# Edit the text file on the server, make sure the samples are there, name it
# On server: 




# We chose primers that were supposed to hit bacteria, so we are taking out everything not classified to bacteria before even looking at the results of the analysis. 
# INTERESTING IDEA: sequence all ribosomal transcripts 


###### When you are done, run everything again to get a single log file of everything and delete all of the other log files 






c############################## Working Directory ###############################

setwd("/Users/dilerhaji/void/simon-lab/microbiome/nz-cicada-microbiome/batch2-projectK/analysis/NZprojK")


################################ Libraries ###############################
 
 
#source("https://bioconductor.org/biocLite.R")
#biocLite("qiime2R")
library(qiime2R)
library(tidyr)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggtree)
library(phyloseq)
library(yaml)
library(biomformat)
library(ape)
library(stringr)
library(scales)
library(RColorBrewer)
library(vegan)
library(data)
library(data.table)
library(ggpubr)
library(indicspecies)
library(MASS)
library(pscl)
library(DESeq)
library(microbiome)
library(cluster)
library(ggrepel)
library(lme4)
library(vegan3d)
library(plotly)
library(RColorBrewer)
library(psadd)
library(dendextend)
library(DivNet)
library(tidyverse)
library(reshape2)
library(DESeq2)
library(structSSI)
library(Biostrings)

################################### read.qza ###################################

read_qza <- function(file, tmp, rm) {

if(missing(tmp)){tmp <- tempdir()}
if(missing(file)){stop("Path to artifact (.qza) not provided")}
if(missing(rm)){rm=TRUE} #remove the decompressed object from tmp

unzip(file, exdir=tmp)
unpacked<-unzip(file, exdir=tmp, list=TRUE)

artifact<-read_yaml(paste0(tmp,"/", paste0(gsub("/..+","", unpacked$Name[1]),"/metadata.yaml"))) #start by loading in the metadata not assuming it will be first file listed
artifact$contents<-data.frame(files=unpacked)
artifact$contents$size=sapply(paste0(tmp, "/", artifact$contents$files), file.size)
artifact$version=read.table(paste0(tmp,"/",artifact$uuid, "/VERSION"))

#if(sum(artifact$version$V2==c("2","4","2018.4.0"))!=3){warning("Artifact was not generated with Qiime2 2018.4, if data is not successfully imported, please report here github.com/jbisanz/qiime2R/issues")}#check version and throw warning if new format

 #get data dependent on format
if(grepl("BIOMV", artifact$format)){
 suppressWarnings(artifact$data<-as(biom_data(read_biom(paste0(tmp, "/", artifact$uui,"/data/feature-table.biom"))),"matrix")) #suppressing warning about \n characters
} else if (artifact$format=="NewickDirectoryFormat"){
 artifact$data<-read.tree(paste0(tmp,"/",artifact$uuid,"/data/tree.nwk"))
} else if (artifact$format=="DistanceMatrixDirectoryFormat") {
 artifact$data<-as.dist(read.table(paste0(tmp,"/", artifact$uuid, "/data/distance-matrix.tsv"), header=TRUE, row.names=1))
} else if (grepl("StatsDirFmt", artifact$format)) {
 if(paste0(artifact$uuid, "/data/stats.csv") %in% artifact$contents$files.Name){artifact$data<-read.csv(paste0(tmp,"/", artifact$uuid, "/data/stats.csv"), header=TRUE, row.names=1)}
 if(paste0(artifact$uuid, "/data/stats.tsv") %in% artifact$contents$files.Name){artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/stats.tsv"), header=TRUE, row.names=1, sep='\t')} #can be tsv or csv
} else if (artifact$format=="TSVTaxonomyDirectoryFormat"){
 artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/taxonomy.tsv"), sep='\t', header=TRUE, quote="")
} else if (artifact$format=="OrdinationDirectoryFormat"){

 linesplit<-suppressWarnings(readLines(paste0(tmp,"/", artifact$uuid, "/data/ordination.txt")))
 linesplit<-linesplit[sapply(linesplit, function(x) x!="")]
 
 for (i in 1:length(linesplit)){
   if(grepl("^Eigvals\\t|^Proportion explained\\t|^Species\\t|^Site\\t|^Biplot\\t|^Site constraints\\t", linesplit[i])){
	 curfile=strsplit(linesplit[i],"\t")[[1]][1]
   } else {
	write(linesplit[i], paste0(tmp,"/", artifact$uuid, "/data/",curfile,".tmp"), append=TRUE)
   }
 }

 for (outs in list.files(paste0(tmp,"/", artifact$uuid,"/data"), full.names = TRUE, pattern = "\\.tmp")){
   NewLab<-gsub(" ", "", toTitleCase(gsub("\\.tmp", "", basename(outs))))
   artifact$data[[NewLab]]<-read.table(outs,sep='\t', header=FALSE)
   if(NewLab %in% c("Eigvals","ProportionExplained")){colnames(artifact$data[[NewLab]])<-paste0("PC",1:ncol(artifact$data[[NewLab]]))}
   if(NewLab %in% c("Site","SiteConstraints")){colnames(artifact$data[[NewLab]])<-c("SampleID", paste0("PC",1:(ncol(artifact$data[[NewLab]])-1)))}
   if(NewLab %in% c("Species")){colnames(artifact$data[[NewLab]])<-c("FeatureID", paste0("PC",1:(ncol(artifact$data[[NewLab]])-1)))}
 }
 
 artifact$data$Vectors<-artifact$data$Site #Rename Site to Vectors so this matches up with the syntax used in the tutorials
 artifact$data$Site<-NULL
 
} else if (artifact$format=="DNASequencesDirectoryFormat") {
 artifact$data<-readDNAStringSet(paste0(tmp,"/",artifact$uuid,"/data/dna-sequences.fasta"))
} else if (artifact$format=="AlignedDNASequencesDirectoryFormat") {
 artifact$data<-readDNAMultipleAlignment(paste0(tmp,"/",artifact$uuid,"/data/aligned-dna-sequences.fasta"))
} else if (grepl("EMPPairedEndDirFmt|EMPSingleEndDirFmt|FastqGzFormat|MultiplexedPairedEndBarcodeInSequenceDirFmt|MultiplexedSingleEndBarcodeInSequenceDirFmt|PairedDNASequencesDirectoryFormat|SingleLanePerSamplePairedEndFastqDirFmt|SingleLanePerSampleSingleEndFastqDirFmt", artifact$format)) {
 artifact$data<-data.frame(files=list.files(paste0(tmp,"/", artifact$uuid,"/data")))
 artifact$data$size<-format(sapply(artifact$data$files, function(x){file.size(paste0(tmp,"/",artifact$uuid,"/data/",x))}, simplify = TRUE))
} else if (artifact$format=="AlphaDiversityDirectoryFormat") {
 artifact$data<-read.table(paste0(tmp, "/", artifact$uuid, "/data/alpha-diversity.tsv"))
} else {
 message("Format not supported, only a list of internal files and provenance is being imported.")
 artifact$data<-list.files(paste0(tmp,"/",artifact$uuid, "/data"))
}

pfiles<-paste0(tmp,"/", grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE))
artifact$provenance<-lapply(pfiles, read_yaml)
names(artifact$provenance)<-grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE)
if(rm==TRUE){unlink(paste0(tmp,"/", artifact$uuid), recursive=TRUE)}
return(artifact)
}

qza_to_phyloseq<-function(features,tree,taxonomy,metadata, tmp){
  
  if(missing(features) & missing(tree) & missing(taxonomy) & missing(metadata)){
    stop("At least one required artifact is needed (features/tree/taxonomy/) or the metadata.")
  }
  
  if(missing(tmp)){tmp="/tmp/"}
  
  argstring<-""
  
  if(!missing(features)){
    features<-read_qza(features, tmp=tmp)$data
    argstring<-paste(argstring, "otu_table(features, taxa_are_rows=T),")
  }
  
  if(!missing(taxonomy)){
    taxonomy<-read_qza(taxonomy, tmp=tmp)$data
    taxt<-suppressWarnings(do.call(rbind, strsplit(as.character(taxonomy$Taxon),"\\; ")))
    rownames(taxt)<-taxonomy$Feature.ID
    argstring<-paste(argstring, "tax_table(taxt),")
  }
  
  if(!missing(tree)){
    tree<-read_qza(tree, tmp=tmp)$data
    argstring<-paste(argstring, "phy_tree(tree),")
  }
  
  if(!missing(metadata)){
    metadata<-read.table(metadata, row.names=1, sep='\t', quote="", header=TRUE)
    argstring<-paste(argstring, "sample_data(metadata),")
    sample_data(metadata)
  }
  
  argstring<-gsub(",$","", argstring) #remove trailing ","
  
  physeq<-eval(parse(text=paste0("phyloseq(",argstring,")")))
  
  return(physeq)
}


################################### taxonomy gza to taxonomy phyloseq function (uses read.qza) ###################################

tax_gza_to_phyloseq <- function(path_tax_gza) {
  taxonomy <- read_qza(path_tax_gza)
  taxonomy <- data.frame(taxonomy$data)
  kingdom <- sapply(str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 1), "__"), "[", 2)
  phylum <- sapply(str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 2), "__"), "[", 2)
  class <- sapply(str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 3), "__"), "[", 2)
  order <- sapply(str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 4), "__"), "[", 2)
  family <- sapply(str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 5), "__"), "[", 2)
  genus <- sapply(str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 6), "__"), "[", 2)
  feature_id <- taxonomy$Feature.ID
  taxonomy <- data.frame(kingdom = kingdom, phylum = phylum, class = class, order = order, family = family, genus = genus)
  rownames(taxonomy) <- feature_id
  return(taxonomy) }
 
 
 
############################### making phyloseq ###################################

tree <- read_qza("../rooted-tree.qza")
tree <- tree$data
taxonomy <- as.matrix(tax_gza_to_phyloseq("../taxonomy.qza"))
table <- read_qza("../table.qza")
table <- table$data
metadata <- read.csv("metadata.csv")
metadata <- sample_data(metadata)
metadata$elevation <- as.numeric(as.character(metadata$elevation))
rownames(metadata) <- metadata$X.SampleID
phy <- phyloseq(otu_table(table, taxa_are_rows = TRUE), tax_table(taxonomy), phy_tree(tree), sample_data(metadata))
colnames(tax_table(phy)) <- c("Kingdom", "Phylum","Class","Order","Family","Genus")
sample_data(phy)$depth_original <- as.numeric(sample_sums(phy))



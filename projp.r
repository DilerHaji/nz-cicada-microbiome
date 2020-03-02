setwd("/Users/dilerhaji/Desktop/")
load(".RData")

################################ Libraries ###############################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#install.packages("yaml") 
library(yaml)
#install.packages("Biostrings") 
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("ape") 
library(ape)
#BiocManager::install("phyloseq")
library(phyloseq)
#BiocManager::install("biomformat")
library(biomformat)
#install.packages("ggplot2") 
library(ggplot2)
#install.packages("stringr") 
library(stringr)
#install.packages("vegan") 
library(vegan)
#install.packages("dplyr") 
library(dplyr)
#install.packages("ggrepel") 
library(ggrepel)
#install.packages("devtools") 
library(devtools)
#install.packages("ggpubr") 
library(ggpubr)
#install.packages("tidyr") 
library(tidyr)
#BiocManager::install("decontam") 
library(decontam)


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
  artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/taxonomy.tsv"), sep='\t', header=TRUE)
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
    taxt<-strsplit(as.character(taxonomy$Taxon),"\\;")
    taxt<-lapply(taxt, function(x){length(x)=7;return(x)})
    taxt<-do.call(rbind, taxt)
    taxt<-apply(taxt,2, function(x) replace(x, grepl("^[kpcofgs]__$", x), NA))
    rownames(taxt)<-taxonomy$Feature.ID
    colnames(taxt)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    argstring<-paste(argstring, "tax_table(taxt),")
  }

  if(!missing(tree)){
    tree<-read_qza(tree, tmp=tmp)$data
    argstring<-paste(argstring, "phy_tree(tree),")
  }

  if(!missing(metadata)){
    
    defline<-suppressWarnings(readLines(metadata)[2])
    if(grepl("^#q2:types", defline)){
      metadata<-read_q2metadata(metadata)
      rownames(metadata)<-metadata$SampleID
      metadata$SampleID<-NULL
    } else{
      metadata<-read.table(metadata, row.names=1, sep='\t', quote="", header=TRUE)
    }
    argstring<-paste(argstring, "sample_data(metadata),")
    sample_data(metadata)
  }

  argstring<-gsub(",$","", argstring) #remove trailing ","

  physeq<-eval(parse(text=paste0("phyloseq(",argstring,")")))

  return(physeq)
}

read.nex <- function(x){
	
	x <- scan(x, what = "c", quiet = TRUE)
		
	## eliminate comments
	## ------------------
	left <- grep("\\[", x)
	right <- grep("\\]", x)
	if ( length(left) > 0 ){
	  m <- cbind(left, right)
	  x <- x[-unlist(apply(m, 1, function(x) x[1]:x[2]))]
	}
	
	x <- x[x != ""]
	
  ## getting number of taxa
  ## ----------------------
	ntax <- x[grep("ntax", x, ignore.case = TRUE)]
	ntax <- gsub("[[:alpha:]]|[[:punct:]]", "", ntax )
	nb <- ntax <- as.numeric(unique(ntax))
		
	## getting number of characters	
  ## ----------------------------
	ncha <- x[grep("nchar", x, ignore.case = TRUE)]
	ncha <- gsub("[[:alpha:]]|[[:punct:]]", "", ncha )
	ncha <- as.numeric(unique(ncha))
	
	## get beginning and end of matrix
  ## -------------------------------
	start <- grep("^\t?matrix$", x, ignore.case = TRUE)
	end <- grep(";", x)
	end <- min(end[end > start])
	M <- x[(start + 1):(end - 1)]
	
	# assemble DNAbin object:
	# -----------------------
	nblock <- ceiling(ncha / nchar(M[2]))
	id <- seq(1, 2 * ntax, by = 2)
	nam <- M[id]
	fuse <- function(s, M, nblock, ntax){
	  paste(M[seq(s, length.out = nblock, by = ntax * 2)], collapse = "")
	}
	seq <- lapply(id + 1, fuse, M = M, nblock = nblock, ntax = ntax)
	obj <- list(nb = ntax, seq = seq, nam = nam, com = NA)
	class(obj) <- "alignment"
	as.DNAbin(obj)
}


phylo <- qza_to_phyloseq(features="nzbiome-dada-table.qza", tree="nzbiome-dada-filtered-alignment-rooted-tree.qza", taxonomy="nzbiome-dada-rep-seqs-taxonomy.qza", metadata="projP-final-metadata.tsv")

phylo <- subset_samples(phylo, !is.na(type) 
	& !type %in% c("uncertain", "Unknown9","bacteriome")
	& !ID %in% c("C1", "C2", "T3", "T5", "T7", "T23")
	)

meta <- sample_data(phylo)






######## FUNCTIONS #######

#### functions for mean relative abundance 
relphy <- function(phy){apply(otu_table(phy), 2, function(x){x/sum(x)})}
relphy <- function(phy){apply(otu_table(phy), 2, function(x){x/sum(x)})}

mrb <- function(otu){apply(apply(otu_table(otu), 2, function(x){x/sum(x)}), 1, function(x){mean(x[x > 0])})[is.finite(apply(apply(otu_table(otu), 2, function(x){x/sum(x)}), 1, function(x){mean(x[x > 0])}))]}
mrb2 <- function(otu){mrb(otu)[mrb(otu) >= quantile(mrb(otu))[3]]}
mrb3 <- function(otu){mrb(otu)[mrb(otu) <= quantile(mrb(otu))[5]]}

pa_phy <- function(phy){
	comm <- as(object = phyloseq::otu_table(phy), Class = "matrix")
	comm_std <- vegan::decostand(comm, "pa")
	phyloseq::otu_table(phy) <- phyloseq::otu_table(comm_std, taxa_are_rows = TRUE)
	return(phy)
}

















#### Ophio DNA

mito <- phylo
#p <- transform_sample_counts(mito, function(x)(x/sum(x, na.rm = TRUE)))
#mito <- subset_taxa(p, Family %in% c("D_4__Mitochondria"))
#write.table( taxa_names(mito), "mito_asv.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

mito <- subset_samples(mito, sample_sums(mito) > 0)
mito <- subset_samples(mito, !type %in% "control")
mito <- subset_taxa(mito, taxa_sums(mito) > 0)
mito_ophio <- subset_taxa(mito, taxa_names(mito) %in% c(
"e1b3e502889eda99faec882decc61084",
"b4cc0fe078d755f6dc41dd4bdb4ec7f0",
"b6b811e7b71729058ad58cf4173a8f38",
"bd7c33818f8aec75f7161d4e7185badb",
"7198d35d03014e52e939411f719a695b",
"adbb584c88e2f713dc2e6215ef878a02",
"776cb3ae1cc229b9999dbd9c7b800fc4",
"bd637283f4c956d6eb3f60aba645062d",
'8a31918806e315293752775581aa6227',
"a597872fcf87b49f099e6525ce8e0bdd",
"dd4a893a3874a25e1e1f9cdca607eb8c",
"3ad79c231b558392426ed778987cfd28",
"dd69e8131ec9521228e03aa69b7685cc",
"fd3fa68f678cb377e70ae99a348b458b",
"59e8f6fedbcceff4ca6c982cc46d1562",
"9fa8c30cf03d22a90f75fcc93cb22b1f",
"5f3a1e0d63fa00f3ecd7c4a43759576e",
"5eb79a1bc542e2bf20251039b5a3892b",
"755666aedf0d0425062f31bb42c1fc8a",
"f347ba578d6e5fc132aafadf7b467e41",
"7c15c08a32970d2123f81c79a841dcda",
"a10c98cca957ba6b03dc4129da14e75a",
"8b7826d9e54b2dfc747d851087520271",
"3b4a02054a191da163cdbb3402ab1adb",
"2f19bacb43636b5d04a5af15efe85bcd",
"f87631331b247a657230e771a3e67f99",
"427f6aaa7865aeb2611e969886dd4b4c",
"9db6287327eae687af1b9ca9b89b7894",
"b8b9f752ace6d5c1558c53881c1a3c24",
"cf56e57f332051fbc8076be97342c1d1",
"532042554bb7b83f4f87972af256d7c6",
"f16746b4c811ba0ddbcc063c94858ae4",
"b66f62979ed505d71a9050f44d03eafc",
"fc1168706670fa562f282a2cc64db9d1",
"0986e084d515a7c589486d57f012e290",
"7fd0409aa73324573665cedccca947f7",
"d95af54ccf6d283644fbb2ebf448b3fe",
"790ccb092610629a233eba365be48fa7",
"71efa5bd59720222774a0e29a2963cc6",
"0a2757312143a76b0e38b780ba067875",
"2e1081dbfe0c6c8ffdb56c1167c72bf6",
"326fda0ea300b8bc64e930f5112f15e9",
"40ef5a4fd881601fdaaf9f35a658325f",
"c0cce88070e16e87cb1e8a0410ded5db",
"f49f87ee110aaf95f83f1e476d3a32b7",
"a75e5c8106e8135d50e2ba9327cefd73",
"84419617fd969dbddba9901ff4ed16d9",
"646f5e07f65b733c72cace35a704f30d",
"6d925e800f29f71a5a1ddaf29b641dc1",
"ab72d8d94a01bf5f6a8990c0686a5d04",
"ab2383c7be692b1f25d72487da6df733",
"f440d0c74ede6a04004b11fd859b1ac0",
"d7526f29c9eb1a0c8d1ab3748edab74f",
"5abaea3f41fe71f043900d5106be7165",
"c31ec8d17f0a762dcf5094edae97c306",
"0d8331a0106b98f8682c7c8627831afa"))


#plot(phy_tree(tip_glom(mito_ophio, h = 0.03)))

mito_ophio2 <- tip_glom(mito_ophio, h = 0.035)
#mito_ophio3 <- merge_samples(mito_ophio2, "species")

mito_ophio3 <- mito_ophio2
tax_table(mito_ophio3)[,6] <- seq(1, length(tax_table(mito_ophio3)[,6]), 1)

tax_table(mito_ophio3)[,6] <- rownames(tax_table(mito_ophio3))


#mito_ophio3 <- transform_sample_counts(mito_ophio3, function(x){x/sum(x)})

plot_bar(mito_ophio3, fill = "Genus") +
    #scale_x_discrete(limits = rownames(sample_data(mito_ophio3)[order(sample_data(mito_ophio3)$genus),])) +
    facet_wrap(~species, scales = "free", shrink = TRUE, drop = TRUE) +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    	  axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    	  legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())





otr <- ggtree(phy_tree(mito_ophio2)) + geom_tiplab(size = 1)

ggsave("otr.png", otr, scale = 0.7)

















#### Relative Sulcia, hodgkinia, mitochondria, archea, and other bacteria 

phy <- phylo

a <- taxa_names(subset_taxa(phy, Genus %in% c("D_5__Candidatus Sulcia")))
b <- taxa_names(subset_taxa(phy, Genus %in% c("D_5__Candidatus Hodgkinia")))
c <- taxa_names(subset_taxa(phy, Family %in% c("D_4__Mitochondria") & !taxa_names(phy) %in% taxa_names(mito_ophio)))
d <- taxa_names(subset_taxa(phy, Order %in% c("D_3__Chloroplast")))
e <- taxa_names(subset_taxa(phy, Kingdom %in% c("D_0__Archaea")))
f <- taxa_names(subset_taxa(phy, !taxa_names(phy) %in% c(a,b,c,d,e) & Kingdom %in% "D_0__Bacteria"))
g <- taxa_names(mito_ophio)

pd <- data.frame(tax_table(phy))
pd[,6] <- as.character(pd[,6])

pd[rownames(pd) %in% f, colnames(pd) %in% "Genus"] <- "Bacteria"
pd[rownames(pd) %in% a, colnames(pd) %in% "Genus"] <- "Candidatus Sulcia"
pd[rownames(pd) %in% b, colnames(pd) %in% "Genus"] <- "Candidatus Hodgkinia"
pd[rownames(pd) %in% c, colnames(pd) %in% "Genus"] <- "Mitochondria"
pd[rownames(pd) %in% d, colnames(pd) %in% "Genus"] <- "Chloroplast"
pd[rownames(pd) %in% e, colnames(pd) %in% "Genus"] <- "Archaea"
pd[rownames(pd) %in% g, colnames(pd) %in% "Genus"] <- "Ophiocordyceps"

pd[,6] <- as.factor(pd[,6])

pd2 <- tax_table(pd)
taxa_names(pd2) <- rownames(pd)
colnames(pd2) <- colnames(tax_table(phy))

tax_table(phy) <- pd2
phy <- subset_taxa(phy, is.finite(Genus))
#phy <- tax_glom(phy, taxrank = rank_names(phy)[6])

phy <- merge_taxa(phy, taxa_names(subset_taxa(phy, Genus %in% "Bacteria")))
tax_table(phy)[,6][is.na(tax_table(phy)[,6]),] <- "Bacteria"
phy <- merge_taxa(phy, taxa_names(subset_taxa(phy, Genus %in% "Chloroplast")))
tax_table(phy)[,6][is.na(tax_table(phy)[,6]),] <- "Chloroplast"
phy <- merge_taxa(phy, taxa_names(subset_taxa(phy, Genus %in% "Archaea")))
tax_table(phy)[,6][is.na(tax_table(phy)[,6]),] <- "Archaea"
phy <- merge_taxa(phy, taxa_names(subset_taxa(phy, Genus %in% "Mitochondria")))
tax_table(phy)[,6][is.na(tax_table(phy)[,6]),] <- "Mitochondria"
phy <- merge_taxa(phy, taxa_names(subset_taxa(phy, Genus %in% "Candidatus Hodgkinia")))
tax_table(phy)[,6][is.na(tax_table(phy)[,6]),] <- "Candidatus Hodgkinia"
phy <- merge_taxa(phy, taxa_names(subset_taxa(phy, Genus %in% "Candidatus Sulcia")))
tax_table(phy)[,6][is.na(tax_table(phy)[,6]),] <- "Candidatus Sulcia"
phy <- merge_taxa(phy, taxa_names(subset_taxa(phy, Genus %in% "Ophiocordyceps")))
tax_table(phy)[,6][is.na(tax_table(phy)[,6]),] <- "Ophiocordyceps"


phy_first <- phy
phy_first_rel <- transform_sample_counts(phy,  function(x){x/sum(x)})


plot_bar(phy_first, fill = "Genus") + 
           geom_bar(stat = "identity") +
           scale_fill_brewer(palette = "Paired") +
           scale_x_discrete(limits = rownames(meta[order(meta$type),])) +
           theme_bw() +
           theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))


phy_b1 <- subset_samples(phy_first_rel, batch %in% "B1")
phy_b2 <- subset_samples(phy_first_rel, batch %in% "B2")
phy_b3 <- subset_samples(phy_first_rel, batch %in% "B3")
phy_b4 <- subset_samples(phy_first_rel, batch %in% "B4")

b1_plot <- plot_bar(phy_b1, , fill = "Genus") + 
           geom_bar(stat = "identity") +
           scale_fill_brewer(palette = "Paired") +
           scale_x_discrete(limits = rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$species),])) +
           theme_bw() +
           theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
           

b2_plot <- plot_bar(phy_b2, , fill = "Genus") + 
           geom_bar(stat = "identity") +
           scale_fill_brewer(palette = "Paired") +
           scale_x_discrete(limits = rownames(sample_data(phy_b2)[order(sample_data(phy_b2)$type, sample_data(phy_b2)$genus, sample_data(phy_b2)$species),])) +
           theme_bw() +
           theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

b3_plot <- plot_bar(phy_b3, , fill = "Genus") + 
           geom_bar(stat = "identity") +
           scale_fill_brewer(palette = "Paired") +
           scale_x_discrete(limits = rownames(sample_data(phy_b3)[order(sample_data(phy_b3)$type, sample_data(phy_b3)$genus, sample_data(phy_b3)$species),])) +
           theme_bw() +
           theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

b4_plot <- plot_bar(phy_b4, , fill = "Genus") + 
           geom_bar(stat = "identity") +
           scale_fill_brewer(palette = "Paired") +
           scale_x_discrete(limits = rownames(sample_data(phy_b4)[order(sample_data(phy_b4)$type, sample_data(phy_b4)$genus, sample_data(phy_b4)$species),])) +
           theme_bw() +
           theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


gar <- ggarrange(b1_plot, b2_plot, b3_plot, b4_plot, common.legend = TRUE)

ggsave("ggarrange.png", gar)

ggsave("b1.png", b1_plot)
ggsave("b2.png", b2_plot)
ggsave("b3.png", b3_plot)
ggsave("b4.png", b4_plot)































### Finding the best tip agglomeration 





































### Host tree for all data 

unique(meta[meta$genus %in% c(x"Maoricicada"), colnames(meta) %in% "species"]$species)

meta[meta$species %in% "mangu",5]



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


####### FUN #######

library(cluster)

otu <- data.frame(otu_table(phy1))

dd <- as.dist(cophenetic.phylo(phy_tree(phy1)))
psclust = cutree(as.hclust(hclust(dd)), h = 0.1)
cliques = levels(factor(psclust))[tapply(psclust, factor(psclust))]

notu <- c()
rnames <- c()
for(i in unique(cliques)) {
	a <- otu[rownames(otu) %in% names(psclust[psclust == i]),]
	rnames <- c(rnames, rownames(a)[1])
	notu <- rbind(notu, apply(a, 2, sum))
}

notu <- data.frame(notu)
rownames(notu) <- rnames


notu$taxa <- rownames(notu)
colnames(notu) <- colnames(otu_table(phy1))

otu2 <- gather(notu, key = sample, value = abund, 1:(dim(notu)[2]-1))
#otu2$sample <- gsub("\\.", "", otu2$sample)

id <- sample_data(phy1)$ID
id <- as.character(id)
met <- data.frame(sample_data(phy1))

otu2 <- cbind(otu2, met[match(otu2$sample, id),])
colnames(otu2)
str(otu2)

ggplot(aes(x = type, y = log(abund)), data = otu2) +
	geom_boxplot() +
	theme_bw()


ag <- aggregate(otu2$abund, by = list(otu2$ID, otu2$type, otu2$batch), function(x){log(length(x[x > 0]))})
ggplot(ag, aes(x = Group.1, y = x, col = Group.2, shape = Group.3)) +
	geom_point(size = 10, alpha = 0.8) +
	geom_hline(yintercept = mean(ag[ag[,2] %in% "control",3])) +
	scale_color_brewer(palette = "Paired") +
    scale_x_discrete(limits = as.character(ag[order(ag$Group.2),1])) +
    theme_bw()
    


























#### Controls

projpc <- subset_samples(phy1, type %in% "control")
projpc_transfer <- subset_samples(projpc, ID %in% c("T1", "T2","T4","T6", "T22", "T24", "T25"))
projpc_diss <- subset_samples(projpc, ID %in% c("T11", "T12", "T13", "T14", "T15", "T16"))
projpc_wash <- subset_samples(projpc, ID %in% c("T17", "T18", "T19", "T20", "T21"))
projpc_extr <- subset_samples(projpc, ID %in% c("T8", "T9", "T10", "PlateControlPC1"))
projpc_pcr <- subset_samples(projpc, ID %in% c("PcrcontrolPCR1", "PcrcontrolPCR2","PcrcontrolPCR3","pcrcontrolPCR4", "PlateControlPC1"))
projpc_marsprojk <- subset_samples(projpc, ID %in% c("MARSblank1"))

proj_control <- list(projpc_transfer, projpc_diss, projpc_wash, projpc_extr, projpc_pcr, projpc_marsprojk)
pcontrol_names <- c()
for(i in 1:6) {
	phy <- proj_control[[i]] 
	otu <- apply(otu_table(phy), 2, function(x){x/sum(x)})
	otu2 <- apply(otu, 1, mean)
	otu2 <- otu2[otu2 > 0]
	lower <- quantile(otu2, probs = seq(0,1,0.05))[2]
	pcontrol_names <- c(pcontrol_names, names(otu2[otu2 > lower]))	
}

#tnames_projpc_transfer <- taxa_names(projpcphy)
#tnames_projpc_diss <- taxa_names(projpcphy)
#tnames_projpc_wash <- taxa_names(projpcphy)
#tnames_projpc_extr <- taxa_names(projpcphy)
#tnames_projpc_pcr <- taxa_names(projpcphy)
#tnames_controls <- c(tnames_projpc_transfer, tnames_projpc_diss, tnames_projpc_wash, tnames_projpc_extr, tnames_projpc_pcr)
#table(taxa_names(projptphy) %in% tnames_controls)
#table(taxa_names(projptphy) %in% taxa_names(projpcphy))























################## Filtering ##############


#### Removing taxonomic units 

phy1 <- phylo	#6871 ASVs
phy1 <- tip_glom(phy1, h = 0.03)	#3544
save.image("start.RData")
phy1 <- subset_taxa(phy1, !Genus %in% c("D_5__Candidatus Sulcia"))	#3520 ASVs
phy1 <- subset_taxa(phy1, !Kingdom %in% c("D_0__Eukaryota"))	#3519
phy1 <- subset_taxa(phy1, !Kingdom %in% c("D_0__Archaea"))	#3402
phy1 <- subset_taxa(phy1, !Family %in% c("D_4__Mitochondria"))	#3360
phy1 <- subset_taxa(phy1, !Order %in% c("D_3__Chloroplast"))	#3346
phy1 <- subset_taxa(phy1, !Genus %in% c("D_5__Candidatus Hodgkinia"))	#3329
phy1 <- subset_taxa(phy1, !is.na(Kingdom))	#3250
phy1 <- subset_taxa(phy1, !Kingdom %in% "Unassigned")	#2378
phy1 <- subset_taxa(phy1, taxa_sums(phy1) > 0)
phy1 <- subset_samples(phy1, sample_sums(phy1) > 0)



#### Subsetting just gut samples 
 
phy1 <- subset_samples(phy1, type %in% c("gut", "control"))
































#### Control quantile by sample quantile plot

projp_mrb_qant <- list()
for(i in 1:500){
	projp <- subset_samples(phy1, batch %in% "B4" & !type %in% "control")
	projp_same_names <- sample_names(projp)
	projp_same_names_pick <- projp_same_names[round(sample(1:length(projp_same_names), 2, replace=FALSE))]
	projp <- subset_samples(projp, sample_names(projp) %in% projp_same_names_pick)
	projp <- subset_samples(projp, sample_sums(projp) > 0)
	projp <- subset_taxa(projp, taxa_sums(projp) > 0)
	projp_mrb <- mrb(projp)
	projp_mrb_qant_temp <- c()
	for(j in 1:length(projp_mrb)){	
		taxon <- names(projp_mrb[j])
		quant <- ecdf(as.numeric(projp_mrb))(as.numeric(projp_mrb[j]))
		projp_mrb_qant_temp <- rbind(projp_mrb_qant_temp, data.frame(taxon, as.numeric(quant)))
		}
		projp_mrb_qant[[i]] <- projp_mrb_qant_temp
	} 

projpc_mrb_qant <- list()
for(i in 1:500){
	projp <- subset_samples(phy1, ID %in% c("T1", "T2","T4","T6", "T22", "T24", "T25", "T11", "T12", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T8", "T9", "T10", "PlateControlPC1", "PcrcontrolPCR1", "PcrcontrolPCR2","PcrcontrolPCR3","pcrcontrolPCR4", "PlateControlPC1"))
	projp_same_names <- sample_names(projp)
	projp_same_names_pick <- projp_same_names[round(sample(1:length(projp_same_names), 2, replace=FALSE))]
	projp <- subset_samples(projp, sample_names(projp) %in% projp_same_names_pick)
	projp <- subset_samples(projp, sample_sums(projp) > 0)
	projp <- subset_taxa(projp, taxa_sums(projp) > 0)
	projp_mrb <- mrb(projp)
	projp_mrb_qant_temp <- c()
	for(j in 1:length(projp_mrb)){	
		taxon <- names(projp_mrb[j])
		quant <- ecdf(as.numeric(projp_mrb))(as.numeric(projp_mrb[j]))
		projp_mrb_qant_temp <- rbind(projp_mrb_qant_temp, data.frame(taxon, as.numeric(quant)))
		
		}
		projpc_mrb_qant[[i]] <- projp_mrb_qant_temp
	} 

projp_mrb_qantSIM <- list()
for(i in 1:500){
	projp <- subset_samples(phy1, batch %in% "B4" & !type %in% "control")
	projp_same_names <- sample_names(projp)
	projp_same_names_pick <- projp_same_names[round(sample(1:length(projp_same_names), 2, replace=FALSE))]
	projp <- subset_samples(projp, sample_names(projp) %in% projp_same_names_pick)
	projp <- subset_samples(projp, sample_sums(projp) > 0)
	projp <- subset_taxa(projp, taxa_sums(projp) > 0)
	projp_mrb <- mrb(projp)
	projp_mrb_qant_temp <- c()
	for(j in 1:length(projp_mrb)){	
		taxon <- names(projp_mrb[j])
		quant <- ecdf(as.numeric(projp_mrb))(as.numeric(projp_mrb[j]))	
		projp_mrb_qant_temp <- rbind(projp_mrb_qant_temp, data.frame(taxon, as.numeric(quant)))
		
		}
		projp_mrb_qantSIM[[i]] <- projp_mrb_qant_temp
	} 

projpc_mrb_qantSIM <- list()
for(i in 1:500){
	projp <- subset_samples(phy1, batch %in% "B4" & !type %in% "control")
	projp_same_names <- sample_names(projp)
	projp_same_names_pick <- projp_same_names[round(sample(1:length(projp_same_names), 2, replace=FALSE))]
	projp <- subset_samples(projp, sample_names(projp) %in% projp_same_names_pick)
	projp <- subset_samples(projp, sample_sums(projp) > 0)
	projp <- subset_taxa(projp, taxa_sums(projp) > 0)
	projp_mrb <- mrb(projp)
	projp_mrb_qant_temp <- c()
	for(j in 1:length(projp_mrb)){	
		taxon <- names(projp_mrb[j])
		quant <- ecdf(as.numeric(projp_mrb))(as.numeric(projp_mrb[j]))
		projp_mrb_qant_temp <- rbind(projp_mrb_qant_temp, data.frame(taxon, as.numeric(quant)))
		}
		projpc_mrb_qantSIM[[i]] <- projp_mrb_qant_temp
	} 


projp_mrb_qant2 <- do.call("rbind", projp_mrb_qant)
projpc_mrb_qant2 <- do.call("rbind", projpc_mrb_qant)
projp_mrb_qantSIM2 <- do.call("rbind", projp_mrb_qantSIM)
projpc_mrb_qantSIM2 <- do.call("rbind", projpc_mrb_qantSIM)


projp_mrb_qant3 <- aggregate(projp_mrb_qant2[,2], by = list(projp_mrb_qant2[,1]), function(x){mean(x[x > 0])})
projpc_mrb_qant3 <- aggregate(projpc_mrb_qant2[,2], by = list(projpc_mrb_qant2[,1]), function(x){mean(x[x > 0])})
projp_mrb_qantSIM3 <- aggregate(projp_mrb_qantSIM2[,2], by = list(projp_mrb_qantSIM2[,1]), function(x){mean(x[x > 0])})
projpc_mrb_qantSIM3 <- aggregate(projpc_mrb_qantSIM2[,2], by = list(projpc_mrb_qantSIM2[,1]), function(x){mean(x[x > 0])})


projp_mrb_qant3_upper <- aggregate(projp_mrb_qant2[,2], by = list(projp_mrb_qant2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[20]  })
projpc_mrb_qant3_upper <- aggregate(projpc_mrb_qant2[,2], by = list(projpc_mrb_qant2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[20]  })
projp_mrb_qantSIM3_upper <- aggregate(projp_mrb_qantSIM2[,2], by = list(projp_mrb_qantSIM2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[20]  })
projpc_mrb_qantSIM3_upper <- aggregate(projpc_mrb_qantSIM2[,2], by = list(projpc_mrb_qantSIM2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[20]  })

projp_mrb_qant3_lower <- aggregate(projp_mrb_qant2[,2], by = list(projp_mrb_qant2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[2]   })
projpc_mrb_qant3_lower <- aggregate(projpc_mrb_qant2[,2], by = list(projpc_mrb_qant2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[2]  })
projp_mrb_qantSIM3_lower <- aggregate(projp_mrb_qantSIM2[,2], by = list(projp_mrb_qantSIM2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[2]   })
projpc_mrb_qantSIM3_lower <- aggregate(projpc_mrb_qantSIM2[,2], by = list(projpc_mrb_qantSIM2[,1]), function(x){  quantile(x, probs = seq(0,1,0.05))[2]   })



qf <- projp_mrb_qant3[match(projpc_mrb_qant3[,1], projp_mrb_qant3[,1]),]

qfupper <- projp_mrb_qant3_upper[match(projpc_mrb_qant3[,1], projp_mrb_qant3_upper[,1]),]
qfupper <- qfupper[is.finite(qfupper[,1]),]
qflower <- projp_mrb_qant3_lower[match(projpc_mrb_qant3[,1], projp_mrb_qant3_lower[,1]),]
qflower <- qflower[is.finite(qflower[,1]),]

qf <- projp_mrb_qant3[match(projpc_mrb_qant3[,1], projp_mrb_qant3[,1]),]

qf2upper <- projpc_mrb_qant3_upper[match(projpc_mrb_qant3[,1], projpc_mrb_qant3_upper[,1]),]
qf2upper <- qf2upper[is.finite(qf[,1]),]

qf2lower <- projpc_mrb_qant3_lower[match(projpc_mrb_qant3[,1], projpc_mrb_qant3_lower[,1]),]
qf2lower <- qf2lower[is.finite(qf[,1]),]

qf2 <- projpc_mrb_qant3[is.finite(qf[,1]),]
qf <- qf[is.finite(qf[,1]),]


qf3 <- data.frame(name = qf[,1], sample = qf[,2], control = qf2[,2], sample_upper = qfupper[,2], sample_lower = qflower[,2], control_lower = qf2lower[,2], control_upper = qf2upper[,2])



qfSIM <- projp_mrb_qantSIM3[match(projpc_mrb_qantSIM3[,1], projp_mrb_qantSIM3[,1]),]
qf2SIM <- projpc_mrb_qantSIM3[is.finite(qfSIM[,1]),]
qfSIM <- qfSIM[is.finite(qfSIM[,1]),]
qf3SIM <- data.frame(name = qfSIM[,1], sample = qfSIM[,2], control = qf2SIM[,2])


ttab <- data.frame(tax_table(phy1))
ttab2 <- ttab[rownames(ttab) %in% qf3[,1],]
qf3$genus <- as.character(ttab2[,6])
qf3$phylum <- as.character(ttab2[,2])

ttab <- data.frame(tax_table(phy1))
ttab2 <- ttab[rownames(ttab) %in% qf3SIM[,1],]
qf3SIM$genus <- as.character(ttab2[,6])
qf3SIM$phylum <- as.character(ttab2[,2])


qf3$genus <- unlist(lapply(str_split(qf3$genus, "__"), "[", 2))
qf3$phylum <- unlist(lapply(str_split(qf3$phylum, "__"), "[", 2))

qf3SIM$genus <- unlist(lapply(str_split(qf3SIM$genus, "__"), "[", 2))
qf3SIM$phylum <- unlist(lapply(str_split(qf3SIM$phylum, "__"), "[", 2))



qf3$sample_bar <- qf3$sample_upper - qf3$sample_lower
qf3$control_bar <- qf3$control_upper - qf3$control_lower

hist(c(qf3$sample_bar, qf3$control_bar), breaks = 20)

cuts <- cut(sort(c(qf3$sample_bar, qf3$control_bar)), breaks = 5)
uno <- range(sort(c(qf3$sample_bar, qf3$control_bar))[cuts == unique(cuts)[1]])
dos <- range(sort(c(qf3$sample_bar, qf3$control_bar))[cuts == unique(cuts)[2]])
tres <- range(sort(c(qf3$sample_bar, qf3$control_bar))[cuts == unique(cuts)[3]])
quat <- range(sort(c(qf3$sample_bar, qf3$control_bar))[cuts == unique(cuts)[4]])
cinc <- range(sort(c(qf3$sample_bar, qf3$control_bar))[cuts == unique(cuts)[5]])



cquant_plot <- ggplot() +	
	#geom_point(data = qf3SIM, aes(x = control, y = sample), shape = 1, color = "grey") +
	geom_errorbar(data = qf3[qf3$sample_bar <= uno[2] & qf3$sample_bar >= uno[1],], aes(x = control, ymin=sample_lower, ymax=sample_upper), alpha = 0.5, width=.01, position=position_dodge(.5), col = "black") +
	geom_errorbar(data = qf3[qf3$sample_bar <= dos[2] & qf3$sample_bar >= dos[1],], aes(x = control, ymin=sample_lower, ymax=sample_upper), alpha = 0.5, width=.01, position=position_dodge(.5), col = "grey60") +
	geom_errorbar(data = qf3[qf3$sample_bar <= tres[2] & qf3$sample_bar >= tres[1],], aes(x = control, ymin=sample_lower, ymax=sample_upper), alpha = 0.5, width=.01, position=position_dodge(.5), col = "grey70") +
	geom_errorbar(data = qf3[qf3$sample_bar <= quat[2] & qf3$sample_bar >= quat[1],], aes(x = control, ymin=sample_lower, ymax=sample_upper), alpha = 0.5, width=.01, position=position_dodge(.5), col = "grey80") +
	geom_errorbar(data = qf3[qf3$sample_bar <= cinc[2] & qf3$sample_bar >= cinc[1],], aes(x = control, ymin=sample_lower, ymax=sample_upper), alpha = 0.5, width=.01, position=position_dodge(.5), col = "grey99") +
	geom_errorbarh(data = qf3[qf3$control_bar <= uno[2] & qf3$control_bar >= uno[1],], aes(y = sample, xmin=control_lower, xmax=control_upper), alpha = 0.5, height=.01, position=position_dodge(.5), col = "black") +
	geom_errorbarh(data = qf3[qf3$control_bar <= dos[2] & qf3$control_bar >= dos[1],], aes(y = sample, xmin=control_lower, xmax=control_upper), alpha = 0.5, height=.01, position=position_dodge(.5), col = "grey60") + 
	geom_errorbarh(data = qf3[qf3$control_bar <= tres[2] & qf3$control_bar >= tres[1],], aes(y = sample, xmin=control_lower, xmax=control_upper), alpha = 0.5, height=.01, position=position_dodge(.5), col = "grey70") +
	geom_errorbarh(data = qf3[qf3$control_bar <= quat[2] & qf3$control_bar >= quat[1],], aes(y = sample, xmin=control_lower, xmax=control_upper), alpha = 0.5, height=.01, position=position_dodge(.5), col = "grey80") +
	geom_errorbarh(data = qf3[qf3$control_bar <= cinc[2] & qf3$control_bar >= cinc[1],], aes(y = sample, xmin=control_lower, xmax=control_upper), alpha = 0.5, height=.01, position=position_dodge(.5), col = "grey99") +
	geom_point(data = qf3, aes(x = control, y = sample, col = phylum), shape = 16, size = 1) +
	geom_text_repel(data = qf3, aes(x = control, y = sample, label = genus, col = phylum), hjust = 0, vjust = 0, size = 2) +
	xlim(0,1) +
	ylim(0,1) + 
	xlab("Control Quantile") +
	ylab("Non-Control Quantile") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggarrange(cquant_plot, cquant_plot2, common.legend = TRUE)

ggsave("cquant_plot.png", cquant_plot, scale = 0.9)



plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 2))
lines(polygon(density(qf3SIM$control), col = "lightgrey", border = NA))
lines(polygon(density(qf3SIM$sample), col = rgb(1, 0, 0,0.1), border = NA))
lines(density(qf3$control), col = "black")
lines(density(qf3$sample), col = "red")



















#### Took out 23 taxa

# 26 control samples 

controls <- subset_samples(phy1, ID %in% c("T1", "T2","T4","T6", "T22", "T24", "T25", "T11", "T12", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T8", "T9", "T10", "PlateControlPC1", "PcrcontrolPCR1", "PcrcontrolPCR2","PcrcontrolPCR3","pcrcontrolPCR4", "PlateControlPC1"))

c2 <- apply(otu_table(controls), 1, sum)
nam <- names(c2[c2 > 0])
phy1 <- subset_taxa(phy1, !taxa_names(phy1) %in% nam)
phy1 <- subset_samples(phy1, sample_sums(phy1) > 0)
















### removing 5 projk control asvs 

projpc_marsprojk <- subset_samples(phy1, ID %in% c("MARSblank1"))
nam <- names(taxa_sums(projpc_marsprojk)[taxa_sums(projpc_marsprojk) > 0])

phy1 <- subset_taxa(phy1, !taxa_names(phy1) %in% nam)
phy1 <- subset_samples(phy1, sample_sums(phy1) > 0)
phy1 <- subset_samples(phy1, !type %in% "control")













### Removing samples with low total abundance


ss <- log(sample_sums(phy1))
nam <- names(ss[exp(ss) < 200])
phy1 <- subset_samples(phy1, !sample_names(phy1) %in% nam)
sample_data(phy1)$depth <- sample_sums(phy1)













### 

tsum <- taxa_sums(phy1)
hist(log(tsum[tsum > 0]), breaks = 100)

b1 <-  subset_samples(phy1, batch %in% c("B1"))
b2 <- subset_samples(phy1, batch %in% c("B2"))
b3 <- subset_samples(phy1, batch %in% "B3")
b4 <- subset_samples(phy1, batch %in% "B4")

tsum <- taxa_sums(b1b2); hist(log(tsum[tsum > 0]), breaks = 100)

tsum <- taxa_sums(b3); hist(log(tsum[tsum > 0]), breaks = 100)

tsum <- taxa_sums(b4); hist(log(tsum[tsum > 0]), breaks = 100)


# 15th percentil 
tsum <- taxa_sums(b1); hist(log(tsum[tsum > 0]), breaks = 100)
tsumfilt <- as.numeric(quantile(log(tsum[tsum > 0]), prob = seq(0,1,0.05))[4])
b1low <- names(tsum[log(tsum) > tsumfilt])

# 15th percentile 
tsum <- taxa_sums(b2); hist(log(tsum[tsum > 0]), breaks = 100)
tsumfilt <- as.numeric(quantile(log(tsum[tsum > 0]), prob = seq(0,1,0.05))[4])
b2low <- names(tsum[log(tsum) > tsumfilt])

# 10th percentile 
tsum <- taxa_sums(b3); hist(log(tsum[tsum > 0]), breaks = 100)
tsumfilt <- as.numeric(quantile(log(tsum[tsum > 0]), prob = seq(0,1,0.05))[3])
b3low <- names(tsum[log(tsum) > tsumfilt])

# 10th percentile 
tsum <- taxa_sums(b4); hist(log(tsum[tsum > 0]), breaks = 100)
tsumfilt <- as.numeric(quantile(log(tsum[tsum > 0]), prob = seq(0,1,0.05))[3])
b4low <- names(tsum[log(tsum) > tsumfilt])



allbatch <- phy1
tsum <- taxa_sums(allbatch); hist(log(tsum[tsum > 0]), breaks = 100)
tsumfilt <- as.numeric(quantile(log(tsum[tsum > 0]), prob = seq(0,1,0.01))[6]); abline(v = tsumfilt)
allbatchlow <- names(tsum[log(tsum) <= tsumfilt])



phy1 <- subset_taxa(phy1, !taxa_names(phy1) %in% allbatchlow)
phy1 <- subset_samples(phy1, sample_sums(phy1) > 0)

save.image()


















#### Abund prev by batch 

otu_b1 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B1"), taxrank = "Family") ))
otu_b2 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B2"), taxrank = "Family")   ))
otu_b3 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B3"), taxrank = "Family") ))
otu_b4 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B4"), taxrank = "Family")  ))




a <- base::apply(otu_b1, 1, function(x){ sum(x) })
b <- base::apply(otu_b2, 1, function(x){  sum(x) })
c <- base::apply(otu_b3, 1, function(x){  sum(x) })
d <- base::apply(otu_b4, 1, function(x){   sum(x) })


a2 <- base::apply(otu_b1, 1, function(x){ length(x[x > 0]) })
b2 <- base::apply(otu_b2, 1, function(x){ length(x[x > 0]) })
c2 <- base::apply(otu_b3, 1, function(x){ length(x[x > 0]) })
d2 <- base::apply(otu_b4, 1, function(x){ length(x[x > 0]) })

ref <- tax_glom(subset_samples(phy1, batch %in% "B4"), taxrank = "Family")
tax <- tax_table(ref)


adata <- data.frame(abund = as.numeric(a/sum(a)), prev = as.numeric(a2)/dim(otu_b1)[2], batch = rep("B1", length(a)), tax = as.character(tax[,5]), class = as.character(tax[,3]))
bdata <- data.frame(abund = as.numeric(b/sum(b)), prev = as.numeric(b2)/dim(otu_b1)[2], batch = rep("B2", length(a)), tax = as.character(tax[,5]), class = as.character(tax[,3]) )
cdata <- data.frame(abund = as.numeric(c/sum(c)), prev = as.numeric(c2)/dim(otu_b1)[2], batch = rep("B3", length(a)), tax = as.character(tax[,5]), class = as.character(tax[,3]) )
ddata <- data.frame(abund = as.numeric(d/sum(d)), prev = as.numeric(d2)/dim(otu_b1)[2], batch = rep("B4", length(a)), tax = as.character(tax[,5]), class = as.character(tax[,3]) )

adata$tax <- unlist(lapply(str_split(adata$tax, "__"), "[", 2))
adata$class <- unlist(lapply(str_split(adata$class, "__"), "[", 2))

bdata$tax <- unlist(lapply(str_split(bdata$tax, "__"), "[", 2))
bdata$class <- unlist(lapply(str_split(bdata$class, "__"), "[", 2))

cdata$tax <- unlist(lapply(str_split(cdata$tax, "__"), "[", 2))
cdata$class <- unlist(lapply(str_split(cdata$class, "__"), "[", 2))

ddata$tax <- unlist(lapply(str_split(ddata$tax, "__"), "[", 2))
ddata$class <- unlist(lapply(str_split(ddata$class, "__"), "[", 2))




dd <- data.frame(rbind(adata, bdata, cdata, ddata))


abundprev_plot <- ggplot(dd, aes(x = prev, y = abund, shape = batch)) +
	geom_text_repel(data = dd[dd$prev > 0.5 | dd$abund > 0.02,], aes(label = tax), size = 2) +
	geom_point(size = 3) +
	geom_point(data = dd[dd$prev > 0.5 | dd$abund > 0.02,], aes(x = prev, y = abund, col = class), size = 3) +
	scale_color_brewer(palette = "Paired") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
ggsave("abundprev_plot.png", abundprev_plot, scale = 0.7)


tax[data.frame(tax[,3])[,1] == "D_2__vadinHA49",]


































#### Filtered dataset ordination and highlighting batch effects after taking control taxa out 

#### Total ordination showing technical differences 


phy <- phylo
phy <- subset_taxa(phy, !Genus %in% c("D_5__Candidatus Sulcia"))
phy <- subset_taxa(phy, !Kingdom %in% c("D_0__Eukaryota"))
phy <- subset_taxa(phy, !Kingdom %in% c("D_0__Archaea"))
phy <- subset_taxa(phy, !Family %in% c("D_4__Mitochondria"))
phy <- subset_taxa(phy, !Order %in% c("D_3__Chloroplast"))
phy <- subset_taxa(phy, !Genus %in% c("D_5__Candidatus Hodgkinia"))
phy <- subset_taxa(phy, !Kingdom %in% c("Unassigned"))
phy <- subset_samples(phy, sample_sums(phy) > 100)


ord <- ordinate(pa_phy(phy), method = "PCoA", distance = "jaccard")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy)$batch)
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = batch)) + 
	geom_point(alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	theme_bw() +
	stat_ellipse(type = "t") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op1 <- op



ord <- ordinate(pa_phy(phy), method = "PCoA", distance = "unifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy)$batch)
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = batch)) + 
	geom_point(alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op2 <- op



ord <- ordinate(pa_phy(phy1), method = "PCoA", distance = "unifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy1)$batch)
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = batch)) + 
	geom_point(alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op3 <- op


ord <- ordinate(phy1, method = "PCoA", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy1)$batch)
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = batch)) + 
	geom_point(alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op4 <- op


ordvec <- data.frame(ord$vectors)[,1:2]
ordvec$batch <- as.character(data.frame(meta[match(rownames(ordvec), meta$ID), colnames(meta) == "batch"])$batch)
ordvec <- data.frame(ordvec)
ordvecg <- gather(ordvec, key = axis, value = val, 1:2)

ggplot(ordvecg, aes(x = axis, y = val, col = batch)) +
	geom_boxplot()




ggsave("s.png", ggarrange(op1, op2, op3, op4, common.legend = TRUE))







## Permanova for batch effects 
###### NOT GOING TO WORK BECAUSE THE DISPERSIONS ACROSS BATCH GROUPS ARE NOT HOMOGENIOUS


disJ <- phyloseq::distance(pa_phy(phy), method = "jaccard")
hisU2 <- phyloseq::distance(pa_phy(phy1), method = "unifrac")
hisU3 <- phyloseq::distance(phy1, method = "wunifrac")

adonis(disJ ~ batch + type + species, data = data.frame(sample_data(phy)))
permutest(betadisper(disJ,  data.frame(sample_data(phy))$batch))

adonis(hisU2 ~ batch, data = data.frame(sample_data(phy1)))
permutest(betadisper(hisU2,  data.frame(sample_data(phy1))$batch))

adonis(hisU3 ~ batch, data = data.frame(sample_data(phy1)))
permutest(betadisper(hisU3,  data.frame(sample_data(phy1))$batch))



















##### Are there microbial communities and how complex are they? 

#### Pairwise sample comparisons

phy <- phy_tip_glommed
phy <- subset_taxa(phy, !is.na(Kingdom))
phy <- subset_taxa(phy, !Kingdom %in% "Unassigned")
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- subset_samples(phy, sample_sums(phy) > 0) 
phy <- subset_samples(phy, type %in% c("gut", "control"))
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
allbatch <- phy
tsumfilt <- as.numeric(quantile(log(tsum[tsum > 0]), prob = seq(0,1,0.01))[6]); abline(v = tsumfilt)
allbatchlow <- names(tsum[log(tsum) <= tsumfilt])
ss <- log(sample_sums(phy))
nam <- names(ss[exp(ss) < 200])
phy <- subset_samples(phy, !sample_names(phy) %in% nam)
phy <- subset_taxa(phy, !taxa_names(phy) %in% allbatchlow)
phy <- subset_samples(phy, sample_sums(phy) > 0)

phyotu <- otu_table(phy)

lap <- matrix(nrow = length(sample_names(phy)), ncol = length(sample_names(phy)))
for(i in 1:length(sample_names(phy))) {
	one <- as.numeric(phyotu[,i])/sum(as.numeric(phyotu[,i]))
	for(j in 1:length(sample_names(phy))){
		two <- as.numeric(phyotu[,j])/sum(as.numeric(phyotu[,j]))
		lap[i,j] <- as.numeric(table(one >= 0.01 & two >= 0.01)[2])
		}
	}
save.image()

ids <- colnames(phyotu)
types <- meta[meta$ID %in% ids, colnames(meta) == "type"]
batch <- meta[meta$ID %in% ids, colnames(meta) == "batch"]

lapp <- lap
lapp[lower.tri(lapp, diag = FALSE)] <- NA


lapp <- data.frame(lapp)
a <- as.numeric(as.matrix(lapp[types == "gut" & batch == "B4", types == "gut" & batch == "B4"]))
b <- as.numeric(as.matrix(lapp[types == "control" & batch == "B4", types == "control" & batch == "B4"]))
c <- as.numeric(as.matrix(lapp[types == "gut" & batch == "B4", types == "control" & batch == "B4"]))

abc <- rbind(data.frame(comp = a, temp = rep("Non-control/Non-control", length(a))), data.frame(comp = b, temp = rep("Control/Control", length(b))), data.frame(comp = c, temp = rep("Control/Non-control", length(c))))

consam <- ggplot(abc, aes(x = temp, y = comp)) +
	geom_jitter() +
	scale_x_discrete(limits = c("Non-control/Non-control", "Control/Control", "Control/Non-control")) +
	xlab("") +
	ylab("Number of overlapping ASVs") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("consam.png", consam, scale = 0.6)




library(MASS)
a2 <- dpois(seq(1, max(c(a,b,c), na.rm = TRUE),1), lambda = as.numeric(fitdistr(a[!is.na(a)], "Poisson")[[1]]))
b2 <- dpois(seq(1, max(c(a,b,c), na.rm = TRUE),1), lambda = as.numeric(fitdistr(b[!is.na(b)], "Poisson")[[1]]))
c2 <- dpois(seq(1, max(c(a,b,c), na.rm = TRUE),1), lambda = as.numeric(fitdistr(c[!is.na(c)], "Poisson")[[1]]))

dpois_consam <- ggplot() +
	geom_line(aes(x = seq(1, max(c(a,b,c), na.rm = TRUE),1), y = a2), col = "red", size = 1) +
	geom_line(aes(x = seq(1, max(c(a,b,c), na.rm = TRUE),1), y = b2), col = "lightblue", size = 1) +
	geom_line(aes(x = seq(1, max(c(a,b,c), na.rm = TRUE),1), y = c2), col = "black", size = 1) +
	xlab("Number of overlapping ASVs") +
	ylab("Log Probability") +
	xlim(c(1,5)) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("dpois_consam.png", dpois_consam, scale = 0.4)



a4 <- rpois(n = 1000, lambda = as.numeric(fitdistr(a[!is.na(a)], "Poisson")[[1]]))
a4 <- jitter(a4, factor = 10)

b4 <- rpois(n = 1000, lambda = as.numeric(fitdistr(b[!is.na(b)], "Poisson")[[1]]))
b4 <- jitter(b4, factor = 10)

c4 <- rpois(n = 1000, lambda = as.numeric(fitdistr(c[!is.na(c)], "Poisson")[[1]]))
c4 <- jitter(c4, factor = 10)

ks.test(a4, b4)
ks.test(a4, c4)
ks.test(b4, c4)








phyloseq::distance(phy, method = "unifrac", type = "samples")








### 

library(cluster)
library(ggtree)


otu <- data.frame(otu_table(phy1))
tax <- data.frame(tax_table(phy1))
met <- data.frame(sample_data(phy1))

abund <- apply(otu, 1, sum)
prev <- apply(otu, 1, function(x){ length(x[(x/sum(x)) > 0.01])})
phylum <- unlist(lapply(str_split(tax[, colnames(tax) == "Phylum"], "__"), "[", 2))


pa <- data.frame(abund = log(abund), prev = prev, phylum = phylum)
pa <- gather(pa, key = measure, value = value, 1:2)

pa2 <- pa[pa$measure == "prev",]
pam <- aggregate(pa2$value, by = list(pa2$phylum, pa2$measure), FUN =  max)
pam <- pam[order(-pam[,3]), 1]


ggplot(pa, aes(x = phylum, y = value, col = measure)) +
	geom_violin() +
	theme_bw() +
	scale_x_discrete(limits = pam) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 


pphylum <- tax_glom(phy1, "Phylum")
p <- phy_tree(pphylum)
p$tip.label <- unlist(lapply(str_split(data.frame(tax_table(pphylum))$Phylum, "__"), "[", 2))
gg <- cbind(taxa = p$tip.label, data.frame(tax_table(pphylum)), data.frame(otu_table(pphylum), stringsAsFactors = FALSE, check.names = FALSE))
(ggtree(p)) %<+% gg +
	geom_tiplab(size = 4) +
	theme(legend.position = c(0.8,0.2), legend.title = element_blank(), legend.key = element_blank())

	
	
	
	
	
	
shan <- vegan::diversity(t(otu), "shannon")	
sim <- vegan::diversity(t(otu), "simpson")	

div3 <- cbind(shan = shan, sim = sim, met)
div3$spge <- paste(div3$species, div3$batch)


shanp <- ggplot(div3, aes(x = batch, y = shan)) + 
	geom_boxplot() +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 

simp <- ggplot(div3, aes(x = batch, y = sim)) + 
	geom_boxplot() +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 

batchdiv <- ggarrange(shanp, simp)
	
ggsave("batchdiv.png", batchdiv, scale = 0.6)
	
	


ggplot(div3, aes(x = spge, y = shan, col = genus)) + 
	geom_boxplot() +
	scale_x_discrete(limits = unique(div3[order(div3$batch, div3$genus), colnames(div3) == "spge"])) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
### Bacterial tree and heatmap 
	
nam <- paste( "D_1__", as.character(pam[1:5]), sep = "")
phy2 <- subset_taxa(phy1, Phylum %in% nam)
phy2 <- tip_glom(phy2, h = 0.4)
p <- phy_tree(phy2)
p$tip.label <- paste(seq(1,length(p$tip.label),1), unlist(lapply(str_split(data.frame(tax_table(phy2))$Class, "__"), "[", 2)), sep = "_")

dd <- cbind(taxa = p$tip.label, data.frame(tax_table(phy2)), data.frame(otu_table(phy2), stringsAsFactors = FALSE, check.names = FALSE))
rownames(dd) <- p$tip.label
#dd <- as.matrix(dd)
#dd[is.na(dd)] <- "Unassigned"
#dd <- data.frame(dd)

dd2 <- dd[, 9:dim(dd)[2]]
dd3 <- c()
for(i in unique(paste(met$genus, met$species, met$batch))) {
	 a <- dd2[, colnames(dd2) %in% as.character(met[paste(met$genus, met$species, met$batch) == i, colnames(met) == "ID"])]
	 if(!is.null(dim(a))){
	 	dd3 <- cbind(dd3, apply(a, 1, sum))
	 } else {
	 	dd3 <- cbind(dd3, a)
	 }
}

dd3 <- data.frame(dd3)
colnames(dd3) <- unique(paste(met$genus, met$species, met$batch))

dd3[dd3 > 0] <- log(dd3[dd3 >0])

p <- drop.tip(p, c("26_NA", "27_NA"))

tree <- (ggtree(p)) %<+% dd +
	geom_tiplab(aes(color = Phylum), size = 2) +
	theme(legend.position = c(0.8,0.2), legend.title = element_blank(), legend.key = element_blank())

bat <- unlist(lapply(str_split(colnames(dd3), " "), "[", 3))
gen <- unlist(lapply(str_split(colnames(dd3), " "), "[", 1))
spe <- unlist(lapply(str_split(colnames(dd3), " "), "[", 2))

dd3 <- dd3[, order(bat, gen, spe)]

tree2 <- gheatmap(tree, dd3, offset = 0.4, width=0.5, colnames=F, low = "white", high = "black")

df = get_heatmap_column_position(tree2, by="bottom")
df$batch <- unlist(lapply(str_split(df$label, " "), "[", 3))
ln <- aggregate(df$x, by = list(df$batch), function(x){ max(x)})

treeall <- tree2 +	geom_text(data=df, aes(x, y, label=label), nudge_y=-10, angle=90, size = 2) +
	geom_vline(xintercept = ln$x, col = "grey") +
	scale_x_continuous(breaks = ln$x, labels = ln[,1]) 















### A zoomed in look 

zoom <- c(
		"104_Bacilli", 
		"105_Bacilli", 
		"51_Alphaproteobacteria", 
		"53_Alphaproteobacteria", 
		"54_Alphaproteobacteria",
		"50_Alphaproteobacteria",
		"49_Alphaproteobacteria",
		"34_Deltaproteobacteria",
		"43_NA",
		"40_NA",
		"38_Gammaproteobacteria",
		"37_Gammaproteobacteria",
		"48_Gammaproteobacteria",
		"13_Actinobacteria", 
		"93_Bacteroidia")

ztax <- phy_tree(phy2)$tip.lab[match(zoom, p$tip.lab)]

dd = as.dist(cophenetic.phylo(phy_tree(subset_taxa(phy1, Phylum %in% nam))))
psclust = cutree(as.hclust(agnes(dd)), h = 0.4)

ps2 <- psclust[names(psclust) %in% ztax]

ps3 <- c()
for(i in unique(as.character(ps2))){
	ps3 <- c(ps3, names(psclust[as.character(psclust) == i]))
}


phy3 <- prune_taxa(taxa = ps3, phy1)
phy3 <- tip_glom(phy3, h = 0.2)
p <- phy_tree(phy3)
p$tip.label <- paste(seq(1,length(p$tip.label),1), unlist(lapply(str_split(data.frame(tax_table(phy3))$Class, "__"), "[", 2)), unlist(lapply(str_split(data.frame(tax_table(phy3))$Family, "__"), "[", 2)), sep = "_")

dd <- cbind(taxa = p$tip.label, data.frame(tax_table(phy3)), data.frame(otu_table(phy3), stringsAsFactors = FALSE, check.names = FALSE))
rownames(dd) <- p$tip.label
#dd <- as.matrix(dd)
#dd[is.na(dd)] <- "Unassigned"
#dd <- data.frame(dd)

dd2 <- dd[, 9:dim(dd)[2]]
dd3 <- c()
for(i in unique(paste(met$genus, met$species, met$batch))) {
	 a <- dd2[, colnames(dd2) %in% as.character(met[paste(met$genus, met$species, met$batch) == i, colnames(met) == "ID"])]
	 if(!is.null(dim(a))){
	 	dd3 <- cbind(dd3, apply(a, 1, function(x){sum(x)}))
	 } else {
	 	dd3 <- cbind(dd3, a)
	 }
}

dd3 <- data.frame(dd3)
colnames(dd3) <- unique(paste(met$genus, met$species, met$batch))

dd3[dd3 > 0] <- log(dd3[dd3 >0])

p2 <- root(p, node = MRCA(p, c("62_Bacteroidia_Marinifilaceae", "63_Bacteroidia_Marinifilaceae", "64_Deltaproteobacteria_bacteriap25")))
p2 <- drop.tip(p2, c("24_NA_NA"))

tree <- (ggtree(p2)) %<+% dd +
	geom_tiplab(aes(color = Phylum), size = 2) +
	theme(legend.position = c(0.8,0.2), legend.title = element_blank(), legend.key = element_blank())

bat <- unlist(lapply(str_split(colnames(dd3), " "), "[", 3))
gen <- unlist(lapply(str_split(colnames(dd3), " "), "[", 1))
spe <- unlist(lapply(str_split(colnames(dd3), " "), "[", 2))

dd3 <- dd3[, order(bat, gen, spe)]

tree2 <- gheatmap(tree, dd3, offset = 0.4, width=0.5, colnames=F, low = "white", high = "black")

df = get_heatmap_column_position(tree2, by="bottom")

treezoom <- tree2 +	geom_text(data=df, aes(x, y, label=label), nudge_y=-10, angle=90, size = 2)











#### Plotting and saving trees 

save.image()

ggsave("treeall.png", treeall)
ggsave("treezoom.png", treezoom)

 
















#### Richness of each phylum in quartiles of distribution of each sample's relative abundances 

otu <- data.frame(otu_table(phy1))
tax <- data.frame(tax_table(phy1))

otur <- apply(otu, 2, function(x){x/sum(x)})

oturq <- apply(otur, 2, function(x){quantile(x[x > 0], prob = seq(0,1,0.02))})

dv2 <- c()
for(j in 1:(length(seq(0,1,0.02))-1)) {
	dv <- list()
	for(i in 1:dim(otur)[2]) {
		a <- otur[,i]
		b <- names(a[a <= oturq[j+1,i] & a >= oturq[j,i]])
		c <- tax[rownames(tax) %in% b,]
		#c2 <- c[c$Phylum == "Proteobacteria",]
		
		if(dim(c)[1] > 0) {
			d <- aggregate(c$Class, by = list(c$Class), function(x){length(x)})
			dv[[i]] <- data.frame(id = rep(i, dim(d)[1]), d)
		} else {
			dv[[i]] <- NULL
		}
		}
	e <- do.call("rbind", dv)
	f <- cbind(quant = rep(j, dim(e)[1]), e)
	dv2 <- rbind(dv2, f)
	}

colnames(dv2) <- c("quant", "id", "phylum", "div")
ag <- aggregate(dv2$div, by = list(dv2$phylum), sum)
ag <- ag[order(-ag$x),]
dv3 <- dv2[dv2$phylum %in% as.character(ag[1:7, 1]),]
dv3$quant <- as.factor(dv3$quant)
levels(dv3$quant) <- seq(1,50,1)
dv3$quant2 <- as.numeric(dv3$quant)
dv3$phylum <- unlist(lapply(str_split(dv3$phylum, "__"), "[", 2))

divplot <- ggplot(dv3[as.numeric(dv3$div) > 1,], aes(x = phylum, y = div, shape = quant, color = quant2, fill = quant2)) +
	geom_boxplot() +
	#geom_jitter(alpha = 0.8, size = 3) +
	scale_color_gradientn(colors = rev(heat.colors(50))) +
	scale_fill_gradientn(colors = rev(heat.colors(50))) +
	scale_x_discrete(limits = c("Bacteroidia", "Actinobacteria", "Bacilli", "Clostridia", "Alphaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria")) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(dv3[as.numeric(dv3$div) > 1,], aes(x = phylum, y = log(div))) +
	geom_boxplot() +
	#geom_jitter(alpha = 0.8, size = 3) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 



ggsave("divplot.png", divplot, scale = 0.6)
















###### Co-occurance

#install.packages("cooccur")
library(cooccur)


otu <- data.frame(otu_table(phy1))

p2 <- tax_glom(phy1, "Phylum")

otu2 <- data.frame(otu_table(p2))

otu3 <- as.matrix(otu[10:11, 1:200])
otu3[otu3 > 0] <- 1
otu3 <- data.frame(otu3)

cooccur.otu <- cooccur(otu3, type = "spp_site", thresh = TRUE, spp_names = TRUE)


library(igraph)




phyotu <- data.frame(otu_table(phy1))

cooc <- matrix(nrow = length(taxa_names(phy1)), ncol = length(taxa_names(phy1)))
for(i in 1:length(taxa_names(phy1))) {
	one <- as.numeric(phyotu[i,])/sum(as.numeric(phyotu[i,]))
	for(j in 1:length(taxa_names(phy1))){
		two <- as.numeric(phyotu[j,])/sum(as.numeric(phyotu[j,]))
		if(!is.na(as.numeric(table(one >= 0.01 & two >= 0.01)[2]))){
			cooc[i,j] <- as.numeric(table(one >= 0.01 & two >= 0.01)[2])
		} else { 
			cooc[i,j] <- 0
		}
		}
	}
save.image("cooc_backup.RData")




tax <- data.frame(tax_table(phy1))

ids <- unlist(lapply(str_split(tax[,4], "__"), "[", 2))
types <- meta[meta$ID %in% ids, colnames(meta) == "type"]
batch <- meta[meta$ID %in% ids, colnames(meta) == "batch"]

lapp <- cooc
lapp[lower.tri(lapp, diag = TRUE)] <- NA

lapp <- data.frame(lapp)
a <- as.matrix(lapp)

a[is.na(a) | a < 1] <- 0

rownames(a) <- colnames(a) <- ids

g <- graph.adjacency(a, weighted = TRUE, diag = FALSE, mode = "upper")
g <- simplify(g)
g<- cluster_edge_betweenness(g, weights = E(g)$weight)

save.image()

V(g)$phylum <- unlist(lapply(str_split(tax[,2], "__"), "[", 2))

V(g)$label.color <- "black"
V(g)$label <- NA
V(net)$color <- colors[V(g)$phylum]

E(g)$width <- E(g)$weight
E(g)$edge.color <- "gray80"
E(g)$width <- 1+E(g)$weight/12


plot(g, vertex.shape="none", vertex.label=V(g)$phylum)


plot(g)

E(g)




plot(g)
ed <- get.edgelist(g)
edremove <- ed[ed[,1] == ed[,2],]

g


g2 <- plot(g)

delete_edges(g, edremove)

























############# PHilR

#source("https://bioconductor.org/biocLite.R")
#biocLite("philr")
library(philr)


phy <- phy1

phy <- transform_sample_counts(phy, function(x) x+1)

phy_tree(phy) <- makeNodeLabel(phy_tree(phy), method="number", prefix='n')

name.balance(phy_tree(phy), tax_table(phy), 'n1')


otu.table <- t(otu_table(phy))
tree <- phy_tree(phy)
metadata <- sample_data(phy)
tax <- tax_table(phy)


gp.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')

gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(phy, 'PCoA', distance=gp.dist)
plot_ordination(phy, gp.pcoa, color='genus') + geom_point(size=4)


gp <- data.frame(gp.philr)

colnames(gp)

gp2 <- gather(gp)
gp2$node <- seq(1, dim(gp2)[1],1)


plot(gp2$value, gp2$node)





gpm <- apply(gp, 2, mean)
gpm <- data.frame(mean = gpm, node = seq(1, length(gpm),1))

allbal_plot <- ggplot(gpm, aes(x = node, y = mean)) + 
	geom_point() +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text( hjust = 1)) 

ggsave("allbal_plot.png", allbal_plot, scale = .5)






















#### PGLMM

library(pez)

htree <- read.beast("host.tre")

ggtree(htree) + geom_tiplab()

met <- data.frame(sample_data(phy1), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALS)
otu <- data.frame(otu_table(phy1), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)

sp <- unique(met[, colnames(met) == "species.tree"])

notu <- c()
nspe <- c()
for(i in sp[!is.na(sp)]){
	a <- data.frame(otu[, !is.na(match(met$species.tree, i))])
	nspe <- c(nspe, dim(a)[2])
	b <- apply(a, 1, sum)
	notu <- cbind(notu, b)
}
colnames(notu) <- sp[!is.na(sp)]

nphy <- phyloseq(otu_table(notu, taxa_are_rows = TRUE), tax_table(phy1), phy_tree(phy1))

ord <- ordinate(nphy, method = "PCoA", distance = "unifrac")

ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], nspe = nspe, sp = rownames(ord$vectors))

nspe_plot <- ggplot(ordDF, aes(x = PC2, y = nspe)) +
	geom_point(size = 4) +
	theme_bw() +
	xlab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text( hjust = 1)) 

ggsave("nspe_plot_pc1.png", nspe_plot, scale = 0.4)
ggsave("nspe_plot_pc2.png", nspe_plot, scale = 0.4)


ord <- ordinate(phy1, method = "PCoA", distance = "unifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], genus = sample_data(phy1)$genus)
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = genus)) + 
	stat_ellipse(type = "t", alpha = 0.3) +
	geom_point(alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

opuni <- op
opwuni <- op

ggsave("opuni.png", opuni, scale = 0.6)
ggsave("opwuni.png", opwuni, scale = 0.6)



ord <- ordinate(nphy, method = "PCoA", distance = "unifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], genus = sample_data(phy1)[match(sample_names(nphy), sample_data(phy1)$species.tree), colnames(sample_data(phy1)) == "genus"]$genus  )
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = genus)) + 
	#stat_ellipse(type = "t", alpha = 0.3) +
	geom_point(size = 4, alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

opuni <- op
opwuni <- op

ggsave("nphyopuni.png", opuni, scale = 0.6)
ggsave("nphyopwuni.png", opwuni, scale = 0.6)















#### True tree

ordDF <- data.frame(ord$vectors, nspe = nspe, sp = rownames(ord$vectors))

t <- read.nexus("host.tre")
t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist(lapply(str_split( unlist(lapply(str_split(t$tip.lab, "_"), "[", 1)), "'" ), "[", 2))

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))

tr <- compute.brlen(t, method = "Grafen", power = 0.5)
vtr <- vcv(tr)
vtr <- vtr/(det(vtr)^(1/length(tr$tip.lab)))
vtr <- t(chol(vtr))

ordDF <- ordDF[match(t$tip.label, rownames(ordDF)),]

nsite = length(t$tip.label)
nspp = 5

X <- matrix(ordDF$nspe, nrow=1, ncol=nsite)
XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow = nspp * nsite, ncol = 1)
site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol = 1)), nrow = nspp * nsite, ncol = 1)
sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp), nrow = nspp * nsite, ncol = 1)

Y <- t(ordDF[,1:5])
Y <- matrix(Y, nrow = nspp, ncol = nsite)
rownames(Y) <- 1:nspp; colnames(Y) <- 1:nsite

YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))

rownames(vtr) <- colnames(vtr) <- levels(dat$site)

re.site <- list(1, site = dat$site, covar = vtr)

model <- communityPGLMM(Y ~ X, data = dat, family = "gaussian", sp = dat$sp, site = dat$site, random.effects = list(re.site), REML = TRUE, verbose = FALSE)
summary(model)








### Random tips 

ordDF <- data.frame(ord$vectors, nspe = nspe, sp = rownames(ord$vectors))

t <- read.nexus("host.tre")
t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist(lapply(str_split( unlist(lapply(str_split(t$tip.lab, "_"), "[", 1)), "'" ), "[", 2))

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))

t$tip.label <- t$tip.label[runif(length(t$tip.label), min = 1, max = length(t$tip.label))]

tr <- compute.brlen(t, method = "Grafen", power = 0.5)
vtr <- vcv(tr)
vtr <- vtr/(det(vtr)^(1/length(tr$tip.lab)))
vtr <- t(chol(vtr))

ordDF <- ordDF[match(t$tip.label, rownames(ordDF)),]

nsite = length(t$tip.label)
nspp = 5

X <- matrix(ordDF$nspe, nrow=1, ncol=nsite)
XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow = nspp * nsite, ncol = 1)
site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol = 1)), nrow = nspp * nsite, ncol = 1)
sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp), nrow = nspp * nsite, ncol = 1)

Y <- t(ordDF[,1:5])
Y <- matrix(Y, nrow = nspp, ncol = nsite)
rownames(Y) <- 1:nspp; colnames(Y) <- 1:nsite

YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))

rownames(vtr) <- colnames(vtr) <- levels(dat$site)

re.site <- list(1, site = dat$site, covar = vtr)

model2 <- communityPGLMM(Y ~ X, data = dat, family = "gaussian", sp = dat$sp, site = dat$site, random.effects = list(re.site), REML = TRUE, verbose = FALSE)

summary(model)
summary(model2)






library(phytools)

lambda <- c()
k <-c()
for(i in 1:40) { 

	ctrait <- ordDF[,i]
	names(ctrait) <- rownames(ordDF)
	k <- c(k, phylosig(t, ctrait, method="K", test=FALSE, nsim=1000, se=NULL, start=NULL, control=list()))

	ctrait <- ordDF[,i]
	names(ctrait) <- rownames(ordDF)
	lambda <- c(lambda, phylosig(t, ctrait, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL, control=list())$lambda)

}

psig <- cbind(lambda, k, PC = seq(1, 40, 1))
psig <- data.frame(psig)

ggplot(psig, aes(x = PC, y = lambda)) +
	geom_point()


ctrait <- ctrait[ match(t$tip.label, names(ctrait), nomatch = 0)[match(t$tip.label, names(ctrait),  nomatch = 0) > 0]    ]

phylosig(t, ctrait, method="K", test=TRUE, nsim=1000)











### Host tree and bacterial treeheatmap

htree <- read.beast("host.tre")

t <- htree

t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t@phylo$tip.label <-  unlist(lapply(str_split( unlist(lapply(str_split(t@phylo$tip.lab, "_"), "[", 1)), "'" ), "[", 2))

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))

t@phylo$tip.label <- t@phylo$tip.label[runif(length(t@phylo$tip.label), min = 1, max = length(t@phylo$tip.label))]


ggtree(t)



phy <- subset_samples(phy1, batch %in% c("B1", "B2", "B3"))

phy <- subset_sample(phy, taxa_sums(phy) < exp(4))

met <- data.frame(sample_data(phy), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALS)
otu <- data.frame(otu_table(phy), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)

sp <- unique(met[, colnames(met) == "species.tree"])

notu <- c()
nspe <- c()
for(i in sp[!is.na(sp)]){
	a <- data.frame(otu[, !is.na(match(met$species.tree, i))])
	nspe <- c(nspe, dim(a)[2])
	b <- apply(a, 1, function(x){ length(x[x > 1]) })
	notu <- cbind(notu, b)
}
colnames(notu) <- sp[!is.na(sp)]


nphy <- phyloseq(otu_table(notu, taxa_are_rows = TRUE), tax_table(phy), phy_tree(phy))
nphy <- subset_taxa(nphy, Phylum %in% c("D_1__Proteobacteria", "D_1__Firmicutes", "D_1__Bacteroidetes", "D_1__Actinobacteria", "D_1__Acidobacteria"))

nphyg <- tip_glom(nphy, h = 0.3)
nphyg <- transform_sample_counts(nphyg, function(x){x/sum(x)})
dd <- data.frame(t(data.frame(otu_table(nphyg), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)
#dd <- data.frame(taxa = rownames(dd), dd, check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)

t <- read.nexus("host.tre")
t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist(lapply(str_split( unlist(lapply(str_split(t$tip.lab, "_"), "[", 1)), "'" ), "[", 2))
t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))

ghtree <- (ggtree(t) + geom_tiplab(size = 4)) %<+% dd +
	theme(legend.position = "none", legend.title = element_blank(), legend.key = element_blank())


#dd[dd > 0] <- log(dd[dd >0])
colnames(dd) <- paste(unlist(lapply(str_split(tax_table(nphyg)[,2], "__"), "[", 2)),   unlist(lapply(str_split(tax_table(nphyg)[,3], "__"), "[", 2)), sep = "_")

dd <- dd[, order(colnames(dd))]
colnames(dd) <- paste(seq(1,dim(dd)[2],1), colnames(dd))


ghtree2b123 <- gheatmap(ghtree, dd, offset = 0.09, width=2, colnames=T, low = "white", high = "black", colnames_angle = 90, colnames_offset_y = -2, font.size = 0.8) +
	theme(legend.position = "none", legend.title = element_blank(), legend.key = element_blank())



ggsave("ghtree2b4.png", ghtree2b4, scale = 0.9)
ggsave("ghtree2b123.png", ghtree2b123, scale = 0.9)

ggsave("ghtree.png", ghtree2b123, scale = 0.9)






ord <- ordinate(nphy, method = "PCoA", distance = "wunifrac")



















## Map?  

library(rgdal"")






























### Distance boxplots for genus and species 

pdis_wunifrac <- phyloseq::distance(phy1, method="wunifrac", type="samples")
pdis_unifrac <- phyloseq::distance(phy1, method="unifrac", type="samples")

pdis <- pdis_unifrac
pdis <- as.matrix(pdis)
pdis[lower.tri(pdis, diag = TRUE)] <- NA

sp <- as.character(unique(sample_data(phy1)$species))
sam <- sample_data(phy1)
bat <- as.character(unique(sample_data(phy1)$batch))

co <- c()
for(i in sp) {

	combinedtmp <- c()
	
	for(j in bat) { 
	
	same <- pdis[match(as.character(sam$ID[sam$species == i & sam$batch == j]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species == i & sam$batch == j])]
	diff <- pdis[match(as.character(sam$ID[sam$species == i & sam$batch == j]), rownames(pdis)),  !colnames(pdis) %in% as.character(sam$ID[sam$species == i & sam$batch == j])]
	same2 <- data.frame(dis = as.numeric(same[!is.na(same)]), comp = rep("same", length(as.numeric(same[!is.na(same)]))))
	diff2 <- data.frame(dis = as.numeric(diff[!is.na(diff)]), comp = rep("diff", length(as.numeric(diff[!is.na(diff)]))))
	combined <- rbind(same2, diff2)
	combined1 <- data.frame(combined, batch = rep(j, dim(combined)[1]))
	combinedtmp <- rbind(combinedtmp, combined1)
	
	}
	
	combined2 <- data.frame(combinedtmp, species = rep(i, dim(combinedtmp)[1]))
	co <- rbind(co, combined2)
	
	}


co1 <- aggregate(co$dis, by = list(co$species, co$comp, co$batch), mean)
co2 <- aggregate(co$dis, by = list(co$species, co$comp, co$batch), function(x) sd(x)/sqrt(length(x)))

co1 <- spread(co1, 2, 4)
co1$d <- co1$same - co1$diff
co2 <-  spread(co2, 2, 4)
co1$sam_sterror <- co2$same
co1$dif_sterror <- co2$diff

co3 <- gather(co1, key = comp, value = dis, 3:4)
co3 <- gather(co3, key = comp_error, value = error, 4:5)
#co3[is.na(co3$error), 6] <- max(co3[, 6], na.rm = TRUE)
colnames(co3) <- c("species", "batch", "dif", "comp", "dis", "comp_error", "error")

gplot <- ggplot(co3, aes(x = species, y = dis, col = comp)) +
	#geom_errorbar(aes(ymin=dis-error, ymax=dis+error), width=0.5)+
	geom_boxplot(data = co, aes(x = species, y = dis, fill = comp)) +
	facet_wrap(~batch, ncol=1, strip.position = "left") +
    scale_x_discrete(limits = unique(co3[order(co3$dis), 1])) +
    theme_bw() +
	theme(legend.title = element_blank(), legend.key = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("sp_unifrac.png", gplot, scale = 1.1)






## all combined

pdis <- pdis_unifrac
pdis <- as.matrix(pdis)
pdis[lower.tri(pdis, diag = TRUE)] <- NA

sp <- as.character(unique(sample_data(phy1)$species))
sam <- sample_data(phy1)
bat <- as.character(unique(sample_data(phy1)$batch))

co <- c()
for(i in sp) {
	same <- pdis[match(as.character(sam$ID[sam$species == i]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species == i])]
	diff <- pdis[match(as.character(sam$ID[sam$species == i]), rownames(pdis)),  !colnames(pdis) %in% as.character(sam$ID[sam$species == i])]
	same2 <- data.frame(dis = as.numeric(same[!is.na(same)]), comp = rep("same", length(as.numeric(same[!is.na(same)]))))
	diff2 <- data.frame(dis = as.numeric(diff[!is.na(diff)]), comp = rep("diff", length(as.numeric(diff[!is.na(diff)]))))
	combined <- rbind(same2, diff2)
	combined1 <- data.frame(combined, species = rep(i, dim(combined)[1]))
	co <- rbind(co, combined1)
	}

co1 <- aggregate(co$dis, by = list(co$species, co$comp), mean)
co2 <- aggregate(co$dis, by = list(co$species, co$comp), function(x) sd(x)/sqrt(length(x)))

co1 <- spread(co1, 2, 3)
co1$d <- co1$same - co1$diff
co2 <-  spread(co2, 2, 3)
co1$sam_sterror <- co2$same
co1$dif_sterror <- co2$diff

co3 <- gather(co1, key = comp, value = dis, 2:3)
co3 <- gather(co3, key = comp_error, value = error, 3:4)
#co3[is.na(co3$error), 6] <- max(co3[, 6], na.rm = TRUE)
colnames(co3) <- c("species", "dif", "comp", "dis", "comp_error", "error")

gplot <- ggplot(co3, aes(x = species, y = dis, col = comp)) +
	#geom_errorbar(aes(ymin=dis-error, ymax=dis+error), width=0.5)+
	geom_boxplot(data = co, aes(x = species, y = dis, fill = comp)) +
    scale_x_discrete(limits = unique(co3[order(co3$dis), 1])) +
    theme_bw() +
	theme(legend.title = element_blank(), legend.key = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("sp_unifrac.png", gplot, scale = 0.8)












### All combined but by genus 

pdis <- pdis_wunifrac
pdis <- as.matrix(pdis)
pdis[lower.tri(pdis, diag = TRUE)] <- NA

sp <- as.character(unique(sample_data(phy1)$genus))
sam <- sample_data(phy1)
bat <- as.character(unique(sample_data(phy1)$batch))

co <- c()
for(i in sp) {
	same <- pdis[match(as.character(sam$ID[sam$genus == i]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$genus == i])]
	diff <- pdis[match(as.character(sam$ID[sam$genus == i]), rownames(pdis)),  !colnames(pdis) %in% as.character(sam$ID[sam$genus == i])]
	same2 <- data.frame(dis = as.numeric(same[!is.na(same)]), comp = rep("same", length(as.numeric(same[!is.na(same)]))))
	diff2 <- data.frame(dis = as.numeric(diff[!is.na(diff)]), comp = rep("diff", length(as.numeric(diff[!is.na(diff)]))))
	combined <- rbind(same2, diff2)
	combined1 <- data.frame(combined, genus = rep(i, dim(combined)[1]))
	co <- rbind(co, combined1)
	}

co1 <- aggregate(co$dis, by = list(co$genus, co$comp), mean)
co2 <- aggregate(co$dis, by = list(co$genus, co$comp), function(x) sd(x)/sqrt(length(x)))

co1 <- spread(co1, 2, 3)
co1$d <- co1$same - co1$diff
co2 <-  spread(co2, 2, 3)
co1$sam_sterror <- co2$same
co1$dif_sterror <- co2$diff

co3 <- gather(co1, key = comp, value = dis, 2:3)
co3 <- gather(co3, key = comp_error, value = error, 3:4)
#co3[is.na(co3$error), 6] <- max(co3[, 6], na.rm = TRUE)
colnames(co3) <- c("genus", "dif", "comp", "dis", "comp_error", "error")

gplot <- ggplot(co, aes(x = genus, y = dis, col = comp)) +
	#geom_errorbar(aes(ymin=dis-error, ymax=dis+error), width=0.5)+
	geom_boxplot(data = co, aes(x = genus, y = dis, fill = comp)) +
    scale_x_discrete(limits = unique(co3[order(co3$dis), 1])) +
    theme_bw() +
	theme(legend.title = element_blank(), legend.key = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

pvals <- c()
for(i in sp) {
	mw <- co[co$genus == i,]
	dis1 <- mw[mw$comp == "same", 1]
	dis2 <- mw[mw$comp == "diff", 1]
	test <- wilcox.test(dis1, dis2) 
	pvals <- c(pvals, test$p.value)
}
names(pvals) <- sp
pvalsadj <- p.adjust(pvals, method = "bonferroni", n = length(pvals))



ggsave("ge_wunifrac.png", gplot, scale = 0.8)
















#### MAntel tests 

t <- read.nexus("host.tre")
t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist(lapply(str_split( unlist(lapply(str_split(t$tip.lab, "_"), "[", 1)), "'" ), "[", 2))

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))






phy <- subset_sample(phy, taxa_sums(phy) < exp(4))

met <- data.frame(sample_data(phy), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALS)
otu <- data.frame(otu_table(phy), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)

sp <- unique(met[, colnames(met) == "species.tree"])

notu <- c()
nspe <- c()
for(i in sp[!is.na(sp)]){
	a <- data.frame(otu[, !is.na(match(met$species.tree, i))])
	nspe <- c(nspe, dim(a)[2])
	b <- apply(a, 1, sum)
	notu <- cbind(notu, b)
}
colnames(notu) <- sp[!is.na(sp)]


nphy <- phyloseq(otu_table(notu, taxa_are_rows = TRUE), tax_table(phy), phy_tree(phy))
nphy <- subset_taxa(nphy, Phylum %in% c("D_1__Proteobacteria", "D_1__Firmicutes", "D_1__Bacteroidetes", "D_1__Actinobacteria", "D_1__Acidobacteria"))

nphyg <- tip_glom(nphy, h = 0.3)
nphyg <- transform_sample_counts(nphyg, function(x){x/sum(x)})




tdis <- as.matrix(cophenetic(t))

np <- merge_samples(phy1, "species.tree")
otu_table(np) <- otu_table(np)[  match(rownames(tdis), sample_names(np), nomatch = 0)[match(rownames(tdis), sample_names(np), nomatch = 0) > 0], ]

pdis_wunifrac <- phyloseq::distance(np, method="wunifrac", type="samples")

pdis_wunifrac <- as.matrix(pdis_wunifrac)
rownames(pdis_wunifrac) == rownames(tdis)


mantel.test(tdis, pdis_wunifrac)
vegan::mantel(tdis, pdis_wunifrac, method="pearson", permutations=999)






















#### Hybrids


pdis <- as.matrix(pdis_unifrac)

muta_muta <- pdis[match(as.character(sam$ID[sam$species == "muta"]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species == "muta"])]
muta_hybrid <- pdis[match(as.character(sam$ID[sam$species == "muta" ]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "muta_tuta"])]

tuta_tuta <- pdis[match(as.character(sam$ID[sam$species == "tuta"]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species == "tuta"])]
tuta_hybrid <- pdis[match(as.character(sam$ID[sam$species == "tuta"]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "muta_tuta"])]

tuta_muta_hybrid <- pdis[match(as.character(sam$ID[sam$species_other == "muta_tuta"]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "muta_tuta"])]

westn_westn <- pdis[match(as.character(sam$ID[sam$species_other == "north" ]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "north" ])]
westn_hybrid <- pdis[colnames(pdis) %in% as.character(sam$ID[sam$ID[sam$species_other == "northwestlandica_southwestlandica"]]), match(as.character(sam$ID[sam$species_other == "north" ]), rownames(pdis))]

wests_wests <- pdis[match(as.character(sam$ID[sam$species_other == "south" ]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "south" ])]
wests_hybrid <- pdis[colnames(pdis) %in% as.character(sam$ID[sam$species_other == "northwestlandica_southwestlandica"]), match(as.character(sam$ID[sam$species_other == "south" ]), rownames(pdis))]

wests_westn_hybrid <- pdis[match(as.character(sam$ID[sam$species_other == "northwestlandica_southwestlandica"]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "northwestlandica_southwestlandica"])]


nels <-  pdis[match(as.character(sam$ID[sam$species == "nelsonensis" ]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species == "nelsonensis" ])]
nels_hybrid <-  pdis[match(as.character(sam$ID[sam$species == "nelsonensis" ]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "nelsonensis/tuta" | sam$species_other == "nelsonensis-tuta" ])]
nels_muta_hybrid <-  pdis[match(as.character(sam$ID[sam$species_other == "nelsonensis/tuta" | sam$species_other == "nelsonensis-tuta" ]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$species_other == "nelsonensis/tuta" | sam$species_other == "nelsonensis-tuta" ])]



muta_muta <- data.frame(dis = as.numeric(unique(muta_muta)), comp = rep("muta_muta", length(as.numeric(unique(muta_muta)))))
tuta_tuta <- data.frame(dis = as.numeric(unique(tuta_tuta)), comp = rep("tuta_tuta", length(as.numeric(unique(tuta_tuta)))))
muta_hybrid <- data.frame(dis = as.numeric(unique(muta_hybrid)), comp = rep("muta_hybrid", length(as.numeric(unique(muta_hybrid)))))
tuta_hybrid <- data.frame(dis = as.numeric(unique(tuta_hybrid)), comp = rep("tuta_hybrid", length(as.numeric(unique(tuta_hybrid)))))
tuta_muta_hybrid <- data.frame(dis = as.numeric(unique(tuta_muta_hybrid)), comp = rep("tuta_muta_hybrid", length(as.numeric(unique(tuta_muta_hybrid)))))


westn_westn <- data.frame(dis = as.numeric(unique(westn_westn)), comp = rep("westn_westn", length(as.numeric(unique(westn_westn)))))
westn_hybrid <- data.frame(dis = as.numeric(unique(westn_hybrid)), comp = rep("westn_hybrid", length(as.numeric(unique(westn_hybrid)))))
wests_wests <- data.frame(dis = as.numeric(unique(wests_wests)), comp = rep("wests_wests", length(as.numeric(unique(wests_wests)))))
wests_hybrid <- data.frame(dis = as.numeric(unique(wests_hybrid)), comp = rep("wests_hybrid", length(as.numeric(unique(wests_hybrid)))))
wests_westn_hybrid <- data.frame(dis = as.numeric(unique(wests_westn_hybrid)), comp = rep("wests_westn_hybrid", length(as.numeric(unique(wests_westn_hybrid)))))


nels <- data.frame(dis = as.numeric(unique(nels)), comp = rep("nels", length(as.numeric(unique(nels)))))
nels_hybrid <- data.frame(dis = as.numeric(unique(nels_hybrid)), comp = rep("nels_hybrid", length(as.numeric(unique(nels_hybrid)))))
nels_muta_hybrid <- data.frame(dis = as.numeric(unique(nels_muta_hybrid)), comp = rep("nels_muta_hybrid", length(as.numeric(unique(nels_muta_hybrid)))))



hybrid <- rbind(muta_muta, muta_hybrid, tuta_tuta, tuta_hybrid, tuta_muta_hybrid , westn_westn, westn_hybrid, wests_wests, wests_hybrid, wests_westn_hybrid, nels, nels_hybrid, nels_muta_hybrid)

hybrid$hybrid <- unlist(lapply(str_split(hybrid$comp, "_"), "[", 3))


hybrid_plot <- ggplot(hybrid[hybrid$dis > 0,], aes(x = comp, y = dis, col = hybrid)) +
	geom_boxplot() +
	geom_jitter(width = 0.2) + 
	#scale_x_discrete(limits = unique(co3[order(co3$dis), 1])) +
    theme_bw() +
	theme(legend.title = element_blank(), legend.key = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))



	
ggsave("hybrid2_plot.png", hybrid_plot, scale = 0.8)


	
t.test(hybrid[hybrid$dis > 0 & hybrid$com == "muta_muta", 1], hybrid[hybrid$dis > 0 & hybrid$com == "muta_hybrid",1])
t.test(hybrid[hybrid$dis > 0 & hybrid$com == "tuta_tuta", 1], hybrid[hybrid$dis > 0 & hybrid$com == "tuta_hybrid",1])
t.test(hybrid[hybrid$dis > 0 & hybrid$com == "westn_westn", 1], hybrid[hybrid$dis > 0 & hybrid$com == "westn_hybrid",1])
t.test(hybrid[hybrid$dis > 0 & hybrid$com == "wests_wests", 1], hybrid[hybrid$dis > 0 & hybrid$com == "wests_hybrid",1])

	
	
	
	
	
	
	
	
	































### Biome distance by geographic distance 

library(geosphere)
library(Imap)

sam <- sample_data(phy1)

sam[sam$lat_sign == "neg", colnames(sam) == "lat"] <- -sam[sam$lat_sign == "neg", colnames(sam) == "lat"]
sam[sam$lon_sign == "neg", colnames(sam) == "lon"] <- -sam[sam$lon_sign == "neg", colnames(sam) == "lon"]



ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
   # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
   # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

   if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
   if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
   else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
   else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
   m[tri] <- t(m)[tri]
   return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
   # Returns a matrix (M) of distances between geographic points.
   # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
   # (df.geopoints$lat[j], df.geopoints$lon[j]).
   # The row and column names are given by df.geopoints$name.

   GeoDistanceInMetres <- function(g1, g2){
      # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
      # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
      # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
      # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
      # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
      DistM <- function(g1, g2){
         require("Imap")
         return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
      }
      return(mapply(DistM, g1, g2))
   }

   n.geopoints <- nrow(df.geopoints)

   # The index column is used to ensure we only do calculations for the upper triangle of points
   df.geopoints$index <- 1:n.geopoints

   # Create a list of lists
   list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})

   # Get a matrix of distances (in metres)
   mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

   # Set the row and column names
   rownames(mat.distances) <- df.geopoints$name
   colnames(mat.distances) <- df.geopoints$name

   return(mat.distances)
}



sam2 <- data.frame(id = sam$ID, lat = sam$lat, lon =sam$lon, row.names = rownames(sam))

geodis <- GeoDistanceInMetresMatrix(sam2)
geodis <- as.matrix(geodis)
#geodis[lower.tri(geodis, diag = TRUE)] <- NA
geodis <- data.frame(geodis)
rownames(geodis) <- colnames(geodis) <- rownames(sam2)

pdis <- pdis_wunifrac
pdis <- as.matrix(pdis)
#pdis[lower.tri(pdis, diag = TRUE)] <- NA


sp <- as.character(unique(sample_data(phy1)$genus))


sp <- c("Kikihia", "Rhodopsalta", "Maoricicada", "Amphipsalta")


sam <- sam[sam$genus %in% sp & sam$ID != "P10",]

geodis <- geodis[rownames(geodis) %in% rownames(sam), colnames(geodis) %in% rownames(sam)]
pdis <- pdis[rownames(pdis) %in% rownames(sam), colnames(pdis) %in% rownames(sam)]


co <- c()
for(i in sp) {
	
	sameg <- geodis[match(as.character(sam$ID[sam$genus == i]), rownames(geodis)),  colnames(geodis) %in% as.character(sam$ID[sam$genus == i])]
	samec <- pdis[match(as.character(sam$ID[sam$genus == i]), rownames(pdis)),  colnames(pdis) %in% as.character(sam$ID[sam$genus == i])]
	samegg2 <- gather(data.frame(sameg), key = id, value = distance, 1:(dim(sameg)[2]))
	samecpdis2 <- gather(data.frame(samec), key = id, value = distance, 1:(dim(samec)[2]))
	
	same <- data.frame(id = samegg2[,1], gdis = samegg2[,2], comdis = samecpdis2[,2], com = rep("same", length(samegg2[,1])))
	#same <- gather(data.frame(same), key = distype, value = dis, 2:3)

	diffg <- geodis[match(as.character(sam$ID[sam$genus == i]), rownames(geodis)),  !colnames(geodis) %in% as.character(sam$ID[sam$genus == i])]
	diffc <- pdis[match(as.character(sam$ID[sam$genus == i]), rownames(pdis)),  !colnames(pdis) %in% as.character(sam$ID[sam$genus == i])]
	diffgg2 <- gather(data.frame(diffg), key = id, value = distance, 1:(dim(diffg)[2]))
	diffcpdis2 <- gather(data.frame(diffc), key = id, value = distance, 1:(dim(diffc)[2]))
	
	diff <- data.frame(id = diffgg2[,1], gdis = diffgg2[,2], comdis = diffcpdis2[,2], com = rep("diff", length(diffgg2[,1])))
	#diff <- gather(data.frame(diff), key = distype, value = dis, 2:3)

	samediff <- rbind(same, diff)
	samediff <- data.frame(samediff, taxon = rep(i, dim(samediff)[1]))
	
	co <- rbind(co, samediff)
	
	}


library(betareg)
library(lmtest)

co2 <- co[co$comdis > 0 & co$comdis < 1,]
co2 <- data.frame(co2[as.character(co2$com) == "same",])

betaco <- betareg(comdis ~ gdis + taxon, data = co2)
#betaco <- lm(comdis ~ gdis + taxon, data = co2)


pred <- plogis(coef(betaco)[1] + 
	seq(1, max(co2$gdis, na.rm = TRUE), 1000)*coef(betaco)[2])
	coef(betaco)[3] +
	coef(betaco)[4] +
	coef(betaco)[5] 


#diffpred <- plogis(
	coef(betaco)[1] + 
	seq(1, max(co2$gdis, na.rm = TRUE), 1000)*coef(betaco)[2] +
	coef(betaco)[3] +
	coef(betaco)[4] +
	coef(betaco)[5] +
	coef(betaco)[6] +
	coef(betaco)[7] +
	coef(betaco)[8] +
	coef(betaco)[9] +
	coef(betaco)[10]) 
	
#diffpred <-
	coef(betaco)[1] + 
	seq(1, max(co2$gdis, na.rm = TRUE), 1000)*coef(betaco)[2] +
	coef(betaco)[3] +
	coef(betaco)[4] +
	coef(betaco)[5] 
	
#sampred <-  plogis(
	coef(betaco)[1] + 
	seq(1, max(co2$gdis, na.rm = TRUE), 1000)*coef(betaco)[2] +
	coef(betaco)[4] +
	coef(betaco)[5] +
	
	
#sampred <-  
	coef(betaco)[1] + 
	seq(1, max(co2$gdis, na.rm = TRUE), 1000)*coef(betaco)[2] +
	coef(betaco)[4] +
	coef(betaco)[5] +
	coef(betaco)[6] +
	coef(betaco)[7] +
	coef(betaco)[8] +
	coef(betaco)[9] +
	coef(betaco)[10]
	
	
#newp <- cbind(gdis = seq(1, max(co2$gdis, na.rm = TRUE), 1000), diffpred, sampred)
newp <- cbind(gdis = seq(1, max(co2$gdis, na.rm = TRUE), 1000), pred)

newp <- data.frame(newp)

gcom_plot <- ggplot(co2, aes(x = gdis, y = comdis)) +
	geom_point(alpha = 0.3) +
	geom_line(data = newp, aes(x = gdis, y = pred), col = "red") +
	#geom_line(data = newp, aes(x = gdis, y = sampred), col = "blue") +
	xlab("Geographic Distance") +
	ylab("Community Distance") + 
	theme_bw() +
	theme(legend.title = element_blank(), legend.key = element_blank(), axis.text.x = element_text(hjust = 1))

#gcom_sw <- gcom_plot
#gcom_dw <- gcom_plot
#gcom_s <- gcom_plot
#gcom_d <- gcom_plot


gcom_all <- ggarrange(gcom_sw, gcom_dw, gcom_s, gcom_d)
ggsave("gcom_all.png", gcom_all, scale = 0.8)



est <- coef(summary(betaco))$mean[,"Estimate"]
#est <- coef(summary(betaco))[,"Estimate"]

ster <- coef(summary(betaco))$mean[,"Std. Error"]
#ster <- coef(summary(betaco))[,"Std. Error"]


modp <- data.frame(est = est, ster = ster, comp = names(est))
modp$est2 <- plogis(modp[1,1]+modp$est)
#modp$est2 <- modp[1,1]+modp$est

modp[1, "est2"] <- plogis(modp[1,1])
#modp[1, "est2"] <- modp[1,1]

#modp <- modp[-1,]

breg_plot <- ggplot(modp, aes(x = comp, y = est2)) +
	geom_point(size = 3) +
	geom_errorbar(aes(ymin=est2-ster, ymax=est2+ster), width=0.2) + 
	scale_x_discrete(limits = modp[order(-modp$est2),"comp"]) +
	geom_hline(yintercept = plogis(est[1]), col = "red") +
	#geom_hline(yintercept = est[1], col = "red") +
	theme_bw() +
	theme(legend.title = element_blank(), legend.key = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))



#same_unifrac <- data.frame(modp, t = rep("same_unifrac", dim(modp)[1]))
#diff_unifrac <- data.frame(modp, t = rep("diff_unifrac", dim(modp)[1]))

#same_wunifrac <- data.frame(modp, t = rep("same_wunifrac", dim(modp)[1]))
#diff_wunifrac <- data.frame(modp, t = rep("diff_wunifrac", dim(modp)[1]))


modp <- rbind(same_unifrac, diff_unifrac, same_wunifrac, diff_wunifrac)
modp$samdiff <- unlist(lapply(str_split(modp$t, "_"), "[", 1))
modp$distance <- unlist(lapply(str_split(modp$t, "_"), "[", 2))

breg_plot <- ggplot(modp, aes(x = comp, y = est2)) +
	geom_point(size = 3) +
	geom_errorbar(aes(ymin=est2-ster, ymax=est2+ster), width=0.2) + 
	#scale_x_discrete(limits = modp[order(-modp$est2),"comp"]) +
	#geom_hline(yintercept = plogis(est[1]), col = "red") +
	#geom_hline(yintercept = est[1], col = "red") +
	facet_wrap(~samdiff+distance, scales = "free_y", strip.position = "left") +
	theme_bw() +
	theme(legend.title = element_blank(), legend.key = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

modp2 <- modp[modp$com %in% c("(Intercept)", "gdis"),]

ggplot(modp2, aes(x = comp, y = est2)) +
	geom_point(size = 3) +
	geom_errorbar(aes(ymin=est2-ster, ymax=est2+ster), width=0.2) + 
	facet_wrap(~samdiff+distance, strip.position = "left") 





#ggsave("breg_plot.png", breg_plot, scale = 0.5)

#ggsave("gcom_plots.png", gcom_plot, scale = 0.5)
#ggsave("gcom_plotd.png", gcom_plot, scale = 0.5)
#ggsave("gcomw_plots.png", gcom_plot, scale = 0.5)
#ggsave("gcomw_plotd.png", gcom_plot, scale = 0.5)


































### batch effect ordination for species shared among the three batches (post taking control taxa out as above)


### Merged individuals within species within batch for species samples in 2 or more batches 
### Calculated new counts by taking the mean relative abundance of each ASV across individuals and multiplying by the mean depth of all individuals within that species

phy <- phy1

a <- unique(sample_data(subset_samples(phy, batch %in% "B4"))[, colnames(sample_data(phy)) == "species"])$species
b <- unique(sample_data(subset_samples(phy, batch %in% "B3"))[, colnames(sample_data(phy)) == "species"])$species

sample_data(phy)$species[sample_data(phy)$species %in% c("southwestlandica", "northwestlandica")] <- "westlandica"

b2 <- unique(sample_data(subset_samples(phy, batch %in% c("B1","B2")))[, colnames(sample_data(phy)) == "species"])$species

sp <- as.character(a[a %in% c(as.character(b), as.character(b2))])

bphy <- subset_samples(phy, species %in% sp & batch %in% c("B1", "B2", "B3", "B4"))
bphy_temp <- bphy

bphy <- transform_sample_counts(bphy, function(x){x/sum(x)})

#bphy <- tax_glom(bphy, taxrank = "Genus")

new_merge <- function(phy_object, fun) {
	new_otu <- c()
	for(i in as.character(unique(sample_data(phy_object)$species))) {
	
		df <- data.frame(sample_data(phy_object))
		ids <- as.character(rownames(df[df$species == i,]))
		otu <- as(otu_table(phy_object), "matrix")
		otu2 <- otu[, colnames(otu) %in% ids]
		
		
		if(!is.null(dim(otu2)[2])){
			msum <- round(mean(apply(otu2, 2, sum)), 0)
			otu3 <- apply(otu2, 2, function(x){x/sum(x)})
			new_otu <- cbind(new_otu, data.frame(apply(otu3, 1, function(x){ round( mean(x)*msum, digits =  0)}))[,1])
		} else {
			new_otu <- cbind(new_otu, data.frame(otu2)[,1])
		}
		
		}
		
		new_otu <- data.frame(new_otu)
		colnames(new_otu) <- as.character(unique(sample_data(phy_object)$species))
		rownames(new_otu) <- rownames(otu)
		df2 <- data.frame(species = as.character(unique(sample_data(phy_object)$species)), row.names = as.character(unique(sample_data(phy_object)$species)))
		new_phy <- phyloseq(otu_table(new_otu, taxa_are_rows = TRUE), sample_data(df2), tax_table(tax_table(phy_object)), phy_tree(phy_object))		
	
		return(new_phy)
	}

bphyall <- new_merge(bphy, "species")
bphyprojp <- new_merge(subset_samples(bphy, batch %in% "B4"), "species")
bphyprojk <- new_merge(subset_samples(bphy, batch %in% "B3"), "species")
bphyjays <- new_merge(subset_samples(bphy, batch %in% c("B1","B2")), "species")

bphyall_otu <- otu_table(bphyall)
colnames(bphyall_otu) <- paste(colnames(otu_table(bphyall)), "all", sep = "_")
bphyall_sam <- sample_data(bphyall)
rownames(bphyall_sam) <- paste(rownames(sample_data(bphyall)), "all", sep = "_")

bphyprojp_otu <- otu_table(bphyprojp)
colnames(bphyprojp_otu) <- paste(colnames(otu_table(bphyprojp)), "projp", sep = "_")
bphyprojp_sam <- sample_data(bphyprojp)
rownames(bphyprojp_sam) <- paste(rownames(sample_data(bphyprojp)), "projp", sep = "_")

bphyprojk_otu <- otu_table(bphyprojk)
colnames(bphyprojk_otu) <- paste(colnames(otu_table(bphyprojk)), "projk", sep = "_")
bphyprojk_sam <- sample_data(bphyprojk)
rownames(bphyprojk_sam) <- paste(rownames(sample_data(bphyprojk)), "projk", sep = "_")

bphyjays_otu <- otu_table(bphyjays)
colnames(bphyjays_otu) <- paste(colnames(otu_table(bphyjays)), "jays", sep = "_")
bphyjays_sam <- sample_data(bphyjays)
rownames(bphyjays_sam) <- paste(rownames(sample_data(bphyjays)), "jays", sep = "_")

newotu <- cbind(bphyall_otu, bphyprojp_otu, bphyprojk_otu, bphyjays_otu)
newsam <- rbind(bphyall_sam, bphyprojp_sam, bphyprojk_sam, bphyjays_sam)

newphy <- phyloseq(otu_table(newotu, taxa_are_rows = TRUE), sample_data(newsam), tax_table(tax_table(phy)), phy_tree(phy))
sample_data(newphy)$species <- as.character(lapply(str_split(rownames(sample_data(newphy)), "_"), "[", 1))
sample_data(newphy)$batch <- as.character(lapply(str_split(rownames(sample_data(newphy)), "_"), "[", 2))
newphy <- subset_taxa(newphy, taxa_sums(newphy) > 0)

sample_data(newphy)$depth <- sample_sums(newphy)

logtaxc <- log(taxa_sums(newphy))
newphy_temp <- newphy
newphy <- newphy_temp
newphy <- subset_taxa(newphy, taxa_names(newphy) %in% names(logtaxc[logtaxc < 2]))
newphy <- subset_samples(newphy, sample_sums(newphy) > 0)


phy <- subset_samples(newphy, !batch %in% "all")
ord <- ordinate(phy, method = "PCoA", distance = "unifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy)$batch, species = sample_data(phy)$species, depth = sample_data(phy)$depth)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	geom_point(aes(color = batch, size = depth), alpha = 0.8) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op1 <- op


all_merge <- new_merge(subset_samples(newphy, !batch %in% "all"))
all_merge_otu <- otu_table(all_merge)
colnames(all_merge_otu) <- paste(colnames(otu_table(all_merge)), "all_merge", sep = "_")
all_merge_sam <- sample_data(all_merge)
rownames(all_merge_sam) <- paste(rownames(sample_data(all_merge)), "all_merge", sep = "_")
all_merge_sam$species <- as.character(all_merge_sam$species)
all_merge2 <- phyloseq(otu_table(all_merge_otu, taxa_are_rows = TRUE), sample_data(all_merge_sam), tax_table(tax_table(all_merge)), phy_tree(all_merge))
sample_data(all_merge2)$species <- as.character(lapply(str_split(rownames(sample_data(all_merge2)), "_"), "[", 1))
sample_data(all_merge2)$batch <- as.character(lapply(str_split(rownames(sample_data(all_merge2)), "_"), "[", 2))
all_merge2 <- subset_taxa(all_merge2, taxa_sums(all_merge2) > 0)
phy <- all_merge2
ord <- ordinate(all_merge2, method = "PCoA", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], species = sample_data(phy)$species)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	geom_point(alpha = 0.8, size = 4) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op2 <- op


phy <- subset_samples(newphy, !batch %in% "all")
ord <- ordinate(phy, method = "PCoA", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy)$batch, species = sample_data(phy)$species, depth = sample_data(phy)$depth)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	geom_point(aes(color = batch, size = depth), alpha = 0.8) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op2 <- op




opr <- ggarrange(op1, op2, common.legend = TRUE)

ggsave("share_species_ord_unifrac_wunifrac.png", opr, scale = 0.8)
























#### Ordination of species within each batch separetely

phy <- subset_samples(newphy, batch %in% "projk")
phy <- subset_samples(phy, species %in% sample_data(subset_samples(newphy, batch %in% "projp"))$species)
ord1 <- ordinate(pa_phy(phy), method = "NMDS", distance = "wunifrac")
ord <- ord1
ordDF <- data.frame(PC1 = ord$points[,1], PC2 = ord$points[,2], species = sample_data(phy)$species)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	geom_point(alpha = 0.8, size = 4) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op1 <- op


phy <- subset_samples(newphy, batch %in% "projp")
phy <- subset_samples(phy, species %in% sample_data(subset_samples(newphy, batch %in% "projk"))$species)
ord2 <- ordinate(pa_phy(phy), method = "NMDS", distance = "wunifrac")
ord <- ord1
ordDF <- data.frame(PC1 = ord$points[,1], PC2 = ord$points[,2], species = sample_data(phy)$species)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	geom_point( alpha = 0.8, size = 4) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op2 <- op


o1 <- ord1$points[,c(1,2)]
rownames(o1) <- unlist(lapply(str_split(rownames(o1), "_"), "[", 1))

o2 <- ord2$points[,c(1,2)]
rownames(o2) <- unlist(lapply(str_split(rownames(o2), "_"), "[", 1))

o2 <- o2[match(rownames(o1), rownames(o2)),]


pro <- procrustes(o1, o2)
ordDF <- data.frame(species = rownames(pro$Yrot), PC1 = pro$X[,1], PC2 = pro$X[,2], PC11 = pro$Yrot[,1], PC22 = pro$Yrot[,2])

q <- ordDF[, c(1,2,3)]
colnames(q) <- c("species", "PC1", "PC2")
q$species <- paste(q$species, "projk", sep = "_")
q2 <- ordDF[, c(1,4,5)]
colnames(q2) <- c("species", "PC1", "PC2")
q2$species <- paste(q2$species, "projp", sep = "_")

ordDF <- rbind(q, q2)
ordDF$ord <- unlist(lapply(str_split(ordDF$species, "_"), "[", 2))
ordDF$species <- unlist(lapply(str_split(ordDF$species, "_"), "[", 1))



op <- ggplot(data = ordDF, aes(x = PC1, y = PC2)) + 
	geom_text_repel(aes(label = species), size = 3, segment.colour = "grey", segment.size = 1) +
	geom_line(aes(group = species)) +
	geom_point(aes(color = ord), size = 3) +
	xlab("NMDS1") + 
	ylab("NMDS2") + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op2 <- op

ggsave("projk_projp_procrustes.png", op2, scale = 0.8)





























### shared species batch effects species identificaion

# First, what are the unique ASVs in the common species between batches


phy <- subset_samples(newphy, !batch %in% "all")
ord <- ordinate(pa_phy(phy), method = "PCoA", distance = "jaccard")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy)$batch, species = sample_data(phy)$species)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	#geom_line(aes(group = species), col = "grey") +
	geom_point(aes(color = batch), alpha = 0.8, size = 4) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op1 <- op

usp <- c()
for(i in unique(sample_data(newphy)$species)){
	sp <- subset_samples(newphy, species %in% i)
	sp_otu <- as(otu_table(sp), "matrix")
	sp_otu[sp_otu > 0] <- 1
	sp_occ <- apply(sp_otu, 1, function(x){sum(x)/length(x)})
	usp <- c(usp, names(sp_occ[sp_occ > 1/3 & sp_occ < 2/3]))
}

n <- subset_taxa(newphy, taxa_names(newphy) %in% usp)
n2 <- subset_samples(n, sample_sums(n) > 0)
n3 <- subset_samples(n2, !batch %in% "all")

phy <- n3
ord <- ordinate(pa_phy(phy), method = "PCoA", distance = "jaccard")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy)$batch, species = sample_data(phy)$species)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	#geom_line(aes(group = species), col = "grey") +
	geom_point(aes(color = batch), alpha = 0.8, size = 4) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op2 <- op



usp <- c()
for(i in unique(sample_data(newphy)$species)){
	sp <- subset_samples(newphy, species %in% i)
	sp_otu <- as(otu_table(sp), "matrix")
	sp_otu[sp_otu > 0] <- 1
	sp_occ <- apply(sp_otu, 1, function(x){sum(x)/length(x)})
	usp <- c(usp, names(sp_occ[sp_occ == 1]))
}

n <- subset_taxa(newphy, taxa_names(newphy) %in% usp)
n2 <- subset_samples(n, sample_sums(n) > 0)
n3 <- subset_samples(n2, !batch %in% "all")

phy <- n3
ord <- ordinate(pa_phy(phy), method = "PCoA", distance = "jaccard")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], batch = sample_data(phy)$batch, species = sample_data(phy)$species)
op <- ggplot(ordDF, aes(x = PC1, y = PC2)) + 
	#geom_line(aes(group = species), col = "grey") +
	geom_point(aes(color = batch), alpha = 0.8, size = 4) +
	geom_text_repel(aes(label = species), size = 3) +
	xlab(paste("PC1", round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	ylab(paste("PC2", round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), sep = ": ")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op3 <- op




ggarrange(op1, op2, op3, common.legend = TRUE)
























### DESeq


sample_data(newphy_counts)$species <- as.factor(sample_data(newphy_counts)$species)

ds <- phyloseq_to_deseq2(subset_samples(newphy, !batch %in% "all") , ~species)


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
diagdds = estimateSizeFactors(ds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(newphy)[rownames(sigtab), ], "matrix"))
head(sigtab)


names <- rownames(sigtab)

phyord <- prune_taxa(taxa = names, bphy)
phyord <- subset_samples(phyord, sample_sums(phyord) > 0)
phyord <- merge_samples(phyord, "species")
sample_data(phyord)$species <- rownames(sample_data(phyord))
ordgq <- ordinate(phyord, method = "PCoA", distance = "bray")
ord_filt <- plot_ordination(phyord, ordgq, color = "species", axes = c(1,2)) +
  geom_text(aes(label = species),  size = 3, vjust = 1.5) +
  geom_point(aes(alpha = 0.5, size = 4), alpha = 0.1) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt




























#### RCM ordination with batch effects 

phy_rcm <- RCM(subset_samples(newphy, batch %in% "all"))

plot(phy_rcm)



















####### Differential abundance using MicrobiomeDDA (depends on functions that need to be defined) 


















### top most abundent

phy <- projp

otu <- apply(otu_table(phy), 2, function(x){x/sum(x)})

hist(log(apply(otu, 1, sum)), breaks = 1000)

sums <- log(apply(otu, 1, sum))
n <- names(sums[sums > -10 & sums < -7])
n <- names(sums[sums > 0])

phy <- subset_taxa(phy, taxa_names(phy) %in% n)
phy <- subset_samples(phy, sample_sums(phy) > 0)

meta <- data.frame(sample_data(phy))

plot_bar(relphy(phy), fill = "Phylum") + 
           geom_bar(stat = "identity") +
           #scale_fill_brewer(palette = "Paired") +
           scale_x_discrete(limits = rownames(meta[order(meta$type),])) +
           theme_bw() +
           theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

ordgq <- ordinate(phy, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "species", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt









### Top most relatively abundant in ProjP vs. the others 

meta[,2]
colnames(meta)

projp <- subset_samples(phy1, batch %in% "projp")

proj <- subset_samples(phy1, !batch %in% "projp")











#### Barplot without 30 controls samples 

projptphy2 <- subset_taxa(projptphy, !taxa_names(projptphy) %in% tnames_controls)
projptphy2 <- subset_samples(projptphy2, sample_sums(projptphy2) > 0)

proj <- subset_taxa(phy1, !taxa_names(projptphy) %in% tnames_controls)
proj <- subset_samples(proj, sample_sums(projptphy2) > 0)


phy <- projptphy2

phy <- proj

ordgq <- ordinate(phy, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "batch", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt








#### Top most relative abundant taxa in each batch after control filtering 

max = 100

#projp
p <- sort(apply(apply(otu_table(projptphy2), 2, function(x){x/sum(x)}), 1, mean), decreasing = TRUE)[1:max]

#jay
j <- sort(apply(apply(otu_table(subset_samples(proj, batch %in% "jays")), 2, function(x){x/sum(x)}), 1, mean), decreasing = TRUE)[1:max]

#projk 
k <- sort(apply(apply(otu_table(subset_samples(proj, batch %in% "projk")), 2, function(x){x/sum(x)}), 1, mean), decreasing = TRUE)[1:max]

phy <- subset_taxa(phy1, taxa_names(phy1) %in% unique(names(c(p, j, k))))
phy <- subset_samples(phy, sample_sums(phy) > 0)


meta <- data.frame(sample_data(phy))

plot_bar(relphy(phy), fill = "Genus") + 
           geom_bar(stat = "identity") +
           facet_wrap(~batch, scale = "free_x") +
           scale_fill_brewer(palette = "Paired") +
        	theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
           
        

ordgq <- ordinate(phy, method = "PCoA", distance = "wunifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "batch", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt



























#### Abundance prevelance plots


phyglom <- tax_glom(newphy, "Genus")

sam <- sample_data(phy1)

phy2 <- transform_sample_counts(phy1, function(x){x/sum(x)})
sam <- sam[,colnames(sample_data(phy2)) %in% c("genus", "species")]
sam <- paste(sam$genus, sam$species, sep = "_")
tax <- as.character(tax_table(phy2)[,6])
wmean <- list()
for(i in unique(tax)){	
	samp <- c()
	for(j in unique(sam)){
		p <- subset_samples(phy2, paste(sample_data(phy2)$genus, sample_data(phy2)$species, sep = "_") %in% j)
		tsum <- sum(sample_sums(p))
		p2 <- subset_taxa(p, Genus %in% i)
		ttaxsum <- taxa_sums(p2)
		samp <- c(samp, ttaxsum/tsum)
	}
	wmean[[i]] <- samp
}
	
wm <- unlist(lapply(wmean, mean))
names(wm[wm > 0.01])

par(mfrow = c(1,2))

hist(wm, breaks = 100)
hist(apply(otu_table(relphy(phy2)), 1, mean), breaks = 100)

relwm <- apply(otu_table(relphy(phy2)), 1, mean)

prevx = as.numeric(apply(otu_table(phy2), 1, function(x){length(x[x > 0])}))

abund = as.numeric(apply(data.frame(otu_table(phy2)), 1, FUN = function(x){mean(x, na.rm = TRUE)}))

p <- data.frame(tax_table(phy2), prev = prevx, abund = as.numeric(as.character(abund)), stringsAsFactors = TRUE)

n <- length(sample_names(phy2))

k <- aggregate(p$prev, list(paste(p[,1], p[,2], p[,3], p[,4], p[,5], p[,6], sep = "_")), function(x){mean(x/n)})

k$abund <- aggregate(p$abund, list(paste(p[,1], p[,2], p[,3], p[,4], p[,5], p[,6], sep = "_")), function(x){max(x)})[,2]

k$Phylum <- as.character(lapply(str_split(k[,1], "__"), "[", 3))

k$Genus <- as.character(lapply(str_split(k[,1], "__"), "[", 7))

k2 <- k[log(k$abund) > 1 & k$x > 0.006,]

ggplot(data = k, aes(x = x, y = log(abund), col = Phylum)) +
	geom_point()


n <- taxa_sums(subset_samples(phy2, genus %in% "Kikihia"))
n2 <- taxa_sums(subset_samples(phy2, !Genus %in% "Kikihia"))

n3 <- data.frame(n = n, n2 = n2)

name <- rownames(n3[n3$n > 0 & n3$n2 == 0, ])

n4 <- subset_taxa(phy2, taxa_names(phy2) %in% name)


tx <- as.character(k2[,5])
tx <- paste("D_5__", tx, sep = "")

phy2 <- subset_taxa(phy1, Genus %in% tx)
phy2glom <- tax_glom(phy2, "Genus")

hist(apply(otu_table(phy2glom), 1, function(x){length(x[x > 0])}), breaks = 1000)
p <- apply(otu_table(phy2glom), 1, function(x){length(x[x > 0])})
p2 <- as.character(tax_table(subset_taxa(phy2glom, taxa_names(phy2glom) %in% names(p[p > 10])))[,6])

phy3 <- subset_taxa(phy2, Genus %in% p2)
phy3glom <- tax_glom(phy3, "Genus")

hist(taxa_sums(phy3glom), breaks = 100)
p <- names(taxa_sums(phy3glom)[taxa_sums(phy3glom) > 20000])
p2 <- as.character(tax_table(subset_taxa(phy3glom, taxa_names(phy3glom) %in% p))[,6])

phy4 <- subset_taxa(phy3, Family %in% p2)
phy4glom <- tax_glom(phy4, "Family")


phy <- phy4glom
ordgq <- ordinate(phy, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "batch", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt

str(ordgq)
vec <- cbind(row.names(ordgq$vectors), ordgq$vectors[,c(1,2)])
n <- vec[vec[,2] > 0.2 & ordgq$vectors[,3] < 0.2, 1]

















#### Controls 

taxa_sums(phy2)[taxa_sums(phy2) >= 0]

phycontrol <- subset_samples(phy2, type %in% "control")
sample_data(phycontrol)$species
phycontrol <- subset_samples(phy2, species %in% c("dneasy", "pcr", "powersoil"))
phycontrol <- subset_samples(phy2, species %in% c("dneasy", "pcr"))
phycontrol <- subset_samples(phy2, species %in% c("dneasy", "pcr", "powersoil"))

phycontrol <- merge_samples(phycontrol, "type", fun = sum)

n = 100
t <- apply(otu_table(phy2), 2, function(x){length(x[x > n & taxa_sums(phycontrol) > 0])})

t2 <- t[names(t) %in% sample_names(subset_samples(phy2, !type %in% "control"))]

hist(t, breaks = 100)
hist(t2, breaks = 100)


phy <- subset_samples(phy2, sample_names(phy2) %in% names(t2[t2 < 6]))


ordgq <- ordinate(phy, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "batch", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt













#### Microbes unique to each batch 

batch <- merge_samples(phy2, "batch", fun = sum)
phy <- batch

length(otu_table(phy)[1, otu_table(phy)[2,] == 0 & otu_table(phy)[3,] == 0])
length(otu_table(phy)[2, otu_table(phy)[1,] == 0 & otu_table(phy)[3,] == 0])
length(otu_table(phy)[3, otu_table(phy)[1,] == 0 & otu_table(phy)[2,] == 0])

jays <- taxa_names(otu_table(phy)[1, otu_table(phy)[2,] == 0 & otu_table(phy)[3,] == 0])
projk <- taxa_names(otu_table(phy)[2, otu_table(phy)[1,] == 0 & otu_table(phy)[3,] == 0])
projp <- taxa_names(otu_table(phy)[3, otu_table(phy)[1,] == 0 & otu_table(phy)[2,] == 0])

phy <- subset_taxa(phy2, !taxa_names(phy2) %in% c(jays, projk, projp))
phy <- subset_samples(phy, sample_sums(phy) > 0)


ordgq <- ordinate(phy, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "batch", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt



ordgq <- ordinate(phy2, method = "PCoA", distance = "wunifrac")
ord_filt <- plot_ordination(phy2, ordgq, color = "batch", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt

















#### Basic diverity metrics



phy <- phy2
phy <- subset_samples(phy2, type %in% "gut")

sim <- diversity(otu_table(phy), index = "simpson", MARGIN = 2)
shan <- diversity(otu_table(phy), index = "shannon", MARGIN = 2)
isim <- diversity(otu_table(phy), index = "invsimpson", MARGIN = 2)

even <- apply(otu_table(phy), 2, FUN = function(x){diversity(x)/log(specnumber(x))})

sums <- sample_sums(phy)


ggplot(data = data.frame(cbind(sim, shan, isim, even, sums)), aes(x = sums)) +
	geom_point(aes(y = sim, col = "blue")) +
	geom_point(aes(y = shan, col = "red")) +
	geom_point(aes(y = isim)) 













#### PCA post filtering 

phy <- phy2

ordgq <- ordinate(phy, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "type", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt










#### merging and ordinating by genus

phy2 <- merge_samples(phy1, "genus2", fun = sum)

sample_data(phy2)$genus2 <- row.names(sample_data(phy2))


ordgq <- ordinate(phy2, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy2, ordgq, label = "genus2", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt





#### Merging and ordinating by type 

phy2 <- merge_samples(phy1, "type", fun = sum)

sample_data(phy2)$type <- row.names(sample_data(phy2))


ordgq <- ordinate(phy2, method = "PCoA", distance = "unifrac")
ord_filt <- plot_ordination(phy2, ordgq, label = "type", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt





### Just Kikihia 

phy2 <- subset_samples(phy1, genus %in% "Kikihia")
phy2 <- merge_samples(phy2, "species", fun = sum)

sample_data(phy2)$species <- row.names(sample_data(phy2))


ordgq <- ordinate(phy2, method = "NMDS", distance = "unifrac")
ord_filt <- plot_ordination(phy2, ordgq, label = "species", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt








#### Just Maoricicada

phy2 <- subset_samples(phy1, genus %in% "Maoricicada")
phy2 <- merge_samples(phy2, "species", fun = sum)

sample_data(phy2)$species <- row.names(sample_data(phy2))


ordgq <- ordinate(phy2, method = "NMDS", distance = "unifrac")
ord_filt <- plot_ordination(phy2, ordgq, label = "species", axes = c(1,2)) +
  #geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  #geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt











#### Mapping samples 

library(rgdal)
library(maptools)
sample_data(phylo)[sample_data(phylo)$lat_sign == "neg" & is.finite(sample_data(phylo)$lat), colnames(sample_data(phylo)) == "lat"] <- -(sample_data(phylo)[sample_data(phylo)$lat_sign == "neg" & is.finite(sample_data(phylo)$lat), colnames(sample_data(phylo)) == "lat"])
sample_data(phylo)[sample_data(phylo)$lon_sign == "neg" & is.finite(sample_data(phylo)$lon), colnames(sample_data(phylo)) == "lon"] <- -(sample_data(phylo)[sample_data(phylo)$lon_sign == "neg" & is.finite(sample_data(phylo)$lon), colnames(sample_data(phylo)) == "lon"])


lnd <- readOGR(dsn = "nz-coast/", "nz-coastlines-topo-150k")
elev <- readOGR(dsn = "nz-coast/", "nz-height-points-topo-150k")

sdata <- sample_data(phylo)[sample_data(phylo)$island %in% c("South Island", "North Island", "Chatham Islands"),]

coord <- cbind(sdata$lon, sdata$lat)
coord <- coord[is.finite(coord[,1]),]
pts <- SpatialPoints(coord, proj4string = CRS(proj4string(lnd)))


points(pts)


library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggspatial)
library(sf)
library(raster)
  
plot(lnd)
plot(elev, add = TRUE)
plot(pts, add = TRUE, col = "white")
  
da <- data.frame(elev)[, c(3,4,2)]
colnames(da) <- c("x", "y", "z")
rasterFromXYZ(da)


e <- extent(elev[,1:2])
e <- e + 1000

r <- raster(e, ncol=10, nrow=2)

x <- rasterize(da[, 1:2], r, da[,3], fun=mean)

plot(x)

over(pts, x)
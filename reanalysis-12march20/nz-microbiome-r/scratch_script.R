---
  title: "NZM"
output: html_notebook
---

#### Importing data from Qiime2 data artifacts ####
  

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R") 
#install.packages("ggmsa")
#install.packages("seqmagick")
library(qiime2R)
library(phyloseq)
library(stringr)
library(seqinr)
library(ggplot2)

# Loading packages. 


tax_gza_to_phyloseq <- function(path_tax_gza) {
  taxonomy <- read_qza(path_tax_gza)
  taxonomy <- data.frame(taxonomy$data)
  kingdom <- sapply(stringr::str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 1), "__"), "[", 2)
  phylum <- sapply(stringr::str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 2), "__"), "[", 2)
  class <- sapply(stringr::str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 3), "__"), "[", 2)
  order <- sapply(stringr::str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 4), "__"), "[", 2)
  family <- sapply(stringr::str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 5), "__"), "[", 2)
  genus <- sapply(stringr::str_split(sapply(str_split(as.character(taxonomy$Taxon),";D_"), "[", 6), "__"), "[", 2)
  feature_id <- taxonomy$Feature.ID
  confidence <- taxonomy$Confidence
  taxonomy <- data.frame(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Confidence = confidence)
  rownames(taxonomy) <- feature_id
  return(taxonomy) 
}

b1_table <- otu_table(read_qza("../jk_table.qza")$data, taxa_are_rows = TRUE) 
b2_table <- otu_table(read_qza("../jt_table.qza")$data, taxa_are_rows = TRUE) 
b3_table <- otu_table(read_qza("../projk_table.qza")$data, taxa_are_rows = TRUE) 
b4_table <- otu_table(read_qza("../projp_table.qza")$data, taxa_are_rows = TRUE) 

b1_classification <- tax_table(as.matrix(tax_gza_to_phyloseq("../jk_classified.qza")))
b2_classification <- tax_table(as.matrix(tax_gza_to_phyloseq("../jt_classified.qza")))
b3_classification <- tax_table(as.matrix(tax_gza_to_phyloseq("../projk_classified.qza")))
b4_classification <- tax_table(as.matrix(tax_gza_to_phyloseq("../projp_classified.qza")))

b1_phylogeny <- phy_tree(read_qza("../jk_tree/rooted_tree.qza")$data)
b2_phylogeny <- phy_tree(read_qza("../jt_tree/rooted_tree.qza")$data)
b3_phylogeny <- phy_tree(read_qza("../projk_tree/rooted_tree.qza")$data)
b4_phylogeny <- phy_tree(read_qza("../projp_tree/rooted_tree.qza")$data)

metadata <- read.csv("../metadata.csv", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$SampleID

b1_metadata <- sample_data(metadata[metadata$batch == "B1", ])
b2_metadata <- sample_data(metadata[metadata$batch == "B2", ])
b3_metadata <- sample_data(metadata[metadata$batch == "B3", ])
b4_metadata <- sample_data(metadata[metadata$batch == "B4", ])

b1 <- phyloseq(b1_table, b1_classification, b1_phylogeny, b1_metadata)
b2 <- phyloseq(b2_table, b2_classification, b2_phylogeny, b2_metadata)
b3 <- phyloseq(b3_table, b3_classification, b3_phylogeny, b3_metadata)
b4 <- phyloseq(b4_table, b4_classification, b4_phylogeny, b4_metadata)

b1
b2
b3
b4

# A function is defined to reformat taxonomic classifications for compatability with phyloseq package. ASV tables, taxonomys, phylogenys, and metadata are imported and then combined into phyloseq objects associated with each dataset.  



#### Non-target features ####

for(i in c("b1", "b2", "b3", "b4")){
  rep_seqs <- read.fasta(paste("../", i, "_seqs.fasta", sep = ""), as.string = "TRUE", forceDNAtolower = FALSE)
  arch <- taxa_names(subset_taxa(get(i), Kingdom %in% c("Archaea")))
  chlor <- taxa_names(subset_taxa(get(i), Order %in% c("Chloroplast")))
  mito <- taxa_names(subset_taxa(get(i), Family %in% c("Mitochondria")))
  hodg <- taxa_names(subset_taxa(get(i), Genus %in% c("Candidatus Sulcia")))
  sulc <- taxa_names(subset_taxa(get(i), Genus %in% c("Candidatus Hodgkinia")))
  for(j in c('arch', 'chlor', 'mito', 'hodg', 'sulc')){
    rep_seqs_subset <- rep_seqs[names(rep_seqs) %in% get(j)]
    names(rep_seqs_subset) <- paste(i, names(rep_seqs_subset), sep = "_")
    write.fasta(sequences = rep_seqs_subset, names = names(rep_seqs_subset), file.out = paste(paste(i, j, sep = "_"), ".fasta", sep = ""))
    rm(j)
  }
  rm(rep_seqs)
}


sulc <- Biostrings::readDNAStringSet("b4_mito.fasta")
sulc_aligned <- msa::msaClustalW(sulc)
Biostrings::writeXStringSet(sulc_aligned@unmasked, "test.fasta")
t <- phangorn::phyDat(ape::read.dna("test.fasta", format = "fasta"))
tm <- phangorn::dist.ml(t, model="F81")
tm2 <- phangorn::NJ(tm)
plot(tm2)
ggmsa::ggmsa("test.fasta", font = NULL, color = "Chemistry_NT")

print(sulc_aligned)

# Sequences of ASVs classified to each endosymbiont (either Hodgkinia or Sulcia), mitochondria, chloroplast, and archaea are exported into respective fasta files. Mitochondrial ASVs are then aligned to a reference set of insect and Ophiocordyceps (symbiont and pathogen) 16s sequences extracted from published datasets. This alignment is used to construct a neighbor joining tree showing that high-abundance mitochondrial ASVs group with known Ophiocordyceps symbionts of cicadas. 


#### Relative abundance of all classifications ####

phy <- b1

a <- taxa_names(subset_taxa(phy, Genus %in% c("Candidatus Sulcia")))
b <- taxa_names(subset_taxa(phy, Genus %in% c("Candidatus Hodgkinia")))
c <- taxa_names(subset_taxa(phy, Family %in% c("Mitochondria") & !taxa_names(phy) %in% taxa_names(mito_ophio)))
d <- taxa_names(subset_taxa(phy, Order %in% c("Chloroplast")))
e <- taxa_names(subset_taxa(phy, Kingdom %in% c("Archaea")))
f <- taxa_names(subset_taxa(phy, !taxa_names(phy) %in% c(a,b,c,d,e) & Kingdom %in% "Bacteria"))
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
bac_rel <- subset_taxa(phy_first_rel, Genus %in% "Bacteria")
sample_data(phy_first_rel)$bac_rel <- sample_sums(bac_rel)
sample_data(phy_first_rel)$depth <- sample_sums(phylo)

plot_bar(phy_first_rel, fill = "Genus") + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rownames(sample_data(phy_first)[order(sample_data(phy_first)$type),])) +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))


phy_b1 <- subset_samples(phy_first_rel, batch %in% "B1")
b1_plot <- plot_bar(phy_b1, , fill = "Genus") + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), 
                   labels = paste(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree), "species.tree"]$species.tree, rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), sep = " ")
  ) +
  xlab("") +		
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.08, size = 8))

phy_b1 <- subset_samples(phy_first_rel, batch %in% "B2")
b2_plot <- plot_bar(phy_b1, , fill = "Genus") + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), 
                   labels = paste(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree), "species.tree"]$species.tree, rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), sep = " ")
  ) +
  xlab("") +		
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.08, size = 8))

phy_b1 <- subset_samples(phy_first_rel, batch %in% "B3")
b3_plot <- plot_bar(phy_b1, , fill = "Genus") + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), 
                   labels = paste(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree), "species.tree"]$species.tree, rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), sep = " ")
  ) +
  xlab("") +		
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.08, size = 8))

phy_b1 <- subset_samples(phy_first_rel, batch %in% "B4")
b4_plot <- plot_bar(phy_b1, , fill = "Genus") + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), 
                   labels = paste(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree), "species.tree"]$species.tree, rownames(sample_data(phy_b1)[order(sample_data(phy_b1)$type, sample_data(phy_b1)$genus, sample_data(phy_b1)$species.tree),]), sep = " ")
  ) +
  xlab("") +		
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.08, size = 8))


gar <- ggarrange(b1_plot, b2_plot, b3_plot, b4_plot, common.legend = TRUE)

ggsave("ggarrange.png", gar)

ggsave("b1.png", b1_plot)
ggsave("b2.png", b2_plot)
ggsave("b3.png", b3_plot)
ggsave("b4.png", b4_plot)


limits <- rownames(sample_data(phy_first)[order(sample_data(phy_first)$batch, sample_data(phy_first)$type, sample_data(phy_first)$genus, sample_data(phy_first)$species),])
labels <- paste(sample_data(phy_first)[order(sample_data(phy_first)$batch, sample_data(phy_first)$type, sample_data(phy_first)$genus, sample_data(phy_first)$species), "species"]$species, limits, sep = " ")

otu <- data.frame(otu_table(phy_first_rel), check.rows = FALSE,
                  check.names = FALSE, fix.empty.names = FALSE,
                  stringsAsFactors = FALSE)
otu$genus <- as.character(tax_table(phy_first_rel)[,6])
otu <- gather(otu, key = sample , value = abund, 1:(dim(otu)[2]-1))

label_data= data.frame( 
  xpos = seq(1, length(labels)),
  ID = labels,
  ypos = 1
)

number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$xpos-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle) 

batch <- sample_data(phy_first)[order(sample_data(phy_first)$batch, sample_data(phy_first)$type, sample_data(phy_first)$genus, sample_data(phy_first)$species), "batch"]$batch
type <- sample_data(phy_first)[order(sample_data(phy_first)$batch, sample_data(phy_first)$type, sample_data(phy_first)$genus, sample_data(phy_first)$species), "type"]$type


g <- ggplot(data = otu, aes(x = sample, y = abund, fill = genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = limits, 
                   labels = labels
  ) +
  xlab("") +	
  ylim(-1,ypos+0.04) +
  coord_polar(start = 0) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),  plot.margin = unit(rep(0.6,4), "cm"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_segment(aes(x = 1, y = -0.25, xend = 40, yend = -0.25), colour = "grey80", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_text(aes(x = 1+20, y = -0.45, label="B1"), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  geom_segment(aes(x = 41, y = -0.25, xend = 76, yend = -0.25), colour = "grey46", alpha=0.8, size=4, inherit.aes = FALSE )  +
  geom_text(aes(x = 41+17, y = -0.45, label="B2"), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  geom_segment(aes(x = 77, y = -0.25, xend = 145, yend = -0.25), colour = "grey80", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_text(aes(x = 77+34, y = -0.45, label="B3"), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  geom_segment(aes(x = 146, y = -0.25, xend = 292, yend = -0.25), colour = "grey46", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_text(aes(x = 146+73, y = -0.45, label="B4"), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  
  geom_segment(aes(x = 1, y = -0.1, xend = 76, yend = -0.1), colour = "grey80", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_segment(aes(x = 79, y = -0.1, xend = 131, yend = -0.1), colour = "grey80", alpha=0.8, size=4, inherit.aes = FALSE )  +
  geom_segment(aes(x = 183, y = -0.1, xend = 292, yend = -0.1), colour = "grey80", alpha=0.8, size=4, inherit.aes = FALSE )  +
  
  geom_segment(aes(x = 132, y = -0.1, xend = 145, yend = -0.1), colour = "grey46", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_segment(aes(x = 146, y = -0.1, xend = 173, yend = -0.1), colour = "grey80", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_segment(aes(x = 77, y = -0.1, xend = 78, yend = -0.1), colour = "grey46", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_segment(aes(x = 174, y = -0.1, xend = 182, yend = -0.1), colour = "grey46", alpha=0.8, size=4 , inherit.aes = FALSE )  +
  geom_text(data=label_data, aes(x=xpos, y=ypos+0.04, label=ID, hjust=hjust), color="black", size=2, angle= label_data$angle, inherit.aes = FALSE )

# Plots of relative abundance of bacteria, each endosymbiont, chloroplast, and archeae across samples and datasets.


#### Exploring control samples ####

table(sample_data(b1)$type)
table(sample_data(b2)$type)
table(sample_data(b3)$type)
table(sample_data(b4)$type)

# Each row represents a count of of the number of samples in each datset for the associated tissue type. The first two datasets are gut-only datasets. However, dataset B2 was generated from homogenized gut and reproductive tissue. The B3 and B4 datasets include increased tissue sampling and controls, with 4 controls in dataset B1 and 31 controls in dataset B4.  


b3_controls <- subset_samples(b3, type %in% "control")
b3_controls <- subset_taxa(b3_controls, taxa_sums(b3_controls) > 0)
b4_controls <- subset_samples(b4, type %in% "control")
b4_controls <- subset_taxa(b4_controls, taxa_sums(b3_controls) > 0)



psmelt(b3_controls)


plot_bar(b3_controls)



```






# Filtering feature tables 

```{r}

```


# Removing contaminants with Decontam



```{r}

qpcr <- read.csv("../qpcr.csv")
qpcr_mean <- aggregate(qpcr$Cq, by = list(qpcr$Target), mean)

nzm_filtered_decontam_plots <- list()
nzm_filtered_decontam_genera <- list()

decontam_threshold <- 0.5

```
Importing qPCR data and 

```{r}

# B1 dataset
index <- 1
input <- nzm_filtered[[index]]
sample_data(input)$qpcr <- qpcr_mean[,2][match(sample_data(input)$ID, qpcr_mean[,1])]
input <- subset_samples(input,qpcr > 0)
freq <- decontam::isContaminant(input, method="frequency", conc="qpcr", threshold = 0.5)

contaminants <- subset_taxa(input, taxa_names(input) %in% taxa_names(input)[c(which(freq$contaminant))])
nzm_filtered[[index]] <- subset_taxa(input, !taxa_names(input) %in% taxa_names(input)[c(which(freq$contaminant))])
nzm_filtered_decontam_plots[[index]] <- plot_frequency(input, taxa_names(input)[c(which(freq$contaminant))], conc="qpcr")
nzm_filtered_decontam_genera[[index]] <- dtable(contaminants)$Genus

# B3 dataset
index <- 3
input <- nzm_filtered[[index]]
sample_data(input)$qpcr <- qpcr_mean[,2][match(sample_data(input)$ID, qpcr_mean[,1])]
input <- subset_samples(input,qpcr > 0)
freq <- decontam::isContaminant(input, method="frequency", conc="qpcr", threshold = 0.3)
contaminants <- subset_taxa(input, taxa_names(input) %in% taxa_names(input)[c(which(freq$contaminant))])
nzm_filtered[[index]] <- subset_taxa(input, !taxa_names(input) %in% taxa_names(input)[c(which(freq$contaminant))])
nzm_filtered_decontam_plots[[index]] <- plot_frequency(input, taxa_names(input)[c(which(freq$contaminant))], conc="qpcr")
nzm_filtered_decontam_genera[[index]] <- dtable(contaminants)$Genus





# B4 dataset
index <- 4
input <- nzm_filtered[[index]]
sample_data(input)$qpcr <- qpcr_mean[,2][match(sample_data(input)$ID, qpcr_mean[,1])]
input <- subset_samples(input,qpcr > 0)
freq <- decontam::isContaminant(input, method="frequency", conc="qpcr", threshold = 0.5)
contaminants <- subset_taxa(input, taxa_names(input) %in% taxa_names(input)[c(which(freq$contaminant))])
nzm_filtered[[index]] <- subset_taxa(input, !taxa_names(input) %in% taxa_names(input)[c(which(freq$contaminant))])
nzm_filtered_decontam_plots[[index]] <- plot_frequency(input, taxa_names(input)[c(which(freq$contaminant))], conc="qpcr")
nzm_filtered_decontam_genera[[index]] <- dtable(contaminants)$Genus






```

#### Relative abundance stuff ####
major_table_new <- setNames(aggregate(major_table$abundance, by = list(major_table$sample, major_table$Genus), sum), c("Sample", "Genus", "Absolute_Abundance"))
major_table_new <- spread(major_table_new, key = Sample, value = Absolute_Abundance)
rel_abund <- data.frame(apply(major_table_new[,-1], 2, function(x){ x / sum(x)}))
rel_abund <- cbind(Genus = major_table_new$Genus, rel_abund)
major_table_new <- gather(rel_abund, key = "Sample", value = "Relative_Abundance", 2:dim(rel_abund)[2])
major_table_new$Sample <- sub("[.]", "-", major_table_new$Sample)
major_table_new$Sample <- sub("X", "", major_table_new$Sample)







##### Permanova stuff #####

dis <- distance(phy, method = j)
permanova <- adonis2(dis ~ genus + species + depth + island + elevation + type, data = ordDF, permutations = 1000, method = j)

genus_sig <- paste("Genus", paste(paste("R2", round(permanova$R2, digits = 3)[1], sep = ":"), paste("Significance", round(permanova$`Pr(>F)`[1], digits = 3), sep = ":"), sep = ", "), sep = " = ") 

species_sig <- paste("Species", paste(paste("R2", round(permanova$R2, digits = 3)[2], sep = ":"), paste("Significance", round(permanova$`Pr(>F)`[2], digits = 3), sep = ":"), sep = ", "), sep = " = ")

depth_sig <- paste("Species", paste(paste("R2", round(permanova$R2, digits = 3)[3], sep = ":"), paste("Significance", round(permanova$`Pr(>F)`[3], digits = 3), sep = ":"), sep = ", "), sep = " = ")

island_sig <- paste("Genus", paste(paste("R2", round(permanova$R2, digits = 3)[4], sep = ":"), paste("Significance", round(permanova$`Pr(>F)`[4], digits = 3), sep = ":"), sep = ", "), sep = " = ") 

elevation_sig <- paste("Species", paste(paste("R2", round(permanova$R2, digits = 3)[5], sep = ":"), paste("Significance", round(permanova$`Pr(>F)`[5], digits = 3), sep = ":"), sep = ", "), sep = " = ")

type_sig <- paste("Species", paste(paste("R2", round(permanova$R2, digits = 3)[6], sep = ":"), paste("Significance", round(permanova$`Pr(>F)`[6], digits = 3), sep = ":"), sep = ", "), sep = " = ")

sig[[j]][[i]] <- paste(genus_sig, species_sig, depth_sig, island_sig, elevation_sig, type_sig, sep = "\n")







##### Combinde library size plot #####

```{r}
library(grid)

df_plot <- data.frame()
for(i in 1:3){
  x <- nzm_filtered2[[i]]
  df <- as.data.frame(sample_data(x)) 
  df$LibrarySize <- sample_sums(x)
  df <- df[!is.na(df$conc) & !is.na(df$LibrarySize), ]
  df <- df[df$type %in% c("bacteriome", "control", "egg", "gut", "reproductive"), ]
  df_plot <- rbind(df_plot, data.frame(df))
}

libsize_combined <- ggplot(data = df_plot, aes(x = log(conc+1), y = log(LibrarySize), col = type, shape = batch)) +
  geom_point(alpha = 1, size = 2) +
  labs(x = "log(Post-PCR DNA Concentration (ng/ul))", y = "log(Total ASV Abundance)", col = "", shape = "") +
  scale_color_manual(values = tissue_type_colors) +
  scale_shape_manual(values = dataset_shapes) +
  guides(shape = FALSE) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    #legend.position = c(1,0),
    legend.position = "top",
    legend.justification = c("left"), 
    legend.background = element_blank(),
    #legend.box.background = element_rect(colour = "grey50"),
    legend.spacing.x = unit(0, 'mm'),
    legend.spacing.y = unit(0, 'mm'),
    legend.key.height = unit(1, "mm"),
    legend.key.width = unit(4, "mm"),
    legend.text=element_text(size=6), 
    legend.margin = margin(0, 0, 0, 0))


libsize_combined

```


##### Maping and elevation #####


library(sf)
library(units)
library(tmap)
library(tmaptools)
library(OpenStreetMap)
library(raster)
library(ggmap)

metadata <- rbind(b1_metadata, b2_metadata, b3_metadata)
metadata <- metadata[!is.na(metadata$lon) & !is.na(metadata$lat), ]
metadata$code2 <- paste(stringsplit(metadata$code, "[.]", 3), stringsplit(metadata$code, "[.]", 4), sep = ".")
metadata <- data.frame(unique(paste(metadata$lat_sign, metadata$lon_sign, metadata$code2)))
metadata$lat_sign = stringsplit(metadata[,1], " ", 1)
metadata$lon_sign = stringsplit(metadata[,1], " ", 2)
metadata$code2 = stringsplit(metadata[,1], " ", 3)

contours = st_read(dsn = "../lds-nz-contours-topo-1500k-SHP", layer = 'nz-contours-topo-1500k')
data("World")
World <-  st_transform(World, CRS("+init=epsg:4326") )
nz <- World[World$sovereignt == "New Zealand", ]

bbox <- st_bbox(nz)
names(bbox) <- c("left", "bottom", "right", "top")
bbox <- bbox + c(-1,-1,1,1)
nz_map <- get_stamenmap(bbox = bbox, zoom = 7, maptype = "terrain-background")

ggmap(nz_map) +
  theme_new()


tm_basemap("Staman") +
  tm_shape(nz) +
  tm_borders() +
  tm_tiles("Staman")


contours_lonlat <- st_transform(contours, crs = 4326)

nz_lonlat <- st_transform(nz, crs = 4326)

grid <- st_make_grid(nz_lonlat, n = c(100, 100))
grid <- st_bind_cols(grid, cell = seq(1, length(grid), 1))
contours_within_grid <- st_within(contours_lonlat, grid)

contours_lonlat2 <- cbind(contours_lonlat, cell = as.numeric(contours_within_grid))
contours_lonlat2 <- aggregate(contours_lonlat2["elevation"], by = list(as.numeric(contours_within_grid)), mean)

elevation <- contours_lonlat2$elevation[match(grid$cell, contours_lonlat2$Group.1)]
elevation[is.na(elevation)] <- 0

grid2 <- st_bind_cols(grid, "Elevation" = elevation)

points <- st_as_sf(metadata, coords = c("lon_sign", "lat_sign"), crs = 4326) 

map <- tm_shape(grid2) +
  tm_polygons(col = "Elevation", border.alpha = 0, alpha = 0.5) +
  #tm_shape(nz_lonlat) +
  #  tm_borders() +
  tm_shape(points) +
  tm_dots() +
  tm_text(text = "code2", size = 0.2, auto.placement = TRUE)

tmap_save(map, "map.png", width = 4, height = 5, units = "in")
map




phylo <- qza_to_phyloseq(features="nzbiome-dada-table.qza", tree="nzbiome-dada-filtered-alignment-rooted-tree.qza", taxonomy="nzbiome-dada-rep-seqs-taxonomy.qza", metadata="projP-final-metadata.tsv")

phylo <- subset_samples(phylo, !is.na(type) 
	& !type %in% c("uncertain", "Unknown9","bacteriome")
	& !ID %in% c("C1", "C2", "T3", "T5", "T7", "T23")
	)

phylo <- prune_taxa(taxa_sums(phylo) > 0, phylo)
for(i in 1:dim(tax_table(phylo))[2]) {
	tax_table(phylo)[,i] <- unlist(lapply(str_split(tax_table(phylo)[,i], "__"), "[", 	2))
}


##### General data characteristics #########

mito_ophio <- prune_taxa(taxa_names(phylo) %in% c(
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
"0d8331a0106b98f8682c7c8627831afa"), phylo)





## Clustered ophio ASVs 

mito_ophio

cluster1 <- c(
"e1b3e502889eda99faec882decc61084",
"9fa8c30cf03d22a90f75fcc93cb22b1f",
"a597872fcf87b49f099e6525ce8e0bdd",
"5eb79a1bc542e2bf20251039b5a3892b",
"a10c98cca957ba6b03dc4129da14e75a",
"326fda0ea300b8bc64e930f5112f15e9",
"3b4a02054a191da163cdbb3402ab1adb",
"b66f62979ed505d71a9050f44d03eafc",
"0986e084d515a7c589486d57f012e290",
"f16746b4c811ba0ddbcc063c94858ae4",
"5abaea3f41fe71f043900d5106be7165"
)


cluster2 <- c(
"b6b811e7b71729058ad58cf4173a8f38",
"bd7c33818f8aec75f7161d4e7185badb",
"dd69e8131ec9521228e03aa69b7685cc",
"5f3a1e0d63fa00f3ecd7c4a43759576e",
"f347ba578d6e5fc132aafadf7b467e41",
"f87631331b247a657230e771a3e67f99",
"427f6aaa7865aeb2611e969886dd4b4c",
'9db6287327eae687af1b9ca9b89b7894',
"b8b9f752ace6d5c1558c53881c1a3c24",
"532042554bb7b83f4f87972af256d7c6",
"790ccb092610629a233eba365be48fa7",
"0a2757312143a76b0e38b780ba067875",
"2e1081dbfe0c6c8ffdb56c1167c72bf6",
"c0cce88070e16e87cb1e8a0410ded5db",
"646f5e07f65b733c72cace35a704f30d",
'6d925e800f29f71a5a1ddaf29b641dc1',
"ab72d8d94a01bf5f6a8990c0686a5d04",
"f440d0c74ede6a04004b11fd859b1ac0",
"0d8331a0106b98f8682c7c8627831afa",
"adbb584c88e2f713dc2e6215ef878a02",
"3ad79c231b558392426ed778987cfd28",
"cf56e57f332051fbc8076be97342c1d1",
"2f19bacb43636b5d04a5af15efe85bcd",
"59e8f6fedbcceff4ca6c982cc46d1562",
"c31ec8d17f0a762dcf5094edae97c306",
"7198d35d03014e52e939411f719a695b",
"fd3fa68f678cb377e70ae99a348b458b",
"7fd0409aa73324573665cedccca947f7",
"84419617fd969dbddba9901ff4ed16d9",
"7c15c08a32970d2123f81c79a841dcda",
"fc1168706670fa562f282a2cc64db9d1"
)

cluster3 <- c(
"b4cc0fe078d755f6dc41dd4bdb4ec7f0",
"dd4a893a3874a25e1e1f9cdca607eb8c",
"755666aedf0d0425062f31bb42c1fc8a",
"8b7826d9e54b2dfc747d851087520271",
"d95af54ccf6d283644fbb2ebf448b3fe",
"a75e5c8106e8135d50e2ba9327cefd73",
'bd637283f4c956d6eb3f60aba645062d',
"f49f87ee110aaf95f83f1e476d3a32b7",
"ab2383c7be692b1f25d72487da6df733"
)

cluster4 <- c(
"776cb3ae1cc229b9999dbd9c7b800fc4",
"8a31918806e315293752775581aa6227",
"71efa5bd59720222774a0e29a2963cc6",
"40ef5a4fd881601fdaaf9f35a658325f",
"d7526f29c9eb1a0c8d1ab3748edab74f"
)

clusters <- list(cluster1, cluster2, cluster3, cluster4)

new_otu <- c()
for(i in 1:4) { 
	a <- otu_table(prune_taxa(clusters[[i]], mito_ophio))
	b <- apply(a, 2, sum)
	new_otu <- rbind(new_otu, b)
}
rownames(new_otu) <- c("cluster1", "cluster2", "cluster3", "cluster4")
mito_ophio2 <- phyloseq(otu_table(new_otu, taxa_are_rows = TRUE), sample_data(mito_ophio))


library(phytools)
library(ggtree)
library(pheatmap)
library(ggstance)

htree <- read.tree("raxml.tre")

htree <- reroot(htree, node = 70)

t <- htree

t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist( lapply( str_split(t$tip.lab, "_"), "[", 1) )

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))

phy <- mito_ophio2
phy <- subset_samples(phy, sample_sums(phy) > 1000 & !type %in% "control" & !batch %in% c("B1", "B2"))
#phy <- subset_samples(phy, !is.na(species.tree))
t2 <- drop.tip(t, t$tip.label[!t$tip.label %in% sample_data(phy)$species.tree])
#phy <- merge_samples(phy, "species")


g <- (ggtree(t2) + geom_tiplab(size = 2)) %<+% data.frame(t(otu_table(phy))) +
	#geom_tiplab(align = TRUE, size = 3) +
	geom_nodelab( nudge_x = 0.005, nudge_y = 0.001, size = 2.5) + 
	theme(legend.position = "none", legend.title = element_blank(), legend.key = element_blank())

ggsave("host_mito_tree.png", g)

gdata <- data.frame(g$data)
gdata <- gdata[order(gdata$y, decreasing = TRUE), ]
tiporder <- gdata[is.na(as.numeric(gdata$label)), "label"]
tiporder <- tiporder[tiporder != ""]

sd <- sample_data(phy)
#sd[sd$species.tree == "muta-NI", "species.tree"] <- "muta-SI"
sd$order <- match(sd$species.tree, tiporder)
limits = sd[order(sd$genus, sd$species, sd$species.tree, decreasing = TRUE), "ID"]$ID
labels = paste(sd[order(sd$genus, sd$species, sd$species.tree, decreasing = TRUE), "species.tree"]$species.tree, sd[order(sd$genus, sd$species, sd$species.tree,  decreasing = TRUE), "code"]$code, sep = " ; ")

otu <- data.frame(otu_table(phy))
otu$taxa <- rownames(otu)
otu <- gather(otu, key = species, value = abund, 1:(dim(otu)[2]-1))


g <- ggplot(data = otu, aes(x = species, y = sqrt(abund), fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = limits, labels = labels) +
    xlab("") +
    ylab("") +
    labs(fill = "Bacterial Family") +
    coord_flip() +
    theme_minimal() +
    theme(plot.margin = unit(c(1,10,0.5,1), "cm"), 
    	legend.position = "right", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 3.5, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 3.5)) 

ggsave("fungal_bar.png", g)




library(adegenet)

snps <- fasta2DNAbin("mito_seqs.fasta", snpOnly = T)
pca <- prcomp(dist.dna(snps, model = "N"), scale=T)
pca <- data.frame(pca$x)
pca$cluster <- rep(0, dim(pca)[1],1)
pca$asv <- rownames(pca)

pca[!is.na(match(rownames(pca), cluster1)), "cluster"] <- "cluster1"
pca[!is.na(match(rownames(pca), cluster2)), "cluster"] <- "cluster2"
pca[!is.na(match(rownames(pca), cluster3)), "cluster"] <- "cluster3"
pca[!is.na(match(rownames(pca), cluster4)), "cluster"] <- "cluster4"

asv_abund <- c()
for(i in 1:dim(pca)[1]) {
	a <- rownames(pca)[i]
	b <- prune_taxa(a, mito_ophio)
	c <- sum(otu_table(b))
	asv_abund <- c(asv_abund, c)
	}

pca$abund <- as.numeric(asv_abund)

g <- ggplot(pca, aes(x = PC1, y = PC2, col = cluster)) +
	geom_point(aes(size = sqrt(abund))) +
	scale_color_manual(values = col_vector) +
	theme_bw() +
    theme(
    	legend.position = "right", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 5)) 

ggsave("fungal_cluster.png", g, scale = 0.7)















phy <- phylo

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

ggsave("overview.png", g, scale = 2)    



phy_first_rel
apply(otu_table(phy_first), 1, function(x) length(x[x > 100]))/dim(otu_table(phy_first))[2]*100
tax_table(phy_first)

sum(sample_sums(subset_samples(phy_first, type %in% "control")))
sum(sample_sums(subset_taxa(subset_samples(phy_first, type %in% "control"), Genus %in% "Bacteria")))
142925/165095
sum(sample_sums(subset_taxa(subset_samples(phy_first, type %in% "control"), Genus %in% "Archaea")))
485/165095
sum(sample_sums(subset_taxa(subset_samples(phy_first, type %in% "control"), Genus %in% "Ophiocordyceps")))
845/165095

sum(sample_sums(subset_taxa(subset_samples(phy_first, type %in% "control"), Genus %in% "Candidatus Hodgkinia")))
0/165095
sum(sample_sums(subset_taxa(subset_samples(phy_first, type %in% "control"), Genus %in% "Candidatus Sulcia")))
1874/165095
sum(sample_sums(subset_taxa(subset_samples(phy_first, type %in% "control"), Genus %in% "Chloroplast")))
18019/165095


d <- data.frame(s = sample_sums(subset_taxa(subset_samples(phy_first, !type %in% "control"), Genus %in% "Ophiocordyceps")))
d$oph <- d$s > 10
sample_data(phy_first)[sample_data(phy_first)$ID %in% rownames(d[d$oph == FALSE,]), "species"]






################## Filtering ##############

phylo	#6676

phy1 <- subset_taxa(phylo, !Genus %in% c("Candidatus Sulcia"))	#125
phy1 <- subset_taxa(phy1, !Kingdom %in% c("Eukaryota"))	#2
phy1 <- subset_taxa(phy1, !Kingdom %in% c("Archaea"))	#167
phy1 <- subset_taxa(phy1, !Family %in% c("Mitochondria"))	#145
phy1 <- subset_taxa(phy1, !Order %in% c("Chloroplast"))	#47
phy1 <- subset_taxa(phy1, !Genus %in% c("Candidatus Hodgkinia"))	#47
phy1 <- subset_taxa(phy1, !is.na(Kingdom))	#1595
phy1 <- subset_taxa(phy1, !is.na(Phylum))	#763
phy1 <- subset_samples(phy1, sample_sums(phy1) > 500) #259 samples (from 292)
phy1 <- subset_taxa(phy1, taxa_sums(phy1) > 0)	#32

phy1 #3753

phy1 <- tip_glom(phy1, h = 0.03) #1865
phy1_topglom_backup <- phy1

phy1 #1888

phy1 <- subset_samples(phy1, !type %in% "egg" )

phy1_genus <- tax_glom(phy1, taxrank = "Genus") #657 genera 




#### DECONTAM 
library(decontam)


## Prev based for B4
## We performed a series of Decontam (R package) filtering steps. First, we used Decontam's prevalence based filtering to remove bacterial taxa that were prevalent across controls (threshold = 0.6) 
## Based on the plot, increase prevelence of bacteria acroos gut samples showed a qualitative increase with prevalence across control samples, suggesting that contaminants ar widespread in gut samples. 
 

phy <- subset_samples(phy1, batch %in% c("B4"))
phy <- subset_taxa(phy, taxa_sums(phy) > 0) #1102

sample_data(phy)$is.neg <- sample_data(phy)$type == "control"
b4_prev <- isContaminant(phy, method="prevalence", neg="is.neg", threshold=0.68)
table(b4_prev$contaminant)
ps.pa <- transform_sample_counts(phy, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$type == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$type == "gut", ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=b4_prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(alpha = 0.5) +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

decontam1 <- rownames(b4_prev[which(b4_prev$contaminant),])

unique(tax_table(prune_taxa(taxa_names(phy1) %in% decontam1, phy1))[,6])

phy1_decontam <- prune_taxa(!taxa_names(phy1) %in% decontam1, phy1) #-36







## Freq-based decontam for B1-B2

phy <- subset_samples(phy1, batch %in% c("B1","B2") & conc > 0)
phy <- subset_samples(phy, conc > 0)
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
b1b2_freq <- isContaminant(phy, method="frequency", conc="conc", threshold = 0.15)
table(b1b2_freq$contaminant)
hist(b1b2_freq$p, breaks = 100)
head(which(b1b2_freq$contaminant))
#plot_frequency(phy, taxa_names(phy)[c(which(b1b2_freq$contaminant))], conc="conc")

decontam2 <- rownames(b1b2_freq[which(b1b2_freq$contaminant),])
unique(tax_table(prune_taxa(taxa_names(phy1) %in% decontam2, phy1))[,6])
phy1_decontam <- prune_taxa(!taxa_names(phy1_decontam) %in% decontam2, phy1_decontam) #-59




## Freq-based decontam for B3


phy <- subset_samples(phy1, batch %in% c("B3") & conc > 0)
phy <- subset_samples(phy, conc > 0)
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
b1b2_freq <- isContaminant(phy, method="frequency", conc="conc", threshold = 0.15)
table(b1b2_freq$contaminant)
hist(b1b2_freq$p, breaks = 200)
head(which(b1b2_freq$contaminant))
#plot_frequency(phy, taxa_names(phy)[c(which(b1b2_freq$contaminant))], conc="conc")

decontam3 <- rownames(b1b2_freq[which(b1b2_freq$contaminant),])
unique(tax_table(prune_taxa(taxa_names(phy1) %in% decontam3, phy1))[,6])

phy1_decontam <- prune_taxa(!taxa_names(phy1_decontam) %in% decontam3, phy1_decontam) #-44




## Freq-based decontam for B4

phy <- subset_samples(phy1, batch %in% c("B4") & conc > 0)
phy <- subset_samples(phy, !type %in% "control")
phy <- subset_samples(phy, conc > 0)
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
b1b2_freq <- isContaminant(phy, method="frequency", conc="conc", threshold = 0.15)
table(b1b2_freq$contaminant)
hist(b1b2_freq$p, breaks = 100)
head(which(b1b2_freq$contaminant))
#plot_frequency(phy, taxa_names(phy)[c(which(b1b2_freq$contaminant))], conc="conc")

decontam4 <- rownames(b1b2_freq[which(b1b2_freq$contaminant),])
unique(tax_table(prune_taxa(taxa_names(phy1) %in% decontam4, phy1))[,6])
phy1_decontam <- prune_taxa(!taxa_names(phy1_decontam) %in% decontam4, phy1_decontam) #-19


phy1_decontam <- subset_samples(phy1_decontam, sample_sums(phy1_decontam) > 0)


table(sample_data(phy1)$batch) - table(sample_data(phy1_decontam)$batch)














## Plot of library sizes showing that controls have lower library sizes, with many of the gut samples with similar library sizes as controls. Eggs have larger library sizes and do not overlap the controls, but are also dominated by sulcia, which we know must be vertically transmitted very generation. The relative abundance of Sulcia in the eggs was near 100% for all egg samples, suggesting that the amount of bacteria in the eggs is negligible to non-existent compared to the amount of Sulcia. Given that the gut sample library sizes overlap library sizes of the controls, it is unclear to what extent the gut bacterial diversity we uncovered reflects the true diversity and abundances of bacterria in the gut or the influence of contaminant bacteria. In other words, the most abundant bacteria in gut samples may be contamninant and differences in library sizes are driven by variable bacterial biomass in cicada guts.   

phy <- phy1
df <- as.data.frame(sample_data(phy)) 
df$LibrarySize <- sample_sums(phy)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

df <- data.frame(df)
meanc <- max(df[df$type == "control" & df$batch == "B4", "LibrarySize"], na.rm = TRUE)

df <- df[df$batch != "B2",]

g <- ggplot(data = df, aes(x = conc, y = log(LibrarySize), col = type)) +
	scale_color_brewer(palette = "Set1") +
	geom_point(alpha = 0.8, size = 2) +
	xlab("Library Concentration") +
	ylab("Log Bacterial Abundance") +
	labs(col = "Sample Type", shape = "Sample Set") +
	geom_hline(yintercept = log(1500), col = "red") +
	theme_bw() +
    theme(panel.grid = element_blank(),  plot.margin = unit(rep(0.6,4), "cm")) +
    facet_wrap(~batch)

ggsave("bacteria_abundance_predecontam.png", g)
ggsave("bacteria_abundance_postdecontam.png", g)

m <- mean(df[df$type == "control" & df$batch == "B4", "LibrarySize"])
m2 <-  mean(df[df$type == "control" & df$batch == "B4", "conc"])
s <- df[df$type != "control" & df$batch == "B4", "LibrarySize"]
s2 <- df[df$type != "control" & df$batch == "B4", "conc"]

length(s[s > m])/length(s)
length(s2[s2 > m2])/length(s2)

m <- mean(df[df$type == "control" & df$batch == "B3", "LibrarySize"])
m2 <-  mean(df[df$type == "control" & df$batch == "B3", "conc"])
s <- df[df$type != "control" & df$batch == "B3", "LibrarySize"]
s2 <- df[df$type != "control" & df$batch == "B3", "conc"]

length(s[s > m])/length(s)
length(s2[s2 > m2])/length(s2)













#### Removing samples with abundance < 1000 after decontam filtering

sums <- sample_sums(phy1_decontam)
hist(log(sums), breaks = 100)


phy1_decontam2 <- subset_samples(phy1_decontam, sample_sums(phy1_decontam) > 1500)

table(sample_data(phy1_decontam)$batch) - table(sample_data(phy1_decontam2)$batch)













## Taking out samples in B4 that cluster near controls 

controls <- subset_samples(phy1, type %in% "control" & batch %in% "B4")
controls_otu <- otu_table(controls)

phy <- subset_samples(phy1, batch %in% "B4" & type %in% "gut")
for(i in 1:length(rownames(otu_table(phy)))) {
	if(rownames( otu_table(phy)[i,] ) %in% rownames(otu_table(phy1_decontam2))) {
		otu_table(phy)[i,] <- otu_table(phy)[i,]
	} else if ( !rownames( otu_table(phy)[i,] ) %in% rownames(otu_table(phy1_decontam2)) ) {
		otu_table(phy)[i,] <- rep(0, length(otu_table(phy)[i,]))
	}
}

colnames(controls_otu) <- paste(colnames(controls_otu), "pre-decontam", sep = "_")
otu <- cbind(otu_table(phy), controls_otu)

otu <- otu[apply(otu, 1, sum) > 0,]
otu <- otu[, apply(otu, 2, sum) > 0 ]

type <- unlist(lapply(str_split(colnames(otu), "_"), "[", 2))
type[is.na(type)] <- "post-decontam"

phy2 <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), phy_tree(phy1))

ord <- ordinate(phy2, method = "NMDS", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$points[,1], PC2 = ord$points[,2], type = type )
xaverage <- mean(ordDF[ordDF$type == "pre-decontam", 1])
yaverage <- mean(ordDF[ordDF$type == "pre-decontam", 2])
ggplot(ordDF, aes(x = PC1, y = PC2, col = type)) + 
	geom_point(alpha = 0.8, size = 2) +
	geom_vline(xintercept = -0.11) +
	geom_vline(xintercept = 0.04) +
	geom_hline(yintercept = 0.045) +
	geom_hline(yintercept = -0.025) +
	geom_point(aes(x = xaverage, y = yaverage), col = "red") +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	stat_ellipse() +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pd <- ordDF[ordDF$type == "post-decontam", c(1,2)]

udist <- function(x){ sqrt( ((x[1] - xaverage) ^ 2) + ((x[2] - yaverage) ^ 2) )}

hist(apply(pd, 1, udist), breaks = 100)
distances <- log(apply(pd, 1, udist))
hist(distances, breaks = 100)
samples_remove <- names(distances[distances <= quantile(distances, probs = seq(0,1,0.01))[26]])

sample_data(subset_samples(phy1, ID %in% samples_remove))


ordDF$excluded <- !is.na(match(rownames(ordDF), samples_remove))
g <- ggplot(ordDF, aes(x = PC1, y = PC2, col = type)) + 
	geom_point(alpha = 0.8, size = 2) +
	geom_text(data = ordDF[ordDF$excluded == TRUE,], aes(x = PC1, y = PC2, label = "x"), col = "black") +
	geom_point(aes(x = xaverage, y = yaverage), col = "red") +
	xlab("NMDS1") +
	ylab("NMDS2") +
	scale_color_brewer(palette = "Paired") +
	#stat_ellipse() +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("B4_control_ordination.png", g, scale = 0.8)

phy1_decontam_reduced <- subset_samples(phy1_decontam2, !sample_names(phy1_decontam2) %in% samples_remove)







### FINAL FILTERING OF JUST GUT 

phy2 <- subset_samples(phy1, type %in% c("gut", "reproductive"))
phy1_decontam2 <- subset_samples(phy1_decontam2, type %in% c("gut", "reproductive"))
phy1_decontam_reduced <- subset_samples(phy1_decontam_reduced, type %in%  c("gut", "reproductive"))


table(sample_data(phy1_decontam2)$batch)
table(sample_data(phy1_decontam_reduced)$batch)

met <-  read.csv("projP-final-metadata.csv")
rownames(met) <- met$ID
sample_data(phy1) <- met
sample_data(phy2) <- met
sample_data(phy1_decontam2) <- met
sample_data(phy1_decontam_reduced) <- met




save.image("start.RData")









#### What's in the controls?

library(RColorBrewer)
n <- length(unique(tax_table(controls)[,4]))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))



controls <- subset_samples(phy1, type %in% "control")
controls <- subset_taxa(controls, taxa_sums(controls) > 0)

unique(tax_table(controls)[,3])

otu <- data.frame(otu_table(controls), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(controls)[,4]), as.character(tax_table(controls)[,3]))
otu$taxa <- paste(as.character(tax_table(controls)[,5]), as.character(tax_table(controls)[,6]), sep = "; ")

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]


limits <- sample_data(controls)[ order(sample_data(controls)[, "species"]$species), "ID"]$ID
labels <- paste(sample_data(controls)[ order(sample_data(controls)[, "species"]$species), "species"]$species, sample_data(controls)[ order(sample_data(controls)[, "species"]$species), "batch"]$batch)

g <- ggplot(data = data.frame(otu2), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    ylab("Abundance") +
    xlab("") +
    theme_bw() +
    scale_x_discrete(limits = limits, labels =  labels) +
    theme(panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(angle = 45, hjust = 1)
    	)
    	
ggsave("control__bar.png", g)




# at the family level 
tax1 <- unique(unlist(lapply(str_split(otu2$taxa, " "), "[" , 1)))
tax2 <- unique(unlist(lapply(str_split(otu2$taxa, " "), "[" , 2)))

guts <- subset_samples(phy1_decontam_reduced, type %in% "gut")
guts <- subset_taxa(guts, Family %in% tax1 )
guts <- subset_taxa(guts, taxa_sums(guts) > 0)

guts2 <- subset_samples(phy1_decontam_reduced, type %in% "gut")
guts2 <- subset_taxa(guts2, !Family %in% tax1 )
guts2 <- subset_taxa(guts2, taxa_sums(guts2) > 0)

mean(taxa_sums(guts))/mean(taxa_sums(guts2))
sum(taxa_sums(guts))/sum(taxa_sums(guts2))


#at the ASV level 
guts <- subset_samples(phy1_decontam_reduced, type %in% "gut")
guts <- subset_taxa(guts, Family %in% tax1 & Genus %in% tax2)
guts <- subset_taxa(guts, taxa_sums(guts) > 0)

guts2 <- subset_samples(phy1_decontam_reduced, type %in% "gut")
guts2 <- subset_taxa(guts2, !Family %in% tax1 & !Genus %in% tax2 | is.na(Family) | is.na(Genus))
guts2 <- subset_taxa(guts2, taxa_sums(guts2) > 0)

mean(taxa_sums(guts))/mean(taxa_sums(guts2))
sum(taxa_sums(guts))/sum(taxa_sums(guts2))

guts <- subset_samples(phy1_decontam_reduced, type %in% "gut")
guts <- subset_taxa(guts, Genus %in% unique(tax_table(controls)[,6]))
guts <- subset_taxa(guts, taxa_sums(guts) > 0)

a <- otu_table(subset_samples(phy1, type %in% "control"))
b <- otu_table(guts)
a <- a[rownames(a) %in% rownames(b),]
c <- cbind(a, b)
d <- phyloseq(otu_table(c, taxa_are_rows = TRUE), tax_table(phy1), sample_data(phy1))
d <- tax_glom(d, taxrank = "Genus")
d <- transform_sample_counts(d, function(x) x/sum(x))
d <- subset_samples(d, sample_sums(d) > 0)

e <- subset_taxa(d, taxa_names(d) %in% rownames(otu_table(d)[apply(otu_table(d), 1, function(x) length(x[x > 0 ])) > 5,]))
e <- subset_samples(e, sample_sums(e) > 0)


heatmap(otu_table(e), Colv = NA, scale="column", col = colorRampPalette(brewer.pal(8, "Blues"))(10), labRow = paste(tax_table(e)[,3], tax_table(e)[,5], tax_table(e)[,6], sep = "; "), cexRow = 0.5)

otu <- data.frame(otu_table(e), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
otu$taxa <- paste(as.character(tax_table(e)[,5]), as.character(tax_table(e)[,6]))
otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]


ggplot(data = data.frame(otu2), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = rownames(sample_data(e)[order(sample_data(e)$type),])) +
    theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))













#### Hybrids !!!!

phy1_decontam_reduced_rel <- transform_sample_counts(phy1_decontam2, function(x) x/sum(x))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B1", "B2"))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,5])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) mean(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.00115, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
order <- order(sd$batch, sd$species.tree)
limits = rownames(sd[order,])
labels = paste(sd[order, "species.tree"]$species.tree, sd[order, "code"]$code, sep = " ; ")

g <- ggplot(data = data.frame(otu3), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = limits, labels = labels) +
    xlab("") +
    ylab("") +
    labs(fill = "Bacterial Family") +
    coord_flip() +
    theme_minimal() +
    theme(plot.margin = unit(c(1,10,0.5,1), "cm"), 
    	legend.position = "right", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 5)) 
    	
ggsave("hybrid_relative_abundance.png", g)







### Process date barplots

phy1_decontam_reduced_rel <- transform_sample_counts(phy1_decontam2, function(x) x/sum(x))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B1"))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,6])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) mean(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.0025, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
order <- order(sd$process_date, sd$batch, sd$species.tree)
limits = rownames(sd[order,])
labels = paste(sd[order, "species.tree"]$species.tree, sd[order, "code"]$code, sd[order, "process_date"]$process_date,  sep = " ; ")

g <- ggplot(data = data.frame(otu3), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = limits, labels = labels) +
    xlab("") +
    ylab("") +
    labs(fill = "Bacterial Family") +
    coord_flip() +
    theme_minimal() +
    theme(plot.margin = unit(c(1,10,0.5,1), "cm"), 
    	legend.position = "right", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 5)) 
    	
ggsave("B1_process_barplot_genus.png", g)


phy1_decontam_reduced_rel <- transform_sample_counts(phy1_decontam2, function(x) x/sum(x))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B1"))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,5])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) mean(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.0015, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
order <- order(sd$process_date, sd$batch, sd$species.tree)
limits = rownames(sd[order,])
labels = paste(sd[order, "species.tree"]$species.tree, sd[order, "code"]$code, sd[order, "process_date"]$process_date,  sep = " ; ")

g <- ggplot(data = data.frame(otu3), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = limits, labels = labels) +
    xlab("") +
    ylab("") +
    labs(fill = "Bacterial Family") +
    coord_flip() +
    theme_minimal() +
    theme(plot.margin = unit(c(1,10,0.5,1), "cm"), 
    	legend.position = "right", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 5)) 
    	
ggsave("B1_process_barplot_family.png", g)






	 ## Ordinations
	 
muta <- subset_samples(phy1_decontam2, batch %in% "B1")
#muta <- transform_sample_counts(muta, function(x) log(x+1))
meta <- sample_data(muta)$species.tree
meta2 <- sample_data(muta)$code
meta3 <-  sample_sums(muta)
meta4 <- sample_data(muta)$process_date
ord <- ordinate(muta, method = "NMDS", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$points[,1], PC2 = ord$points[,2], meta = meta, meta2 = meta2, meta3 = meta3, meta4 = meta4)
g <- ggplot(ordDF, aes(x = PC1, y = PC2, col = meta, shape = meta4)) + 
	geom_point(alpha = 0.8, size = log(order(meta3))) +
	geom_text_repel(aes(label = meta2), size = 2) +
	xlab("NMDS1") +
	ylab("NMDS2") +
	labs(col = "") +
	scale_color_brewer(palette = "Paired") +
	#stat_ellipse(data = ordDF, aes(x = PC1, y = PC2, shape = meta4), inherit.aes = FALSE) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("hybrid_muta_ordination.png", g, scale = 0.6)



westlandica <- subset_samples(phy1_decontam2, batch %in% "B2")
meta <- as.character(sample_data(westlandica)$species.tree)
meta[is.na(meta)] <- "westlandica-south-north-hybrid"
meta2 <- sample_data(westlandica)$code
meta3 <-  sample_sums(westlandica)
ord <- ordinate(westlandica, method = "NMDS", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$points[,1], PC2 = ord$points[,2], meta = meta, meta2 = meta2, meta3 = meta3)
g <- ggplot(ordDF, aes(x = PC1, y = PC2, col = meta)) + 
	geom_point(alpha = 0.8, size = log(order(meta3))) +
	geom_text_repel(aes(label = meta2), size = 2) +
	xlab("NMDS1") +
	ylab("NMDS2") +
	labs(col = "") +
	scale_color_brewer(palette = "Set1") +
	#stat_ellipse() +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("hybrid_westlandica_ordination.png", g, scale = 0.6)





	# Distance boxplots


phy <- subset_samples(phy1_decontam2, batch %in% "B1")
sam <- sample_data(phy)

pdis_unifrac <- phyloseq::distance(phy, method="wunifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])
dis <- pdis2
dis$elev_diff <- abs(as.numeric(unlist(lapply(str_split(dis$ID, "_"), "[", 2))) - as.numeric(unlist(lapply(str_split(dis$id, "_"), "[", 2))))
dis$habitat_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 3)) == unlist(lapply(str_split(dis$id, "_"), "[", 3)), yes = 0, no = 1)
dis$island_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 4)) == unlist(lapply(str_split(dis$id, "_"), "[", 4)), yes = 0, no = 1)
dis$batch_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 5)) == unlist(lapply(str_split(dis$id, "_"), "[", 5)), yes = 0, no = 1)
dis$speciesID <- unlist(lapply(str_split(dis$ID, "_"), "[", 6))
dis$speciesid <- unlist(lapply(str_split(dis$id, "_"), "[", 6))
head(dis)
dis$hybrid <- rep(NA, dim(dis)[1])
#dis2 <- dis[!dis$speciesid %in% "muta-tuta-hybrid",]
dis2 <- dis
dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "muta-SI", "hybrid"] <- "MutaSI / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "muta-SI", "hybrid"] <- "MutaSI / Hybrid"

dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "muta-NI", "hybrid"] <- "MutaNI / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "muta-NI", "hybrid"] <- "MutaNI / Hybrid"

dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "tuta-clade1-northwestSI", "hybrid"] <- "TutaClade1 / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "tuta-clade1-northwestSI", "hybrid"] <- "TutaClade1 / Hybrid"

dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "tuta-clade2-northeastSI", "hybrid"] <- "TutaClade2 / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "tuta-clade2-northeastSI", "hybrid"] <- "TutaClade2 / Hybrid"

dis2[dis2$speciesID %in% c("muta-SI") & dis2$speciesid %in% c("muta-SI"), "hybrid"] <- "MutaSI/MutaSI"
dis2[dis2$speciesID %in% c("tuta-clade1-northwestSI") & dis2$speciesid %in% c("tuta-clade1-northwestSI"), "hybrid"] <- "TutaClade1 / TutaClade1"
dis2[dis2$speciesID %in% c("tuta-clade2-northeastSI") & dis2$speciesid %in% c("tuta-clade2-northeastSI"), "hybrid"] <- "TutaClade2 / TutaClade2"
dis2$pdis <- as.numeric(dis2$pdis)
dis2 <- dis2[!is.na(dis2$hybrid),]
dis_wunifrac_muta <- dis2



pdis_unifrac <- phyloseq::distance(phy, method="unifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])
dis <- pdis2
dis$elev_diff <- abs(as.numeric(unlist(lapply(str_split(dis$ID, "_"), "[", 2))) - as.numeric(unlist(lapply(str_split(dis$id, "_"), "[", 2))))
dis$habitat_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 3)) == unlist(lapply(str_split(dis$id, "_"), "[", 3)), yes = 0, no = 1)
dis$island_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 4)) == unlist(lapply(str_split(dis$id, "_"), "[", 4)), yes = 0, no = 1)
dis$batch_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 5)) == unlist(lapply(str_split(dis$id, "_"), "[", 5)), yes = 0, no = 1)
dis$speciesID <- unlist(lapply(str_split(dis$ID, "_"), "[", 6))
dis$speciesid <- unlist(lapply(str_split(dis$id, "_"), "[", 6))
head(dis)
dis$hybrid <- rep(NA, dim(dis)[1])
#dis2 <- dis[!dis$speciesid %in% "muta-tuta-hybrid",]
dis2 <- dis
dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "muta-SI", "hybrid"] <- "MutaSI / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "muta-SI", "hybrid"] <- "MutaSI / Hybrid"

dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "muta-NI", "hybrid"] <- "MutaNI / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "muta-NI", "hybrid"] <- "MutaNI / Hybrid"

dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "tuta-clade1-northwestSI", "hybrid"] <- "TutaClade1 / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "tuta-clade1-northwestSI", "hybrid"] <- "TutaClade1 / Hybrid"

dis2[dis2$speciesID == "muta-tuta-hybrid" & dis2$speciesid == "tuta-clade2-northeastSI", "hybrid"] <- "TutaClade2 / Hybrid"
dis2[dis2$speciesid == "muta-tuta-hybrid" & dis2$speciesID == "tuta-clade2-northeastSI", "hybrid"] <- "TutaClade2 / Hybrid"

dis2[dis2$speciesID %in% c("muta-SI") & dis2$speciesid %in% c("muta-SI"), "hybrid"] <- "MutaSI/MutaSI"
dis2[dis2$speciesID %in% c("tuta-clade1-northwestSI") & dis2$speciesid %in% c("tuta-clade1-northwestSI"), "hybrid"] <- "TutaClade1 / TutaClade1"
dis2[dis2$speciesID %in% c("tuta-clade2-northeastSI") & dis2$speciesid %in% c("tuta-clade2-northeastSI"), "hybrid"] <- "TutaClade2 / TutaClade2"
dis2$pdis <- as.numeric(dis2$pdis)
dis2 <- dis2[!is.na(dis2$hybrid),]
dis_unifrac_muta <- dis2



phy <- subset_samples(phy1_decontam2, batch %in% "B2")
sam <- sample_data(phy)
sam$species.tree <- as.character(sam$species.tree)
sam$species.tree[is.na(sam$species.tree)] <- "westlandica-north-south-hybrid"

pdis_unifrac <- phyloseq::distance(phy, method="wunifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])
dis <- pdis2
dis$elev_diff <- abs(as.numeric(unlist(lapply(str_split(dis$ID, "_"), "[", 2))) - as.numeric(unlist(lapply(str_split(dis$id, "_"), "[", 2))))
dis$habitat_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 3)) == unlist(lapply(str_split(dis$id, "_"), "[", 3)), yes = 0, no = 1)
dis$island_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 4)) == unlist(lapply(str_split(dis$id, "_"), "[", 4)), yes = 0, no = 1)
dis$batch_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 5)) == unlist(lapply(str_split(dis$id, "_"), "[", 5)), yes = 0, no = 1)
dis$speciesID <- unlist(lapply(str_split(dis$ID, "_"), "[", 6))
dis$speciesid <- unlist(lapply(str_split(dis$id, "_"), "[", 6))
head(dis)
dis$hybrid <- rep(NA, dim(dis)[1])
#dis2 <- dis[!dis$speciesid %in% "muta-tuta-hybrid",]
dis2 <- dis
dis2[dis2$speciesID == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesid == "westlandica-south", "hybrid"] <- "WestlandicaSouth / Hybrid"
dis2[dis2$speciesid == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesID == "westlandica-south", "hybrid"] <- "WestlandicaSouth / Hybrid"

dis2[dis2$speciesID == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesid == "westlandica-north", "hybrid"] <- "WestlandicaNorth / Hybrid"
dis2[dis2$speciesid == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesID == "westlandica-north", "hybrid"] <- "WestlandicaNorth / Hybrid"

dis2[dis2$speciesID == "westlandica-south" & dis2$speciesid == "westlandica-south", "hybrid"] <- "WestlandicaSouth / WestlandicaSouth"
dis2[dis2$speciesid == "westlandica-south" & dis2$speciesID == "westlandica-south", "hybrid"] <- "WestlandicaSouth / WestlandicaSouth"

dis2[dis2$speciesID == "westlandica-north" & dis2$speciesid == "westlandica-north", "hybrid"] <- "WestlandicaNorth / WestlandicaNorth"
dis2[dis2$speciesid == "westlandica-north" & dis2$speciesID == "westlandica-north", "hybrid"] <- "WestlandicaNorth / WestlandicaNorth"

dis2$pdis <- as.numeric(dis2$pdis)
dis2 <- dis2[!is.na(dis2$hybrid),]
dis_wunifrac_westlandica <- dis2

phy <- subset_samples(phy1_decontam2, batch %in% "B2")
sam <- sample_data(phy)
sam$species.tree <- as.character(sam$species.tree)
sam$species.tree[is.na(sam$species.tree)] <- "westlandica-north-south-hybrid"



pdis_unifrac <- phyloseq::distance(phy, method="unifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])
dis <- pdis2
dis$elev_diff <- abs(as.numeric(unlist(lapply(str_split(dis$ID, "_"), "[", 2))) - as.numeric(unlist(lapply(str_split(dis$id, "_"), "[", 2))))
dis$habitat_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 3)) == unlist(lapply(str_split(dis$id, "_"), "[", 3)), yes = 0, no = 1)
dis$island_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 4)) == unlist(lapply(str_split(dis$id, "_"), "[", 4)), yes = 0, no = 1)
dis$batch_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 5)) == unlist(lapply(str_split(dis$id, "_"), "[", 5)), yes = 0, no = 1)
dis$speciesID <- unlist(lapply(str_split(dis$ID, "_"), "[", 6))
dis$speciesid <- unlist(lapply(str_split(dis$id, "_"), "[", 6))
head(dis)
dis$hybrid <- rep(NA, dim(dis)[1])
#dis2 <- dis[!dis$speciesid %in% "muta-tuta-hybrid",]
dis2 <- dis
dis2[dis2$speciesID == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesid == "westlandica-south", "hybrid"] <- "WestlandicaSouth / Hybrid"
dis2[dis2$speciesid == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesID == "westlandica-south", "hybrid"] <- "WestlandicaSouth / Hybrid"

dis2[dis2$speciesID == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesid == "westlandica-north", "hybrid"] <- "WestlandicaNorth / Hybrid"
dis2[dis2$speciesid == "westlandicaSouth-westlandicaNorth-hybrid" & dis2$speciesID == "westlandica-north", "hybrid"] <- "WestlandicaNorth / Hybrid"

dis2[dis2$speciesID == "westlandica-south" & dis2$speciesid == "westlandica-south", "hybrid"] <- "WestlandicaSouth / WestlandicaSouth"
dis2[dis2$speciesid == "westlandica-south" & dis2$speciesID == "westlandica-south", "hybrid"] <- "WestlandicaSouth / WestlandicaSouth"

dis2[dis2$speciesID == "westlandica-north" & dis2$speciesid == "westlandica-north", "hybrid"] <- "WestlandicaNorth / WestlandicaNorth"
dis2[dis2$speciesid == "westlandica-north" & dis2$speciesID == "westlandica-north", "hybrid"] <- "WestlandicaNorth / WestlandicaNorth"

dis2$pdis <- as.numeric(dis2$pdis)
dis2 <- dis2[!is.na(dis2$hybrid),]
dis_unifrac_westlandica <- dis2



wunifrac <- rbind(dis_wunifrac_westlandica, dis_wunifrac_muta)
unifrac <- rbind(dis_unifrac_westlandica, dis_unifrac_muta)

dis <- data.frame(hybrid = wunifrac$hybrid, wunifrac = wunifrac$pdis, unifrac = unifrac$pdis)

unique(dis$hybrid)

mutah <- dis[dis$hybrid %in% c("MutaNI / Hybrid", "TutaClade2 / Hybrid", "TutaClade1 / Hybrid","MutaSI / Hybrid"), "wunifrac"]
mutah2 <- dis[dis$hybrid %in% c("TutaClade2 / TutaClade2", "TutaClade1 / TutaClade1", "MutaSI/MutaSI"), "wunifrac"]

westh <- dis[dis$hybrid %in% c("WestlandicaSouth / Hybrid", "WestlandicaNorth / Hybrid"), "wunifrac"]
westh2 <- dis[dis$hybrid %in% c("WestlandicaSouth / WestlandicaSouth", "WestlandicaNorth / WestlandicaNorth"), "wunifrac"]

t.test(mutah, mutah2)
t.test(westh, westh2)
wilcox.test(mutah, mutah2)
wilcox.test(westh, westh2)



qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


g <- ggplot(data = dis, aes(x = hybrid, y = -unifrac)) +
	geom_boxplot(data = dis, aes(x = hybrid, y = -wunifrac), alpha = 0.1, col = col_vector[11], fill = "white", outlier.colour = NULL) +
	geom_boxplot(data = dis, aes(x = hybrid, y = -unifrac), alpha = 0.1, col =  col_vector[22], fill = "white", outlier.colour = NULL) +
	xlab("") +
	ylab("") +
	geom_jitter(data = dis, aes(x = hybrid, y = -wunifrac), width = 0.1, shape = 1, col = col_vector[11]) +
	geom_jitter(data = dis, aes(x = hybrid, y = -unifrac), width = 0.1, shape = 1, col = col_vector[22]) +
	theme_bw() + 
    theme(panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    	)


ggsave("hybrid_pairwise.png", g, scale = 0.8)











## Hybrid PhilR

#BiocManager::install("philr")
library(philr)

phy <- subset_samples(phy1_decontam2, batch %in% "B1")
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
plot_ordination(phy, gp.pcoa, color='species') + geom_point(size=4)
gp <- data.frame(gp.philr)
gp$species <- sample_data(metadata)$species
gp2 <- gather(gp, key = node, value = balance, 1:(dim(gp)[2]-1))
gp3 <- gp2[gp2$balance < -1 | gp2$balance > 1,]
balname <- apply(gp3, 1, function(x) {  name.balance(phy_tree(phy), tax_table(phy), x["node"])  })
order <- order(as.numeric(unlist(lapply(str_split(gp3$node, "n"), "[", 2))))
gp3$balname <- balname
g <- ggplot(data = gp3, aes(x = balname, y = balance, col = species)) +
	geom_point(alpha = 0.7, size = 3) +
	geom_point(shape = 1, size = 3, col = "black") +
	coord_flip() +
	xlab("") +
	ylab("") +
	scale_x_discrete(limits = unique(gp3$balname)) +
	scale_color_manual(values = cols) +
	theme_bw() + 
    theme(
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(),
    	axis.text.y = element_text(size = 6)
    	)

ggsave("hybrid_B2_balance.png", g)














phy <- subset_samples(phy1_decontam2, batch %in% "B1")
phy <- transform_sample_counts(phy, function(x) x+1)
phy_tree(phy) <- makeNodeLabel(phy_tree(phy), method="number", prefix='n')
name.balance(phy_tree(phy), tax_table(phy), 'n1')
otu.table <- t(otu_table(phy))
tree <- phy_tree(phy)
metadata <- sample_data(phy)
tax <- tax_table(phy)
gp.philr <- philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(phy, 'PCoA', distance=gp.dist)
plot_ordination(phy, gp.pcoa, color='species') + geom_point(size=4)
gp <- data.frame(gp.philr)
gp$species <- sample_data(metadata)$species
gp2 <- gather(gp, key = node, value = balance, 1:(dim(gp)[2]-1))
gp3 <- gp2[gp2$balance < -0.5 | gp2$balance > 0.5,]
balname <- apply(gp3, 1, function(x) {  name.balance(phy_tree(phy), tax_table(phy), x["node"])  })
order <- order(as.numeric(unlist(lapply(str_split(gp3$node, "n"), "[", 2))))
gp3$balname <- balname
a <- aggregate(gp3$balance, by = list(gp3$balname), function(x) length(x))
hist(a[,2], breaks = 50)
gp4 <- gp3[gp3$balname %in% a[a[,2] > 10, 1],]
qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) == 'Paired',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
g1 <- ggplot(data = gp4, aes(x = balname, y = balance, col = species)) +
	geom_point(alpha = 0.7, size = 3) +
	#geom_point(shape = 1, size = 3, col = "black") +
	coord_flip() +
	xlab("") +
	ylab("") +
	scale_x_discrete(limits = unique(gp4$balname)) +
	scale_color_manual(values = col_vector[c(1,5,2), ]) +
	theme_bw() + 
    theme(
        legend.position = "bottom",
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(),
    	axis.text.y = element_text(size = 6)
    	)


phy <- subset_samples(phy1_decontam2, batch %in% "B2")
phy <- transform_sample_counts(phy, function(x) x+1)
phy_tree(phy) <- makeNodeLabel(phy_tree(phy), method="number", prefix='n')
name.balance(phy_tree(phy), tax_table(phy), 'n1')
otu.table <- t(otu_table(phy))
tree <- phy_tree(phy)
metadata <- sample_data(phy)
tax <- tax_table(phy)
gp.philr <- philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(phy, 'PCoA', distance=gp.dist)
plot_ordination(phy, gp.pcoa, color='species') + geom_point(size=4)
gp <- data.frame(gp.philr)
gp$species <- sample_data(metadata)$species
gp2 <- gather(gp, key = node, value = balance, 1:(dim(gp)[2]-1))
gp3 <- gp2[gp2$balance < -0.5 | gp2$balance > 0.5,]
balname <- apply(gp3, 1, function(x) {  name.balance(phy_tree(phy), tax_table(phy), x["node"])  })
order <- order(as.numeric(unlist(lapply(str_split(gp3$node, "n"), "[", 2))))
gp3$balname <- balname
a <- aggregate(gp3$balance, by = list(gp3$balname), function(x) length(x))
hist(a[,2], breaks = 50)
gp4 <- gp3[gp3$balname %in% a[a[,2] > 17, 1],]
qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) == 'Paired',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
g2 <- ggplot(data = gp4, aes(x = balname, y = balance, col = species)) +
	geom_point(alpha = 0.7, size = 3) +
	#geom_point(shape = 1, size = 3, col = "black") +
	coord_flip() +
	xlab("") +
	ylab("") +
	scale_x_discrete(limits = unique(gp4$balname)) +
	scale_color_manual(values = col_vector[c(1,2,5), ]) +
	theme_bw() + 
    theme(
    	legend.position = "bottom",
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(),
    	axis.text.y = element_text(size = 6)
    	)
    	
    	
ggsave("hybrid_B1_reduced.png", g1, scale = 0.8)
ggsave("hybrid_B2_reduced.png", g2, scale = 0.8)














### Map of hybrids 


#install.packages("ggmap")
library(ggmap)
#install.packages("maps")
library(maps)
#install.packages("mapdata")
library(mapdata)
#install.packages("mapproj")
library(mapproj)
#install.packages("maptools")
library(maptools)
#install.packages("rgdal")
library(rgdal)
library(raster)

sdata <- sample_data(phy1_decontam_reduced)[sample_data(phy1_decontam_reduced)$island %in% c("South Island", "North Island"),]
sdata <- sdata[sdata$batch %in% c("B1", "B2"), ]
sdata <- data.frame(sdata)
sdata$lat_sign <- as.numeric(sdata$lat_sign)
sdata$lon_sign <- as.numeric(sdata$lon_sign)
sdata$species <- as.character(sdata$species)
sdata$code2 <- unlist(lapply(lapply(str_split(sdata$code, "[.]"), "[", 3:4), function(x){paste(x[1], x[2], sep = ".")}))

sdata[!is.na(sdata$species.tree) & sdata$species == "hybrid", "species"] <- 
sdata[is.na(sdata$species.tree), ] <- "westlandica_hybrid"

chi_bb <- c(left = min(sdata$lon_sign-0.5) ,
            bottom = min(sdata$lat_sign-0.5)  ,
            right = max(sdata$lon_sign+0.5) ,
            top = max(sdata$lat_sign)+0.5)

mp <- get_stamenmap(bbox = chi_bb, zoom = 8, maptype = "terrain-background")

g <- ggmap(mp) +
	geom_point(data = sdata, aes(x = lon_sign, y = lat_sign, col = species), size = 3) +
	geom_point(data = sdata, aes(x = lon_sign, y = lat_sign), shape = 1, size = 3, col = "black") +
	theme_minimal() +
	scale_color_brewer(palette = "Set2") +
	geom_text_repel(data = sdata[!duplicated(sdata$code2), ], aes(x = lon_sign, y = lat_sign, label = code2), size = 3, segment.size = 0.4, segment.alpha = 0.8)

ggsave("B1_hybrid_map.png", g)

























####### Map of B3 with barplots
library(RColorBrewer)
n <- 17
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
col_vector <- sample(col_vector, n)


phy1_decontam_reduced_rel <- transform_sample_counts(subset_samples(phy1_decontam2, type %in% "gut"), function(x) x/sum(x))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B3"))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,5])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) mean(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.0015, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
order <- order(sd$batch, sd$genus, sd$species, sd$species.tree)
limits = rownames(sd[order,])
labels = paste(sd[order, "species.tree"]$species.tree, sd[order, "code"]$code, sep = " ; ")

g <- ggplot(data = data.frame(otu3), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = limits, labels = labels) +
    xlab("") +
    ylab("") +
    labs(fill = "Bacterial Family") +
    coord_flip() +
    theme_minimal() +
    theme(plot.margin = unit(c(1,10,0.5,1), "cm"), 
    	legend.position = "right", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 6.5,  hjust = 0.5, vjust = 0.5)) 
    	
ggsave("B3_relative_abundance_family.png", g)





phy1_decontam_reduced_rel <- transform_sample_counts(subset_samples(phy1_decontam2, type %in% "gut"), function(x) x/sum(x))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B3"))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,6])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) max(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.1936, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
order <- order(sd$batch, sd$genus, sd$species, sd$species.tree)
limits = rownames(sd[order,])
labels = paste(sd[order, "species.tree"]$species.tree, sd[order, "code"]$code, sep = " ; ")

g <- ggplot(data = data.frame(otu3), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = limits, labels = labels) +
    xlab("") +
    ylab("") +
    labs(fill = "Bacterial Genus") +
    coord_flip() +
    theme_minimal() +
    theme(plot.margin = unit(c(1,10,0.5,1), "cm"), 
    	legend.position = "left", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 6.5,  hjust = 0.5, vjust = 0.5)) 
    	
ggsave("B3_relative_abundance_genus.png", g)














### Hybrids modeling 

library(pscl)
library(MASS)
library(boot)


phy1_decontam_reduced_rel <- subset_samples(phy1_decontam2, type %in% "gut")
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B3"))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,5])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) mean(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.0015, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
order <- order(sd$batch, sd$genus, sd$species, sd$species.tree)
limits = rownames(sd[order,])
labels = paste(sd[order, "species.tree"]$species.tree, sd[order, "code"]$code, sep = " ; ")



mod <- glm.nb(abund ~ sample, data = otu3)
summary(mod)





















## B3 With phylogeny, heatmap

library(phytools)
library(ggtree)
library(pheatmap)
library(ggstance)

htree <- read.tree("raxml.tre")

htree <- reroot(htree, node = 70)

t <- htree

t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist( lapply( str_split(t$tip.lab, "_"), "[", 1) )

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))

phy <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B3"))

t2 <- drop.tip(t, t$tip.label[!t$tip.label %in% sample_data(phy)$species.tree])

g <- ggtree(t2) + 
	#geom_tiplab(align = TRUE, size = 3) +
	geom_nodelab( nudge_x = 0.01, nudge_y = 0.001, size = 2.5) + 
	xlim(0,0.5)

ggsave("B3_tree.png", g)


gdata <- data.frame(g$data)
gdata <- gdata[order(gdata$y, decreasing = TRUE), ]
tiporder <- gdata[is.na(as.numeric(gdata$label)), "label"]
tiporder <- tiporder[tiporder != ""]







## Class barplot

phy1_decontam_reduced_rel <- transform_sample_counts(subset_samples(phy1_decontam2, type %in% "gut"), function(x) x/sum(x))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B3"))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,3])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) mean(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.00028, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
sd[sd$species.tree == "muta-NI", "species.tree"] <- "muta-SI"
sd$order <- match(sd$species.tree, tiporder)
limits = sd[order(sd$order, decreasing = TRUE), "ID"]$ID
labels = paste(sd[order(sd$order, decreasing = TRUE), "species.tree"]$species.tree, sd[order(sd$order,  decreasing = TRUE), "code"]$code, sep = " ; ")

g <- ggplot(data = data.frame(otu3), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = limits, labels = labels) +
    xlab("") +
    ylab("") +
    labs(fill = "Bacterial Class") +
    coord_flip() +
    theme_minimal() +
    theme(plot.margin = unit(c(1,10,0.5,1), "cm"), 
    	legend.position = "right", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    	axis.text.y = element_text(size = 6.5,  hjust = 1, vjust = 0.5)) 
    	
ggsave("B3_relative_abundance_class.png", g)











## B3 Map

#install.packages("ggmap")
library(ggmap)
#install.packages("maps")
library(maps)
#install.packages("mapdata")
library(mapdata)
#install.packages("mapproj")
library(mapproj)
#install.packages("maptools")
library(maptools)
#install.packages("rgdal")
library(rgdal)
library(raster)

sdata <- sample_data(phy1_decontam_reduced)[sample_data(phy1_decontam_reduced)$island %in% c("South Island", "North Island"),]
sdata <- sdata[sdata$batch %in% c("B3"), ]
sdata <- data.frame(sdata)
sdata$lat_sign <- as.numeric(sdata$lat_sign)
sdata$lon_sign <- as.numeric(sdata$lon_sign)
sdata$species <- as.character(sdata$species)
sdata$code2 <- unlist(lapply(lapply(str_split(sdata$code, "[.]"), "[", 3:4), function(x){paste(x[1], x[2], sep = ".")}))

chi_bb <- c(left = min(sdata$lon_sign[!is.na(sdata$lon_sign)]-0.5) ,
            bottom = min(sdata$lat_sign[!is.na(sdata$lat_sign)]-0.5)  ,
            right = max(sdata$lon_sign[!is.na(sdata$lon_sign)]+0.5) ,
            top = max(sdata$lat_sign[!is.na(sdata$lat_sign)]+0.5))

mp <- get_stamenmap(bbox = chi_bb, zoom = 8, maptype = "terrain-background")

library(RColorBrewer)
n <- 21
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
col_vector2 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector2 <- sample(col_vector2, n)


g <- ggmap(mp) +
	geom_point(data = sdata, aes(x = lon_sign, y = lat_sign, col = genus), size = 3) +
	geom_point(data = sdata, aes(x = lon_sign, y = lat_sign), shape = 1, size = 3, col = "grey10") +
	theme_minimal() +
	geom_text_repel(data = sdata[!duplicated(sdata$code2), ], aes(x = lon_sign, y = lat_sign, label = code2), size = 2.5, segment.size = 0.4, segment.alpha = 0.8) +
	scale_color_brewer(palette = "Set2") 
	
ggsave("B3_map.png", g)




























#### B3 heatmap

library(pheatmap)

pheat <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B3"))
pheat <- transform_sample_counts(pheat, function(x) x/sum(x))
#pheat <- merge_samples(pheat, "species.tree")
asv_prev <- apply(otu_table(pheat), 1, function(x) length(x[x > 0]))
asv <- names(asv_prev[asv_prev > 1])
pheat <- prune_taxa(taxa = asv, pheat)


pheatotu <- otu_table(pheat)
abund_ranks <- apply(pheatotu, 1, rank)
hist(abund_ranks, breaks = 100)
abund_ranks <- abund_ranks - 26
abund_ranks[abund_ranks < 1] <- 1
abund_ranks <- data.frame(abund_ranks)
#abund_ranks <- abund_ranks[, apply(abund_ranks, 2, function(x) min(x[x > 1])) > 10]
otu_table(pheat) <- otu_table(abund_ranks, taxa_are_rows = FALSE)


heat_sam <- sample_data(pheat)
heat <- otu_table(pheat)
colnames(heat) <- paste(tax_table(pheat)[,2], paste(tax_table(pheat)[,5], seq(1, length(tax_table(pheat)[,2])), sep = "."), sep = " - ")
rownames(heat) <- paste(heat_sam$ID, heat_sam$code, heat_sam$species)
heat <- heat[order(as.numeric(unlist(lapply(str_split(lapply(str_split(heat_sam$ID, "K"), "[", 2), "g"), "[", 1)))),]

heat_clust <- hclust(dist(heat), method = "complete")
heat_clust_cut <- cutree(tree = heat_clust, k = 7)
bac_clust <-  hclust(dist(t(heat)), method = "complete")
bac_clust_cut <- cutree(tree = bac_clust, k = 2)

heat_data <- data.frame(cluster = heat_clust_cut)
bac_data <- data.frame(cluster = bac_clust_cut)

heat_data <- as.numeric(heat_sam$elevation)
names(heat_data) <- rownames(heat)
heat_data <- data.frame(heat_data)

g <- pheatmap(heat, fontsize_col = 6, cluster_rows = FALSE, fontsize_row = 6, annotation_legend = FALSE)

ggsave("B3_heatmap.png", g)

















	 ## Ordinations
	 
phy <- subset_samples(phy1_decontam2, batch %in% "B3")
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- tax_glom(phy, taxrank = "Genus")
sam <- sample_data(phy)

meta <- sample_data(phy)$genus
meta2 <- paste(sample_data(phy)$species.tree, paste(lapply(str_split(sample_data(phy)$code, "[[.]]"), "[", 3), lapply(str_split(sample_data(phy)$code, "[[.]]"), "[", 4), lapply(str_split(sample_data(phy)$code, "[[.]]"), "[", 5), sep = "."), sep = " ")
meta3 <-  sample_sums(phy)
ord <- ordinate(phy, method = "PCoA", distance = "unifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], meta = meta, meta2 = meta2, meta3 = meta3)
sam$ord <- rep(1, dim(sam)[1])
sam[sam$ID %in% rownames(ordDF[ordDF$PC1 < -0.05,]), "ord"] <- "group1"
sam[sam$ID %in% rownames(ordDF[ordDF$PC1 > -0.05,]), "ord"] <- "group2"
sample_data(phy) <- sam
meta4 <- sample_data(phy)$ord
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], meta = meta, meta2 = meta2, meta3 = meta3, meta4 = meta4)
g <- ggplot(ordDF, aes(x = PC1, y = PC2, col = meta)) + 
	geom_point(alpha = 0.8, size = log(order(meta3))) +
	geom_text_repel(aes(label = meta2), size = 2, segment.size = 0.3, segment.alpha = 0.5) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	labs(col = "") +
	scale_color_brewer(palette = "Set2") +
	stat_ellipse(data = ordDF, aes(x = PC1, y = PC2, shape = meta4), inherit.aes = FALSE, alpha = 0.3) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("B3_ordination.png", g, scale = 0.8)




#save.image("start2.RData")
#BiocManager::install("DESeq2")
library(DESeq2)

phy <- subset_samples(phy1_decontam2, batch %in% "B3")
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- tax_glom(phy, taxrank = "Family")
sam <- sample_data(phy)
sam$ord <- rep(1, dim(sam)[1])
sam[sam$ID %in% rownames(ordDF[ordDF$PC1 < -0.05,]), "ord"] <- "group1"
sam[sam$ID %in% rownames(ordDF[ordDF$PC1 > -0.05,]), "ord"] <- "group2"
sample_data(phy) <- sam

ds <- phyloseq_to_deseq2(phy, ~ord)

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
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
head(sigtab)

phy2 <- transform_sample_counts(phy, function(x) { x/sum(x) })
phy2 <- subset_taxa(phy2, Family %in% sigtab$Family)
phy2 <- subset_taxa(phy2, taxa_sums(phy2) > 0)
phy2 <- subset_samples(phy2, sample_sums(phy2) > 0)
#phy2 <- merge_samples(phy2, group = "ord")
sam <- sample_data(phy2)
limits <- sam[order(sam$ord), "ID"]$ID
labels <-  paste(sam[order(sam$ord), "species.tree"]$species.tree, sam[order(sam$ord), "ord"]$ord, sep = " ")

g <- plot_bar(phy2, fill = "Family") +
	scale_fill_brewer(palette = "Set3") +
	scale_x_discrete(limits = limits, labels = labels) +
	theme_minimal() +
	coord_flip() +
	theme(
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1))

ggsave("B3_ord_deseq.png", g, scale = 0.8)
























#### Cophenetic distance 

phy <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B3"))
phy <- subset_samples(phy, !is.na(species.tree))
sam <- data.frame(sample_data(phy))
sam[sam$species.tree == "muta-NI", "species.tree"] <- "muta-SI" 
sample_data(phy) <- sam
sam <- sample_data(phy)


library(phytools)
library(ggtree)
library(pheatmap)
library(ggstance)

htree <- read.tree("raxml.tre")

htree <- reroot(htree, node = 70)

t <- htree

t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist( lapply( str_split(t$tip.lab, "_"), "[", 1) )

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))

t2 <- drop.tip(t, t$tip.label[!t$tip.label %in% sample_data(phy)$species.tree])







## Mantel tests 

phy <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B4"))
phy <- merge_samples(phy, group = "species.tree")
phy <- tax_glom(phy, taxrank = "Genus")
sam <- data.frame(sample_data(phy))
sam[sam$species.tree == "muta-NI", "species.tree"] <- "muta-SI" 
sample_data(phy) <- sam
sam <- sample_data(phy)
phy <- subset_samples(phy, rownames(sample_data(phy)) %in% t2$tip.label)


coph <- as.matrix(cophenetic(t2))
pdis_unifrac <- phyloseq::distance(phy, method="wunifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)

coph <- coph[match(rownames(pdis_unifrac), rownames(coph)),]
coph <- coph[, match(colnames(pdis_unifrac), colnames(coph))]
mantel(coph, pdis_unifrac)


ord <- ordinate(phy, method = "PCoA", distance = "wunifrac")

phylosig(t2, ord$vectors[,2], method = "K", test = TRUE)












phy <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B3"))
phy <- subset_samples(phy, !is.na(species.tree))
sam <- data.frame(sample_data(phy))
sam[sam$species.tree == "muta-NI", "species.tree"] <- "muta-SI" 
sample_data(phy) <- sam
sam <- sample_data(phy)


pdis_unifrac <- phyloseq::distance(phy, method="wunifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sam$code, sep = "_")

pdis_unifrac2 <- phyloseq::distance(phy, method="unifrac", type="samples")
pdis_unifrac2 <- as.matrix(pdis_unifrac2)
pdis_unifrac2[lower.tri(pdis_unifrac2, diag = TRUE)] <- "diag"
colnames(pdis_unifrac2) <- rownames(pdis_unifrac2) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sam$code, sep = "_")

ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
ptemp2 <- data.frame(cbind(pdis_unifrac2, ID = rownames(pdis_unifrac2)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)

pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])
pdis22 <- gather(ptemp2, key = id, value = pdis, 1:dim(sam)[1])

head(pdis2)
head(pdis22)

dis <- data.frame(pdis2, pdis2 = as.numeric(pdis22$pdis))

dis$elev_diff <- abs(as.numeric(unlist(lapply(str_split(dis$ID, "_"), "[", 2))) - as.numeric(unlist(lapply(str_split(dis$id, "_"), "[", 2))))
dis$habitat_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 3)) == unlist(lapply(str_split(dis$id, "_"), "[", 3)), yes = 0, no = 1)
dis$island_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 4)) == unlist(lapply(str_split(dis$id, "_"), "[", 4)), yes = 0, no = 1)
dis$batch_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 5)) == unlist(lapply(str_split(dis$id, "_"), "[", 5)), yes = 0, no = 1)
dis$speciesID <- unlist(lapply(str_split(dis$ID, "_"), "[", 6))
dis$speciesid <- unlist(lapply(str_split(dis$id, "_"), "[", 6))
dis$codeID <- unlist(lapply(str_split(dis$ID, "_"), "[", 7))
dis$codeid <- unlist(lapply(str_split(dis$id, "_"), "[", 7))
dis$species_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 6)) == unlist(lapply(str_split(dis$id, "_"), "[", 6)), yes = 0, no = 1)



coph <- cophenetic(t2)
coph[lower.tri(coph, diag = FALSE)] <- "diag"
coph <- data.frame(coph, speciesID = rownames(coph), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
coph <- gather(coph, key = speciesid, value = coph, 1:(dim(coph)[2]-1))

copho_vec <- c()
for(i in 1:dim(dis)[1]){
	
	dis[i, "speciesID"]
	dis[i, "speciesid"]
	
	a <- coph[coph$speciesID == dis[i, "speciesID"] & coph$speciesid == dis[i, "speciesid"],]
	
	if(dim(a)[1] > 0) {
		copho_vec <- c(copho_vec, a$coph)
	} else (
		copho_vec <- c(copho_vec, "NA")
	)
}

dis$copho <- copho_vec

dis <- dis[!dis$pdis %in% "diag",]
dis <- dis[!dis$copho %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis <- dis[dis$copho != "NA", ]
dis$copho <- as.numeric(dis$copho)


dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


temp_gen <- c()
for(i in 1:dim(dis)[1]) { 

	a <- as.character(unique(sam[sam$species.tree == dis[i, "speciesid"], "genus"]$genus))
	b <- as.character(unique(sam[sam$species.tree == dis[i, "speciesID"], "genus"]$genus))
	
	if(a == b) {
		a <- as.character(unique(sam[sam$species.tree == dis[i, "speciesid"], "genus"]$genus))
		temp_gen <- c(temp_gen, a)
	} else { 
		a <- paste(a, b, sep = " / ")
		temp_gen <- c(temp_gen, a)
	}
	
}

dis$genus <- temp_gen
#dis <- dis[!is.na(dis$copho2),]

dis[dis$genus %in% c("Kikihia / Rhodopsalta", "Rhodopsalta / Kikihia"), "genus"] <- "Kikihia / Rhodopsalta"
dis[dis$genus %in% c("Kikihia / Maoricicada", "Maoricicada / Kikihia"), "genus"] <- "Kikihia / Maoricicada"
dis[dis$genus %in% c("Rhodopsalta / Maoricicada", "Maoricicada / Rhodopsalta"), "genus"] <- "Rhodopsalta / Maoricicada"


dis$copho2 <- cut(dis$copho, breaks = quantile(dis$copho, prob = seq(0, 1, 0.1)), include.lowest = TRUE)







library(betareg)
library(lmtest)


summary, (betareg(pdis ~ copho , dis[dis$genus %in% c("Kikihia"), ]))

g1 <- ggplot(data = dis[dis$genus %in% c("Kikihia"), ], aes(x = copho2, y = pdis)) +
	#geom_boxplot(alpha = 0) + 
	geom_jitter(width = 0.05, alpha = 0.8) + 
	#geom_smooth(data =  dis[dis$genus %in% c("Kikihia", "Rhodopsalta"), ], aes(x = copho, y = pdis), method = "lm", se = FALSE, col = "grey", alpha = 0.5) +
	scale_color_brewer(palette = "Set2") +
	xlab("") +
	ylab("") +
	theme_bw() +
    theme(plot.margin = unit(c(1,5,1,5), "cm"), 
    	legend.position = "bottom", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    	axis.text.y = element_text(size = 10))


summary(betareg(pdis ~ copho , dis[!dis$genus %in% c("Kikihia", "Rhodopsalta", "Maoricicada"), ]))

g2 <- ggplot(data = dis[!dis$genus %in% c("Kikihia", "Rhodopsalta", "Maoricicada"), ], aes(x = copho2, y = pdis, col = genus)) +
	#geom_boxplot(data = dis[!dis$genus %in% c("Kikihia", "Rhodopsalta", "Maoricicada"), ], aes(x = copho2, y = pdis), alpha = 0) + 
	geom_jitter(width = 0.05, alpha = 0.8) + 
	#geom_smooth(data = dis[!dis$genus %in% c("Kikihia", "Rhodopsalta", "Maoricicada"), ], aes(x = copho, y = pdis), method = "lm", se = FALSE, col = "grey", alpha = 0.5) +
	scale_color_brewer(palette = "Set2") +
	xlab("") +
	ylab("") +
	theme_bw() +
    theme(plot.margin = unit(c(1,5,1,5), "cm"), 
    	legend.position = "bottom", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    	axis.text.y = element_text(size = 10))



ggsave("B3_copho1_dis.png", g1, scale = 0.8)
ggsave("B3_copho2_dis.png", g2, scale = 0.8)






## Habitat diff

temp_gen <- c()
for(i in 1:dim(dis)[1]) { 

	a <- as.character(unique(sam[sam$species.tree == dis[i, "speciesid"], "habitat"]$habitat))
	b <- as.character(unique(sam[sam$species.tree == dis[i, "speciesID"], "habitat"]$habitat))
	
	if(is.na(a) | is.na(b)) {
	
		temp_gen <- c(temp_gen, NA)
	
	} else {  
	
		if(a == b) {
			temp_gen <- c(temp_gen, a)
		} else { 
			a <- paste(a, b, sep = " / ")
			temp_gen <- c(temp_gen, a)
	}
	}
}

dis$habitat <- temp_gen
#dis <- dis[!is.na(dis$habitat),]




dis[dis$habitat %in% c("forest / grass", "grass / forest"), "habitat"] <- "forest / grass"
dis[dis$habitat %in% c("forest / shrub", "shrub / forest"), "habitat"] <- "forest / shrub"
dis[dis$habitat %in% c("grass / shrub", "shrub / grass"), "habitat"] <- "grass / shrub"

dis <- gather(dis, key = distance, value = pdis, 3:4)

kikidis <- dis[dis$genus %in% c("Kikihia"), ]

g <- ggplot(data = kikidis, aes(x = habitat, y = pdis,  col = distance)) +
	geom_boxplot(data = kikidis[kikidis$distance == "pdis", ],alpha = 0, col = "grey") + 
	geom_boxplot(data = kikidis[kikidis$distance == "pdis2", ], alpha = 0, col = "grey") + 
	geom_jitter(width = 0.1) + 
	#scale_color_brewer(palette = "Set2") +
	scale_x_discrete(limits = c("forest", "grass", "shrub", "forest / grass", "forest / shrub", "grass / shrub")) +
	xlab("") +
	ylab("") +
	coord_flip() +
	theme_bw() +
    theme(plot.margin = unit(c(1,5,1,5), "cm"), 
    	legend.position = "bottom", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    	axis.text.y = element_text(size = 10))

ggsave("B3_kikihia_habitat.png", g, scale = 0.8)





a <- dis[!dis$habitat %in% c("forest", "grass", "shrub") & dis$distance == "pdis", "pdis"]
b <- dis[dis$habitat %in% c("forest", "grass", "shrub") & dis$distance == "pdis", "pdis"]
c <- dis[!dis$habitat %in% c("forest", "grass", "shrub") & dis$distance == "pdis2", "pdis"]
d <- dis[dis$habitat %in% c("forest", "grass", "shrub") & dis$distance == "pdis2", "pdis"]

wilcox.test(a, b)
wilcox.test(c, d)
t.test(a, b)
t.test(c, d)


g <- ggplot() + 
	geom_density(aes(x = a[a > 0]), col = col_vector[2], size = 1) +
	geom_density(aes(x = b[b > 0]), col = col_vector[3], size = 1) +
	geom_density(aes(x = c[c > 0]), col = col_vector[2], size = 1) +
	geom_density(aes(x = d[d > 0]), col = col_vector[3], size = 1) +
	xlab("") +
	ylab("") +
	theme_bw() +
    theme( 
    	legend.position = "bottom", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    	axis.text.y = element_text(size = 10))

ggsave("B3_habitat_density.png", g, scale = 0.5)








### Philr for habitat 

#BiocManager::install("philr")
library(philr)

phy <- subset_samples(phy1_decontam2, batch %in% "B3")
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
plot_ordination(phy, gp.pcoa, color='species') + geom_point(size=4)
gp <- data.frame(gp.philr)
gp$habitat <- sample_data(metadata)$habitat
gp <- gp[!is.na(gp$habitat),]
gp2 <- gather(gp, key = node, value = balance, 1:(dim(gp)[2]-1))
gp3 <- gp2[gp2$balance < - 0.5 | gp2$balance > 0.5,]
balname <- apply(gp3, 1, function(x) {  name.balance(phy_tree(phy), tax_table(phy), x["node"])  })
order <- order(as.numeric(unlist(lapply(str_split(gp3$node, "n"), "[", 2))))
gp3$balname <- balname
g <- ggplot(data = gp3, aes(x = balname, y = balance, col = habitat)) +
	geom_point(alpha = 0.7, size = 3) +
	geom_point(shape = 1, size = 3, col = "black") +
	coord_flip() +
	xlab("") +
	ylab("") +
	scale_x_discrete(limits = unique(gp3$balname)) +
	#scale_color_manual(values = cols) +
	theme_bw() + 
    theme(
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(),
    	axis.text.y = element_text(size = 3)
    	)
    	
ggsave("B3_habitat_balance.png", g, scale = 1.2)










## DeSeq for habitat
library(DESeq2)

phy <- subset_samples(phy1_decontam2, batch %in% "B3")
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- tax_glom(phy, taxrank = "Family")
phy <- subset_samples(phy, !is.na(habitat))
sam <- sample_data(phy)

ds <- phyloseq_to_deseq2(phy, ~habitat)

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
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
head(sigtab)

phy2 <- transform_sample_counts(phy, function(x) { x/sum(x) })
phy2 <- subset_taxa(phy2, Family %in% sigtab$Family)
phy2 <- subset_taxa(phy2, taxa_sums(phy2) > 0)
phy2 <- subset_samples(phy2, sample_sums(phy2) > 0)
#phy2 <- merge_samples(phy2, group = "ord")
sam <- sample_data(phy2)
limits <- sam[order(sam$habitat), "ID"]$ID
labels <-  paste(sam[order(sam$habitat), "species.tree"]$species.tree, sam[order(sam$habitat), "habitat"]$habitat, sep = " ")

g <- plot_bar(phy2, fill = "Family") +
	scale_fill_brewer(palette = "Set3") +
	scale_x_discrete(limits = limits, labels = labels) +
	theme_minimal() +
	coord_flip() +
	theme(
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1))

ggsave("B3_habitat_deseq.png", g, scale = 0.8)















# Elevation difference 

diselev <- dis[dis$genus %in% c("Kikihia"), ]
diselev$elev_diff2 <- cut(diselev$elev_diff, breaks = quantile(diselev$elev_diff, prob = seq(0, 1, 0.2)), include.lowest = TRUE)

summary(betareg(pdis ~ elev_diff , data = diselev[diselev$distance == "pdis", ]))
summary(betareg(pdis ~ elev_diff , data = diselev[diselev$distance == "pdis2", ]))

g <- ggplot(data = diselev, aes(x = elev_diff2, y = pdis, col = distance)) +
	#geom_smooth(method = "lm", se = FALSE, col = "grey", alpha = 0.5) +
	#geom_text(aes(x = diselev$elev_diff, y = diselev$pdis, label = paste(diselev$speciesID, diselev$speciesid, sep = ":")), size = 1.5) +
	#geom_boxplot(alpha = 0) + 
	geom_jitter(width = 0.1) + 
	geom_smooth(method = "lm", se = FALSE, col = "black", alpha = 0.5, size = 1) +
	#scale_color_brewer(palette = "Set2") +
	xlab("") +
	ylab("") +
	theme_bw() +
    theme( 
    	legend.position = "bottom", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    	axis.text.y = element_text(size = 10))

ggsave("B3_elevation.png", g, scale = 0.5)











### Processing date and distances 

disproc <- dis

disproc$procID <- as.numeric( unlist(lapply(str_split(unlist(lapply(str_split(unlist(lapply(str_split(dis$ID, "_"), "[", 1)), "K"), "[", 2)), "g"), "[", 1))  )
disproc$procid <- as.numeric( unlist(lapply(str_split(unlist(lapply(str_split(unlist(lapply(str_split(dis$id, "_"), "[", 1)), "K"), "[", 2)), "g"), "[", 1))  )
disproc$disproc <- abs(disproc$procID - disproc$procid)
disproc <- disproc[!is.na(disproc$disproc),]
disproc$disproc2 <- cut(disproc$disproc, breaks = quantile(disproc$disproc, prob = seq(0, 1, 0.1)), include.lowest = TRUE)

summary(betareg(pdis ~ elev_diff , data = diselev[diselev$distance == "pdis", ]))


disproc <- disproc[disproc$genus %in% c("Kikihia"), ]

ggplot(data = disproc[disproc$distance == "pdis", ], aes(x = disproc, y = copho, col = distance)) +
	#geom_smooth(method = "lm", se = FALSE, col = "grey", alpha = 0.5) +
	#geom_text(aes(x = diselev$elev_diff, y = diselev$pdis, label = paste(diselev$speciesID, diselev$speciesid, sep = ":")), size = 1.5) +
	#geom_boxplot(alpha = 0) + 
	geom_point() + 
	#geom_smooth(method = "lm", se = FALSE, col = "black", alpha = 0.5, size = 1) +
	#scale_color_brewer(palette = "Set2") +
	xlab("") +
	ylab("") +
	theme_bw() +
    theme( 
    	legend.position = "bottom", 
    	panel.grid.major = element_blank(), 
    	panel.grid.minor = element_blank(), 
    	axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    	axis.text.y = element_text(size = 10))


unique(disproc$id)







### Are species more likely to be in different habitats 

sam$samcode <- unlist(lapply(str_split(sam$code, "[[.]]"), "[", 3))

stsp <- aggregate(sam$species.tree, by  = list(sam$samcode), function(x) { length( unique(x) )  })
hist(stsp[,2], breaks = 10)

































#### Abund prev by batch 

otu_b1 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B1"), taxrank = "Genus") ))
otu_b2 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B2"), taxrank = "Genus")   ))
otu_b3 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B3"), taxrank = "Genus") ))
otu_b4 <- data.frame(otu_table(  tax_glom(subset_samples(phy1, batch %in% "B4"), taxrank = "Genus")  ))

a <- base::apply(otu_b1, 1, function(x){ sum(x) })
b <- base::apply(otu_b2, 1, function(x){  sum(x) })
c <- base::apply(otu_b3, 1, function(x){  sum(x) })
d <- base::apply(otu_b4, 1, function(x){   sum(x) })

a2 <- base::apply(otu_b1, 1, function(x){ length(x[x > 0]) })
b2 <- base::apply(otu_b2, 1, function(x){ length(x[x > 0]) })
c2 <- base::apply(otu_b3, 1, function(x){ length(x[x > 0]) })
d2 <- base::apply(otu_b4, 1, function(x){ length(x[x > 0]) })

ref <- tax_glom(subset_samples(phy1, batch %in% "B4"), taxrank = "Genus")
tax <- tax_table(ref)

adata <- data.frame(abund = as.numeric(a/sum(a)), prev = as.numeric(a2)/dim(otu_b1)[2], batch = rep("B1", length(a)), tax = as.character(tax[,6]), class = as.character(tax[,3]))
bdata <- data.frame(abund = as.numeric(b/sum(b)), prev = as.numeric(b2)/dim(otu_b1)[2], batch = rep("B2", length(a)), tax = as.character(tax[,6]), class = as.character(tax[,3]) )
cdata <- data.frame(abund = as.numeric(c/sum(c)), prev = as.numeric(c2)/dim(otu_b1)[2], batch = rep("B3", length(a)), tax = as.character(tax[,6]), class = as.character(tax[,3]) )
ddata <- data.frame(abund = as.numeric(d/sum(d)), prev = as.numeric(d2)/dim(otu_b1)[2], batch = rep("B4", length(a)), tax = as.character(tax[,6]), class = as.character(tax[,3]) )

adata$tax <- unlist(lapply(str_split(adata$tax, "__"), "[", 2))
adata$class <- unlist(lapply(str_split(adata$class, "__"), "[", 2))

bdata$tax <- unlist(lapply(str_split(bdata$tax, "__"), "[", 2))
bdata$class <- unlist(lapply(str_split(bdata$class, "__"), "[", 2))

cdata$tax <- unlist(lapply(str_split(cdata$tax, "__"), "[", 2))
cdata$class <- unlist(lapply(str_split(cdata$class, "__"), "[", 2))

ddata$tax <- unlist(lapply(str_split(ddata$tax, "__"), "[", 2))
ddata$class <- unlist(lapply(str_split(ddata$class, "__"), "[", 2))

dd <- data.frame(rbind(adata, bdata, cdata, ddata))

abundprev_plot <- ggplot(dd, aes(x = prev, y = abund)) +
	geom_text_repel(data = dd[dd$prev > 0.25 | dd$abund > 0.015,], aes(label = tax), size = 2) +
	geom_point(size = 3) +
	geom_point(data = dd[dd$prev > 0.25 | dd$abund > 0.015,], aes(x = prev, y = abund, col = class), size = 3) +
	scale_color_brewer(palette = "Paired") +
	#facet_wrap(~batch, scale = "free") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
ggsave("abundprev_plot.png", abundprev_plot, scale = 0.7)

#tax[data.frame(tax[,3])[,1] == "D_2__vadinHA49",]








## looking by species 

phyg <- tax_glom(phy1t, taxrank = "Order")
phyb <- subset_samples(phyg, type %in% "gut")

ar2 <- c()

for(j in c("B1", "B2", "B3", "B4")) {

	phy <- subset_samples(phyb, batch %in% j)
	tax <- tax_table(phy)
	ar <- list()
	species <- unique(sample_data(phy)$species.tree)[!is.na(unique(sample_data(phy)$species.tree))]
	species <- species[species %in% names(table(sample_data(phy)$species.tree)[table(sample_data(phy)$species.tree) > 2])]

	for(i in species) {
		id <- i
		ar[[id]] <- data.frame(otu_table( subset_samples(phy, species.tree %in% i)))
		a <- base::apply(ar[[id]], 1, function(x){ sum(x) })
		a2 <- base::apply(ar[[id]], 1, function(x){ length(x[x > 0]) })
		adata <- data.frame(abund = as.numeric(a/sum(a)), prev = as.numeric(a2)/dim(ar[[id]])[2], batch = rep(j, length(a)), species = rep(i, length(a)), tax = as.character(tax[,4]), class = as.character(tax[,3]))
		adata$tax <- unlist(lapply(str_split(adata$tax, "__"), "[", 2))
		adata$class <- unlist(lapply(str_split(adata$class, "__"), "[", 2))
		ar[[id]] <- adata
	}

ar2 <- rbind(ar2, do.call("rbind", ar))

}

ar3 <- ar2


library(RColorBrewer)
n <- 12
qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) == 'Paired',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

cols <- c("Actinobacteria" = col_vector[1],
		  "Alphaproteobacteria" = col_vector[2],
		  "Bacilli" = col_vector[3],
		  "Babeliae" = col_vector[4],
		  "Bacteroidia" = col_vector[5],
		  "Clostridia" = col_vector[6],
		  "Deltaproteobacteria" = col_vector[7],
		  "Gammaproteobacteria" = col_vector[8],
		  "Mollicutes" = col_vector[9],
		  "Gemmatimonadetes" = col_vector[1],
		  "Negativicutes" = col_vector[11],
		  "Oxyphotobacteria" = col_vector[12]
		  )


ar2 <- ar3[ar3$batch == "B3",]

g1 <- ggplot(ar2, aes(x = prev, y = abund, shape = batch)) +
	geom_text_repel(data = ar2[ar2$prev > 0.45 | ar2$abund > 0.1,], aes(label = species, col = class), size = 2) +
	#geom_line(data = ar2[ar2$prev > 0.55 | ar2$abund > 0.5,], aes(group = tax, col = class)) +
	geom_point(size = 3, alpha = 0.5) +
	geom_point(data = ar2[ar2$prev > 0.45 | ar2$abund > 0.1,], aes(x = prev, y = abund, col = class), size = 3) +
	scale_color_manual(values = cols) +
	#facet_wrap(~batch, scale = "free") +
	theme_bw() +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
ar2 <- ar3[ar3$batch == "B4",]

g2 <- ggplot(ar2, aes(x = prev, y = abund, shape = batch)) +
	geom_text_repel(data = ar2[ar2$prev > 0.45 | ar2$abund > 0.1,], aes(label = species, col = class), size = 2) +
	#geom_line(data = ar2[ar2$prev > 0.55 | ar2$abund > 0.5,], aes(group = tax, col = class)) +
	geom_point(size = 3, alpha = 0.5) +
	geom_point(data = ar2[ar2$prev > 0.45 | ar2$abund > 0.1,], aes(x = prev, y = abund, col = class), size = 3) +
	scale_color_manual(values = cols) +
	#facet_wrap(~batch, scale = "free") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

abundprev_species <- ggarrange(g1, g2)

ggsave("abundprev_species.png", abundprev_species, scale = 0.8)
ggsave("abundprev_species_legend.png", g2, scale = 0.8)














# Ordinations

phy <- tax_glom(phy1t, taxrank = "Genus")
phy <- subset_samples(phy, type %in% "gut" & genus %in% c("Kikihia", "Maoricicada", "Rhodopsalta"))
phy <- subset_samples(phy, sample_sums(phy) > 0)
phy <- subset_samples(phy, batch %in% c("B3", "B4"))
phy <- subset_samples(phy, sample_sums(phy) > 0)

ord <- ordinate(phy, method = "PCoA", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], species = sample_data(phy)$genus, batch = sample_data(phy)$batch)
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = species, shape = batch)) + 
	geom_point(alpha = 0.8, size = 5) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	#stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op4 <- op

ggsave("ord.png", op4, scale = 0.7)


hisU3 <- phyloseq::distance(phy, method = "wunifrac")

adonis(hisU3 ~ species + batch, data = data.frame(sample_data(phy)))
permutest(betadisper(hisU3,  data.frame(sample_data(phy))$batch))
permutest(betadisper(hisU3,  data.frame(sample_data(phy))$species))




met <- data.frame(sample_data(phy), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)
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
nphy <- phyloseq(otu_table(notu, taxa_are_rows = TRUE), tax_table(phy1), phy_tree(phy1))

ord <- ordinate(nphy, method = "PCoA", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], nspe = nspe, sp = rownames(ord$vectors))
op <- ggplot(ordDF, aes(x = PC1, y = PC2, size = nspe)) + 
	geom_point(alpha = 0.8) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	#stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
op4 <- op




nspe_plot <- ggplot(ordDF, aes(x = PC1, y = nspe)) +
	geom_point(size = 4) +
	theme_bw() +
	xlab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text( hjust = 1)) 























#### Host phylogeny by chosen taxa in abund/prev plots
library(phytools)
library(ggtree)

htree <- read.tree("raxml.tre")

htree <- reroot(htree, node = 70)

t <- htree

t <- drop.tip(t, c("Amphipsalta_zelandica", "Kosemia_yezoensis", "Cicadettana_calliope", "phaeoptera"))
t$tip.label <-  unlist( lapply( str_split(t$tip.lab, "_"), "[", 1) )

t <- drop.tip(t, c("phaeoptera", "muta-NI-east", "westlandica-north-coast"))







### Host phylogeny with heatmap (based on abund prev from above)

phyb <- subset_samples(phyg, type %in% "gut")
phyb <- subset_taxa(phyb, taxa_sums(phyb) > 0)
phyb <- transform_sample_counts(phyb, function(x){ log(x+1) })
phyb <- subset_samples(phyb, batch %in% c("B3", "B4"))

table(sample_data(phyb)$species.tree)


dd <- c()
for(i in unique(sample_data(phyb)$species.tree)[!is.na(unique(sample_data(phyb)$species.tree))]) { 
	a <- subset_samples(phyb, species.tree %in% i)
	b <- data.frame(otu_table(a))
	c <- apply(b, 1, function(x){ length(x[x > log(2)])/dim(b)[2] })
	dd <- cbind(dd, c)
}

colnames(dd) <- unique(sample_data(phyb)$species.tree)[!is.na(unique(sample_data(phyb)$species.tree))]

gen <- ar2[ar2$abund > 0.01 | ar2$prev > 0.5, "tax"]
gen <- taxa_names(subset_taxa(phyb, Order %in% paste("D_3__", gen, sep = "")))

dd <- dd[rownames(dd) %in% gen, ]
dd <- data.frame(t(dd), check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE)

#colnames(dd) <- paste(seq(1, dim(dd)[2]), unlist(lapply(str_split(as.character(tax_table(phyb)[rownames(tax_table(phyb)) %in% colnames(dd), 6]), "__"), "[", 2)), sep = "_")
colnames(dd) <- unlist(lapply(str_split(as.character(tax_table(phyb)[rownames(tax_table(phyb)) %in% colnames(dd), 4]), "__"), "[", 2))

#sorted <- sort(apply(dd, 2, function(x){ length( x[is.finite(x) & !is.na(x)] ) }), decreasing = TRUE)
sorted <- sort(apply(dd, 2, function(x){ length(x[x > 0] ) }), decreasing = TRUE)
#sorted2 <- apply(dd, 2, function(x){ mean( x[is.finite(x) & !is.na(x)] ) })
sorted2 <- apply(dd, 2, function(x){ mean( x[ x > 0 ]) })
sorted2[!is.finite(sorted2)] <- 0
sorted2 <- sorted2[match(names(sorted), names(sorted2))]


dd <- dd[, match( names(sorted[sorted > 2]), colnames(dd))]

rm3 <- names(apply(dd, 1, function(x){ mean(x[ x > 0 ])} )[apply(dd, 1, function(x){ mean(x[ x > 0 ])} ) == 1])
rm1 <- names(apply(dd, 1, function(x){ mean(x[is.finite(x) & !is.na(x)]) } )[is.na(apply(dd, 1, function(x){ mean(x[is.finite(x) & !is.na(x)]) } ))])
rm2 <- t$tip.label[!t$tip.label %in% sample_data(phyb)$species.tree]

t2 <- drop.tip(t, c(rm1, rm2, rm3))

plot(t2)
nodelabels()

ghtree <- ggtree(t2) + 
	geom_tiplab(align = TRUE, size = 3) +
	geom_nodelab( nudge_x = 0.01, nudge_y = 0.001, size = 2.5) + 
	ylim(0,39) +
	geom_hilight(node=34, fill="blue", alpha=.1, extendto = 0.31) +
    geom_hilight(node=55, fill="green", alpha=.1, extendto = 0.31) +
    geom_hilight(node=56, fill="yellow", alpha=.1, extendto = 0.31)

#install.packages("pheatmap")
#library(pheatmap)

#heat_b3_logabund <- gheatmap(ghtree, dd, offset = 0.12, width=1.5, colnames=T, low = "white", high = "red", colnames_angle = 45, colnames_offset_x = -0.02, colnames_offset_y = -1, font.size = 2.3)
#heat_b4_logabund <- gheatmap(ghtree, dd, offset = 0.12, width=1.5, colnames=T, low = "white", high = "red", colnames_angle = 45, colnames_offset_x = -0.02, colnames_offset_y = -1, font.size = 2.3)
#heat_b3_prev <- gheatmap(ghtree, dd, offset = 0.12, width=1.5, colnames=T, low = "white", high = "red", colnames_angle = 45, colnames_offset_x = -0.02, colnames_offset_y = -1, font.size = 2.3)
#heat_b4_prev <- gheatmap(ghtree, dd, offset = 0.12, width=1.5, colnames=T, low = "white", high = "red", colnames_angle = 45, colnames_offset_x = -0.02, colnames_offset_y = -1, font.size = 2.3)

heat_all_prev <- gheatmap(ghtree, dd, offset = 0.1, width=1.5, colnames=T, low = "white", high = "red", colnames_angle = 55, colnames_offset_x = 0, colnames_offset_y = 0, font.size = 3, hjust = 0, colnames_position = "top", color = "lightgrey") 

#ggsave("heat_b3_logabund.png", heat_b3_logabund)
#ggsave("heat_b4_logabund.png", heat_b4_logabund)
#ggsave("heat_b3_prev.png", heat_b3_prev)
#ggsave("heat_b4_prev.png", heat_b4_prev)

ggsave("heat_all_prev.png", heat_all_prev, scale = 0.9)
























#### Modeling 

#install.packages("geosphere")
library(geosphere)
#install.packages("Imap")
library(Imap)

phy <- subset_samples(phy1t, type %in% "gut" & batch %in% c("B3", "B4"))
phy <- transform_sample_counts(phy, function(x){ log(x+1) })

sam <- sample_data(phy)
sam2 <- data.frame(id = sam$ID, lat = sam$lat_sign, lon = sam$lon_sign, row.names = rownames(sam))



## GEODIS

geodis <- GeoDistanceInMetresMatrix(sam2)
geodis[lower.tri(geodis, diag = TRUE)] <- "diag"
colnames(geodis) <- rownames(geodis) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")

gtemp <- data.frame(cbind(geodis, ID = rownames(geodis)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)

geodis2 <- gather(gtemp, key = id, value = geodis, 1:dim(sam)[1])



## PDIS

pdis_unifrac <- phyloseq::distance(phy, method="wunifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")

ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)

pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])
head(pdis2)


dis <- data.frame(pdis2, geodis = geodis2$geodis)

dis$elev_diff <- abs(as.numeric(unlist(lapply(str_split(dis$ID, "_"), "[", 2))) - as.numeric(unlist(lapply(str_split(dis$id, "_"), "[", 2))))
dis$habitat_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 3)) == unlist(lapply(str_split(dis$id, "_"), "[", 3)), yes = 0, no = 1)
dis$island_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 4)) == unlist(lapply(str_split(dis$id, "_"), "[", 4)), yes = 0, no = 1)
dis$batch_diff <- ifelse(unlist(lapply(str_split(dis$ID, "_"), "[", 5)) == unlist(lapply(str_split(dis$id, "_"), "[", 5)), yes = 0, no = 1)
dis$speciesID <- unlist(lapply(str_split(dis$ID, "_"), "[", 6))
dis$speciesid <- unlist(lapply(str_split(dis$id, "_"), "[", 6))



## COPH

coph <- cophenetic(t)
coph[lower.tri(coph, diag = TRUE)] <- "diag"
coph <- data.frame(coph, speciesID = rownames(coph), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
coph <- gather(coph, key = speciesid, value = coph, 1:(dim(coph)[2]-1))

copho_vec <- c()
for(i in 1:dim(dis)[1]){
	a <- coph[coph$speciesID == dis[i, "speciesID"] & coph$speciesid == dis[i, "speciesid"],]
	if(dim(a)[1] > 0) {
		copho_vec <- c(copho_vec, a$coph)
	} else (
		copho_vec <- c(copho_vec, "NA")
	)
}

dis$copho <- as.numeric(copho_vec)





dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]
#ggplot(data = dis, aes(x = copho, y = elev_diff)) +
#	geom_point()




sp <- unique(sam[sam$genus == "Kikihia", "species.tree"])$species.tree
dis2 <- dis[dis$speciesID %in% sp & dis$speciesid %in% sp,]


library(rjags)
library(R2jags)

sink("jags.R")
cat("
    model{

    # process model
    for (i in 1:N){
		pdis[i] ~ dbeta(alpha[i], beta[i])
  		alpha[i] <- mu[i] * phi
  		beta[i]  <- (1-mu[i]) * phi
  		logit(mu[i]) <- 
  			a + 
  			b*copho[i] + 
  			c*geodis[i] + 
  			d*elev_diff[i] + 
  			e*habitat_diff[i] + 
  			f*island_diff[i] + 
  			g*island_diff[i]*geodis[i] +
  			#h*elev_diff[i]*copho[i] + 
  			k*island_diff[i]*habitat_diff[i] + 
  			#k*elev_diff[i]*habitat_diff[i] +
  			bj[batch_block[i]]
	}
    
    # random effects 
    for(j in 1:2) {
       bj[j] ~ dnorm(0, taubj)
    }

    # priors
    phi ~ dgamma(.001,.001)
    taubj ~ dgamma(.001,.001)
    a ~ dnorm(0,.001)
	b ~ dnorm(0,.001)
	c ~ dnorm(0,.001)
	d ~ dnorm(0,.001)
	e ~ dnorm(0,.001)
	f ~ dnorm(0,.001)
	g ~ dnorm(0,.001)
	#h ~ dnorm(0,.001)
	k ~ dnorm(0,.001)
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  batch_block = dis2$batch_diff+1,
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff+1,
  island_diff = dis2$island_diff+1
  )

jags.model_kikihia <- jags(data = jags.data,
                   parameters.to.save = c("phi", "taubj", "a", "b", "c", "d","e","f","g","k"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 30000,
                   n.burnin = 1000)

traceplot(jags.model_kikihia)


out <- data.frame(jags.model_kikihia$BUGSoutput$summary)
out$par <- c("Intercept", "Phylogenic Distance", "Geographical Distance", "Elevation Difference", "Deviance", "Habitat Difference", "Island Difference", "Island x Geographical Distance", "Elevation x Phylogenetic Distance", "Variance", "Latent Factor Variance")

out <- out[rownames(out) %in% c("b","c","d","e","f"), ]
out2 <- out[rownames(out) %in% c("c","d"), ]

aplot <- ggplot(data = out, aes(x = par, y = mean)) +
	geom_point() +
	geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=0.1) +
	ylab("Posterior Mean Effect") +
	xlab("Variable") +
	geom_hline(yintercept = 0, col = "red") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


bplot <- ggplot(data = out2, aes(x = par, y = mean)) +
	geom_point() +
	geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=0.1) +
	ylab("Posterior Mean Effect") +
	xlab("Variable") +
	geom_hline(yintercept = 0, col = "red") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("aplot.png", aplot, scale = 0.4)
ggsave("bplot.png", bplot, scale = 0.4)













sp <- unique(sam[sam$genus == "Maoricicada", "species.tree"])$species.tree
dis2 <- dis[dis$speciesID %in% sp & dis$speciesid %in% sp,]


library(rjags)
library(R2jags)

sink("jags.R")
cat("
    model{

    # process model
    for (i in 1:N){
		pdis[i] ~ dbeta(alpha[i], beta[i])
  		alpha[i] <- mu[i] * phi
  		beta[i]  <- (1-mu[i]) * phi
  		logit(mu[i]) <- 
  			a + 
  			b*copho[i] + 
  			c*geodis[i] + 
  			d*elev_diff[i] + 
  			#e*habitat_diff[i] + 
  			f*island_diff[i] + 
  			g*island_diff[i]*geodis[i] +
  			h*elev_diff[i]*copho[i] + 
  			#k*elev_diff[i]*habitat_diff[i] +
  			bj[batch_block[i]]
	}
    
    # random effects 
    for(j in 1:2) {
       bj[j] ~ dnorm(0, taubj)
    }

    # priors
    phi ~ dgamma(.001,.001)
    taubj ~ dgamma(.001,.001)
    a ~ dnorm(0,.001)
	b ~ dnorm(0,.001)
	c ~ dnorm(0,.001)
	d ~ dnorm(0,.001)
	#e ~ dnorm(0,.001)
	f ~ dnorm(0,.001)
	g ~ dnorm(0,.001)
	h ~ dnorm(0,.001)
	#k ~ dnorm(0,.001)
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  batch_block = dis2$batch_diff+1,
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  #habitat_diff = dis2$habitat_diff,
  island_diff = dis2$island_diff
  )

jags.model_maoricicada <- jags(data = jags.data,
                   parameters.to.save = c("phi", "taubj", "a", "b", "c", "d","f","g","h"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 30000,
                   n.burnin = 1000)

traceplot(jags.model_maoricicada)








## diff

library(stringr)

library(RColorBrewer)
n <- 17
qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) == 'Paired',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

cols <- c("Actinobacteria" = col_vector[1],
		  "Alphaproteobacteria" = col_vector[2],
		  "Bacilli" = col_vector[3],
		  "Babeliae" = NA,
		  "Bacteroidia" = col_vector[5],
		  "Clostridia" = col_vector[6],
		  "Deltaproteobacteria" = col_vector[7],
		  "Gammaproteobacteria" = col_vector[8],
		  "Mollicutes" = col_vector[9],
		  "Gemmatimonadetes" = NA,
		  "Negativicutes" = col_vector[11],
		  "Oxyphotobacteria" = NA,
		  "Acidobacteriia" = NA,
		  "Deinococci" = NA,
		  "Fusobacteriia" = NA,
		  "Thermoleophilia" = NA,
		  "Verrucomicrobia" = NA
		  )

newphy <- subset_samples(phy1t, genus %in% "Kikihia" & batch %in% c("B3"))
#newphy <- transform_sample_counts(newphy, function(x){ log(x +1) })
newphy <- tax_glom(newphy, taxrank = "Genus")
newphy <- subset_samples(newphy, sample_sums(newphy) > 0)
newphy <- subset_samples(newphy, !is.na(habitat))
phy <- newphy
#phy <- merge_samples(newphy, "habitat")
for(i in 1:dim(tax_table(phy))[2]) {
	tax_table(phy)[,i] <- unlist(lapply(str_split(tax_table(phy)[,i], "__"), "[", 2))
}
#hist(log(taxa_sums(phy)), breaks = 100)
#nam <- names(log(taxa_sums(phy))[log(taxa_sums(phy)) > 6])
nam <- names(cols)
#phy <- prune_taxa(taxa = nam, phy)
phy <- subset_taxa(phy, Class %in% nam)
phy <- subset_samples(phy, sample_sums(phy) > 0)
phy <- transform_sample_counts(phy, function(x){ x/sum(x) })


b3habitat <- plot_bar(phy, x = "ID", fill = "Class") +
	facet_wrap(~habitat, scale = "free") +
	scale_fill_manual(values = cols) +
	theme_bw() +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))



newphy <- subset_samples(phy1t, genus %in% "Kikihia" & batch %in% c("B4"))
newphy <- tax_glom(newphy, taxrank = "Genus")
newphy <- subset_samples(newphy, sample_sums(newphy) > 0)
newphy <- subset_samples(newphy, !is.na(habitat))
phy <- merge_samples(newphy, "habitat")
hist(log(taxa_sums(phy)), breaks = 100)
nam <- names(log(taxa_sums(phy))[log(taxa_sums(phy)) > 6])
phy <- prune_taxa(taxa = nam, phy)
phy <- subset_samples(phy, sample_sums(phy) > 0)
phy <- transform_sample_counts(phy, function(x){ x/sum(x) })
for(i in 1:dim(tax_table(phy))[2]) {
	tax_table(phy)[,i] <- unlist(lapply(str_split(tax_table(phy)[,i], "__"), "[", 2))
}

b4habitat <- plot_bar(phy, x = "habitat", fill = "Class") +
	#facet_wrap(~habitat, scale = "free") +
	scale_fill_manual(values = cols) +
	theme_bw() +
	theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


ggarrange(b3habitat, b4habitat)

ggsave("b3habitat.png", b3habitat, scale = 0.6)
ggsave("b4habitat.png", b4habitat, scale = 0.6)









newphy <- subset_samples(phy1t, genus %in% "Kikihia" & batch %in% c("B4"))
newphy <- tax_glom(newphy, taxrank = "Genus")
newphy <- subset_samples(newphy, sample_sums(newphy) > 0)
newphy <- subset_samples(newphy, !is.na(habitat))
phy <- newphy
phy <- subset_samples(phy, sample_sums(phy) > 0)
#phy <- transform_sample_counts(phy, function(x){ x/sum(x) })


ord <- ordinate(phy, method = "PCoA", distance = "wunifrac")
ordDF <- data.frame(PC1 = ord$vectors[,1], PC2 = ord$vectors[,2], habitat = sample_data(phy)$habitat)
op <- ggplot(ordDF, aes(x = PC1, y = PC2, color = habitat)) + 
	geom_point(alpha = 0.8, size = 4) +
	xlab(paste("PC1", paste(round(100 * ord$values$Eigenvalues[1]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	ylab(paste("PC2", paste(round(100 * ord$values$Eigenvalues[2]/sum(ord$values$Eigenvalues), 1), "%", sep = " "), sep = ": ")) + 
	scale_color_brewer(palette = "Paired") +
	#stat_ellipse(type = "t") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



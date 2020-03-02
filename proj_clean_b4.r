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

phy <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B4"))

t2 <- drop.tip(t, t$tip.label[!t$tip.label %in% sample_data(phy)$species.tree])

g <- ggtree(t2) + 
	geom_tiplab(align = TRUE, size = 3) +
	geom_nodelab( nudge_x = 0.01, nudge_y = 0.001, size = 2.5) + 
	xlim(0,0.5)

ggsave("B4_tree.png", g)


gdata <- data.frame(g$data)
gdata <- gdata[order(gdata$y, decreasing = TRUE), ]
tiporder <- gdata[is.na(as.numeric(gdata$label)), "label"]
tiporder <- tiporder[tiporder != ""]







## Class barplot

phy1_decontam_reduced_rel <- transform_sample_counts(subset_samples(phy1_decontam2, type %in% "gut"), function(x) x/sum(x))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, batch %in% c("B4"))
phy1_decontam_reduced_rel <- subset_samples(phy1_decontam_reduced_rel, !is.na(species.tree))

otu <- data.frame(otu_table(phy1_decontam_reduced_rel), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE)
#otu$taxa <- paste(as.character(tax_table(phy1_decontam_reduced_rel)[,5]), as.character(tax_table(phy1_decontam_reduced_rel)[,6]))
otu$taxa <- as.character(tax_table(phy1_decontam_reduced_rel)[,6])

otu2 <- gather(otu, key = sample, value = abund, 1:(dim(otu)[2]-1))
a <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) length(x[x > 0]))
b <- a[a[,2] > 1, 1]
otu2 <- otu2[otu2$taxa %in% b, ]

c <- aggregate(otu2$abund, by = list(otu2$taxa), function(x) mean(x))
hist(log(c[,2]), breaks = 100)
d <- c[c[,2] > 0.0038, 1]
otu3 <- otu2[otu2$taxa %in% d, ]
otu3 <- otu3[otu3$abund > 0, ]

sd <- sample_data(phy1_decontam_reduced_rel)
sd[sd$species.tree == "muta", "species.tree"] <- "muta-SI"
sd$order <- match(sd$species.tree, tiporder)
limits = sd[order(sd$order, decreasing = TRUE), "ID"]$ID
labels = paste(sd[order(sd$order, decreasing = TRUE), "species.tree"]$species.tree, sd[order(sd$order,  decreasing = TRUE), "code"]$code, sep = " ; ")

g <- ggplot(data = data.frame(otu3), aes(x = sample, y = abund, fill = taxa)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector[-6]) +
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
    	
ggsave("B4_relative_abundance_genus.png", g)










### DeSeq, New Zealand and non-New Zealand















 
### Diveristy between species vs. among species

phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B4"))
sam <- sample_data(phy)
sp <- names(table(sam$species)[ table(sam$species) >= 2])
phy <- subset_samples(phy, species %in% sp)
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- subset_samples(phy, sample_sums(phy) > 0)
phy <- merge_samples(phy, "species")
phyotu <- otu_table(phy)

spdiv <- apply(phyotu, 1, function(x) { diversity(x, index = "shannon") })

div <- c()
for(i in 1:19) {
	a <- phyotu[sample(1:dim(phyotu)[1], 1),]
	b <- phyotu[sample(1:dim(phyotu)[1], 1),]
	
	d <- a+b
	div <- c(div, diversity(d, "simpson"))
}


divdat <- rbind(cbind(type = rep("within_species", length(spdiv)), div = spdiv), cbind(type = rep("between_species", length(div)), div = div))
divdat <- data.frame(divdat)
divdat$div <- as.numeric(divdat$div)
divdat$sp <- rownames(divdat)

g <- ggplot(divdat, aes(x = type, y = div)) +
	geom_boxplot() +
	#geom_text(data = divdat[divdat$type == "within_species",], aes(x = type, y = div, label = sp), size = 3) +
	theme_bw() 
	
ggsave("B4_div.png",g, scale = 0.5)
	
	
	
	
	
	
	
	
   

## B4 Map

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

sdata <- sample_data(phy1_decontam2)[sample_data(phy1_decontam2)$island %in% c("South Island", "North Island"),]
sdata <- sdata[sdata$batch %in% c("B4"), ]
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


cols <- c(Amphipsalta = "##FEFCFC", Kikihia = "##FEFCFC", Maoricicada = "##FEFCFC", nymph = "##FEFCFC")


g <- ggmap(mp) +
	#geom_text_repel(data = sdata, aes(x = lon_sign, y = lat_sign, label = species.tree, col = genus), size = 2.5, segment.size = 0.4, segment.alpha = 0.8, force = 30) +
	geom_point(data = sdata, aes(x = lon_sign, y = lat_sign), size = 1) +
	theme_minimal() +
	geom_text_repel(data = sdata[!duplicated(sdata$code2), ], aes(x = lon_sign, y = lat_sign, label = code2), size = 2.5, segment.size = 0.4, segment.alpha = 0.8)
		
ggsave("B4_map.png", g)




























#### B4 heatmap

library(pheatmap)

pheat <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B4"))
pheat <- transform_sample_counts(pheat, function(x) x/sum(x))
#pheat <- merge_samples(pheat, "species.tree")
asv_prev <- apply(otu_table(pheat), 1, function(x) length(x[x > 0]))
#asv <- names(asv_prev[asv_prev > 1])
#pheat <- prune_taxa(taxa = asv, pheat)


pheatotu <- otu_table(pheat)
abund_ranks <- apply(pheatotu, 1, rank)
hist(abund_ranks, breaks = 100)
abund_ranks <- abund_ranks - 10
abund_ranks[abund_ranks < 0] <- 1
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

ggsave("B4_heatmap.png", g)





















	 ## Ordinations
	 
phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% "B4")
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- tax_glom(phy, taxrank = "Genus")
phy <- subset_samples(phy, sample_sums(phy) > 0)
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

ggsave("B4_ordination.png", g, scale = 0.8)




#save.image("start2.RData")
#BiocManager::install("DESeq2")
library(DESeq2)

phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% "B4")
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
phy <- tax_glom(phy, taxrank = "Family")
phy <- subset_samples(phy, sample_sums(phy) > 0)
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
alpha = 0.00011
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

ggsave("B4_ord_deseq.png", g, scale = 0.8)
























#### Cophenetic distance 

phy <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B4"))
phy <- subset_samples(phy, !is.na(species.tree))
phy <- subset_samples(phy, sample_sums(phy2) > 0)
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
phy <- subset_samples(phy, rownames(sample_data(phy)) %in% t2$tip.label)
#phy <- tax_glom(phy, taxrank = "Genus")
phy <- subset_samples(phy, sample_sums(phy) > 0)
phy <- subset_taxa(phy, taxa_sums(phy) > 0)
sam <- data.frame(sample_data(phy))
sam[sam$species.tree == "muta-NI", "species.tree"] <- "muta-SI" 
sample_data(phy) <- sam
sam <- sample_data(phy)


coph <- as.matrix(cophenetic(t2))
pdis_unifrac <- phyloseq::distance(phy, method="wunifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)

coph <- coph[match(rownames(pdis_unifrac), rownames(coph)),]
coph <- coph[, match(colnames(pdis_unifrac), colnames(coph))]
mantel(coph, pdis_unifrac)


ord <- ordinate(phy, method = "PCoA", distance = "wunifrac")

phylosig(t2, ord$vectors[,5], method = "K", test = TRUE)












phy <- subset_samples(subset_samples(phy1_decontam2, type %in% "gut"), batch %in% c("B4"))
phy <- subset_samples(phy, !is.na(species.tree))
phy2 <- subset_samples(phy2, sample_sums(phy2) > 0)
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

dis[dis$genus %in% c("Kikihia / Maoricicada", "Maoricicada / Kikihia"), "genus"] <- "Kikihia / Maoricicada"


dis$copho2 <- cut(dis$copho, breaks = quantile(dis$copho, prob = seq(0, 1, 0.1)), include.lowest = TRUE)




library(betareg)
library(lmtest)


summary(betareg(pdis ~ copho , dis[dis$genus %in% c("Kikihia", "Maoricicada"), ]))

g1 <- ggplot(data = dis[dis$genus %in% c("Kikihia", "Maoricicada"), ], aes(x = copho2, y = pdis, col = genus)) +
	#geom_boxplot(alpha = 0) + 
	geom_jitter(width = 0.05, alpha = 0.8) + 
	#geom_smooth(data =  dis[dis$genus %in% c("Kikihia", "Rhodopsalta"), ], aes(x = copho, y = pdis), method = "lm", se = FALSE, col = "grey", alpha = 0.5) +
	scale_color_brewer(palette = "Set1") +
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



ggsave("B4_copho1_dis.png", g1, scale = 0.8)
ggsave("B4_copho2_dis.png", g2, scale = 0.8)






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

dis[dis$habitat %in% c("forest"), "habitat"]

dis <- gather(dis, key = distance, value = pdis, 3:4)

kikidis <- dis[dis$genus %in% c("Kikihia"), ]
kikidis <- dis[!is.na(dis$habitat), ]

g <- ggplot(data = kikidis, aes(x = habitat, y = pdis,  col = distance)) +
	geom_boxplot(data = kikidis[kikidis$distance == "pdis", ],alpha = 0, col = "grey") + 
	geom_boxplot(data = kikidis[kikidis$distance == "pdis2", ], alpha = 0, col = "grey") + 
	geom_jitter(width = 0.1) + 
	#scale_color_brewer(palette = "Set2") +
	scale_x_discrete(limits = c("grass", "shrub", "forest / grass", "forest / shrub", "grass / shrub")) +
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

ggsave("B4_kikihia_habitat.png", g, scale = 0.8)





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

ggsave("B4_habitat_density.png", g, scale = 0.5)






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
    	axis.text.y = element_text(size = 6)
    	)
    	
ggsave("hybrid_B2_balance.png", g)




















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

ggsave("B4_elevation.png", g, scale = 0.5)











### Are species more likely to be in different habitats 

sam$samcode <- unlist(lapply(str_split(sam$code, "[[.]]"), "[", 3))

stsp <- aggregate(sam$species.tree, by  = list(sam$samcode), function(x) { length( unique(x) )  })
hist(stsp[,2], breaks = 10)










### Processing date and distances 

disproc <- dis

disproc$procID <- as.numeric(unlist(lapply(str_split(unlist(lapply(str_split(dis$ID, "_"), "[", 1)), "P"), "[", 2)))
disproc$procid <- as.numeric(unlist(lapply(str_split(unlist(lapply(str_split(dis$id, "_"), "[", 1)), "P"), "[", 2)))

disproc$disproc <- abs(disproc$procID - disproc$procid)
disproc <- disproc[!is.na(disproc$disproc),]
disproc$disproc2 <- cut(disproc$disproc, breaks = quantile(disproc$disproc, prob = seq(0, 1, 0.1)), include.lowest = TRUE)

summary(betareg(pdis ~ elev_diff , data = diselev[diselev$distance == "pdis", ]))


disproc <- disproc[disproc$genus %in% c("Kikihia"), ]

ggplot(data = disproc, aes(x = disproc, y = copho, col = distance)) +
	#geom_smooth(method = "lm", se = FALSE, col = "grey", alpha = 0.5) +
	#geom_text(aes(x = diselev$elev_diff, y = diselev$pdis, label = paste(diselev$speciesID, diselev$speciesid, sep = ":")), size = 1.5) +
	#geom_boxplot(alpha = 0) + 
	geom_point() + 
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


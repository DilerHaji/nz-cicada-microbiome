######### Heatmap of raw everything at Genus level
p <- subset_samples(phy, type %in% "gut")
p <- tax_glom(p, taxrank = "Genus")

v <- apply(data.frame(otu_table(p)), 1, FUN = function(x){sum(x)})
v2 <- sort(v, decreasing = TRUE)
v3 <- names(v2[1:25])

p <- prune_taxa(taxa = v3, p)

p <- prune_samples(samples = names(sample_sums(p)[sample_sums(p) > 1]), p)

plot_heatmap(p, "NMDS", "jaccard", "species", "Genus")






#####################  heatmap tests

s <- sample_sums(subset_taxa(phy, Genus %in% c("Candidatus Sulcia")))
m <- sample_sums(subset_taxa(phy, Family %in% c("Mitochondria")))
c <- sample_sums(subset_taxa(phy, Order %in% c("Chloroplast")))

hist(s, breaks = 20)
sn <- names(s[s < 100000])
phy3 <- subset_samples(phy3, !ID %in% sn)

hist(m)
mn <- names(m[m > 20000])
phy3 <- subset_samples(phy3, !ID %in% mn)

hist(c)
cn <- names(c[c > 2000])
phy3 <- subset_samples(phy3, !ID %in% cn)



p <- subset_samples(phy3, type %in% "gut")
p <- tax_glom(p, taxrank = "Genus")
p <- logphy(p)

pc <- subset_samples(phy3, type %in% "control")
pc <- tax_glom(pc, taxrank = "Genus")
pc <- logphy(pc)



v <- apply(otu_table(p), 1, function(x){length(x[x > 0])})
v2 <- names(v[v > 10])


plot_heatmap(prune_taxa(taxa = v2, p), "PCoA", "jaccard", "species", "Family") + facet_wrap(~state, scale = "free")


plot_heatmap(prune_taxa(taxa = v2, p), "PCoA", "jaccard", "species", "Family") + facet_wrap(~genus, scale = "free_x")









# of all the ASVs across guts, are there more that show in both controls than just one control? 

prevx = as.numeric(apply(otu_table(p), 1, function(x){length(x[x > 0])}))

prevxc = as.numeric(apply(otu_table(pc), 1, function(x){length(x[x > 0])}))

plot(prevxc, prevx)


abund = as.numeric(apply(data.frame(otu_table(pc)), 1, FUN = function(x){sum(x)}))

tax <- tax_table(p)[,6]

cont <- data.frame(prevx, prevxc, abund, tax)
cont$Genus <- as.character(cont$Genus)

ggplot(data = cont, aes(x = prevxc, y = prevx)) +
	geom_text(aes(label = Genus),  size = 3, vjust = 1.5)







#####################  hierarchical testing

p <- phy3
p <-  subset_samples(p, type %in% "gut")
#p <- subset_taxa(p, !Genus %in% "Burkholderia-Caballeronia-Paraburkholderia")
p <- tax_glom(p, taxrank = "Genus")


pd <- phyloseq_to_deseq2(p, design = ~ species + species + island)
geo_mean_protected <- function(x) {
  if (all(x == 0)) {
    return (0)
  }
  exp(mean(log(x[x != 0])))
}
geoMeans <- apply(counts(pd), 1, geo_mean_protected)
pd <- estimateSizeFactors(pd, geoMeans = geoMeans)
pd <- estimateDispersions(pd)
pdab <- getVarianceStabilizedData(pd)
abund_sums <- rbind(data.frame(sum = colSums(pdab),
                               sample = colnames(pdab),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(logphy(p))),
                               sample = rownames(otu_table(logphy(p))),
                               type = "log(1 + x)"))
ggplot(abund_sums) +
  geom_histogram(aes(x = sum)) +
  facet_wrap(~type, scale = "free") +
  xlab("Total abundance within sample")


el <- phy_tree(logphy(p))$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(rownames(pdab), seq_len(phy_tree(logphy(p))$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, pdab, sample_data(logphy(p))$type)

hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)

val <- hfdr_res@p.vals
val <- val[val$adj.significance %in% "***",]

plot(hfdr_res)

rownames(val)

tax_table(prune_taxa(taxa = rownames(val), phy))

plot_bar(relphy(prune_taxa(taxa = rownames(val), phy3)), x = "species", fill = "Genus") + facet_wrap(~type, scales = "free") 









#########  adonis


phy_adonis <- tax_glom(phy3, taxrank = "Genus")
phy_adonis <- phy3

adonis2(t(otu_table(subset_samples(phy_adonis, type %in% "gut"))) ~ genus*species + island*state, data = data.frame(sample_data(subset_samples(phy, type %in% "gut"))), method = "jaccard")

adonis2(t(abs(pdab)) ~ genus*species + island, data = data.frame(sample_data(subset_samples(phy, type %in% "gut"))), method = "bray")

adonis2(t(otu_table(phy1)) ~ type + genus + species + island + locality + elevation2, data = data.frame(sample_data(phy_filtered)), method = "bray")

adonis2(t(otu_table(phy2)) ~ type + genus + species + island + locality + elevation2, data = data.frame(sample_data(phy_filtered)), method = "bray")

adonis2(t(otu_table(phy3)) ~ type  + locality + genus + species, data = data.frame(sample_data(phy)), method = "jaccard")

hist(log(as.numeric(otu_table(prune_samples(samples = "K10g", phy)))), breaks = 100)
hist(log(as.numeric(otu_table(prune_samples(samples = "K66g", phy)))), breaks = 100)








###### Prevelance across phylogeny 

taxon = "Kikihia"
phy_filtered_tglom <- subset_samples(phy3, type %in% "gut")
phy_filtered_tglom <- subset_samples(phy3, genus %in% taxon)
phy_filtered_tglom <- relphy(phy_filtered_tglom)
prevx = as.numeric(apply(otu_table(phy_filtered_tglom), 1, function(x){length(x[x > 0])}))
abund = as.numeric(apply(data.frame(otu_table(relphy(phy_filtered_tglom))), 1, FUN = function(x){mean(x, na.rm = TRUE)}))
p <- data.frame(tax_table(phy_filtered_tglom), prev = prevx, abund = as.numeric(as.character(abund)), stringsAsFactors = TRUE)
p2 <- p
n <- length(sample_names(phy_filtered_tglom))
k <- aggregate(p2$prev, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x/n)})
k$abund <- aggregate(p2$abund, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x)})[,2]
k <- k[k$x > 0,]
k$Phylum <- as.character(lapply(str_split(k[,1], "_"), "[", 1))
k$Genus <- as.character(lapply(str_split(k[,1], "_"), "[", 3))

k2 <- k[k$x > 0.1, 5]
k2 <- k2[-12]



## Taking out prev=1 taxa across guts 
p <- subset_samples(phy3, type %in% "gut")

p_merged <- merge_samples(p, "genus")
p_merged2 <- data.frame(t(otu_table(p_merged)))
p_merged2[p_merged2 > 0] <- 1
p_merged3 <- apply(p_merged2, 1, sum)
unitax <- names(p_merged3[p_merged3 == 1])

p <- prune_taxa(taxa = unitax, p)

temp <- apply(otu_table(p), 1, function(x){length(x[x > 0])})
tempnames <- names(temp[temp > 1])
p <- prune_taxa(taxa = tempnames, p)
spec <- data.frame(unique(sample_data(p)[, colnames(sample_data(p)) == "species"]))


pr <- data.frame(1:dim(tax_table(p))[1])
for(i in spec[,1]) {
a <- subset_samples(p, species %in% i)
b <- as.numeric(apply(otu_table(a), 1, function(x){length(x[x > 0])}))
b <- data.frame(b)
colnames(b) <- i
pr <- cbind(pr, b)
}
pr <- pr[, -1]

rownames(pr) <- paste(as.character(data.frame(tax_table(p)[,6])[,1]), seq(1, dim(tax_table(p))[1],1))

pr <- t(pr)
rownames(pr) <- spec[,1]
pr <- data.frame(pr)

tree_all <- read.nexus("Astral.tre")
to_keep <- c("I6903_02_NZ_HB_OCB_01_LEG_Hemiptera_Cicadidae_Kikihia_muta_east_seq1", "I16244_03NZGBMAR01_Hemiptera_Cicadidae_Rhodopsalta_cruentata_seq1", "I16257_02NZBRMRV02_Hemiptera_Cicadidae_Maoricicada_hamiltoni_seq1", "I27907_18_NZ_KA_KAI_3_LL11_Hemiptera_Cicadidae_Kikihia_paxillulae_seq1", "I16240_01NZTOMAN04_Hemiptera_Cicadidae_Maoricicada_campbelli_seq1", "I6885_08_NZ_SC_IDN_12_LEG_Hemiptera_Cicadidae_Kikihia_angusta_seq1", "I9904_08_NZ_NN_BES_17_Hemiptera_Cicadidae_Kikihia_tuta_seq1", "I6895_11_NZ_NN_LGU_02_LEG_Hemiptera_Cicadidae_Kikihia_flemingi_seq2", "I9972_MM0039_Hemiptera_Cicadidae_Kikihia_rosea_seq1", "I16255_02NZTOTAS03_Hemiptera_Cicadidae_Maoricicada_iolanthe_seq1", "I6893_99_NZ_99_14_LEG_Hemiptera_Cicadidae_Kikihia_cauta_seq1", "I16260_06NZHBPOR12_Hemiptera_Cicadidae_Rhodopsalta_leptomera_seq1", "I6903_02_NZ_HB_OCB_01_LEG_Hemiptera_Cicadidae_Kikihia_muta_east_seq1", "I6914_03_NZ_MC_BPN_02_LEG_Hemiptera_Cicadidae_Kikihia_penisularis_seq1", "I6897_14_NZ_NC_NIG_05_LEG_Hemiptera_Cicadidae_Kikihia_horologium_seq1", "I6901_08_NZ_OL_RSC_01_LEG_Hemiptera_Cicadidae_Kikihia_murihikua_seq1", "I16256_02NZMCLCR02_Hemiptera_Cicadidae_Rhodopsalta_microdora_seq1", "I27902_02_NZ_ND_MAV_09_LL6_Hemiptera_Cicadidae_Kikihia_cutora_seq1", "I6907_14_NZ_NN_JDH_01_LEG_Hemiptera_Cicadidae_Kikihia_nelsonensis_seq1", "I6849_11_NZ_WD_OKT_33_LEG_Hemiptera_Cicadidae_Kikihia_southwestlandica_seq1")
tree_all2 <- keep.tip(tree_all, to_keep)
tree_all2$tip.label <- c("campbelli", "iolanthe", "hamiltoni", "leptomera", "cruentata", "microdora", "muta", "nelsonsensis", "paxillulae", "tuta", "westlandica", "murihikua", "horologium", "rosea", "peninsularis", "angusta", "cutora", "flemingi", "cauta")

g <- ggtree(tree_all2) + geom_tiplab(size = 3) + geom_treescale()
gheatmap(g, pr, width=10, colnames=T, high = "red", low="white", offset = 0.01, colnames_angle = -90, font.size = 1.2) %>% scale_x_ggtree







######### Overlap patterns 


sam_names <- sample_names(subset_samples(phy_adonis, type %in% "gut"))

dis2 <- c()

for(j in sam_names) {
temp <- c()
for(i in sam_names) {
n <- as.numeric(otu_table(prune_samples(samples = j, phy)))
n2 <- as.numeric(otu_table(prune_samples(samples = i, phy)))
n[n > 0] <- 1
n2[n2 > 0] <- 1
n3 <- n + n2
nam <- rownames(otu_table(phy))
temp <- c(temp, as.numeric(length(nam[n3 > 1])))
}
dis2 <- cbind(dis2, temp)
}


colnames(sample_data(phy))
h <- sample_data(phy)[, c(2,7)]
h2 <- h[!is.na(match(h$ID, sam_names))]

d1 <- data.frame(dis2[, h2$genus == "Kikihia"])
d2 <- d1[h2$genus == "Kikihia",]

apply(d2, 2, mean)

d1 <- data.frame(dis2[, h2$genus == "Maoricicada"])
d2 <- d1[h2$genus == "Maoricicada",]

apply(d2, 2, mean)


d1 <- data.frame(dis2[, h2$genus == "Rhodopsalta"])
d2 <- d1[h2$genus == "Rhodopsalta",]

apply(d2, 2, mean)



d1 <- data.frame(dis2[, h2$genus == "Amphipsalta"])
d2 <- d1[h2$genus == "Amphipsalta",]

apply(d2, 2, mean)





####### dendrogram with barplots

phy3 <- phy
phy3 <- subset_taxa(phy3, !Genus %in% c("Candidatus Sulcia"))
phy3 <- subset_taxa(phy3, !Genus %in% c("Candidatus Hodgkinia"))
phy3 <- subset_taxa(phy3, !Family %in% c("Mitochondria"))
phy3 <- subset_taxa(phy3, !Order %in% c("Chloroplast"))
phy3 <- subset_taxa(phy3, !Family %in% c("Neoptera"))
phy3 <- subset_taxa(phy3, !Kingdom %in% c("Eukaryota"))
phy3 <- subset_taxa(phy3, !is.na(Phylum))

p <- subset_samples(phy3, type %in% "gut")
p <- tax_glom(p, "Genus")

nam <- sample_names(p)
top <- c()
for(i in nam){
a <- prune_samples(samples = i, p)
b <- taxa_sums(a)
c <- names(sort(b, decreasing = TRUE)[1:5])
top <- c(top, c)
}
top <- unique(top)


b <- taxa_sums(p)
c <- names(sort(b, decreasing = TRUE)[1:20])
p <- prune_taxa(taxa = c, p)

plot_bar(relphy(p))


n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

top_plot <- plot_bar(relphy(p), x = "ID", fill = "Genus") + 
           geom_bar(stat = "identity", color = "white") +
           scale_x_discrete(limits = row.names(sample_data(p)[order(data.frame(sample_data(p)[, colnames(sample_data(p)) == "genus"]), data.frame(sample_data(p)[, colnames(sample_data(p)) == "species"])), ]), labels = data.frame(sample_data(p)[order(data.frame(sample_data(p)[, colnames(sample_data(p)) == "genus"]), data.frame(sample_data(p)[, colnames(sample_data(p)) == "species"])), ])[, colnames(sample_data(p)) == "species"]) +
           scale_fill_manual(values = unname(col_vector))




ggsave("top_bar_plot.pdf", top_plot, scale = 0.9)

ggsave("astral_nuclear_tree.png", g, scale = 0.9)




















##################### random Forests 

library("randomForest")
library("plyr") 
library("rfUtilities") 
library("caret") 


# Gensus

scaled_ASV <- scale(data.frame(otu_table(subset_samples(phy3, type %in% "gut"))), center = TRUE, scale = TRUE)

scaled_ASV_state <- data.frame(t(scaled_ASV))
scaled_ASV_state$genus <- as.factor(as.character(unlist(metadata[rownames(scaled_ASV_state), "genus"])))
scaled_ASV_state <- scaled_ASV_state[-c(1:4),]

set.seed(151)

RF_state_classify <- randomForest(x=scaled_ASV_state[,1:(ncol(scaled_ASV_state)-1)] , y=scaled_ASV_state[ , ncol(scaled_ASV_state)] , ntree=10000, importance=TRUE, proximities=TRUE, mtry=1, strata = scaled_ASV_state[,1:(ncol(scaled_ASV_state)-1)])


imp <- RF_state_classify$importance
colnames(imp)


ggplot() + geom_histogram(aes(x = imp[,1]), binwidth = 0.00008) + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() + geom_histogram(aes(x = imp[,2]), binwidth = 0.00008) + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() + geom_histogram(aes(x = imp[,3]), binwidth = 0.00008) + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() + geom_histogram(aes(x = imp[,4]), binwidth = 0.00008) + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


imp[,1][imp[,1] > 0]
hist(imp[,2][imp[,2] > 0])
hist(imp[,3][imp[,3] > 0])

as.vector(tax_table(phy3)[which.max(importance(RF_state_classify)), c("Family", "Genus")])

varImpPlot(RF_state_classify, type = 2)
nam <- names(sort(imp[,6], decreasing = TRUE)[3])
nam <- sub("X", "", nam, fixed = TRUE)


p <- prune_taxa(taxa = nam, phy3)

#p <- prune_samples(samples =  names(sample_sums(p)[sample_sums(p) > 0]), p)

plot_bar(logphy(p), x = "ID", fill = "Genus") + facet_wrap(~genus, scale = "free_x")

hist(subset_sample(p, genus %in% "Kikihia")




ordgq <- ordinate(p, method = "PCoA", distance = "wunifrac")
plot_ordination(p, ordgq, color = "genus", axes = c(1,2)) +
  geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  geom_line(aes( group = species), alpha = 0.5) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())







fit_control <- trainControl( method = "LOOCV" )    
RF_state_classify_loocv <- train( scaled_ASV_state[,1:(ncol(scaled_ASV_state)-1)] , y=scaled_ASV_state[, ncol(scaled_ASV_state)] , method="rf", ntree=10000 , tuneGrid=data.frame( mtry=1 ) , trControl=fit_control )







# Species 

p2 <- subset_samples(phy3, !type %in% "control" & genus %in% "Kikihia")

spec <- names(table(sample_data(phy3)$species)[table(sample_data(p2)$species) > 1])[-1]

scaled_ASV <- scale(data.frame(otu_table(subset_samples(phy3, !type %in% "control" & genus %in% "Kikihia" & species %in% spec & !species %in% "hamiltoni"))), center = TRUE, scale = TRUE)

scaled_ASV_state <- data.frame(t(scaled_ASV))
scaled_ASV_state$species <- as.factor(as.character(unlist(metadata[rownames(scaled_ASV_state), "species"])))
scaled_ASV_state <- scaled_ASV_state[-1,]

set.seed(151)

RF_state_classify <- randomForest(x=scaled_ASV_state[,1:(ncol(scaled_ASV_state)-1)] , y=scaled_ASV_state[ , ncol(scaled_ASV_state)] , ntree=1000, importance=TRUE, proximities=TRUE, mtry=1 )


fit_control <- trainControl( method = "LOOCV" )    

RF_state_classify_loocv <- train( scaled_ASV_state[,1:(ncol(scaled_ASV_state)-1)] , y=scaled_ASV_state[, ncol(scaled_ASV_state)] , method="rf", ntree=1000 , tuneGrid=data.frame( mtry=4), trControl=fit_control )




# Type 

scaled_ASV <- scale(data.frame(otu_table(phy2), center = TRUE, scale = TRUE))

scaled_ASV_state <- data.frame(t(scaled_ASV))
scaled_ASV_state$type <- as.factor(as.character(unlist(metadata[rownames(scaled_ASV_state), "type"])))
scaled_ASV_state <- scaled_ASV_state[-4,]

set.seed(151)

RF_state_classify <- randomForest(x=scaled_ASV_state[,1:(ncol(scaled_ASV_state)-1)] , y=scaled_ASV_state[ , ncol(scaled_ASV_state)] , ntree=1000, importance=TRUE, proximities=TRUE, mtry=1 )


fit_control <- trainControl( method = "LOOCV" )    

RF_state_classify_loocv <- train( scaled_ASV_state[,1:(ncol(scaled_ASV_state)-1)] , y=scaled_ASV_state[, ncol(scaled_ASV_state)] , method="rf", ntree=1000 , tuneGrid=data.frame( mtry=4), trControl=fit_control )













##################### Cluster anlaysis on ordination - dendrogram by host tree 

comdend <- hclust(distance(phy, method = "unifrac", type  = "samples"))

m <- as.dist(cophenetic.phylo(phy_tree(tree)))
hostdend <- hclust(m)

comdend <- comdend %>% prune(c("MARSblank1", "17NzAkPkr01g", "17NzBrMat01g", "17NzBrMrv01g", "17NzWnNev21g", "C1", "C2", "K49gr1", "K49gr2")) 

matched_branches <- hostdend$labels[match(comdend$labels, hostdend$labels)]
matched_branches <- matched_branches[!is.na(matched_branches)]
unmatched_branches <- hostdend$labels[!hostdend$labels %in% matched_branches] 

hostdend <- hostdend %>% prune(unmatched_branches)

comdend <- comdend %>% prune(comdend$labels[!comdend$labels %in% hostdend$labels])

meta <- data.frame(sample_data(phy))[data.frame(sample_data(phy))[,1] %in% hostdend$labels,]

hostdend$labels <- paste(hostdend$labels, meta[match(hostdend$labels, meta[,1]), 8], sep = "_")
comdend$labels <- paste(comdend$labels, meta[match(comdend$labels, meta[,1]), 8], sep = "_")


n <- 55
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

str_split(hostdend$labels[hostdend$order], sep = "_")

g <- lapply(str_split(hostdend$labels[hostdend$order], "_"), "[", 2)
g2 <- as.character(meta[,7][match(g, meta[,9])])
g2[g2 %in% "Maoricicada"] <- col_vector[1]
g2[g2 %in% "Kikihia"] <- col_vector[2]
g2[g2 %in% "Rhodopsalta"] <- col_vector[3]
g2[g2 %in% "Amphipsalta"] <- col_vector[4]

dend <- dendlist(as.dendrogram(hostdend), as.dendrogram(comdend))

dend %>% untangle(method = "step1side") %>% tanglegram(lab.cex = 0.5, color_lines = g2, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)

dend %>% untangle(method = "step1side") %>% entanglement	

cor_bakers_gamma(dend)
set.seed(12345)

R <- 100
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- as.dendrogram(comdend)
for(i in 1:R) {
   dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
   cor_bakers_gamma_results[i] <- cor_bakers_gamma(as.dendrogram(comdend), dend_mixed)
}
round(sum(cor_bakers_gamma(dend) < cor_bakers_gamma_results)/ R, 4)















#####################  networks 
library(netassoc)  
library(Corbi)
library(netcom)
library(CINNA)
library(igraph)
library(hierformR)
library(rgr)


plot_network(make_network(phy2, type = "taxa", dist.fun="wunifrac", max.dist=0.9), phy, point_label = "Genus")

plot_net(phy2, type = "taxa", color = "Order", distance="jaccard", maxdist = 0.3, point_label = "Family") #+ scale_color_brewer(palette = "Paired")

plot_net(phy2, type = "taxa", distance="jaccard", maxdist = 0.2, color = "Phylum", point_label = "Genus") #+ scale_color_brewer(palette = "Paired")




# dataset phy1
phy1_m <- subset_samples(phy1, genus %in% c("Maoricicada"))
phy1_m <- prune_taxa(taxa = names(taxa_sums(phy1_m)[taxa_sums(phy1_m) > 0]), phy1_m)

phy1_k <- subset_samples(phy1, genus %in% c("Kikihia"))
phy1_k <- prune_taxa(taxa = names(taxa_sums(phy1_k)[taxa_sums(phy1_k) > 0]), phy1_k)

phy1_r <- subset_samples(phy1, genus %in% c("Rhodopsalta"))
phy1_r <- prune_taxa(taxa = names(taxa_sums(phy1_r)[taxa_sums(phy1_r) > 0]), phy1_r)

phy1_a <- subset_samples(phy1, genus %in% c("Amphipsalta"))
phy1_a <- prune_taxa(taxa = names(taxa_sums(phy1_a)[taxa_sums(phy1_a) > 0]), phy1_a)


phy_net1 <- make_netassoc_network(otu_table(phy1))
save.image("nzprojk-network.RData")
phy_net1_m <- make_netassoc_network(otu_table(phy1_m))
save.image("nzprojk-network.RData")
phy_net1_k <- make_netassoc_network(otu_table(phy1_k))
save.image("nzprojk-network.RData")
phy_net1_r <- make_netassoc_network(otu_table(phy1_r))
save.image("nzprojk-network.RData")
phy_net1_a <- make_netassoc_network(otu_table(phy1_a))
save.image("nzprojk-network.RData")



# dataset phy2
phy2_m <- subset_samples(phy2, genus %in% c("Maoricicada"))
phy2_m <- prune_taxa(taxa = names(taxa_sums(phy2_m)[taxa_sums(phy2_m) > 0]), phy2_m)

phy2_k <- subset_samples(phy2, genus %in% c("Kikihia"))
phy2_k <- prune_taxa(taxa = names(taxa_sums(phy2_k)[taxa_sums(phy2_k) > 0]), phy2_k)

phy2_r <- subset_samples(phy2, genus %in% c("Rhodopsalta"))
phy2_r <- prune_taxa(taxa = names(taxa_sums(phy2_r)[taxa_sums(phy2_r) > 0]), phy2_r)

phy2_a <- subset_samples(phy2, genus %in% c("Amphipsalta"))
phy2_a <- prune_taxa(taxa = names(taxa_sums(phy2_a)[taxa_sums(phy2_a) > 0]), phy2_a)

phy_net2 <- make_netassoc_network(otu_table(phy2))
save.image("nzprojk-network.RData")
phy_net2_m <- make_netassoc_network(otu_table(phy2_m))
save.image("nzprojk-network.RData")
phy_net2_k <- make_netassoc_network(otu_table(phy2_k))
save.image("nzprojk-network.RData")
phy_net2_r <- make_netassoc_network(otu_table(phy2_r))
save.image("nzprojk-network.RData")
phy_net2_a <- make_netassoc_network(otu_table(phy2_a))
save.image("nzprojk-network.RData")



# dataset phy3
phy3_m <- subset_samples(phy3, genus %in% c("Maoricicada"))
phy3_m <- prune_taxa(taxa = names(taxa_sums(phy3_m)[taxa_sums(phy3_m) > 0]), phy3_m)

phy3_k <- subset_samples(phy3, genus %in% c("Kikihia"))
phy3_k <- prune_taxa(taxa = names(taxa_sums(phy3_k)[taxa_sums(phy3_k) > 0]), phy3_k)

phy3_r <- subset_samples(phy3, genus %in% c("Rhodopsalta"))
phy3_r <- prune_taxa(taxa = names(taxa_sums(phy3_r)[taxa_sums(phy3_r) > 0]), phy3_r)

phy3_a <- subset_samples(phy3, genus %in% c("Amphipsalta"))
phy3_a <- prune_taxa(taxa = names(taxa_sums(phy3_a)[taxa_sums(phy3_a) > 0]), phy3_a)

phy_net3 <- make_netassoc_network(otu_table(phy3))
save.image("nzprojk-network.RData")
phy_net3_m <- make_netassoc_network(otu_table(phy3_m))
save.image("nzprojk-network.RData")
phy_net3_k <- make_netassoc_network(otu_table(phy3_k))
save.image("nzprojk-network.RData")
phy_net3_r <- make_netassoc_network(otu_table(phy3_r))
save.image("nzprojk-network.RData")
phy_net3_a <- make_netassoc_network(otu_table(phy3_a))
save.image("nzprojk-network.RData")


# dataset phy1 species 
phy_net_phy1 <- list()
for(i in unique(sample_data(phy1)$species)[-c(6,23)]) {
	p <- subset_samples(phy1, species %in% i)
	p <- prune_taxa(taxa = names(taxa_sums(p)[taxa_sums(p) > 0]), p)
	phy_net_phy1[[i]] <- make_netassoc_network(otu_table(p))
	save.image("nzprojk-network.RData")
}


# dataset phy2 species 
phy_net_phy2 <- list()
for(i in unique(sample_data(phy2)$species)[-c(6,23)]) {
	p <- subset_samples(phy2, species %in% i)
	p <- prune_taxa(taxa = names(taxa_sums(p)[taxa_sums(p) > 0]), p)
	phy_net_phy2[[i]] <- make_netassoc_network(otu_table(p))
	save.image("nzprojk-network.RData")
}


# dataset phy3 species 
phy_net_phy3 <- list()
for(i in unique(sample_data(phy3)$species)[-c(6,23)]) {
	p <- subset_samples(phy3, species %in% i)
	p <- prune_taxa(taxa = names(taxa_sums(p)[taxa_sums(p) > 0]), p)
	phy_net_phy3[[i]] <- make_netassoc_network(otu_table(p))
	save.image("nzprojk-network.RData")
}



plot_network(phy_net2$network_pos)


# Centrality degree 
net <- phy_net2_k$network_pos
e <- centr_degree(net)
s <- as_long_data_frame(net)
e2 <- data.frame(e$res, seq(1, length(e$res),1))
e3 <- e2[e2[,1] > 3, 1]
s2 <- as.character(s[s$from %in% e3, 4])
tax_table(prune_taxa(taxa = s2, phy))

# Centrality degree 
net <- phy_net2_k$network_pos
e <- closeness(net)
s2 <- names(e[order(e, decreasing = TRUE)][1:20])
tax_table(prune_taxa(taxa = s2, phy))



# community structure 
net <- phy_net2$network_pos
net2 <- phy_net1$network_pos


rownames(phy_net2$matrix_spsp_ses_thresholded)
V(phy_net2$network_pos)$name

remove.na(phy_net2$matrix_spsp_ses_thresholded)

m <- phy_net2_k$matrix_spsp_ses_thresholded
m[!rowSums(!is.finite(m)),]
m[!is.finite(m)] <- 0

plot_netassoc_matrix(phy_net2_k$matrix_spsp_ses_thresholded, colors = c("black","blue"))


plot_network(graph_from_adjacency_matrix(phy_net2_k$matrix_spsp_ses_all))


net <- phy_net2_k$network_pos
modularity(net, membership(walktrap.community(as.undirected(net))))
iso <- which(degree(net)==0)
net = delete.vertices(net, iso)
net_plot_k <- plot_network(net)


m <- walktrap.community(as.undirected(net))
member <- membership(m)
plot(m, as.undirected(net), vertex.label = member, vertex.size = 3, arrow.size = 0.1, arrow.width = 0.1, margin = -0.2)

d1 <- row.names(data.frame(tax_table(prune_taxa(taxa = names(member[member == 1]), phy))[,1]))
d2 <- row.names(data.frame(tax_table(prune_taxa(taxa = names(member[member == 11]), phy))[,1]))
d3 <- row.names(data.frame(tax_table(prune_taxa(taxa = names(member[member == 12]), phy))[,1]))
d4 <- row.names(data.frame(tax_table(prune_taxa(taxa = names(member[member == 13]), phy))[,1]))
d5 <- row.names(data.frame(tax_table(prune_taxa(taxa = names(member[member == 14]), phy))[,1]))
d6 <- row.names(data.frame(tax_table(prune_taxa(taxa = names(member[member == 15]), phy))[,1]))

q <- plot_bar(prune_taxa(taxa = d1, phy1_k), x = "clade")
q1 <- plot_bar(prune_taxa(taxa = d2, phy1_k), x = "clade")
q2 <- plot_bar(prune_taxa(taxa = d3, phy1_k), x = "clade")
q3 <- plot_bar(prune_taxa(taxa = d4, phy1_k), x = "clade")
q4 <- plot_bar(prune_taxa(taxa = d5, phy1_k), x = "clade")
q5 <- plot_bar(prune_taxa(taxa = d6, phy1_k), x = "clade")

plot_bar(prune_taxa(taxa = d1, logphy(phy1_k)), x = "ID") + facet_wrap(~species, scale = "free_x")

tax_table(prune_taxa(taxa = d1, logphy(phy1_k)))


ggarrange(q, q1, q2, q3, q4, q5)


net <- phy_net3_m$network_pos
modularity(net, membership(walktrap.community(as.undirected(net))))
iso <- which(degree(net)==0)
net = delete.vertices(net, iso)
net_plot_m <- plot_network(net)

net <- phy_net3_r$network_pos
modularity(net, membership(walktrap.community(as.undirected(net))))
iso <- which(degree(net)==0)
net = delete.vertices(net, iso)
net_plot_r <- plot_network(net)

net <- phy_net3_a$network_pos
modularity(net, membership(walktrap.community(as.undirected(net))))
iso <- which(degree(net)==0)
net = delete.vertices(net, iso)
net_plot_a <- plot_network(net)

net_plot_k
net_plot_m
net_plot_r
net_plot_a











##################### model





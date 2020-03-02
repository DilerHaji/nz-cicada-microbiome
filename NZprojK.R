############### Number of specimens with >10,000 ASVs per species

table(sample_data(subset_samples(samples = names(sample_sums(phy)[sample_sums(phy) > 10000]), phy))[, 9])




################ phy1 (Endosymbionts and orgenelles removed)

phy1 <- phy
phy1 <- subset_taxa(phy1, !Genus %in% c("Candidatus Sulcia"))
phy1 <- subset_taxa(phy1, !Family %in% c("Mitochondria"))
phy1 <- subset_taxa(phy1, !Order %in% c("Chloroplast"))
phy1 <- subset_taxa(phy1, !is.na(Phylum))
sample_data(phy1)$depth <- as.numeric(sample_sums(phy1))







################ phy2

phy_filtered <- phy
	
m <- as.dist(cophenetic.phylo(phy_tree(phy_filtered)))
m2 <- cutree(as.hclust(agnes(m, method = "single")), h = 0.03)
otu_new <- c()
for(i in levels(factor(m2))) {
	m3 <- names(m2[m2 == i])
	m4 <- prune_taxa(taxa = m3, phy_filtered)
	m5 <- t(data.frame(sample_sums(m4)))
	c <- paste(tax_table(m4)[,1], tax_table(m4)[,2], tax_table(m4)[,3], tax_table(m4)[,4], tax_table(m4)[,5], tax_table(m4)[,6], sep = "_")
	m6 <- table(c, useNA = "always")
	m7 <- names(m6[m6 == max(m6)])[1]
	m8 <- rownames(tax_table(m4)[c %in% m7])[1]
	rownames(m5) <- m8
	otu_new <- rbind(otu_new, m5)
}

otu_table(phy_filtered) <- otu_table(otu_new, taxa_are_rows = TRUE)

prevx = apply(X = otu_table(phy_filtered), MARGIN = ifelse(taxa_are_rows(phy_filtered), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
a <- names(prevx[prevx == 0])
prevx = apply(X = otu_table(subset_samples(phy_filtered, type %in% "control")), MARGIN = ifelse(taxa_are_rows(subset_samples(phy_filtered, type %in% "control")), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
a2 <- names(prevx[prevx > 1])
p <- unique(c(a, a2))
phy_filtered <- subset_taxa(phy_filtered, !taxa_names(phy_filtered) %in% p)

phy_filtered <- subset_taxa(phy_filtered, !Kingdom %in% c("Eukaryota"))
phy_filtered <- subset_taxa(phy_filtered, !Genus %in% c("Candidatus Sulcia"))
phy_filtered <- subset_taxa(phy_filtered, !Genus %in% c("Candidatus Hodgkinia"))
phy_filtered <- subset_taxa(phy_filtered, !Family %in% c("Mitochondria"))
phy_filtered <- subset_taxa(phy_filtered, !Order %in% c("Chloroplast"))
phy_filtered <- subset_taxa(phy_filtered, !is.na(Kingdom))
phy_filtered <- subset_taxa(phy_filtered, !is.na(Phylum))
#phy_filtered <- subset_taxa(phy_filtered, !is.na(Class))
#phy_filtered <- subset_taxa(phy_filtered, !is.na(Order))

n <- names(sample_sums(phy_filtered)[sample_sums(phy_filtered) > 99])
phy_filtered <- prune_samples(samples = n, phy_filtered)
sample_data(phy_filtered) <- cbind(sample_data(phy_filtered), depth = as.numeric(sample_sums(phy_filtered)))

phy2 <- phy_filtered

p <- phy2
n <- length(sample_names(p))
prevx = as.numeric(apply(otu_table(p), 1, function(x){length(x[x > 0])}))
abund = as.numeric(apply(data.frame(otu_table(p)), 1, FUN = function(x){sum(x)}))
pa <- data.frame(tax_table(p), prev = prevx, abund = abund, stringsAsFactors = TRUE)
k <- aggregate(pa$prev, list(paste(pa$Phylum, pa$Family, pa$Genus, sep = "_")), function(x){max(x/n)})
k$abund <- aggregate(pa$abund, list(paste(pa$Phylum, pa$Family, pa$Genus, sep = "_")), function(x){max(x)})[,2]
k$Phylum <- as.character(lapply(str_split(k[,1], "_"), "[", 1))
k$Family <- as.character(lapply(str_split(k[,1], "_"), "[", 2))
k$Genus <- as.character(lapply(str_split(k[,1], "_"), "[", 3))
g <- ggplot(data = k, aes(y = log(abund), x = x, col = Phylum)) +
	geom_point() +
	geom_text_repel(data = k[log(k$abund) < 3.7 | log(k$abund) > 9.7 | k$x > 0.25, ], aes(label = Genus, alpha = 0.1), segment.alpha = 0.15) +
	labs(x = "Mean Relative Prevalence", y = "Mean Relative Abundance") +
	scale_y_continuous(breaks = round(seq(min(log(abund)), max(log(abund)), by = 0.5),1)) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("prev_abund.pdf", g, scale = 0.8)



to.exclude <- names(taxa_sums(phy2)[log(taxa_sums(phy2)) < 4.2 ])
phy2 <- subset_taxa(phy2, !taxa_names(phy2) %in% to.exclude)











################ phy2.5 (Reduced phy_filtered data, top X taxa in each sample)

tred <- c()
for(i in sample_names(phy_filtered)){
  x <- prune_samples(samples = i, x = phy_filtered)
  y <- names(sort(taxa_sums(x), decreasing = TRUE)[1:2])
  tred <- c(tred, y)
}

phy2.5 <- prune_taxa(taxa = tred, x = phy2)
  
  
  
  
  
  
  
  
  
  
################ phy3 

phy3 <- phy
phy3 <- subset_taxa(phy3, !Genus %in% c("Candidatus Sulcia"))
phy3 <- subset_taxa(phy3, !Genus %in% c("Candidatus Hodgkinia"))
phy3 <- subset_taxa(phy3, !Family %in% c("Mitochondria"))
phy3 <- subset_taxa(phy3, !Order %in% c("Chloroplast"))
phy3 <- subset_taxa(phy3, !Family %in% c("Neoptera"))
phy3 <- subset_taxa(phy3, !Kingdom %in% c("Eukaryota"))
phy3 <- subset_taxa(phy3, !is.na(Phylum))
sample_data(phy3)$depth <- as.numeric(sample_sums(phy3))

p <- phy3
n <- length(sample_names(p))
prevx = as.numeric(apply(otu_table(p), 1, function(x){length(x[x > 0])}))
abund = as.numeric(apply(data.frame(otu_table(p)), 1, FUN = function(x){sum(x)}))
pa <- data.frame(tax_table(p), prev = prevx, abund = abund, stringsAsFactors = TRUE)
k <- aggregate(pa$prev, list(paste(pa$Phylum, pa$Family, pa$Genus, sep = "_")), function(x){max(x/n)})
k$abund <- aggregate(pa$abund, list(paste(pa$Phylum, pa$Family, pa$Genus, sep = "_")), function(x){max(x)})[,2]
k$Phylum <- as.character(lapply(str_split(k[,1], "_"), "[", 1))
k$Family <- as.character(lapply(str_split(k[,1], "_"), "[", 2))
k$Genus <- as.character(lapply(str_split(k[,1], "_"), "[", 3))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
g <- ggplot(data = k, aes(y = log(abund), x = x, col = Phylum)) +
	geom_point(size = 1) +
	geom_text_repel(data = k[log(k$abund) < -7.5 | log(k$abund) > 9.7 | k$x > 0.25 | (log(k$abund) < 6.7 & k$x > 0.06), ], aes(label = Family, alpha = 0.1), segment.alpha = 0.15) +
	labs(x = "Max Relative Prevalence", y = "Log Max Abundance") +
	scale_y_continuous(breaks = round(seq(min(log(abund)), max(log(abund)), by = 0.5),1)) +
	scale_color_manual(values = getPalette(23)) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


to.exclude <- names(taxa_sums(phy3)[log(taxa_sums(phy3)) < 5.2])
phylum_exclude <- k[k$x > 1/length(sample_names(phy3)),4]
family_exclude <- k[k$x > 1/length(sample_names(phy3)),5]
genera_exclude <- k[k$x > 1/length(sample_names(phy3)),6]

to.include <- unique(paste(phylum_exclude, family_exclude, genera_exclude))

tab <- dtable(phy3)
tab$tab <- paste(tab[,77], tab[,80], tab[,81])
to.include.asv <- rownames(tab[tab$tab %in% to.include, ])

phy3 <- subset_taxa(phy3, !taxa_names(phy3) %in% to.exclude)
phy3 <- prune_taxa(taxa = to.include.asv, phy3)









################ phy3.5 

tred <- c()
for(i in sample_names(phy3)){
  x <- prune_samples(samples = i, x = phy3)
  y <- names(sort(taxa_sums(x), decreasing = TRUE)[1:2])
  tred <- c(tred, y)
}

phy3.5 <- prune_taxa(taxa = tred, x = phy3)
  








################ Heatmap of phy2 

plot_heatmap(logphy(phy3.5), taxa.label = "Family", sample.label = "species", sample.order = "island")

plot_heatmap(logphy(phy2, taxa.label = "Family", sample.label = "species", sample.order = "island")








################ Blanks 

## MARS blank
m <- names(sort(taxa_sums(prune_samples(samples = "MARSblank1", phy)), decreasing = TRUE))
phy_marsblank <- prune_taxa(taxa = m, phy)
phy_marsblank_table <- dtable(phy_marsblank)
write.csv(phy_marsblank_table, "phy_marsblank_table")

## C1 blank
m <- names(sort(taxa_sums(prune_samples(samples = "MARSblank1", phy)), decreasing = TRUE))
phy_c1 <- prune_taxa(taxa = m, phy)
phy_c1_table <- dtable(phy_c1)
write.csv(phy_c1_table, "phy_c1_table")

## C2 blank
m <- names(sort(taxa_sums(prune_samples(samples = "MARSblank1", phy)), decreasing = TRUE))
phy_c2 <- prune_taxa(taxa = m, phy)
phy_c2_table <- dtable(phy_c2)
write.csv(phy_c2_table, "phy_c2_table")








################ Micrococcaceae diversity 

phy_microc <- subset_taxa(phy, Family %in% c("Micrococcaceae"))

tax_table(phy_microc)[,6] <- paste(tax_table(phy_microc)[,6], seq(1, length(tax_table(phy_microc)[,6]),1))

plot_bar(relphy(phy_microc), x = "species", fill = "Genus") +
	facet_wrap(~genus, scale = "free_x")







################  diversity 

phy_brad <- subset_taxa(phy_filtered, Family %in% c("Burkholderiaceae"))

tax_table(phy_brad)[,6] <- paste(tax_table(phy_brad)[,6], seq(1, length(tax_table(phy_brad)[,6]),1))

plot_bar(relphy(phy_brad), x = "species", fill = "Genus") +
	facet_wrap(~genus, scale = "free_x")

system("touch Burkholderia-Caballeronia-Paraburkholderia")
for(i in 1:length(taxa_names(phy_blank1))){
	n <- taxa_names(phy_blank1)[i]
	write.table(n, "n.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
	system("grep -f n.txt sequences.fasta -A 1 >> Burkholderia-Caballeronia-Paraburkholderia")
}






################  diversity 

phy_entero <- subset_taxa(phy, Family %in% c("Enterobacteriaceae"))

tax_table(phy_entero)[,6] <- paste(tax_table(phy_entero)[,6], seq(1, length(tax_table(phy_entero)[,6]),1))

plot_bar(relphy(phy_entero), x = "species", fill = "Genus") +
	facet_wrap(~genus, scale = "free_x")








################## ordinating

# FULL data 

ordgq <- ordinate(phy, method = "NMDS", distance = "wunifrac")
ord_filt <- plot_ordination(phy, ordgq, color = "genus", axes = c(1,2)) +
  geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt


ordgq <- ordinate(phy1, method = "NMDS", distance = "wunifrac")
ord_filt <- plot_ordination(phy1, ordgq, color = "genus", axes = c(1,2)) +
  geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt


ordgq <- ordinate(phy_filtered, method = "NMDS", distance = "wunifrac")
ord_filt <- plot_ordination(phy_filtered, ordgq, color = "genus", axes = c(1,2)) +
  geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt


ordgq <- ordinate(phy2, method = "NMDS", distance = "wunifrac")
ord_filt <- plot_ordination(phy2, ordgq, color = "genus", axes = c(1,2)) +
  geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt


ordgq <- ordinate(phy3, method = "NMDS", distance = "wunifrac")
ord_filt <- plot_ordination(phy3, ordgq, color = "genus", axes = c(1,2)) +
  geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  #stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt





p <- subset_taxa(phy3, Family %in% "Moraxellaceae")
s <- names(sample_sums(p)[sample_sums(p) > 0])
p2 <- prune_samples(samples = s, p)
ordgq <- ordinate(p2, method = "PCoA", distance = "wunifrac")
ord_filt <- plot_ordination(p2, ordgq, color = "species", axes = c(1,2)) +
  geom_text(aes(label = species_unique),  size = 3, vjust = 1.5) +
  geom_point(aes(size = depth_original)) +
  geom_point(aes(alpha = 0.5), alpha = 0.1) +
  stat_ellipse() +
  #labs(col = "Genus") + coord_fixed(sqrt(evals[2] / evals[1])) +
  #scale_color_brewer(palette = "Paired") + 
  theme_minimal() + theme(text = element_text(size = 16))
ord_filt

adonis2(t(otu_table(p2)) ~ type + genus + species + island + locality + elevation2, data = data.frame(sample_data(p2)), method = "bray")















################# barplots

typ = "gut"
tax = "Kikihia"

p_phy <- subset_samples(phy, type %in% typ & genus %in% tax)
p_phy1 <- subset_samples(phy1, type %in% typ & genus %in% tax)
p_phy2 <- subset_samples(phy2, type %in% typ & genus %in% tax)
p_phy3 <- subset_samples(phy3, type %in% typ & genus %in% tax)



# Everything phy 

plot_bar(relphy(p_phy), x = "ID", fill = "Phylum") + 
    geom_bar(stat = "identity") +
    #scale_fill_brewer(palette = "Paired") +
    #scale_x_discrete(limits = data.frame(phy@sam_data[,colnames(phy@sam_data) == "species"])[,1]) +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))


# Everything phy1

plot_bar(relphy(p_phy1), x = "ID", fill = "Phylum") + 
    geom_bar(stat = "identity") +
    #scale_fill_brewer(palette = "Paired") +
    #scale_x_discrete(limits = data.frame(phy@sam_data[,colnames(phy@sam_data) == "genus"])[,1]) +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))


# Everything phy1

plot_bar(relphy(p_phy3), x = "ID", fill = "Phylum") + 
    geom_bar(stat = "identity") +
    #scale_fill_brewer(palette = "Paired") +
    #scale_x_discrete(limits = data.frame(phy@sam_data[,colnames(phy@sam_data) == "genus"])[,1]) +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))




# only highly abundant bacteria 

p_up <- p_phy2

phy_his <- log(taxa_sums(p_up)+1)[log(taxa_sums(p_up)+1) > 0]
up <- quantile(phy_his, prob = seq(0, 1, 0.1), na.rm = TRUE)[10]
up_names <- names(phy_his[phy_his >= up])
p_up <- p_up %>% prune_taxa(taxa = up_names)
plot_bar(relphy(p_up), x = "species_unique", fill = "Family") + 
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired") +
    #scale_x_discrete(limits = data.frame(p_up@sam_data[,colnames(p_up@sam_data) == "species"])[,1]) +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))




# only highly prevelent bactera

p_prev <- p_phy2

prevx <- apply(otu_table(p_prev), 1, function(x){length(x[x > 0])})
prev_up <- quantile(prevx, prob = seq(0, 1, 0.1), na.rm = TRUE)[10]
prev_up_names <- names(prevx[prevx >= prev_up])
p_prev <- p_prev %>% prune_taxa(taxa = prev_up_names)
plot_bar(relphy(p_prev), x = "species_unique", fill = "Family") + 
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Paired") +
    #scale_x_discrete(limits = data.frame(p_prev@sam_data[,colnames(p_prev@sam_data) == "species"])[,1]) +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))



















################# Microbes Unique to sample 

length(sample_names(subset_samples(phy_filtered, genus %in% "Kikihia" & type %in% "gut")))

o <- subset_samples(phy_filtered, genus %in% "Kikihia" & type %in% "gut")

p <- apply(otu_table(o) > 0, 1, function(x){ length(x[x == TRUE])}) 

p[p > 10]

tax_table(prune_taxa(taxa = "ba0b38ba40ff42404964165ce5f0930a", phy_filtered))

k <- subset_samples(phy_filtered, genus %in% "Kikihia" & type %in% "gut")

m <- subset_samples(phy_filtered, genus %in% "Maoricicada" & type %in% "gut")

a <- subset_samples(phy_filtered, genus %in% "Amphipsalta" & type %in% "gut")

r <- subset_samples(phy_filtered, genus %in% "Rhodopsalta" & type %in% "gut")

c <- subset_samples(phy_filtered, type %in% "control")

k1 <- names(taxa_sums(k)[taxa_sums(k) > 0])
m1 <- names(taxa_sums(k)[taxa_sums(m) > 0])
a1 <- names(taxa_sums(k)[taxa_sums(a) > 0])
r1 <- names(taxa_sums(k)[taxa_sums(r) > 0])
c1 <- names(taxa_sums(k)[taxa_sums(c) > 0])

# Unique to Kikihia guts 
j <- k1[match(c(r1, m1, a1), k1)]
j2 <- j[!is.na(j)]
j3 <- k1[!k1 %in% j2]
length(as.character(tax_table(prune_taxa(taxa = j3, phy_filtered))[,6]))

# Unique to Maoricicada guts 
j <- m1[match(c(r1, k1, a1), m1)]
j2 <- j[!is.na(j)]
j3 <- m1[!m1 %in% j2]
length(as.character(tax_table(prune_taxa(taxa = j3, phy_filtered))[,6]))

# Unique to Amphipsalta guts 
j <- a1[match(c(r1, k1, m1), a1)]
j2 <- j[!is.na(j)]
j3 <- a1[!a1 %in% j2]
length(as.character(tax_table(prune_taxa(taxa = j3, phy_filtered))[,6]))

# Unique to Rhodopsalta guts 
j <- r1[match(c(a1, k1, m1), r1)]
j2 <- j[!is.na(j)]
j3 <- r1[!r1 %in% j2]
length(as.character(tax_table(prune_taxa(taxa = j3, phy_filtered))[,6]))

# Unique to controls
j <- c1[match(c(a1, k1, m1, r1), c1)]
j2 <- j[!is.na(j)]
j3 <- c1[!c1 %in% j2]
length(as.character(tax_table(prune_taxa(taxa = j3, phy_filtered))[,6]))






####################### Prevelance rank relationship 

colval <- c("Acidobacteria" = col_values[1],
                 "Actinobacteria" = col_values[2],
                 "Armatimonadetes" = col_values[3],
                 "Bacteroidetes" = col_values[28],
                 "Chlamydiae" = col_values[5],
                 "Chloroflexi"=col_values[6],
                 "Deinococcus-Thermus"=col_values[7],
                 "Euryarchaeota"=col_values[8],
                 "Firmicutes"=col_values[9],
                 "Gemmatimonadetes"=col_values[12],
                 "Patescibacteria"=col_values[11],
                 "Planctomycetes"=col_values[12],
                 "Proteobacteria"=col_values[13],
                 "Tenericutes"=col_values[14],
                 "Verrucomicrobia"=col_values[15],
                 "Cyanobacteria"=col_values[16],
                 "Deferribacteres"=col_values[17],
                 "Dependentiae"=col_values[18],
                 "Elusimicrobia"=col_values[19],
                 "Epsilonbacteraeota"=col_values[20],
                 "Fusobacteria"=col_values[21],
                 "Nitrospinae"=col_values[22],
                 "Patescibacteria"=col_values[23],
                 "Rokubacteria"=col_values[24],
                 "Spirochaetes"=col_values[25],
                 "BRC1"=col_values[26])

taxon = "Kikihia"
phy_filtered_tglom <- subset_samples(phy_filtered, type %in% "gut")
phy_filtered_tglom <- subset_samples(phy_filtered_tglom, genus %in% taxon)
#phy_filtered_tglom <- tax_glom(phy_filtered_tglom, "Genus")
phy_filtered_tglom <- relphy(phy_filtered_tglom)
prevx = as.numeric(apply(otu_table(phy_filtered_tglom), 1, function(x){length(x[x > 0])}))
abund = as.numeric(apply(data.frame(otu_table(relphy(phy_filtered_tglom))), 1, FUN = function(x){mean(x, na.rm = TRUE)}))
p <- data.frame(tax_table(phy_filtered_tglom), prev = prevx, abund = as.numeric(as.character(abund)), stringsAsFactors = TRUE)
#ggplot(aes(x = rank(abund, ties.method = "first"), y = log(abund), size = prev), data = p) +
#	geom_point(alpha = 0.2)  
p2 <- p
n <- length(sample_names(phy_filtered_tglom))
k <- aggregate(p2$prev, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x/n)})
k$abund <- aggregate(p2$abund, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x)})[,2]
k <- k[k$x > 0,]
k$Phylum <- as.character(lapply(str_split(k[,1], "_"), "[", 1))
k$Genus <- as.character(lapply(str_split(k[,1], "_"), "[", 3))
kiki <- ggplot(data = k, aes(y = abund, x = x, col = Phylum)) +
	geom_point() +
	geom_text_repel(data = k[k$abund > 0.01 & k$x > 0.1, ], aes(label = Genus, alpha = 0.1), segment.alpha = 0.15) +
	labs(x = "Mean Relative Prevalence", y = "Mean Relative Abundance") +
	scale_colour_manual(values = colval) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

taxon = "Maoricicada"
phy_filtered_tglom <- subset_samples(phy_filtered, type %in% "gut")
phy_filtered_tglom <- subset_samples(phy_filtered_tglom, genus %in% taxon)
#phy_filtered_tglom <- tax_glom(phy_filtered_tglom, "Genus")
phy_filtered_tglom <- relphy(phy_filtered_tglom)
prevx = as.numeric(apply(otu_table(phy_filtered_tglom), 1, function(x){length(x[x > 0])}))
abund = as.numeric(apply(data.frame(otu_table(relphy(phy_filtered_tglom))), 1, FUN = function(x){mean(x, na.rm = TRUE)}))
p <- data.frame(tax_table(phy_filtered_tglom), prev = prevx, abund = as.numeric(as.character(abund)), stringsAsFactors = TRUE)
#ggplot(aes(x = rank(abund, ties.method = "first"), y = log(abund), size = prev), data = p) +
#	geom_point(alpha = 0.2)  
p2 <- p
n <- length(sample_names(phy_filtered_tglom))
k <- aggregate(p2$prev, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x/n)})
k$abund <- aggregate(p2$abund, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x)})[,2]
k <- k[k$x > 0,]
k$Phylum <- as.character(lapply(str_split(k[,1], "_"), "[", 1))
k$Genus <- as.character(lapply(str_split(k[,1], "_"), "[", 3))
Maori <- ggplot(data = k, aes(y = abund, x = x, col = Phylum)) +
	geom_point() +
	geom_text_repel(data = k[k$abund > 0.02 | k$x > 0.2, ], aes(label = Genus, alpha = 0.1), segment.alpha = 0.15) +
	labs(x = "Mean Relative Prevalence", y = "Mean Relative Abundance") +
	scale_colour_manual(values = colval) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

taxon = "Rhodopsalta"
phy_filtered_tglom <- subset_samples(phy_filtered, type %in% "gut")
phy_filtered_tglom <- subset_samples(phy_filtered_tglom, genus %in% taxon)
#phy_filtered_tglom <- tax_glom(phy_filtered_tglom, "Genus")
phy_filtered_tglom <- relphy(phy_filtered_tglom)
prevx = as.numeric(apply(otu_table(phy_filtered_tglom), 1, function(x){length(x[x > 0])}))
abund = as.numeric(apply(data.frame(otu_table(relphy(phy_filtered_tglom))), 1, FUN = function(x){mean(x, na.rm = TRUE)}))
p <- data.frame(tax_table(phy_filtered_tglom), prev = prevx, abund = as.numeric(as.character(abund)), stringsAsFactors = TRUE)
#ggplot(aes(x = rank(abund, ties.method = "first"), y = log(abund), size = prev), data = p) +
#	geom_point(alpha = 0.2)  
p2 <- p
n <- length(sample_names(phy_filtered_tglom))
k <- aggregate(p2$prev, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x/n)})
k$abund <- aggregate(p2$abund, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x)})[,2]
k <- k[k$x > 0,]
k$Phylum <- as.character(lapply(str_split(k[,1], "_"), "[", 1))
k$Genus <- as.character(lapply(str_split(k[,1], "_"), "[", 3))
Rhodo <- ggplot(data = k, aes(y = abund, x = x, col = Phylum)) +
	geom_point() +
	geom_text_repel(data = k[k$abund > 0.02 | k$x > 0.2, ], aes(label = Genus, alpha = 0.1), segment.alpha = 0.15) +
	labs(x = "Mean Relative Prevalence", y = "Mean Relative Abundance") +
	scale_colour_manual(values = colval) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

taxon = "Amphipsalta"pbtrensim
phy_filtered_tglom <- subset_samples(phy_filtered, type %in% "gut")
phy_filtered_tglom <- subset_samples(phy_filtered_tglom, genus %in% taxon)
#phy_filtered_tglom <- tax_glom(phy_filtered_tglom, "Genus")
phy_filtered_tglom <- relphy(phy_filtered_tglom)
prevx = as.numeric(apply(otu_table(phy_filtered_tglom), 1, function(x){length(x[x > 0])}))
abund = as.numeric(apply(data.frame(otu_table(relphy(phy_filtered_tglom))), 1, FUN = function(x){mean(x, na.rm = TRUE)}))
p <- data.frame(tax_table(phy_filtered_tglom), prev = prevx, abund = as.numeric(as.character(abund)), stringsAsFactors = TRUE)
#ggplot(aes(x = rank(abund, ties.method = "first"), y = log(abund), size = prev), data = p) +
#	geom_point(alpha = 0.2)  
p2 <- p
n <- length(sample_names(phy_filtered_tglom))
k <- aggregate(p2$prev, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x/n)})
k$abund <- aggregate(p2$abund, list(paste(p2$Phylum, p2$Family, p2$Genus, sep = "_")), function(x){mean(x)})[,2]
k <- k[k$x > 0,]
k$Phylum <- as.character(lapply(str_split(k[,1], "_"), "[", 1))
k$Genus <- as.character(lapply(str_split(k[,1], "_"), "[", 3))
Amphi <- ggplot(data = k, aes(y = abund, x = x, col = Phylum)) +
	geom_point() +
	geom_text_repel(data = k[k$abund > 0.02 | k$x > 0.2, ], aes(label = Genus, alpha = 0.1), segment.alpha = 0.15) +
	labs(x = "Mean Relative Prevalence", y = "Mean Relative Abundance") +
	scale_colour_manual(values = colval) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
					
ggsave("kiki.pdf", kiki, scale = 0.8)
ggsave("maori.pdf", Maori, scale = 0.8)
ggsave("rhodo.pdf", Rhodo, scale = 0.8)
ggsave("amphi.pdf", Amphi, scale = 0.8)






################## Host genetic divergence and prev/abund 

tree <- read.nexus("rooted_co1_tree")

tree$tip.label <- as.character(lapply(str_split(tree$tip.label, "_"), "[", 1))

tree$tip.label <- c("Maoricicada1", "Maoricicada2", "Maoricicada3", "Maoricicada4", "Maoricicada5", "Maoricicada6",
"Maoricicada7", "Maoricicada8", "Maoricicada9", "Maoricicada10",
"Maoricicada11", "Maoricicada12", "Maoricicada13", "Maoricicada14", "Maoricicada15",
"Maoricicada16", "Maoricicada17", "K19g", "Maoricicada18", "Maoricicada19",
"K11g", "K34g", "K10g", "Maoricicada21", "Maoricicada22",
"Maoricicada23", "Maoricicada24", "Maoricicada25", "Maoricicada26", "Maoricicada27",
"Maoricicada28", "Maoricicada29", "Maoricicada30", "Maoricicada31", "Maoricicada32",
"Maoricicada33", "Maoricicada34", "Maoricicada35", "K43r", "K44r",        
"K44g", "Maoricicada36", "K43g", "Maoricicada37", "Maoricicada38",
"Maoricicada39", "Maoricicada40", "Maoricicada41", "Maoricicada42", "Maoricicada43",
"Maoricicada44", "K46g", "K46r", "Notopsalta1",  "Amphipsalta1",
"K21g", "K56g", "K57g", "Kikihia1", "Kikihia2",    
"K47g", "K57r", "K20g", "Kikihia3", "K64g1",        
"Kikihia4", "K62g0", "K48r", "Kikihia5", "Kikihia6",    
"K62r", "K64g2", "K62g", "K64g", "Kikihia7",    
"Kikihia8", "Kikihia9", "Kikihia10", "Kikihia11", "K31g",      
"K45r", "K49g2", "K45g", "K56r", "K13g",      
"Kikihia12", "Kikihia13", "Kikihia14", "K58g1", "K53r",      
"K51g0", "K53g", "K52r", "K52g", "K51g",      
"Kikihia15", "Kikihia16", "Kikihia17", "Kikihia18", "K59g1",     
"K58g", "K59r", "K17g", "K17r", "Kikihia19",   
"K60r", "K60g1", "K60g", "K59g", "K55r",      
"K50g", "K30g", "K49g", "Kikihia20", "K66g",      
"K48g", "K28g", "K22g", "K41g", "K41r",      
"Kikihia21", "Kikihia22", "Kikihia23", "K36g", "K39g",      
"K6g", "K65r", "K27g", "K65g1", "K26g",      
"Kikihia24", "K65g", "Kikihia25", "Kikihia26", "K16g",        
"Kikihia27", "K54g", "Kikihia28", "Kikihia29", "K18g",       
"Kikihia30", "K35g", "Kikihia31", "K32g", "K33g",      
"Rhodopsalta1", "Rhodopsalta2", "Rhodopsalta3", "K42g", "Rhodopsalta4",
"K55g", "K50r", "K61r", "K61g", "K23g",        
"K37g", "K38g", "K15g", "K54r", "K40g",       
"K9g", "K24g")

m <- as.dist(cophenetic.phylo(phy_tree(tree)))
o <- pcoa(m)
e <- sample_data(phy_filtered)$genus[match(names(o$vector[,2]), sample_data(phy_filtered)$ID)]
e2 <- sample_data(phy_filtered)$species[match(names(o$vector[,2]), sample_data(phy_filtered)$ID)]
oo <- data.frame(x = o$vectors[,1], y = o$vectors[,2], n = e, n2 = e2, id = names(o$vectors[,2]))

sample_data(phy_filtered)$dis <- oo$id[match(sample_names(phy_filtered), oo$id)]






#################### Phylogeny and heatmap of specific microbes 

phy_filtered_tglom2 <- relphy(phy)
phy_filtered_tglom2 <- subset_samples(phy_filtered_tglom2, genus %in% "Kikihia")
phy_filtered_tglom2 <- subset_taxa(phy_filtered_tglom2, type %in% "gut")
phy_filtered_tglom2 <- subset_taxa(phy_filtered_tglom2, Genus %in% "Pseudomonas")

s <- data.frame(otu_table(phy_filtered_tglom2), asv = rownames(otu_table(phy_filtered_tglom2)))

s2 <- gather(s, key = sample, value = rel_abund, 1:49)

s <- t(otu_table(phy_filtered_tglom2))
rownames(s) <- data.frame(sample_data(phy_filtered_tglom2)[, colnames(sample_data(phy_filtered_tglom2)) == "species_unique"])[,1]

p <- ggtree(phy_tree(phy_filtered_tglom2), ladderize = FALSE) +
  geom_tiplab(size = 0, linesize = 1) +
  theme(legend.position="right") + 
  scale_color_brewer(palette = "Paired") +
  geom_treescale()

heat <- phy_filtered_tglom2@otu_table[, c(order(sample_data(phy_filtered_tglom2)$clade))]
colnames(heat) <- paste(seq(1,49,1), sort(sample_data(phy_filtered_tglom2)$clade), sep = "_")

gheatmap(p, log(heat), width=10,
         colnames=T, high = "red", low="white", offset = 0.01, colnames_angle = -45, font.size = 1.2) %>% scale_x_ggtree















##################### Sample divergence by diversity 

p <- phy
m <- as.dist(cophenetic.phylo(phy_tree(tree)))
o <- metaMDS(m)
e <- sample_data(p)$genus[match(names(o$points[,2]), sample_data(p)$ID)]
e2 <- sample_data(p)$species[match(names(o$points[,1]), sample_data(p)$ID)]
oo <- data.frame(x = o$points[,1], y = o$points[,2], n = e, n2 = e2, id = names(o$points[,2]))
sample_data(p)$dis1 <- oo$x[match(sample_names(p), oo$id)]
sample_data(p)$dis2 <- oo$y[match(sample_names(p), oo$id)]
sample_data(p)$div <- as.numeric(apply(otu_table(p), 2, function(x){length(x[x > 0])}))
plot_ly(data.frame(sample_data(p)), x = ~dis1, y = ~dis2, z = ~div) %>%
  add_markers(color = ~genus) 



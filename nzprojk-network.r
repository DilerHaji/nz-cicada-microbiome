library(phyloseq)
library(netassoc)
library(igraph)
library(infotheo)

load("nzprojk.RData")


# dataset phy1
phy1_m <- subset_samples(phy1, genus %in% c("Maoricicada"))
phy1_k <- subset_samples(phy1, genus %in% c("Kikihia"))
phy1_r <- subset_samples(phy1, genus %in% c("Rhodopsalta"))
phy1_a <- subset_samples(phy1, genus %in% c("Amphipsalta"))

phy_net1 <- make_netassoc_network(otu_table(phy1))
phy_net1_m <- make_netassoc_network(otu_table(phy1_m))
phy_net1_k <- make_netassoc_network(otu_table(phy1_k))
phy_net1_r <- make_netassoc_network(otu_table(phy1_r))
phy_net1_a <- make_netassoc_network(otu_table(phy1_a))



# dataset phy2
phy2_m <- subset_samples(phy2, genus %in% c("Maoricicada"))
phy2_k <- subset_samples(phy2, genus %in% c("Kikihia"))
phy2_r <- subset_samples(phy2, genus %in% c("Rhodopsalta"))
phy2_a <- subset_samples(phy2, genus %in% c("Amphipsalta"))

phy_net2 <- make_netassoc_network(otu_table(phy2))
phy_net2_m <- make_netassoc_network(otu_table(phy2_m))
phy_net2_k <- make_netassoc_network(otu_table(phy2_k))
phy_net2_r <- make_netassoc_network(otu_table(phy2_r))
phy_net2_a <- make_netassoc_network(otu_table(phy2_a))



# dataset phy3
phy3_m <- subset_samples(phy3, genus %in% c("Maoricicada"))
phy3_k <- subset_samples(phy3, genus %in% c("Kikihia"))
phy3_r <- subset_samples(phy3, genus %in% c("Rhodopsalta"))
phy3_a <- subset_samples(phy3, genus %in% c("Amphipsalta"))

phy_net3 <- make_netassoc_network(otu_table(phy3))
phy_net3_m <- make_netassoc_network(otu_table(phy3_m))
phy_net3_k <- make_netassoc_network(otu_table(phy3_k))
phy_net3_r <- make_netassoc_network(otu_table(phy3_r))
phy_net3_a <- make_netassoc_network(otu_table(phy3_a))



# dataset phy2 species 
phy_net4 <- list()
for(i in unique(sample_data(phy2)$species)[-c(6,23)]) {
	p <- subset_samples(phy2, species %in% i)
	phy_net4[[i]] <- make_netassoc_network(otu_table(p))
}



# dataset phy3 species 
phy_net4 <- list()
for(i in unique(sample_data(phy3)$species)[-c(6,23)]) {
	p <- subset_samples(phy3, species %in% i)
	phy_net4[[i]] <- make_netassoc_network(otu_table(p))
}

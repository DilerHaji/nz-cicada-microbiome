################ exploring mitochondria taxa with blast ####

phy_na_rel <- transform_sample_counts(phy, fun = function(x){x/sum(x)})

tax_means <- function (x) {
  x <- otu_table(x)
  if (taxa_are_rows(x)) {
    apply(x, 1, mean)
  }
  else {
    colSums(x)
  }}

nas <- data.frame(tax_table(phy_na_rel))[,5] %in% c("Mitochondria")
na_names <- rownames(data.frame(tax_table(phy_na_rel))[nas,])

phy_na_rel <- prune_taxa(taxa = na_names, phy_na_rel)

sam_names <- sample_names(phy_na_rel)

sub_names_list <- c()
phy_na_rel_sub <- list(c())
phy_na_rel_sub_name <- list(c())
for(i in 1:length(sam_names)) {
  phy_na_rel_sub[[i]] <- prune_samples(samples = sam_names[i], phy_na_rel)
  phy_na_rel_sub_name[[i]] <- sort(taxa_sums(phy_na_rel_sub[[i]]), decreasing = TRUE)[1:20]
  sub_names_list <- c(sub_names_list, names(sort(taxa_sums(phy_na_rel_sub[[i]]), decreasing = TRUE)[1:3]))
}

unique(sub_names_list)

ophio <- c("6a8c6e7b5575384d84603190f553bbed", "8cbb3d6c3439f9f43ad468ffa7cea5d5",
           "70428f5db572bc38407de883b2d60481", "d551ba53df954e9371989d2e7d435925",
           "51bc85a3f61d8375610d307f87b4c7e6", "0de6ce1622bad3f7ea326eaf061fee18",
           "4870e0cc17e0a7fa57c8438f43e72de9", "12084e97e05b39a9f5a0a5ece30f7a23",
           "04852f6d726489aaf5c7360e9afa6886", "f74cccc440480bb27b15b25ca61e04d8",
           "4c79db3de07d0a9416e64fa05e1ea684", "297d4c52941e236f905079d96ffa11c3",
           "45d24f2938208821534b9fc5a4c2252d", "e3d076519e2a5c3110991b03df94d2fe",
           "29a1d2df2a01d5ddaa30366b3b47bca7", "c1669fe7f9afc4729dd7f36464936b5c",
           "5643f76d22effcf06afb578cf8c8a13a", "be9f5a81c0edd7df470e1c9b24085e68",
           "61116285c0d315b493b3557038794612", "42b99bd6334b35508b5b3d115c539575",
           "0d464d3224311bd64cf5ae2368f0f291")

fungi <- c("d0766aa99f5060dd684eef2f12ffadcf", "bcbff99e4eca9143631d8d5cdc44ce4b",
           "d1279c6f4ce6debf58f4cad02dd1ebe6", "11ed8a174829c04b89f8152f3b357baf",
           "fd4b6e74983a2fafd15663bcd7b4c974", "3fe21970aed637cfb38274f5f89e4acd",
           "9efdc204333fc3b138e4c51bdbef2c64", "f047662694a5ae885af90784d71f9167",
           "b059ea169cee3138470c80e90886fdd1", "a0c2c515d5e95694645d270b8ab05b63",
           "5f3c97585781647cd1ca790529d191da", "60d4c5d766ba0dcf80e047d5820fd39c",
           "b23c5d84713e015a6401e115991312cb", "3d44ca36bee41e0f02c4055225bee252",
           "abbc2c260abb2996184eaf3d2c09d222", "7f115bfccf6ef1dace471f8c906677c0",
           "a413fb0a3efe660811508eeb638719cc", "43a787353d93e2df7831f883aa51370f")

plant <- c("b78ea5a0dbf08f8e748570029046b05b","4cecc105157f1695480389254fc693df",
           "360bdb08d0a39af09a1df2ada7a00f44","396e2a0d8cdd152f635608a7fcd98f5d",
           "c4a851299273e0c0e069224517bbdf25", "66108f4dde01fe4491eb27fbbf478753",
           "15c32b08f96a511dd999c5b61fcfb235", "d47ac03963ba61e7270d352d682d71a0",
           "e52f5e207d4738e85db6c63f1ddbd9fb", "6b41b9ffd9ade8ea7f3e0836aba26407",
           "c680337d325a90de27cbf9355052e1b6")

new <- prune_taxa(taxa = ophio, phy)


ophio_sort <- names(sort(taxa_sums(new), decreasing = TRUE)[1:12])

new2 <- prune_taxa(taxa = ophio_sort, phy_na_rel)
new_ophio_clade <- new

other_fun <- prune_taxa(taxa = fungi, phy)

tax_table(new)[,5] <- seq(1,21,1)
tax_table(new2)[,5] <- seq(1,12,1)
tax_table(new_ophio_clade)[,5] <- c(1,1,1,1,1,1,4,4,1,2,2,2,2,2,2,2,2,3,3,3,3)

tax_table(new2)[,5]


barplot_family(new2, group = "species_unique", upper = 50) +
  facet_wrap(~species, scales = "free_x")



sam <- data.frame(sample_data(new2))
sam_nam_mix <- as.character(sam[sam$batch == c("ok"), 1])



new22 <- prune_samples(samples = sam_nam_mix, new2)

ophio_na <- tax_table(new)[!taxa_names(new) %in% taxa_names(new2),]
ophio_na[,5] <- NA
ophio_good <- tax_table(new)[taxa_names(new) %in% taxa_names(new2),]
oph_na_good <- rbind(ophio_na, ophio_good)

new23 <- new
tax_table(new23) <- oph_na_good

barplot_family(new23, group = "ID", upper = 50) +
  facet_wrap(~species, scales = "free_x")

tax_table(new)


barplot_genus(new2, group = "ID", upper = 50) +
  facet_wrap(~genus, scales = "free_x")


x <- prune_samples(samples = c("17NzWnNev21g"), phy_filtered)
sort(taxa_sums(x), decreasing = TRUE)[1:5]



ls
phynarel_tax_table <- data.frame(tax_table(phy_na_rel))

phynarel_tax_table[rownames(phynarel_tax_table) == "6ae2bc9c8f5ae62cd72209c8d8ed2a45",]


write.table(x = as.factor(row.names(phynarel_tax_table)), file = "na_names.txt", row.names = FALSE)



# trying to look at abundance 

txp <- as.numeric(dim(otu_table(new2))[1])
s <- ctm(phy = new2, taxdown = 1, taxup = txp)

coefs <- s[[2]]
keep_tax <- coefs[coefs$coef > 0, 1]
keep_tax <- as.character(keep_tax[!is.na(keep_tax)])

new3 <- prune_taxa(taxa = keep_tax, new2)
barplot_family(new3, group = "ID", upper = 50) +
  facet_wrap(~species, scales = "free_x")


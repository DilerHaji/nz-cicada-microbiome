library(geosphere)
library(Imap)

phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			(b + bj[library_block[i]])*copho[i] + 
  			(c + bj[library_block[i]])*geodis[i] + 
  			(d + bj[library_block[i]])*elev_diff[i] + 
  			(e + bj[library_block[i]])*habitat_diff[i] +
  			bj[library_block[i]]
		}
    
    # latent factors
    for(j in N2) {
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
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model1 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "taubj", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model1)
























phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
#phy <- transform_sample_counts(phy, function(x){ log(x+1) })
sam <- sample_data(phy)
sam2 <- data.frame(id = sam$ID, lat = sam$lat_sign, lon = sam$lon_sign, row.names = rownames(sam))

## GEODIS
geodis <- GeoDistanceInMetresMatrix(sam2)
geodis[lower.tri(geodis, diag = TRUE)] <- "diag"
colnames(geodis) <- rownames(geodis) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
gtemp <- data.frame(cbind(geodis, ID = rownames(geodis)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
geodis2 <- gather(gtemp, key = id, value = geodis, 1:dim(sam)[1])

## PDIS
pdis_unifrac <- phyloseq::distance(pa_phy(phy), method="unifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])

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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			(b + bj[library_block[i]])*copho[i] + 
  			(c + bj[library_block[i]])*geodis[i] + 
  			(d + bj[library_block[i]])*elev_diff[i] + 
  			(e + bj[library_block[i]])*habitat_diff[i] +
  			bj[library_block[i]]
		}
    
    # latent factors
    for(j in N2) {
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
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model2 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "taubj", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model2)
























phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
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
pdis_unifrac <- phyloseq::distance(phy, method="bray", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])

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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			(b + bj[library_block[i]])*copho[i] + 
  			(c + bj[library_block[i]])*geodis[i] + 
  			(d + bj[library_block[i]])*elev_diff[i] + 
  			(e + bj[library_block[i]])*habitat_diff[i] +
  			bj[library_block[i]]
		}
    
    # latent factors
    for(j in N2) {
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
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model3 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "taubj", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model3)



















phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
#phy <- transform_sample_counts(phy, function(x){ log(x+1) })
sam <- sample_data(phy)
sam2 <- data.frame(id = sam$ID, lat = sam$lat_sign, lon = sam$lon_sign, row.names = rownames(sam))

## GEODIS
geodis <- GeoDistanceInMetresMatrix(sam2)
geodis[lower.tri(geodis, diag = TRUE)] <- "diag"
colnames(geodis) <- rownames(geodis) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
gtemp <- data.frame(cbind(geodis, ID = rownames(geodis)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
geodis2 <- gather(gtemp, key = id, value = geodis, 1:dim(sam)[1])

## PDIS
pdis_unifrac <- phyloseq::distance(pa_phy(phy), method="jaccard", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])

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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			(b + bj[library_block[i]])*copho[i] + 
  			(c + bj[library_block[i]])*geodis[i] + 
  			(d + bj[library_block[i]])*elev_diff[i] + 
  			(e + bj[library_block[i]])*habitat_diff[i] +
  			bj[library_block[i]]
		}
    
    # latent factors
    for(j in N2) {
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
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model4 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "taubj", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model4)









jags1 <- list(jags.model1,
	jags.model2,
	jags.model3,
	jags.model4)

























































phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			e*habitat_diff[i]
		}
    

    # priors
    phi ~ dgamma(.001,.001)
    a ~ dnorm(0,.001)
	b ~ dnorm(0,.001)
	c ~ dnorm(0,.001)
	d ~ dnorm(0,.001)
	e ~ dnorm(0,.001)
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model1.1 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model1.1)
























phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
#phy <- transform_sample_counts(phy, function(x){ log(x+1) })
sam <- sample_data(phy)
sam2 <- data.frame(id = sam$ID, lat = sam$lat_sign, lon = sam$lon_sign, row.names = rownames(sam))

## GEODIS
geodis <- GeoDistanceInMetresMatrix(sam2)
geodis[lower.tri(geodis, diag = TRUE)] <- "diag"
colnames(geodis) <- rownames(geodis) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
gtemp <- data.frame(cbind(geodis, ID = rownames(geodis)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
geodis2 <- gather(gtemp, key = id, value = geodis, 1:dim(sam)[1])

## PDIS
pdis_unifrac <- phyloseq::distance(pa_phy(phy), method="unifrac", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])

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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			e*habitat_diff[i]
		}
 
    # priors
    phi ~ dgamma(.001,.001)
    a ~ dnorm(0,.001)
	b ~ dnorm(0,.001)
	c ~ dnorm(0,.001)
	d ~ dnorm(0,.001)
	e ~ dnorm(0,.001)
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model2.2 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model2.2)
























phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
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
pdis_unifrac <- phyloseq::distance(phy, method="bray", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])

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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			e*habitat_diff[i]
		}
    

    # priors
    phi ~ dgamma(.001,.001)
    a ~ dnorm(0,.001)
	b ~ dnorm(0,.001)
	c ~ dnorm(0,.001)
	d ~ dnorm(0,.001)
	e ~ dnorm(0,.001)
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model3.3 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model3.3)



















phy <- subset_samples(phy1_decontam2, type %in% "gut" & batch %in% c("B3") & genus %in% "Kikihia")
phy <- subset_taxa(phy, taxa_sums(phy) > 100)
#phy <- transform_sample_counts(phy, function(x){ log(x+1) })
sam <- sample_data(phy)
sam2 <- data.frame(id = sam$ID, lat = sam$lat_sign, lon = sam$lon_sign, row.names = rownames(sam))

## GEODIS
geodis <- GeoDistanceInMetresMatrix(sam2)
geodis[lower.tri(geodis, diag = TRUE)] <- "diag"
colnames(geodis) <- rownames(geodis) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
gtemp <- data.frame(cbind(geodis, ID = rownames(geodis)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
geodis2 <- gather(gtemp, key = id, value = geodis, 1:dim(sam)[1])

## PDIS
pdis_unifrac <- phyloseq::distance(pa_phy(phy), method="jaccard", type="samples")
pdis_unifrac <- as.matrix(pdis_unifrac)
pdis_unifrac[lower.tri(pdis_unifrac, diag = TRUE)] <- "diag"
colnames(pdis_unifrac) <- rownames(pdis_unifrac) <- paste(as.character(sam$ID), sam$elevation, sam$habitat, sam$island, sam$batch, sam$species.tree, sep = "_")
ptemp <- data.frame(cbind(pdis_unifrac, ID = rownames(pdis_unifrac)), check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
pdis2 <- gather(ptemp, key = id, value = pdis, 1:dim(sam)[1])

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

## Library size 
sam <- sample_data(phy_first_rel)
library_gen <- c()
for(i in 1:dim(dis)[1]) { 
	a <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "bac_rel"]$bac_rel))
	b <- as.numeric(unique(sam[sam$ID == unlist(lapply(str_split(dis[i, "ID"], "_"), "[", 1)), "depth"]$depth))
	c <- a*b
	library_gen <- c(library_gen, c)
}
dis$library <- library_gen

dis <- dis[!dis$pdis %in% "diag" & !dis$geodis %in% "diag",]
dis$pdis <- as.numeric(dis$pdis)
dis$geodis <- as.numeric(dis$geodis)
dis <- dis[!is.na(dis$copho), ]
dis <- dis[dis$pdis > 0 & dis$pdis < 1 & !is.na(dis$pdis),]


### Within species sampling skew 
table(dis$speciesid)
table(dis$speciesID)
tab <- table(sample_data(phy)$species.tree)

dis$richid <- tab[match(dis$speciesid, names(tab))]
dis$richID <- tab[match(dis$speciesID, names(tab))]
dis$skew <- abs(apply(cbind(dis$richid/(dis$richID + dis$richid), dis$richID/(dis$richID + dis$richid)), 1, min) - 0.5)

dis2 <- dis

dis2$copho <- (dis2$copho - mean(dis2$copho))/sd(dis2$copho)
dis2$geodis <- (dis2$geodis - mean(dis2$geodis))/sd(dis2$geodis)
dis2$elev_diff <- (dis2$elev_diff - mean(dis2$elev_diff))/sd(dis2$elev_diff)
dis2$library <- (dis2$library - mean(dis2$library))/sd(dis2$library)

######## Full, Kikihia #######

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
  			e*habitat_diff[i]
		}
    
    # priors
    phi ~ dgamma(.001,.001)
    a ~ dnorm(0,.001)
	b ~ dnorm(0,.001)
	c ~ dnorm(0,.001)
	d ~ dnorm(0,.001)
	e ~ dnorm(0,.001)
    
  }
  ") ; sink()


jags.data <- list(
  copho = dis2$copho,
  N = dim(dis2)[1],
  N2 = unique(cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE)),
  N3 = as.factor(unique(dis2$skew + 1)),
  skew_block =  as.factor(dis2$skew + 1),
  library_block = cut(dis2$library, breaks =  quantile(dis2$library, probs = seq(0,1,0.25)), include.lowest = TRUE),
  pdis = dis2$pdis,
  geodis = dis2$geodis,
  elev_diff = dis2$elev_diff,
  habitat_diff = dis2$habitat_diff,
  skew = dis2$skew
  )

jags.model4.4 <- jags(data = jags.data,
                   parameters.to.save = c("phi", "a", "b", "c", "d","e"), 
                   model.file = "jags.R",
                   n.chains = 5, 
                   n.iter = 50000,
                   n.burnin = 10000)

#traceplot(jags.model4.4)









jags2 <- list(jags.model1.1,
	jags.model2.2,
	jags.model3.3,
	jags.model4.4)






















jags3 <- c(jags1, jags2)
jags4 <- c()
for(i in 1:length(jags3)){
	out <- data.frame(jags3[[i]]$BUGSoutput$summary)
	out <- out[, colnames(out) %in% c("mean", "X2.5.", "X97.5.")]
	out$par <- rownames(out)
	out$mode <- as.character(i)
	jags4 <- rbind(jags4, out)
	}

jags4 <- jags4[jags4$par %in% c("b", "c", "d", "e"), ]
jags4$model <-  c()
jags4[jags4$mode %in% c("1", "2", "3", "4"), "model"] <- "glm"
jags4[!jags4$mode %in% c("1", "2", "3", "4"), "model"] <- "glmm"
jags4$dis <- c()
jags4[jags4$mode %in% c("1", "5"), "dis"] <- "weighted UniFrac"
jags4[jags4$mode %in% c("2", "6"), "dis"] <- "unweighted UniFrac"
jags4[jags4$mode %in% c("3", "7"), "dis"] <- "Bray"
jags4[jags4$mode %in% c("4", "8"), "dis"] <- "Jaccard"
jags4$par2 <- c()
jags4[jags4$par %in% c("b"), "par2"] <- "Cophenetic Distance"
jags4[jags4$par %in% c("c"), "par2"] <- "Geographical Distance"
jags4[jags4$par %in% c("d"), "par2"] <- "Elevation Difference"
jags4[jags4$par %in% c("e"), "par2"] <- "Habitat Difference"


g <- ggplot(data = jags4, aes(x = paste(par2, dis, sep = ": "), y = mean, col = model)) +
	geom_point() +
	geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=0.1)	+
	geom_hline(yintercept = 0, col = "red") +
	xlab("") +
	ylab("Posterior Estimate") +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 7, angle = 45, hjust = 1))

ggsave("mcmc.png", g, scale = 0.8)









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









library(lme4)
lmer(pdis ~ copho + (copho | library) + geodis + (geodis | library) + elev_diff + (elev_diff | library) + habitat_diff + (habitat_diff | library), dis2)
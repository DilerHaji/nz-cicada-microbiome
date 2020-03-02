###############################

dtable <-  function(ps) {
  if (!is.na(match(colnames(t(ps@otu_table)), row.names(ps@sam_data)))[1]) {
    h <- data.frame(t(ps@otu_table@.Data))
    h$sum <- apply(h, 1, sum)
    tax <- data.frame(ps@tax_table@.Data)
    h2 <- cbind(h, tax)
    h2 <- h2[order(h2$sum, decreasing = TRUE),]
    h2 <- h2[, colnames(h2) != "sum"]
    return(h2)} else {
      h <- data.frame(ps@otu_table@.Data)
      h$sum <- apply(h, 1, sum)
      tax <- data.frame(ps@tax_table@.Data)
      h2 <- cbind(h, tax)
      h2 <- h2[order(h2$sum, decreasing = TRUE),]
      h2 <- h2[, colnames(h2) != "sum"]
    }
}




###############################

ctm <- function(phy, taxup){
  n <- names(sample_sums(phy)[is.na(match(names(sample_sums(phy)),
names(sample_sums(phy)[sample_sums(phy) == 0])))])
  phy <- prune_samples(samples = n, x = phy)
  prevx = apply(X = otu_table(phy), MARGIN = ifelse(taxa_are_rows(phy), yes = 1, no = 2), FUN =
function(x){sum(x > 0)})
  names_to_keep <- names(prevx[prevx > 3])
  phy <- prune_taxa(taxa = names_to_keep, x = phy)
  taxup <- length(unique(rownames(data.frame(otu_table(phy)))))
  s <- data.frame(1:length(sample_names(phy)))
  for(i in 1:taxup){
    t <- names(sort(taxa_sums(phy), decreasing = TRUE)[i])
    p <- prune_taxa(taxa = t, phy)
    t2 <- data.frame(as.numeric(sample_sums(p)))
    colnames(t2) <- t
    s <- cbind(s, t2)
  }
  s <- s[,-1]
  s <- cbind(s, names = rownames(s), sums = as.numeric(sample_sums(phy)))
  f <- taxup-1
  s <- gather(s, key = type, value = abund, 1:f)
  #s <- s[s$abund != 0,]
  x = s
  c <- glm.nb(x$abund ~ x$sums + factor(x$type))
  c.confidence <- 1 - pchisq(summary(c)$deviance, summary(c)$df.residual)
  c.summary <- summary(c)
  strp <- function(x){sapply(strsplit(names(x[,1]), ")"), "[",
2)[!is.na(sapply(strsplit(names(x[,1]), ")"), "[", 2))]}
  c.sig <- strp(coef(summary(c))[coef(summary(c))[,4] < 0.05,])
  c.notsig <- strp(coef(summary(c))[coef(summary(c))[,4] > 0.05,])
  c.sig.keep <- strp(coef(summary(c))[coef(summary(c))[,4] < 0.05 & coef(summary(c))[,1] > 0,])
  c.sig.toss <- strp(coef(summary(c))[coef(summary(c))[,4] < 0.05 & coef(summary(c))[,1] < 0,])
  phy.sig <- prune_taxa(taxa = c.sig, phy)
  phy.notsig <- prune_taxa(taxa = c.notsig, phy)
  phy.sig.keep <- prune_taxa(taxa = c.sig.keep, phy)
  phy.sig.toss <- prune_taxa(taxa = c.sig.toss, phy)
  return(list(phy.sig, phy.notsig, phy.sig.keep, phy.sig.toss))
}




###############################

grand_pre <- function(x, tax_prev) {
  prevx = apply(X = otu_table(x), MARGIN = ifelse(taxa_are_rows(x), yes = 1, no = 2), FUN =
function(x){sum(x > 0)})
  prevx2 = data.frame(Prevalence = prevx, TotalAbundance = taxa_sums(x), tax_table(x))
  y <- plyr::ddply(prevx2, tax_prev, function(df1){cbind(mean_prev = median(df1$Prevalence),
sum_prev = sum(df1$Prevalence), mean_abund = median(df1$TotalAbundance), sum_abund =
sum(df1$TotalAbundance))})
  return(y)
}




###############################

subset_prev <- function(x, prev_thresh, by_abund, tax){
  phypre <- grand_pre(x, tax_prev = tax)
  if(prev_thresh){phypre <- as.character(phypre[,1][phypre$sum_prev <= prev_thresh])}
  if(by_abund){phypre <- as.character(phypre[,1][phypre$mean_abund <= by_abund])}
  return(phypre)
}




###############################

extaxa <- function(phy, taxa){
  u <- tax_table(prune_taxa(taxa = taxa, x = phy))
  return(u)
  }




###############################

toptax <- function(phy, level){ 
  ertab <- c()
  for(i in sample_names(phy)) {
  j <- prune_samples(sample = i, x = phy)
  y <- names(taxa_sums(j)[taxa_sums(j) == max(taxa_sums(j))])
  e <- tax_table(prune_taxa(taxa = y, x = phy))[, level]
  ertab <- rbind(ertab, data.frame(sample = i, taxon = e))
  }
  return(ertab)
  }




###############################

logphy <- function(phy){
  transform_sample_counts(phy, function(x) log(1 + x))
  }




###############################

relphy <- function(phy){
  transform_sample_counts(phy, fun = function(x){x/sum(x)})
}




###############################

rarephy <- function(phy, size){rarefy_even_depth(physeq = phy, rngseed = 1, sample.size = size,
replace = TRUE, trimOTUs = TRUE)}




###############################

barplot_phy <- function(phylo, group, upper){
  
  phy_filtered_type <- merge_samples(phylo, group = group)
  sample_data(phy_filtered_type)[, colnames(sample_data(phy_filtered_type)) == group] <- rownames(sample_data(phy_filtered_type))
  phy_filtered_type <- transform_sample_counts(phy_filtered_type, fun = function(x){x/sum(x)})
  
  prop <- list(c())
  name <- c()
  for(i in 2:length(sample_names(phy_filtered_type))) {
    e <- sample_names(phy_filtered_type)[i]
    prop[[i]] <- prune_samples(samples = e, phy_filtered_type)
    name <- c(name, names(taxa_sums(prop[[i]])[order(taxa_sums(prop[[i]]), decreasing = TRUE)][1:upper]))
  }
  
  phy_filtered_type <- merge_samples(phylo, group = group )
  phy_filtered_type <- prune_taxa(taxa = name, phy_filtered_type)
  sample_data(phy_filtered_type)[, colnames(sample_data(phy_filtered_type)) == group] <- rownames(sample_data(phy_filtered_type))
  phy_filtered_type <- transform_sample_counts(phy_filtered_type, fun = function(x){x/sum(x)})
  
  return(plot_bar(phy_filtered_type, x = group, fill = "Genus") + 
           geom_bar(stat = "identity") +
           #scale_fill_brewer(palette = "Paired") +
           #scale_x_discrete(limits = data.frame(phy@sam_data[,colnames(phy@sam_data) == "type"])[,1]) +
           theme_bw() +
           theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 axis.text.x = element_text(angle = 90, hjust = 1)))
  return(sample_data(phy_filtered_type)[,colnames(sample_data(phy_filtered_type)) == "genus"])
}
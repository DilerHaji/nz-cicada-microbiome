# Just by opening the project, you have already set the working directory ####
getwd()

alpha <- read.table(file="../workshop.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.groups.ave-std.summary",
                    head = TRUE, stringsAsFactors = FALSE)
str(alpha)
alpha2 <- alpha[alpha$method == "ave",]
plot(alpha2$shannon)

alpha[,2]

?sort
?order

library(dplyr)
library(tidyr)

alpha.ave <- filter(alpha, method == "ave")

plot(alpha.ave$coverage, alpha.ave$invsimpson)
plot(alpha.ave$shannon, alpha.ave$invsimpson)

env <- read.table("weir.farm.env.txt", header = TRUE, stringsAsFactors = FALSE)
sample.name <- read.table("sample.txt", header = TRUE, stringsAsFactors = FALSE)

env1 <- right_join(env, sample.name, copy = TRUE)

alpha.env <- right_join(env1, alpha, by = "group")
alpha.env <- filter(alpha.env, method == "ave")

library(ggplot2)
ggplot(data=alpha.env, aes(x = shannon, y= invsimpson)) +
  #geom_smooth(method="lm", col = "gray") +
  geom_point(aes(color=factor(alpha.env$Type), size = sobs), pch = 21) +
  theme_classic() +
  labs(col = "Type", size = "Species observed")
  

# Make model and than this 
geom_line()



############## WHAT YOU SHOULD LOOK AT EVERYTIME YOU OPEN MICROBIOME DATA ###################

# MothurWorkshopDataExploration

#Function to read mothur generated dist into R
parseDistanceDF = function(phylip_file) {
  
  # Read the first line of the phylip file to find out how many sequences/samples it contains
  temp_connection = file(phylip_file, 'r')
  len = readLines(temp_connection, n=1)
  len = as.numeric(len)
  len = len +1
  close(temp_connection)
  
  
  phylip_data = read.table(phylip_file, fill=T, row.names=1, skip=1, col.names=1:len)
  colnames(phylip_data) <- row.names(phylip_data)
  return(phylip_data)
}

# Reading in beta diversity dist matrices

jc <- parseDistanceDF("../workshop.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.jest.0.03.lt.ave.dist")
bc <- parseDistanceDF("../workshop.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.braycurtis.0.03.lt.ave.dist")
tyc <- parseDistanceDF("../workshop.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.thetayc.0.03.lt.ave.dist")


# Run nms on dist

#install.packages("vegan")
library(vegan)

# trymin mean permutation = 50 times 
# wascores - we want to use distance files generated from mothur so we set this to false
# K is the number of dimensions you want
# REMEMBER THAT NMDS is non-metric!!! so don't try to interprate values. ONLY APPLIES TO THIS DATASET AND THAT'S IT 
jc.nms <- metaMDS(as.dist(jc), k = 2, trymin = 50, trymax = 500, wascores = FALSE)
tyc.nms <- metaMDS(as.dist(tyc), k = 2, trymin = 50, trymax = 500, wascores = FALSE)
bc.nms <- metaMDS(as.dist(bc), k = 2, trymin = 50, trymax = 500, wascores = FALSE)

mantel.partial(jc, bc, tyc)




# Stress is the inverse of R-squared. Basically the residuals. The smaller the stress the better.
# Stress over 0.2 (or 20 on a whole number scale) the relationship between the samples and random is the same
# The more dimensions the higher the stress
# Stress is telling you (sort of) about the residuals of regression 
# Run nms a bunch of times 

# The coordinates are from the best iteration from above - THE ONE WITH THE LOWEST STRESS!!!
# view values within a list 
jc.nms$points
# view parts within a list 
jc.nms$points[,1]


# built in plotting function vegan
ordiplot(jc.nms)
ordilabel(jc.nms, display = "sites")

jc.points <- data.frame(jc.nms$points)

#install.packages("RColorBrewer")
library(RColorBrewer)

ggplot(jc.points, aes(x=MDS1, y=MDS2, label=rownames(jc.points))) +
  geom_point(aes(color=factor(alpha.env$Type)), size = 4) +
  theme_bw()
  
ggplot(jc.points, aes(x=MDS1, y=MDS2, label=rownames(jc.points))) +
  geom_point(aes(color=factor(alpha.env$Type), shape = factor(alpha.env$Site)), size = 4) +
  scale_color_brewer(palette = "Accent") +
  theme_bw()

ggsave(file= "jc.nms.jpg")


ggplot(jc.points, aes(x=MDS1, y=MDS2, label=rownames(jc.points))) +
  geom_point(aes(color=factor(alpha.env$Site)), size = 4) +
  theme(axis.title.x  = element_blank(),
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.y  = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank())



# partial mantell - vegan or ecodist 
# Tells you the difference between two matrices 
# Use this to compare between mothur and Qiime 

?mantel.partial


####  Hypothesis testing of total community PERMANOVA ####
# Did my treatements affect the structure of my entire community 
# Permutation based mutliple anova (only catagorical variables)
# Also ADONIS in vegan - works on continuous variables too

# You need to make sure you don't over permute!!!
perm.jc <- adonis(as.dist(jc) ~ alpha.env$Type*alpha.env$Site, perm = 99, rm.na = TRUE)
perm.jc 
# Uses an F model to generate a p-value
# Almost ALWAYS significant. If you increase permutations, you'll decrease p-value 
# R2 is NOT R-squared. Related to the sum of squares. Adds up to one. The amout of variability within your data that is explained by your model
# The first term (TYPE) determins the R2. If the order changes in the formula, the numbers change. Pick what you think will have the biggest contribution to variance first and onward 


### Display alpha diversity with boxplots ####

ggplot(alpha.env, aes(y=invsimpson, x=Type)) +
  geom_boxplot() +
  geom_jitter(col = "red") +
  theme_bw()


#### Indicator Species ####

otu <- read.table(file = "../workshop.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.0.03.subsample.shared", header = TRUE, stringsAsFactors = FALSE, row.names = 2)

library(dplyr)

otu <- select(otu, -label, -numOtus)

# Reading in taxonomy file to replace OTU ID in the OTU names 
# Kendra's code for how to do this!!!! 
# WARNING WARNING : The next line won't work because of the mothur bug - sequences are not all getting 6 taxonomic level assignments 
taxa <- read.table(textConnection(gsub("\\(.+?\\);", "\t", readLines("../workshop.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.0.03.cons.taxonomy"))), col.names=c("OTU", "Size", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), skip=1)

# Subsettting the rows in taxa where the OTU names is the column names
# Filter out the rows in cons.taxonomy that are not in the subsampled otu table 
taxa <- taxa[taxa$OTU %in% names(otu),]

# install.packages("indicspecies")
library(indicspecies)

## Couldn't get taxa to read properly because not all rows had 6 taxonomix levels
# Multipatt allows you to compare multiple groups and multiple combinations of groups 
# 
indic <- multipatt()












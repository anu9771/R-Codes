# Loading required packages
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(BiocManager)
library(microbiome)
library(NetCoMi)
library(WGCNA)

# Setting the working directory
setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology")

#Phyloseq2*****************************************************************************************************************
#Adding .qza files to phyloseq object
physeq=qza_to_phyloseq(
  features = "datafiles/table.qza",
  taxonomy="datafiles/taxonomy.qza",
)

#Reading the meta data
metx=read.delim("datafiles/metadata_table.txt",sep=",",header=T,row.names=sample_names(physeq))
#head(metx)

#Meta data as sample data
metx=metx[,-1]
metx=sample_data(metx)
#metx

#Merging the sample
complete_physeq=merge_phyloseq(physeq,metx)
#complete_physeq
xa=sample_data(complete_physeq)
xb=tax_table(complete_physeq)
xc=otu_table(complete_physeq)

complete_agg <- tax_glom(complete_physeq, 'Genus')
#complete_agg
tax_table(complete_agg)
complete_agg <- aggregate_taxa(complete_agg, 'Genus')
complete_agg

# keep only taxa that were ....................................................
minTotRelAbun <- 5e-5
sums <- taxa_sums(complete_agg)
#sums
keepTaxa <- taxa_names(complete_agg)[which((sums / sum(sums)) > minTotRelAbun)]
#keepTaxa
filt_complete_agg <- prune_taxa(keepTaxa, complete_agg)
filt_complete_agg


saveRDS(filt_complete_agg, "Breed/Complete_SriLanka.rds")


comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)


#SparCC#######################################################################################################

SparCC_data <- readRDS("Breed/Complete_SriLanka.rds")

top_SparCC <- prune_taxa(names(sort(taxa_sums(SparCC_data),TRUE)[1:50]), SparCC_data)
#plot_heatmap(top_JRJRC)
SparCC_data

net_single_SparCC <- netConstruct(SparCC_data,
                                 verbose = 3,
                                 filtTax = "highestFreq",
                                 filtTaxPar = list(highestFreq = 100),
                                 # filtTax = "none",
                                 # filtTaxPar = "totalReads",
                                 filtSamp = "totalReads",
                                 filtSampPar = list(totalReads = 1000),
                                 zeroMethod = "none", normMethod = "none",
                                 measure = "sparcc",
                                 sparsMethod = "threshold", thresh = 0.4,
                                 dissFunc = "signed", 
                                 seed = 123456)

saveRDS(net_single_SparCC, "Breed/Complete_SriLanka/Complete_SriLanka_SparCC.rds")

props_single_SparCC <- netAnalyze(net_single_SparCC, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_SparCC, "Breed/Complete_SriLanka/Complete_SriLanka_SparCC analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_SparCC)
net.summary

p_SparCC <- plot(props_single_SparCC,
                shortenLabels = "none",
                # labelLength = 16,
                # charToRm = "g__",
                labelScale = FALSE,
                rmSingles = "all",
                nodeSize = "eigenvector",
                nodeColor = "cluster",
                hubBorderCol = "blue",
                cexNodes = 1,
                cexLabels = 0.5,
                edgeWidth = 1,
                highlightHubs = TRUE,
                cexHubs = 1.5,
                # cexHubLabels = 2,
                title1 = "SparCC Network on Genus level for Sri Lankan Data.", 
                showTitle = TRUE,
                cexTitle = 1.5)

saveRDS(p_SparCC, "Breed/Complete_SriLanka/Complete_SriLanka_SparCC plot.rds")


from <- p_SparCC$labels$labels1[p_SparCC$q1$Edgelist$from]
to <- p_SparCC$labels$labels1[p_SparCC$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_SparCC$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_SparCC$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_SparCC$q1$Edgelist$from)))
edges$from <- p_SparCC$q1$Edgelist$from
edges$to <- p_SparCC$q1$Edgelist$to
edges$weight <- p_SparCC$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Complete_SriLanka/Complete_SriLanka_SparCC edge data.csv", row.names = FALSE)

hubs <- props_single_SparCC$hubs$hubs1
write(hubs, "Breed/Complete_SriLanka/Complete_SriLanka_SparCC Hubs.txt")

node.lables <- p_SparCC$labels$labels1
clust <- props_single_SparCC$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_SparCC$centralities$degree1[nodes$lable]
evs=props_single_SparCC$centralities$eigenv1[nodes$lable]
betweennesses=props_single_SparCC$centralities$between1[nodes$lable]
closenesses=props_single_SparCC$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/Complete_SriLanka/Complete_SriLanka_SparCC node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Complete_SriLanka/Complete_SriLanka_SparCC hub data.csv", row.names = FALSE)

##Associations########################################################################################################################
Genus_Sparcc_Breed_complete_Top100=data.frame(Node1=net_single_SparCC[["edgelist1"]][["v1"]],
                                           Node2=net_single_SparCC[["edgelist1"]][["v2"]],
                                           Association=net_single_SparCC[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_SparCC[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_SparCC[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Breed_complete_Top100=Genus_Sparcc_Breed_complete_Top100[order(Genus_Sparcc_Breed_complete_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Sparcc_Breed_complete_Top100,"Breed/Complete_SriLanka/Genus_Sparcc_Breed_complete_Top100.csv",row.names = FALSE)



#SPRING#######################################################################################################
SPRING_data <- readRDS("Breed/Complete_SriLanka.rds")

top_SPRING <- prune_taxa(names(sort(taxa_sums(SPRING_data),TRUE)[1:50]), SPRING_data)
#plot_heatmap(top_JRJRC)
SPRING_data

net_single_SPRING <- netConstruct(SPRING_data,
                                  verbose = 3,
                                  filtTax = "highestFreq",
                                  filtTaxPar = list(highestFreq = 100),
                                  # filtTax = "none",
                                  # filtTaxPar = "totalReads",
                                  filtSamp = "totalReads",
                                  filtSampPar = list(totalReads = 1000),
                                  zeroMethod = "none", normMethod = "none",
                                  measure = "spring",
                                  measurePar = list(nlambda = 20, 
                                                    rep.num = 20,
                                                    thresh = 0.1, 
                                                    Rmethod="approx",
                                                    subsample.ratio = 0.7, # 10*sqrt(n)/n for n > 144
                                                    lambda.min.ratio = 0.01,
                                                    lambdaseq = "data-specific",
                                                    ncores=1),
                                  sparsMethod = "none",  dissFunc = "signed", 
                                  seed = 123456)

saveRDS(net_single_SPRING, "Breed/Complete_SriLanka/Complete_SriLanka_SPRING.rds")

props_single_SPRING <- netAnalyze(net_single_SPRING, 
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "eigenvector", 
                                  hubQuant = 0.95)

saveRDS(props_single_SPRING, "Breed/Complete_SriLanka/Complete_SriLanka_SPRING analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_SPRING)
net.summary

p_SPRING <- plot(props_single_SPRING,
                 shortenLabels = "none",
                 # labelLength = 16,
                 # charToRm = "g__",
                 labelScale = FALSE,
                 rmSingles = "all",
                 nodeSize = "eigenvector",
                 nodeColor = "cluster",
                 hubBorderCol = "blue",
                 cexNodes = 1,
                 cexLabels = 0.5,
                 edgeWidth = 1,
                 highlightHubs = TRUE,
                 cexHubs = 1.5,
                 # cexHubLabels = 2,
                 title1 = "SPRING Network on Genus level for Sri Lankan Data.", 
                 showTitle = TRUE,
                 cexTitle = 1.5)

saveRDS(p_SPRING, "Breed/Complete_SriLanka/Complete_SriLanka_SPRING plot.rds")


from <- p_SPRING$labels$labels1[p_SPRING$q1$Edgelist$from]
to <- p_SPRING$labels$labels1[p_SPRING$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_SPRING$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_SPRING$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_SPRING$q1$Edgelist$from)))
edges$from <- p_SPRING$q1$Edgelist$from
edges$to <- p_SPRING$q1$Edgelist$to
edges$weight <- p_SPRING$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Complete_SriLanka/Complete_SriLanka_SPRING edge data.csv", row.names = FALSE)

hubs <- props_single_SPRING$hubs$hubs1
write(hubs, "Breed/Complete_SriLanka/Complete_SriLanka_SPRING Hubs.txt")

node.lables <- p_SPRING$labels$labels1
clust <- props_single_SPRING$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_SPRING$centralities$degree1[nodes$lable]
evs=props_single_SPRING$centralities$eigenv1[nodes$lable]
betweennesses=props_single_SPRING$centralities$between1[nodes$lable]
closenesses=props_single_SPRING$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/Complete_SriLanka/Complete_SriLanka_SPRING node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Complete_SriLanka/Complete_SriLanka_SPRING hub data.csv", row.names = FALSE)


#CCLASSO#######################################################################################################

LASSO_data <- readRDS("Breed/Complete_SriLanka.rds")

top_LASSO <- prune_taxa(names(sort(taxa_sums(LASSO_data),TRUE)[1:50]), LASSO_data)
#plot_heatmap(top_JRJRC)
LASSO_data

net_single_LASSO <- netConstruct(LASSO_data,
                                 verbose = 3,
                                 filtTax = "highestFreq",
                                 filtTaxPar = list(highestFreq = 100),
                                 filtSamp = "totalReads",
                                 filtSampPar = list(totalReads = 1000),
                                 zeroMethod = "none", normMethod = "none",
                                 measure = "cclasso",
                                 sparsMethod = "threshold", thresh = 0.4,
                                 dissFunc = "signed", 
                                  seed = 123456)

saveRDS(net_single_LASSO, "Breed/Complete_SriLanka/Complete_SriLanka_LASSO.rds")

props_single_LASSO <- netAnalyze(net_single_LASSO, 
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "eigenvector", 
                                  hubQuant = 0.95)

saveRDS(props_single_LASSO, "Breed/Complete_SriLanka/Complete_SriLanka_LASSO analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_LASSO)
net.summary

p_LASSO <- plot(props_single_LASSO,
                 shortenLabels = "none",
                 # labelLength = 16,
                 # charToRm = "g__",
                 labelScale = FALSE,
                 rmSingles = "all",
                 nodeSize = "eigenvector",
                 nodeColor = "cluster",
                 hubBorderCol = "blue",
                 cexNodes = 1,
                 cexLabels = 0.5,
                 edgeWidth = 1,
                 highlightHubs = TRUE,
                 cexHubs = 1.5,
                 # cexHubLabels = 2,
                 title1 = "LASSO Network on Genus level for Sri Lankan Data.", 
                 showTitle = TRUE,
                 cexTitle = 1.5)

saveRDS(p_LASSO, "Breed/Complete_SriLanka/Complete_SriLanka_LASSO plot.rds")


from <- p_LASSO$labels$labels1[p_LASSO$q1$Edgelist$from]
to <- p_LASSO$labels$labels1[p_LASSO$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_LASSO$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_LASSO$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_LASSO$q1$Edgelist$from)))
edges$from <- p_LASSO$q1$Edgelist$from
edges$to <- p_LASSO$q1$Edgelist$to
edges$weight <- p_LASSO$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Complete_SriLanka/Complete_SriLanka_LASSO edge data.csv", row.names = FALSE)

hubs <- props_single_LASSO$hubs$hubs1
write(hubs, "Breed/Complete_SriLanka/Complete_SriLanka_LASSO Hubs.txt")

node.lables <- p_LASSO$labels$labels1
clust <- props_single_LASSO$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_LASSO$centralities$degree1[nodes$lable]
evs=props_single_LASSO$centralities$eigenv1[nodes$lable]
betweennesses=props_single_LASSO$centralities$between1[nodes$lable]
closenesses=props_single_LASSO$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/Complete_SriLanka/Complete_SriLanka_LASSO node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Complete_SriLanka/Complete_SriLanka_LASSO hub data.csv", row.names = FALSE)


#CCREPE#######################################################################################################

CCREPE_data <- readRDS("Breed/Complete_SriLanka.rds")

top_CCREPE <- prune_taxa(names(sort(taxa_sums(CCREPE_data),TRUE)[1:50]), CCREPE_data)
#plot_heatmap(top_JRJRC)
CCREPE_data

net_single_CCREPE <- netConstruct(CCREPE_data,
                                  verbose = 3,
                                  filtTax = "highestFreq",
                                  filtTaxPar = list(highestFreq = 100),
                                  filtSamp = "totalReads",
                                  filtSampPar = list(totalReads = 1000),
                                  zeroMethod = "none", normMethod = "none",
                                  measure = "ccrepe",
                                  sparsMethod = "threshold", thresh = 0.4,
                                  dissFunc = "signed", 
                                  seed = 123456)

saveRDS(net_single_CCREPE, "Breed/Complete_SriLanka/Complete_SriLanka_CCREPE.rds")

props_single_CCREPE <- netAnalyze(net_single_CCREPE, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_CCREPE, "Breed/Complete_SriLanka/Complete_SriLanka_CCREPE analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_CCREPE)
net.summary

p_CCREPE <- plot(props_single_CCREPE,
                shortenLabels = "none",
                # labelLength = 16,
                # charToRm = "g__",
                labelScale = FALSE,
                rmSingles = "all",
                nodeSize = "eigenvector",
                nodeColor = "cluster",
                hubBorderCol = "blue",
                cexNodes = 1,
                cexLabels = 0.5,
                edgeWidth = 1,
                highlightHubs = TRUE,
                cexHubs = 1.5,
                # cexHubLabels = 2,
                title1 = "CCREPE Network on Genus level for Sri Lankan Data.", 
                showTitle = TRUE,
                cexTitle = 1.5)

saveRDS(p_CCREPE, "Breed/Complete_SriLanka/Complete_SriLanka_CCREPE plot.rds")


from <- p_CCREPE$labels$labels1[p_CCREPE$q1$Edgelist$from]
to <- p_CCREPE$labels$labels1[p_CCREPE$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_CCREPE$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_CCREPE$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_CCREPE$q1$Edgelist$from)))
edges$from <- p_CCREPE$q1$Edgelist$from
edges$to <- p_CCREPE$q1$Edgelist$to
edges$weight <- p_CCREPE$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Complete_SriLanka/Complete_SriLanka_CCREPE edge data.csv", row.names = FALSE)

hubs <- props_single_CCREPE$hubs$hubs1
write(hubs, "Breed/Complete_SriLanka/Complete_SriLanka_CCREPE Hubs.txt")

node.lables <- p_CCREPE$labels$labels1
clust <- props_single_CCREPE$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_CCREPE$centralities$degree1[nodes$lable]
evs=props_single_CCREPE$centralities$eigenv1[nodes$lable]
betweennesses=props_single_CCREPE$centralities$between1[nodes$lable]
closenesses=props_single_CCREPE$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/Complete_SriLanka/Complete_SriLanka_CCREPE node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Complete_SriLanka/Complete_SriLanka_CCREPE hub data.csv", row.names = FALSE)


#PROPR#######################################################################################################

# PROPR_data <- readRDS("Breed/Complete_SriLanka.rds")
# 
# top_PROPR <- prune_taxa(names(sort(taxa_sums(PROPR_data),TRUE)[1:50]), PROPR_data)
# #plot_heatmap(top_JRJRC)
# PROPR_data
# 
# net_single_PROPR <- netConstruct(PROPR_data,
#                                  verbose = 3,
#                                  filtTax = "highestFreq",
#                                  filtTaxPar = list(highestFreq = 100),
#                                  filtSamp = "totalReads",
#                                  filtSampPar = list(totalReads = 1000),
#                                  zeroMethod = "pseudo", normMethod = "none",
#                                  measure = "propr",
#                                  sparsMethod = "threshold", thresh = 0.2,
#                                  dissFunc = "signed",  
#                                  seed = 123456)
# 
# saveRDS(net_single_PROPR, "Breed/Complete_SriLanka/Complete_SriLanka_PROPR.rds")
# 
# props_single_PROPR <- netAnalyze(net_single_PROPR, 
#                                   clustMethod = "cluster_fast_greedy",
#                                   hubPar = "eigenvector", 
#                                   hubQuant = 0.95)
# 
# saveRDS(props_single_PROPR, "Breed/Complete_SriLanka/Complete_SriLanka_PROPR analysis.rds")
# 
# #?summary.microNetProps
# net.summary <- summary(props_single_PROPR)
# net.summary
# 
# p_PROPR <- plot(props_single_PROPR,
#                  shortenLabels = "none",
#                  # labelLength = 16,
#                  # charToRm = "g__",
#                  labelScale = FALSE,
#                  rmSingles = "all",
#                  nodeSize = "eigenvector",
#                  nodeColor = "cluster",
#                  hubBorderCol = "blue",
#                  cexNodes = 1,
#                  cexLabels = 0.5,
#                  edgeWidth = 1,
#                  highlightHubs = TRUE,
#                  cexHubs = 1.5,
#                  # cexHubLabels = 2,
#                  title1 = "PROPR Network on Genus level for Sri Lankan Data.", 
#                  showTitle = TRUE,
#                  cexTitle = 1.5)
# 
# saveRDS(p_PROPR, "Breed/Complete_SriLanka/Complete_SriLanka_PROPR plot.rds")
# 
# 
# from <- p_PROPR$labels$labels1[p_PROPR$q1$Edgelist$from]
# to <- p_PROPR$labels$labels1[p_PROPR$q1$Edgelist$to]
# direction <- vector(mode = "integer")
# 
# for (x in 1:length(from)){
#   if (net_single_PROPR$assoMat1[from[x],to[x]] < 0){
#     direction[x] <- -1 
#   }
#   else if (net_single_PROPR$assoMat1[from[x],to[x]] > 0){
#     direction[x] <- 1 
#   }
# } 
# 
# edges <- data.frame(row.names = c(1:length(p_PROPR$q1$Edgelist$from)))
# edges$from <- p_PROPR$q1$Edgelist$from
# edges$to <- p_PROPR$q1$Edgelist$to
# edges$weight <- p_PROPR$q1$Edgelist$weight
# edges$association <- direction
# 
# write.csv(edges, file = "Breed/Complete_SriLanka/Complete_SriLanka_PROPR edge data.csv", row.names = FALSE)
# 
# hubs <- props_single_PROPR$hubs$hubs1
# write(hubs, "Breed/Complete_SriLanka/Complete_SriLanka_PROPR Hubs.txt")
# 
# node.lables <- p_PROPR$labels$labels1
# clust <- props_single_PROPR$clustering$clust1[node.lables]
# node.hubs <- integer(length((node.lables)))
# node.hubs[match(hubs,node.lables)] <- 1
# 
# 
# nodes <- data.frame(row.names = 1:length(node.lables))
# nodes$index <- 1:length(node.lables)
# nodes$lable <- node.lables
# nodes$cluster <- clust
# nodes$hubs <- node.hubs
# degrees <- props_single_PROPR$centralities$degree1[nodes$lable]
# evs=props_single_PROPR$centralities$eigenv1[nodes$lable]
# betweennesses=props_single_PROPR$centralities$between1[nodes$lable]
# closenesses=props_single_PROPR$centralities$close1[nodes$lable]
# nodes$degree <- degrees
# nodes$ev=evs
# nodes$between=betweennesses
# nodes$close=closenesses
# 
# 
# write.csv(nodes, file = "Breed/Complete_SriLanka/Complete_SriLanka_PROPR node data.csv", row.names = FALSE)
# 
# hubdata <- subset(nodes, hubs==1, select = -hubs)
# write.csv(hubdata, file = "Breed/Complete_SriLanka/Complete_SriLanka_PROPR hub data.csv", row.names = FALSE)
# 


#SPEC-EASI#######################################################################################################

SPECEASI_data <- readRDS("Breed/Complete_SriLanka.rds")

top_SPECEASI <- prune_taxa(names(sort(taxa_sums(SPECEASI_data),TRUE)[1:50]), SPECEASI_data)
#plot_heatmap(top_JRJRC)
SPECEASI_data

net_single_SPECEASI <- netConstruct(SPECEASI_data,
                                    verbose = 3,
                                    filtTax = "highestFreq",
                                    filtTaxPar = list(highestFreq = 100),
                                    filtSamp = "totalReads",
                                    filtSampPar = list(totalReads = 1000),
                                    zeroMethod = "none", normMethod = "none",
                                    measure = "spieceasi",
                                    measurePar = list(method = "mb", 
                                                      nlambda = 20,
                                                      lambda.min.ratio = 0.01,
                                                      pulsar.params = list(rep.num = 20,
                                                                           thresh = 0.1, 
                                                                           subsample.ratio = 0.7)),
                                    sparsMethod = "none", dissFunc = "signed", 
                                  seed = 123456)

saveRDS(net_single_SPECEASI, "Breed/Complete_SriLanka/Complete_SriLanka_SPECEASI.rds")

props_single_SPECEASI <- netAnalyze(net_single_SPECEASI, 
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "eigenvector", 
                                  hubQuant = 0.95)

saveRDS(props_single_SPECEASI, "Breed/Complete_SriLanka/Complete_SriLanka_SPECEASI analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_SPECEASI)
net.summary

p_SPECEASI <- plot(props_single_SPECEASI,
                 shortenLabels = "none",
                 # labelLength = 16,
                 # charToRm = "g__",
                 labelScale = FALSE,
                 rmSingles = "all",
                 nodeSize = "eigenvector",
                 nodeColor = "cluster",
                 hubBorderCol = "blue",
                 cexNodes = 1,
                 cexLabels = 0.5,
                 edgeWidth = 1,
                 highlightHubs = TRUE,
                 cexHubs = 1.5,
                 # cexHubLabels = 2,
                 title1 = "SPECEASI Network on Genus level for Sri Lankan Data.", 
                 showTitle = TRUE,
                 cexTitle = 1.5)

saveRDS(p_SPECEASI, "Breed/Complete_SriLanka/Complete_SriLanka_SPECEASI plot.rds")


from <- p_SPECEASI$labels$labels1[p_SPECEASI$q1$Edgelist$from]
to <- p_SPECEASI$labels$labels1[p_SPECEASI$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_SPECEASI$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_SPECEASI$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_SPECEASI$q1$Edgelist$from)))
edges$from <- p_SPECEASI$q1$Edgelist$from
edges$to <- p_SPECEASI$q1$Edgelist$to
edges$weight <- p_SPECEASI$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Complete_SriLanka/Complete_SriLanka_SPECEASI edge data.csv", row.names = FALSE)

hubs <- props_single_SPECEASI$hubs$hubs1
write(hubs, "Breed/Complete_SriLanka/Complete_SriLanka_SPECEASI Hubs.txt")

node.lables <- p_SPECEASI$labels$labels1
clust <- props_single_SPECEASI$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_SPECEASI$centralities$degree1[nodes$lable]
evs=props_single_SPECEASI$centralities$eigenv1[nodes$lable]
betweennesses=props_single_SPECEASI$centralities$between1[nodes$lable]
closenesses=props_single_SPECEASI$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/Complete_SriLanka/Complete_SriLanka_SPECEASI node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Complete_SriLanka/Complete_SriLanka_SPECEASI hub data.csv", row.names = FALSE)


#GCODA#######################################################################################################

GCODA_data <- readRDS("Breed/Complete_SriLanka.rds")

top_GCODA <- prune_taxa(names(sort(taxa_sums(GCODA_data),TRUE)[1:50]), GCODA_data)
#plot_heatmap(top_JRJRC)
GCODA_data

net_single_GCODA <- netConstruct(GCODA_data,
                                 verbose = 3,
                                 filtTax = "highestFreq",
                                 filtTaxPar = list(highestFreq = 100),
                                 filtSamp = "totalReads",
                                 filtSampPar = list(totalReads = 1000),
                                 zeroMethod = "none", normMethod = "none",
                                 measure = "gcoda",
                                 measurePar = list(lambda.min.ratio = 0.01,
                                                   nlambda = 20,
                                                   ebic.gamma = 1),
                                 sparsMethod = "none", 
                                 dissFunc = "signed", 
                                    seed = 123456)

saveRDS(net_single_GCODA, "Breed/Complete_SriLanka/Complete_SriLanka_GCODA.rds")

props_single_GCODA <- netAnalyze(net_single_GCODA, 
                                    clustMethod = "cluster_fast_greedy",
                                    hubPar = "eigenvector", 
                                    hubQuant = 0.95)

saveRDS(props_single_GCODA, "Breed/Complete_SriLanka/Complete_SriLanka_GCODA analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_GCODA)
net.summary

p_GCODA <- plot(props_single_GCODA,
                   shortenLabels = "none",
                   # labelLength = 16,
                   # charToRm = "g__",
                   labelScale = FALSE,
                   rmSingles = "all",
                   nodeSize = "eigenvector",
                   nodeColor = "cluster",
                   hubBorderCol = "blue",
                   cexNodes = 1,
                   cexLabels = 0.5,
                   edgeWidth = 1,
                   highlightHubs = TRUE,
                   cexHubs = 1.5,
                   # cexHubLabels = 2,
                   title1 = "GCODA Network on Genus level for Sri Lankan Data.", 
                   showTitle = TRUE,
                   cexTitle = 1.5)

saveRDS(p_GCODA, "Breed/Complete_SriLanka/Complete_SriLanka_GCODA plot.rds")


from <- p_GCODA$labels$labels1[p_GCODA$q1$Edgelist$from]
to <- p_GCODA$labels$labels1[p_GCODA$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_GCODA$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_GCODA$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_GCODA$q1$Edgelist$from)))
edges$from <- p_GCODA$q1$Edgelist$from
edges$to <- p_GCODA$q1$Edgelist$to
edges$weight <- p_GCODA$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Complete_SriLanka/Complete_SriLanka_GCODA edge data.csv", row.names = FALSE)

hubs <- props_single_GCODA$hubs$hubs1
write(hubs, "Breed/Complete_SriLanka/Complete_SriLanka_GCODA Hubs.txt")

node.lables <- p_GCODA$labels$labels1
clust <- props_single_GCODA$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_GCODA$centralities$degree1[nodes$lable]
evs=props_single_GCODA$centralities$eigenv1[nodes$lable]
betweennesses=props_single_GCODA$centralities$between1[nodes$lable]
closenesses=props_single_GCODA$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/Complete_SriLanka/Complete_SriLanka_GCODA node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Complete_SriLanka/Complete_SriLanka_GCODA hub data.csv", row.names = FALSE)













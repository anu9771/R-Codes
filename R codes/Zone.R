#01. SparCC
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

#**************************************************************************************************************************

# # filtering the taxa based on prevalance
# phyla2Filter = c("Synergistetes", "OP8", "TM6",
#                  "SR1","GN02","Fibrobacteres",
#                  "Nitrospirae","WPS-2","Planctomycetes")
# # Filter entries with unidentified Phylum.
# complete_physeq = subset_taxa(complete_physeq, !Phylum %in% phyla2Filter)
# #complete_physeq


#Visualization
# sample_data(complete_physeq)
# otu_table(complete_physeq)
# tax_table(complete_physeq)


#aggregating taxa in the Genus level
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

sample_data(complete_physeq)
tax_table(complete_physeq)
otu_table(complete_physeq)


# select data of Different Agricultural Zone.
Wet_physeq <- subset_samples(filt_complete_agg, AgriculturalZone =="Wet")
Dry_physeq <- subset_samples(filt_complete_agg, AgriculturalZone =="Dry")
Intermediate_physeq <- subset_samples(filt_complete_agg, AgriculturalZone =="Intermediate")



# Saving the phyloseq objects for Agricultural Zones
saveRDS(Wet_physeq, "Zone/Wet Zone.rds")
saveRDS(Dry_physeq, "Zone/Dry Zone.rds")
saveRDS(Intermediate_physeq, "Zone/Intermediate Zone.rds")


#Wet Zone******************************************************************************
Wet_data <- readRDS("Zone/Wet Zone.rds")

top_Wet <- prune_taxa(names(sort(taxa_sums(Wet_data),TRUE)[1:50]), Wet_data)
#plot_heatmap(top_Wet)
Wet_data

#SparCC
net_single_Wet <- netConstruct(Wet_data,
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

saveRDS(net_single_Wet, "Zone/Wet Zone/Wet_Sparcc network.rds")

props_single_Wet <- netAnalyze(net_single_Wet, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_Wet, "Zone/Wet Zone/Wet_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Wet)
net.summary

p_Wet <- plot(props_single_Wet,
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
                title1 = "Wet Zone Network on Genus level with SparCC Method.", 
                showTitle = TRUE,
                cexTitle = 1.5)

saveRDS(p_Wet, "Zone/Wet Zone/Wet_Sparcc plot.rds")

from <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$from]
to <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Wet$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Wet$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Wet$q1$Edgelist$from)))
edges$from <- p_Wet$q1$Edgelist$from
edges$to <- p_Wet$q1$Edgelist$to
edges$weight <- p_Wet$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Wet Zone/Wet edge data.csv", row.names = FALSE)

hubs <- props_single_Wet$hubs$hubs1
write(hubs, "Zone/Wet Zone/Wet Zone Hubs.txt")

node.lables <- p_Wet$labels$labels1
clust <- props_single_Wet$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Wet$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Wet Zone/Wet node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Wet Zone/Wet hub data.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_Sparcc_Zone_Wet_Top100=data.frame(Node1=net_single_Wet[["edgelist1"]][["v1"]],
                                          Node2=net_single_Wet[["edgelist1"]][["v2"]],
                                          Association=net_single_Wet[["edgelist1"]][["asso"]],
                                          Dissimilarity=net_single_Wet[["edgelist1"]][["diss"]],
                                          Adjecency=net_single_Wet[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Zone_Wet_Top100=Genus_Sparcc_Zone_Wet_Top100[order(Genus_Sparcc_Zone_Wet_Top100$Association,
                                                                    decreasing = T),]

write.csv(Genus_Sparcc_Zone_Wet_Top100,"Associations/Genus_Sparcc_Zone_Wet_Top100.csv",row.names = FALSE)

############################


#Dry Zone******************************************************************************
Dry_data <- readRDS("Zone/Dry Zone.rds")

top_Dry <- prune_taxa(names(sort(taxa_sums(Dry_data),TRUE)[1:50]), Dry_data)
#plot_heatmap(top_Dry)
Dry_data

#SparCC
net_single_Dry <- netConstruct(Dry_data,
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

saveRDS(net_single_Dry, "Zone/Dry Zone/Dry_Sparcc network.rds")

props_single_Dry <- netAnalyze(net_single_Dry, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Dry, "Zone/Dry Zone/Dry_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Dry)
net.summary

p_Dry <- plot(props_single_Dry,
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
              title1 = "Dry Zone Network on Genus level with SparCC Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Dry, "Zone/Dry Zone/Dry_Sparcc plot.rds")

from <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$from]
to <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Dry$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Dry$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Dry$q1$Edgelist$from)))
edges$from <- p_Dry$q1$Edgelist$from
edges$to <- p_Dry$q1$Edgelist$to
edges$weight <- p_Dry$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Dry Zone/Dry edge data.csv", row.names = FALSE)

hubs <- props_single_Dry$hubs$hubs1
write(hubs, "Zone/Dry Zone/Dry Zone Hubs.txt")

node.lables <- p_Dry$labels$labels1
clust <- props_single_Dry$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Dry$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Dry Zone/Dry node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Dry Zone/Dry hub data.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_Sparcc_Zone_Dry_Top100=data.frame(Node1=net_single_Dry[["edgelist1"]][["v1"]],
                                        Node2=net_single_Dry[["edgelist1"]][["v2"]],
                                        Association=net_single_Dry[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Dry[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Dry[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Zone_Dry_Top100=Genus_Sparcc_Zone_Dry_Top100[order(Genus_Sparcc_Zone_Dry_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_Sparcc_Zone_Dry_Top100,"Associations/Genus_Sparcc_Zone_Dry_Top100.csv",row.names = FALSE)

############################



#Intermediate Zone******************************************************************************
Intermediate_data <- readRDS("Zone/Intermediate Zone.rds")

top_Intermediate <- prune_taxa(names(sort(taxa_sums(Intermediate_data),TRUE)[1:50]), Intermediate_data)
#plot_heatmap(top_Intermediate)
Intermediate_data

#SparCC
net_single_Intermediate <- netConstruct(Intermediate_data,
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

saveRDS(net_single_Intermediate, "Zone/Intermediate Zone/Intermediate_Sparcc network.rds")

props_single_Intermediate <- netAnalyze(net_single_Intermediate, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Intermediate, "Zone/Intermediate Zone/Intermediate_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Intermediate)
net.summary

p_Intermediate <- plot(props_single_Intermediate,
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
              title1 = "Intermediate Zone Network on Genus level with SparCC Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Intermediate, "Zone/Intermediate Zone/Intermediate_Sparcc plot.rds")

from <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$from]
to <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Intermediate$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Intermediate$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Intermediate$q1$Edgelist$from)))
edges$from <- p_Intermediate$q1$Edgelist$from
edges$to <- p_Intermediate$q1$Edgelist$to
edges$weight <- p_Intermediate$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Intermediate Zone/Intermediate edge data.csv", row.names = FALSE)

hubs <- props_single_Intermediate$hubs$hubs1
write(hubs, "Zone/Intermediate Zone/Intermediate Zone Hubs.txt")

node.lables <- p_Intermediate$labels$labels1
clust <- props_single_Intermediate$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Intermediate$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Intermediate Zone/Intermediate node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Intermediate Zone/Intermediate hub data.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_Sparcc_Zone_Intermediate_Top100=data.frame(Node1=net_single_Intermediate[["edgelist1"]][["v1"]],
                                        Node2=net_single_Intermediate[["edgelist1"]][["v2"]],
                                        Association=net_single_Intermediate[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Intermediate[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Intermediate[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Zone_Intermediate_Top100=Genus_Sparcc_Zone_Intermediate_Top100[order(Genus_Sparcc_Zone_Intermediate_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_Sparcc_Zone_Intermediate_Top100,"Associations/Genus_Sparcc_Zone_Intermediate_Top100.csv",row.names = FALSE)

############################
#02.CCREPE

#Wet Zone******************************************************************************
Wet_data <- readRDS("Zone/Wet Zone.rds")

top_Wet <- prune_taxa(names(sort(taxa_sums(Wet_data),TRUE)[1:50]), Wet_data)
#plot_heatmap(top_Wet)
Wet_data

net_single_Wet <- netConstruct(Wet_data, verbose = 3,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "pseudo", normMethod = "fractions",
                           measure = "ccrepe",
                           sparsMethod = "threshold", thresh = 0.4,
                           dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single_Wet, "Zone/Wet Zone/Wet_ccrepe network.rds")

props_single_Wet <- netAnalyze(net_single_Wet, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Wet, "Zone/Wet Zone/Wet_ccrepe analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Wet)
net.summary

p_Wet <- plot(props_single_Wet,
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
              title1 = "Wet Zone Network on Genus level with ccrepe Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Wet, "Zone/Wet Zone/Wet_ccrepe plot.rds")

from <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$from]
to <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Wet$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Wet$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Wet$q1$Edgelist$from)))
edges$from <- p_Wet$q1$Edgelist$from
edges$to <- p_Wet$q1$Edgelist$to
edges$weight <- p_Wet$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Wet Zone/Wet edge data_ccrepe.csv", row.names = FALSE)

hubs <- props_single_Wet$hubs$hubs1
write(hubs, "Zone/Wet Zone/Wet Zone Hubs_ccrepe.txt")

node.lables <- p_Wet$labels$labels1
clust <- props_single_Wet$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Wet$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Wet Zone/Wet node data_ccrepe.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Wet Zone/Wet hub data_ccrepe.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_ccrepe_Zone_Wet_Top100=data.frame(Node1=net_single_Wet[["edgelist1"]][["v1"]],
                                        Node2=net_single_Wet[["edgelist1"]][["v2"]],
                                        Association=net_single_Wet[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Wet[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Wet[["edgelist1"]][["adja"]]
)

Genus_ccrepe_Zone_Wet_Top100=Genus_ccrepe_Zone_Wet_Top100[order(Genus_ccrepe_Zone_Wet_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_ccrepe_Zone_Wet_Top100,"Associations/Genus_ccrepe_Zone_Wet_Top100.csv",row.names = FALSE)

############################



#Dry Zone******************************************************************************
Dry_data <- readRDS("Zone/Dry Zone.rds")

top_Dry <- prune_taxa(names(sort(taxa_sums(Dry_data),TRUE)[1:50]), Dry_data)
#plot_heatmap(top_Dry)
Dry_data

net_single_Dry <- netConstruct(Dry_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "pseudo", normMethod = "fractions",
                               measure = "ccrepe",
                               sparsMethod = "threshold", thresh = 0.4,
                               dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Dry, "Zone/Dry Zone/Dry_ccrepe network.rds")

props_single_Dry <- netAnalyze(net_single_Dry, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Dry, "Zone/Dry Zone/Dry_ccrepe analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Dry)
net.summary

p_Dry <- plot(props_single_Dry,
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
              title1 = "Dry Zone Network on Genus level with ccrepe Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Dry, "Zone/Dry Zone/Dry_ccrepe plot.rds")

from <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$from]
to <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Dry$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Dry$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Dry$q1$Edgelist$from)))
edges$from <- p_Dry$q1$Edgelist$from
edges$to <- p_Dry$q1$Edgelist$to
edges$weight <- p_Dry$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Dry Zone/Dry edge data_ccrepe.csv", row.names = FALSE)

hubs <- props_single_Dry$hubs$hubs1
write(hubs, "Zone/Dry Zone/Dry Zone Hubs_ccrepe.txt")

node.lables <- p_Dry$labels$labels1
clust <- props_single_Dry$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Dry$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Dry Zone/Dry node data_ccrepe.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Dry Zone/Dry hub data_ccrepe.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_ccrepe_Zone_Dry_Top100=data.frame(Node1=net_single_Dry[["edgelist1"]][["v1"]],
                                        Node2=net_single_Dry[["edgelist1"]][["v2"]],
                                        Association=net_single_Dry[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Dry[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Dry[["edgelist1"]][["adja"]]
)

Genus_ccrepe_Zone_Dry_Top100=Genus_ccrepe_Zone_Dry_Top100[order(Genus_ccrepe_Zone_Dry_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_ccrepe_Zone_Dry_Top100,"Associations/Genus_ccrepe_Zone_Dry_Top100.csv",row.names = FALSE)

############################


#Intermediate Zone******************************************************************************
Intermediate_data <- readRDS("Zone/Intermediate Zone.rds")

top_Intermediate <- prune_taxa(names(sort(taxa_sums(Intermediate_data),TRUE)[1:50]), Intermediate_data)
#plot_heatmap(top_Intermediate)
Intermediate_data

net_single_Intermediate <- netConstruct(Intermediate_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "pseudo", normMethod = "fractions",
                               measure = "ccrepe",
                               sparsMethod = "threshold", thresh = 0.4,
                               dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Intermediate, "Zone/Intermediate Zone/Intermediate_ccrepe network.rds")

props_single_Intermediate <- netAnalyze(net_single_Intermediate, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Intermediate, "Zone/Intermediate Zone/Intermediate_ccrepe analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Intermediate)
net.summary

p_Intermediate <- plot(props_single_Intermediate,
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
              title1 = "Intermediate Zone Network on Genus level with ccrepe Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Intermediate, "Zone/Intermediate Zone/Intermediate_ccrepe plot.rds")

from <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$from]
to <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Intermediate$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Intermediate$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Intermediate$q1$Edgelist$from)))
edges$from <- p_Intermediate$q1$Edgelist$from
edges$to <- p_Intermediate$q1$Edgelist$to
edges$weight <- p_Intermediate$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Intermediate Zone/Intermediate edge data_ccrepe.csv", row.names = FALSE)

hubs <- props_single_Intermediate$hubs$hubs1
write(hubs, "Zone/Intermediate Zone/Intermediate Zone Hubs_ccrepe.txt")

node.lables <- p_Intermediate$labels$labels1
clust <- props_single_Intermediate$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Intermediate$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Intermediate Zone/Intermediate node data_ccrepe.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Intermediate Zone/Intermediate hub data_ccrepe.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_ccrepe_Zone_Intermediate_Top100=data.frame(Node1=net_single_Intermediate[["edgelist1"]][["v1"]],
                                        Node2=net_single_Intermediate[["edgelist1"]][["v2"]],
                                        Association=net_single_Intermediate[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Intermediate[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Intermediate[["edgelist1"]][["adja"]]
)

Genus_ccrepe_Zone_Intermediate_Top100=Genus_ccrepe_Zone_Intermediate_Top100[order(Genus_ccrepe_Zone_Intermediate_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_ccrepe_Zone_Intermediate_Top100,"Associations/Genus_ccrepe_Zone_Intermediate_Top100.csv",row.names = FALSE)

############################


#CCLasso
#Wet Zone******************************************************************************
Wet_data <- readRDS("Zone/Wet Zone.rds")

top_Wet <- prune_taxa(names(sort(taxa_sums(Wet_data),TRUE)[1:50]), Wet_data)
#plot_heatmap(top_Wet)
Wet_data

net_single_Wet <- netConstruct(Wet_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "pseudo", normMethod = "none",
                               measure = "cclasso",
                               sparsMethod = "threshold", thresh = 0.4,
                               dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Wet, "Zone/Wet Zone/Wet_cclasso network.rds")

props_single_Wet <- netAnalyze(net_single_Wet, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Wet, "Zone/Wet Zone/Wet_cclasso analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Wet)
net.summary

p_Wet <- plot(props_single_Wet,
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
              title1 = "Wet Zone Network on Genus level with cclasso Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Wet, "Zone/Wet Zone/Wet_cclasso plot.rds")

from <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$from]
to <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Wet$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Wet$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Wet$q1$Edgelist$from)))
edges$from <- p_Wet$q1$Edgelist$from
edges$to <- p_Wet$q1$Edgelist$to
edges$weight <- p_Wet$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Wet Zone/Wet edge data_cclasso.csv", row.names = FALSE)

hubs <- props_single_Wet$hubs$hubs1
write(hubs, "Zone/Wet Zone/Wet Zone Hubs_cclasso.txt")

node.lables <- p_Wet$labels$labels1
clust <- props_single_Wet$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Wet$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Wet Zone/Wet node data_cclasso.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Wet Zone/Wet hub data_cclasso.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_cclasso_Zone_Wet_Top100=data.frame(Node1=net_single_Wet[["edgelist1"]][["v1"]],
                                        Node2=net_single_Wet[["edgelist1"]][["v2"]],
                                        Association=net_single_Wet[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Wet[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Wet[["edgelist1"]][["adja"]]
)

Genus_cclasso_Zone_Wet_Top100=Genus_cclasso_Zone_Wet_Top100[order(Genus_cclasso_Zone_Wet_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_cclasso_Zone_Wet_Top100,"Associations/Genus_cclasso_Zone_Wet_Top100.csv",row.names = FALSE)

############################



#Dry Zone******************************************************************************
Dry_data <- readRDS("Zone/Dry Zone.rds")

top_Dry <- prune_taxa(names(sort(taxa_sums(Dry_data),TRUE)[1:50]), Dry_data)
#plot_heatmap(top_Dry)
Dry_data

net_single_Dry <- netConstruct(Dry_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "pseudo", normMethod = "none",
                               measure = "cclasso",
                               sparsMethod = "threshold", thresh = 0.4,
                               dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Dry, "Zone/Dry Zone/Dry_cclasso network.rds")

props_single_Dry <- netAnalyze(net_single_Dry, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Dry, "Zone/Dry Zone/Dry_cclasso analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Dry)
net.summary

p_Dry <- plot(props_single_Dry,
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
              title1 = "Dry Zone Network on Genus level with cclasso Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Dry, "Zone/Dry Zone/Dry_cclasso plot.rds")

from <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$from]
to <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Dry$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Dry$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Dry$q1$Edgelist$from)))
edges$from <- p_Dry$q1$Edgelist$from
edges$to <- p_Dry$q1$Edgelist$to
edges$weight <- p_Dry$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Dry Zone/Dry edge data_cclasso.csv", row.names = FALSE)

hubs <- props_single_Dry$hubs$hubs1
write(hubs, "Zone/Dry Zone/Dry Zone Hubs_cclasso.txt")

node.lables <- p_Dry$labels$labels1
clust <- props_single_Dry$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Dry$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Dry Zone/Dry node data_cclasso.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Dry Zone/Dry hub data_cclasso.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_cclasso_Zone_Dry_Top100=data.frame(Node1=net_single_Dry[["edgelist1"]][["v1"]],
                                         Node2=net_single_Dry[["edgelist1"]][["v2"]],
                                         Association=net_single_Dry[["edgelist1"]][["asso"]],
                                         Dissimilarity=net_single_Dry[["edgelist1"]][["diss"]],
                                         Adjecency=net_single_Dry[["edgelist1"]][["adja"]]
)

Genus_cclasso_Zone_Dry_Top100=Genus_cclasso_Zone_Dry_Top100[order(Genus_cclasso_Zone_Dry_Top100$Association,
                                                                  decreasing = T),]

write.csv(Genus_cclasso_Zone_Dry_Top100,"Associations/Genus_cclasso_Zone_Dry_Top100.csv",row.names = FALSE)

############################





#Intermediate Zone******************************************************************************
Intermediate_data <- readRDS("Zone/Intermediate Zone.rds")

top_Intermediate <- prune_taxa(names(sort(taxa_sums(Intermediate_data),TRUE)[1:50]), Intermediate_data)
#plot_heatmap(top_Intermediate)
Intermediate_data

net_single_Intermediate <- netConstruct(Intermediate_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "pseudo", normMethod = "none",
                               measure = "cclasso",
                               sparsMethod = "threshold", thresh = 0.4,
                               dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Intermediate, "Zone/Intermediate Zone/Intermediate_cclasso network.rds")

props_single_Intermediate <- netAnalyze(net_single_Intermediate, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Intermediate, "Zone/Intermediate Zone/Intermediate_cclasso analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Intermediate)
net.summary

p_Intermediate <- plot(props_single_Intermediate,
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
              title1 = "Intermediate Zone Network on Genus level with cclasso Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Intermediate, "Zone/Intermediate Zone/Intermediate_cclasso plot.rds")

from <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$from]
to <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Intermediate$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Intermediate$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Intermediate$q1$Edgelist$from)))
edges$from <- p_Intermediate$q1$Edgelist$from
edges$to <- p_Intermediate$q1$Edgelist$to
edges$weight <- p_Intermediate$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Intermediate Zone/Intermediate edge data_cclasso.csv", row.names = FALSE)

hubs <- props_single_Intermediate$hubs$hubs1
write(hubs, "Zone/Intermediate Zone/Intermediate Zone Hubs_cclasso.txt")

node.lables <- p_Intermediate$labels$labels1
clust <- props_single_Intermediate$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Intermediate$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Intermediate Zone/Intermediate node data_cclasso.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Intermediate Zone/Intermediate hub data_cclasso.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_cclasso_Zone_Intermediate_Top100=data.frame(Node1=net_single_Intermediate[["edgelist1"]][["v1"]],
                                         Node2=net_single_Intermediate[["edgelist1"]][["v2"]],
                                         Association=net_single_Intermediate[["edgelist1"]][["asso"]],
                                         Dissimilarity=net_single_Intermediate[["edgelist1"]][["diss"]],
                                         Adjecency=net_single_Intermediate[["edgelist1"]][["adja"]]
)

Genus_cclasso_Zone_Intermediate_Top100=Genus_cclasso_Zone_Intermediate_Top100[order(Genus_cclasso_Zone_Intermediate_Top100$Association,
                                                                  decreasing = T),]

write.csv(Genus_cclasso_Zone_Intermediate_Top100,"Associations/Genus_cclasso_Zone_Intermediate_Top100.csv",row.names = FALSE)

############################



#Spiec-Easi
#Wet Zone******************************************************************************
Wet_data <- readRDS("Zone/Wet Zone.rds")

top_Wet <- prune_taxa(names(sort(taxa_sums(Wet_data),TRUE)[1:50]), Wet_data)
#plot_heatmap(top_Wet)
Wet_data

net_single_Wet <- netConstruct(Wet_data, verbose = 3,
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
                                                                      subsample.ratio = 0.8)),
                               sparsMethod = "none", dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Wet, "Zone/Wet Zone/Wet_spiec network.rds")

props_single_Wet <- netAnalyze(net_single_Wet, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Wet, "Zone/Wet Zone/Wet_spiec analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Wet)
net.summary

p_Wet <- plot(props_single_Wet,
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
              title1 = "Wet Zone Network on Genus level with spiec Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Wet, "Zone/Wet Zone/Wet_spiec plot.rds")

from <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$from]
to <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Wet$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Wet$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Wet$q1$Edgelist$from)))
edges$from <- p_Wet$q1$Edgelist$from
edges$to <- p_Wet$q1$Edgelist$to
edges$weight <- p_Wet$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Wet Zone/Wet edge data_spiec.csv", row.names = FALSE)

hubs <- props_single_Wet$hubs$hubs1
write(hubs, "Zone/Wet Zone/Wet Zone Hubs_spiec.txt")

node.lables <- p_Wet$labels$labels1
clust <- props_single_Wet$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Wet$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Wet Zone/Wet node data_spiec.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Wet Zone/Wet hub data_spiec.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_spiec_Zone_Wet_Top100=data.frame(Node1=net_single_Wet[["edgelist1"]][["v1"]],
                                         Node2=net_single_Wet[["edgelist1"]][["v2"]],
                                         Association=net_single_Wet[["edgelist1"]][["asso"]],
                                         Dissimilarity=net_single_Wet[["edgelist1"]][["diss"]],
                                         Adjecency=net_single_Wet[["edgelist1"]][["adja"]]
)

Genus_spiec_Zone_Wet_Top100=Genus_spiec_Zone_Wet_Top100[order(Genus_spiec_Zone_Wet_Top100$Association,
                                                                  decreasing = T),]

write.csv(Genus_spiec_Zone_Wet_Top100,"Associations/Genus_spiec_Zone_Wet_Top100.csv",row.names = FALSE)

############################


#Dry Zone******************************************************************************
Dry_data <- readRDS("Zone/Dry Zone.rds")

top_Dry <- prune_taxa(names(sort(taxa_sums(Dry_data),TRUE)[1:50]), Dry_data)
#plot_heatmap(top_Dry)
Dry_data

net_single_Dry <- netConstruct(Dry_data, verbose = 3,
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
                                                                      subsample.ratio = 0.8)),
                               sparsMethod = "none", dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Dry, "Zone/Dry Zone/Dry_spiec network.rds")

props_single_Dry <- netAnalyze(net_single_Dry, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Dry, "Zone/Dry Zone/Dry_spiec analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Dry)
net.summary

p_Dry <- plot(props_single_Dry,
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
              title1 = "Dry Zone Network on Genus level with spiec Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Dry, "Zone/Dry Zone/Dry_spiec plot.rds")

from <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$from]
to <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Dry$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Dry$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Dry$q1$Edgelist$from)))
edges$from <- p_Dry$q1$Edgelist$from
edges$to <- p_Dry$q1$Edgelist$to
edges$weight <- p_Dry$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Dry Zone/Dry edge data_spiec.csv", row.names = FALSE)

hubs <- props_single_Dry$hubs$hubs1
write(hubs, "Zone/Dry Zone/Dry Zone Hubs_spiec.txt")

node.lables <- p_Dry$labels$labels1
clust <- props_single_Dry$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Dry$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Dry Zone/Dry node data_spiec.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Dry Zone/Dry hub data_spiec.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_spiec_Zone_Dry_Top100=data.frame(Node1=net_single_Dry[["edgelist1"]][["v1"]],
                                       Node2=net_single_Dry[["edgelist1"]][["v2"]],
                                       Association=net_single_Dry[["edgelist1"]][["asso"]],
                                       Dissimilarity=net_single_Dry[["edgelist1"]][["diss"]],
                                       Adjecency=net_single_Dry[["edgelist1"]][["adja"]]
)

Genus_spiec_Zone_Dry_Top100=Genus_spiec_Zone_Dry_Top100[order(Genus_spiec_Zone_Dry_Top100$Association,
                                                              decreasing = T),]

write.csv(Genus_spiec_Zone_Dry_Top100,"Associations/Genus_spiec_Zone_Dry_Top100.csv",row.names = FALSE)

############################


#Intermediate Zone******************************************************************************
Intermediate_data <- readRDS("Zone/Intermediate Zone.rds")

top_Intermediate <- prune_taxa(names(sort(taxa_sums(Intermediate_data),TRUE)[1:50]), Intermediate_data)
#plot_heatmap(top_Intermediate)
Intermediate_data

net_single_Intermediate <- netConstruct(Intermediate_data, verbose = 3,
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
                                                                      subsample.ratio = 0.8)),
                               sparsMethod = "none", dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Intermediate, "Zone/Intermediate Zone/Intermediate_spiec network.rds")

props_single_Intermediate <- netAnalyze(net_single_Intermediate, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Intermediate, "Zone/Intermediate Zone/Intermediate_spiec analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Intermediate)
net.summary

p_Intermediate <- plot(props_single_Intermediate,
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
              title1 = "Intermediate Zone Network on Genus level with spiec Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Intermediate, "Zone/Intermediate Zone/Intermediate_spiec plot.rds")

from <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$from]
to <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Intermediate$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Intermediate$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Intermediate$q1$Edgelist$from)))
edges$from <- p_Intermediate$q1$Edgelist$from
edges$to <- p_Intermediate$q1$Edgelist$to
edges$weight <- p_Intermediate$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Intermediate Zone/Intermediate edge data_spiec.csv", row.names = FALSE)

hubs <- props_single_Intermediate$hubs$hubs1
write(hubs, "Zone/Intermediate Zone/Intermediate Zone Hubs_spiec.txt")

node.lables <- p_Intermediate$labels$labels1
clust <- props_single_Intermediate$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Intermediate$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Intermediate Zone/Intermediate node data_spiec.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Intermediate Zone/Intermediate hub data_spiec.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_spiec_Zone_Intermediate_Top100=data.frame(Node1=net_single_Intermediate[["edgelist1"]][["v1"]],
                                       Node2=net_single_Intermediate[["edgelist1"]][["v2"]],
                                       Association=net_single_Intermediate[["edgelist1"]][["asso"]],
                                       Dissimilarity=net_single_Intermediate[["edgelist1"]][["diss"]],
                                       Adjecency=net_single_Intermediate[["edgelist1"]][["adja"]]
)

Genus_spiec_Zone_Intermediate_Top100=Genus_spiec_Zone_Intermediate_Top100[order(Genus_spiec_Zone_Intermediate_Top100$Association,
                                                              decreasing = T),]

write.csv(Genus_spiec_Zone_Intermediate_Top100,"Associations/Genus_spiec_Zone_Intermediate_Top100.csv",row.names = FALSE)

############################



#gCoda
#Wet Zone******************************************************************************
Wet_data <- readRDS("Zone/Wet Zone.rds")

top_Wet <- prune_taxa(names(sort(taxa_sums(Wet_data),TRUE)[1:50]), Wet_data)
#plot_heatmap(top_Wet)
Wet_data

net_single_Wet <- netConstruct(Wet_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "pseudo", normMethod = "none",
                               measure = "gcoda",
                               measurePar = list(lambda.min.ratio = 0.01,
                                                 nlambda = 20,
                                                 ebic.gamma = 1),
                               sparsMethod = "none", 
                               dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Wet, "Zone/Wet Zone/Wet_gcoda network.rds")

props_single_Wet <- netAnalyze(net_single_Wet, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Wet, "Zone/Wet Zone/Wet_gcoda analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Wet)
net.summary

p_Wet <- plot(props_single_Wet,
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
              title1 = "Wet Zone Network on Genus level with gcoda Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Wet, "Zone/Wet Zone/Wet_gcoda plot.rds")

from <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$from]
to <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Wet$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Wet$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Wet$q1$Edgelist$from)))
edges$from <- p_Wet$q1$Edgelist$from
edges$to <- p_Wet$q1$Edgelist$to
edges$weight <- p_Wet$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Wet Zone/Wet edge data_gcoda.csv", row.names = FALSE)

hubs <- props_single_Wet$hubs$hubs1
write(hubs, "Zone/Wet Zone/Wet Zone Hubs_gcoda.txt")

node.lables <- p_Wet$labels$labels1
clust <- props_single_Wet$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Wet$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Wet Zone/Wet node data_gcoda.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Wet Zone/Wet hub data_gcoda.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_gcoda_Zone_Wet_Top100=data.frame(Node1=net_single_Wet[["edgelist1"]][["v1"]],
                                       Node2=net_single_Wet[["edgelist1"]][["v2"]],
                                       Association=net_single_Wet[["edgelist1"]][["asso"]],
                                       Dissimilarity=net_single_Wet[["edgelist1"]][["diss"]],
                                       Adjecency=net_single_Wet[["edgelist1"]][["adja"]]
)

Genus_gcoda_Zone_Wet_Top100=Genus_gcoda_Zone_Wet_Top100[order(Genus_gcoda_Zone_Wet_Top100$Association,
                                                              decreasing = T),]

write.csv(Genus_gcoda_Zone_Wet_Top100,"Associations/Genus_gcoda_Zone_Wet_Top100.csv",row.names = FALSE)

############################


#Dry Zone******************************************************************************
Dry_data <- readRDS("Zone/Dry Zone.rds")

top_Dry <- prune_taxa(names(sort(taxa_sums(Dry_data),TRUE)[1:50]), Dry_data)
#plot_heatmap(top_Dry)
Dry_data

net_single_Dry <- netConstruct(Dry_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "pseudo", normMethod = "none",
                               measure = "gcoda",
                               measurePar = list(lambda.min.ratio = 0.01,
                                                 nlambda = 20,
                                                 ebic.gamma = 1),
                               sparsMethod = "none", 
                               dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Dry, "Zone/Dry Zone/Dry_gcoda network.rds")

props_single_Dry <- netAnalyze(net_single_Dry, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Dry, "Zone/Dry Zone/Dry_gcoda analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Dry)
net.summary

p_Dry <- plot(props_single_Dry,
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
              title1 = "Dry Zone Network on Genus level with gcoda Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Dry, "Zone/Dry Zone/Dry_gcoda plot.rds")

from <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$from]
to <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Dry$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Dry$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Dry$q1$Edgelist$from)))
edges$from <- p_Dry$q1$Edgelist$from
edges$to <- p_Dry$q1$Edgelist$to
edges$weight <- p_Dry$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Dry Zone/Dry edge data_gcoda.csv", row.names = FALSE)

hubs <- props_single_Dry$hubs$hubs1
write(hubs, "Zone/Dry Zone/Dry Zone Hubs_gcoda.txt")

node.lables <- p_Dry$labels$labels1
clust <- props_single_Dry$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Dry$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Dry Zone/Dry node data_gcoda.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Dry Zone/Dry hub data_gcoda.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_gcoda_Zone_Dry_Top100=data.frame(Node1=net_single_Dry[["edgelist1"]][["v1"]],
                                       Node2=net_single_Dry[["edgelist1"]][["v2"]],
                                       Association=net_single_Dry[["edgelist1"]][["asso"]],
                                       Dissimilarity=net_single_Dry[["edgelist1"]][["diss"]],
                                       Adjecency=net_single_Dry[["edgelist1"]][["adja"]]
)

Genus_gcoda_Zone_Dry_Top100=Genus_gcoda_Zone_Dry_Top100[order(Genus_gcoda_Zone_Dry_Top100$Association,
                                                              decreasing = T),]

write.csv(Genus_gcoda_Zone_Dry_Top100,"Associations/Genus_gcoda_Zone_Dry_Top100.csv",row.names = FALSE)

############################


#Intermediate Zone******************************************************************************
Intermediate_data <- readRDS("Zone/Intermediate Zone.rds")

top_Intermediate <- prune_taxa(names(sort(taxa_sums(Intermediate_data),TRUE)[1:50]), Intermediate_data)
#plot_heatmap(top_Intermediate)
Intermediate_data

net_single_Intermediate <- netConstruct(Intermediate_data, verbose = 3,
                                        filtTax = "highestFreq",
                                        filtTaxPar = list(highestFreq = 100),
                                        filtSamp = "totalReads",
                                        filtSampPar = list(totalReads = 1000),
                                        zeroMethod = "pseudo", normMethod = "none",
                                        measure = "gcoda",
                                        measurePar = list(lambda.min.ratio = 0.01,
                                                          nlambda = 20,
                                                          ebic.gamma = 1),
                                        sparsMethod = "none", 
                                        dissFunc = "signed", 
                                        seed = 123456)

saveRDS(net_single_Intermediate, "Zone/Intermediate Zone/Intermediate_gcoda network.rds")

props_single_Intermediate <- netAnalyze(net_single_Intermediate, 
                                        clustMethod = "cluster_fast_greedy",
                                        hubPar = "eigenvector", 
                                        hubQuant = 0.95)

saveRDS(props_single_Intermediate, "Zone/Intermediate Zone/Intermediate_gcoda analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Intermediate)
net.summary

p_Intermediate <- plot(props_single_Intermediate,
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
                       title1 = "Intermediate Zone Network on Genus level with gcoda Method.", 
                       showTitle = TRUE,
                       cexTitle = 1.5)

saveRDS(p_Intermediate, "Zone/Intermediate Zone/Intermediate_gcoda plot.rds")

from <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$from]
to <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Intermediate$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Intermediate$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Intermediate$q1$Edgelist$from)))
edges$from <- p_Intermediate$q1$Edgelist$from
edges$to <- p_Intermediate$q1$Edgelist$to
edges$weight <- p_Intermediate$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Intermediate Zone/Intermediate edge data_gcoda.csv", row.names = FALSE)

hubs <- props_single_Intermediate$hubs$hubs1
write(hubs, "Zone/Intermediate Zone/Intermediate Zone Hubs_gcoda.txt")

node.lables <- p_Intermediate$labels$labels1
clust <- props_single_Intermediate$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Intermediate$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Intermediate Zone/Intermediate node data_gcoda.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Intermediate Zone/Intermediate hub data_gcoda.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_gcoda_Zone_Intermediate_Top100=data.frame(Node1=net_single_Intermediate[["edgelist1"]][["v1"]],
                                                Node2=net_single_Intermediate[["edgelist1"]][["v2"]],
                                                Association=net_single_Intermediate[["edgelist1"]][["asso"]],
                                                Dissimilarity=net_single_Intermediate[["edgelist1"]][["diss"]],
                                                Adjecency=net_single_Intermediate[["edgelist1"]][["adja"]]
)

Genus_gcoda_Zone_Intermediate_Top100=Genus_gcoda_Zone_Intermediate_Top100[order(Genus_gcoda_Zone_Intermediate_Top100$Association,
                                                                                decreasing = T),]

write.csv(Genus_gcoda_Zone_Intermediate_Top100,"Associations/Genus_gcoda_Zone_Intermediate_Top100.csv",row.names = FALSE)

############################






#Spring
#Wet Zone******************************************************************************
Wet_data <- readRDS("Zone/Wet Zone.rds")

top_Wet <- prune_taxa(names(sort(taxa_sums(Wet_data),TRUE)[1:50]), Wet_data)
#plot_heatmap(top_Wet)
Wet_data

net_single_Wet <- netConstruct(Wet_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "none", normMethod = "none",
                               measure = "spring",
                               measurePar = list(nlambda = 20, 
                                                 rep.num = 20,
                                                 thresh = 0.1, 
                                                 subsample.ratio = 0.8, # 10*sqrt(n)/n for n > 144
                                                 lambda.min.ratio = 0.01,
                                                 lambdaseq = "data-specific",
                                                 ncores=1),
                               sparsMethod = "none",  dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Wet, "Zone/Wet Zone/Wet_spring network.rds")

props_single_Wet <- netAnalyze(net_single_Wet, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Wet, "Zone/Wet Zone/Wet_spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Wet)
net.summary

p_Wet <- plot(props_single_Wet,
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
              title1 = "Wet Zone Network on Genus level with spring Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Wet, "Zone/Wet Zone/Wet_spring plot.rds")

from <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$from]
to <- p_Wet$labels$labels1[p_Wet$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Wet$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Wet$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Wet$q1$Edgelist$from)))
edges$from <- p_Wet$q1$Edgelist$from
edges$to <- p_Wet$q1$Edgelist$to
edges$weight <- p_Wet$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Wet Zone/Wet edge data_spring.csv", row.names = FALSE)

hubs <- props_single_Wet$hubs$hubs1
write(hubs, "Zone/Wet Zone/Wet Zone Hubs_spring.txt")

node.lables <- p_Wet$labels$labels1
clust <- props_single_Wet$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Wet$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Wet Zone/Wet node data_spring.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Wet Zone/Wet hub data_spring.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_spring_Zone_Wet_Top100=data.frame(Node1=net_single_Wet[["edgelist1"]][["v1"]],
                                       Node2=net_single_Wet[["edgelist1"]][["v2"]],
                                       Association=net_single_Wet[["edgelist1"]][["asso"]],
                                       Dissimilarity=net_single_Wet[["edgelist1"]][["diss"]],
                                       Adjecency=net_single_Wet[["edgelist1"]][["adja"]]
)

Genus_spring_Zone_Wet_Top100=Genus_spring_Zone_Wet_Top100[order(Genus_spring_Zone_Wet_Top100$Association,
                                                              decreasing = T),]

write.csv(Genus_spring_Zone_Wet_Top100,"Associations/Genus_spring_Zone_Wet_Top100.csv",row.names = FALSE)

############################




#Dry Zone******************************************************************************
Dry_data <- readRDS("Zone/Dry Zone.rds")

top_Dry <- prune_taxa(names(sort(taxa_sums(Dry_data),TRUE)[1:50]), Dry_data)
#plot_heatmap(top_Dry)
Dry_data

net_single_Dry <- netConstruct(Dry_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "none", normMethod = "none",
                               measure = "spring",
                               measurePar = list(nlambda = 20, 
                                                 rep.num = 20,
                                                 thresh = 0.1, 
                                                 subsample.ratio = 0.8, # 10*sqrt(n)/n for n > 144
                                                 lambda.min.ratio = 0.01,
                                                 lambdaseq = "data-specific",
                                                 ncores=1),
                               sparsMethod = "none",  dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Dry, "Zone/Dry Zone/Dry_spring network.rds")

props_single_Dry <- netAnalyze(net_single_Dry, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Dry, "Zone/Dry Zone/Dry_spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Dry)
net.summary

p_Dry <- plot(props_single_Dry,
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
              title1 = "Dry Zone Network on Genus level with spring Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Dry, "Zone/Dry Zone/Dry_spring plot.rds")

from <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$from]
to <- p_Dry$labels$labels1[p_Dry$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Dry$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Dry$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Dry$q1$Edgelist$from)))
edges$from <- p_Dry$q1$Edgelist$from
edges$to <- p_Dry$q1$Edgelist$to
edges$weight <- p_Dry$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Dry Zone/Dry edge data_spring.csv", row.names = FALSE)

hubs <- props_single_Dry$hubs$hubs1
write(hubs, "Zone/Dry Zone/Dry Zone Hubs_spring.txt")

node.lables <- p_Dry$labels$labels1
clust <- props_single_Dry$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Dry$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Dry Zone/Dry node data_spring.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Dry Zone/Dry hub data_spring.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_spring_Zone_Dry_Top100=data.frame(Node1=net_single_Dry[["edgelist1"]][["v1"]],
                                        Node2=net_single_Dry[["edgelist1"]][["v2"]],
                                        Association=net_single_Dry[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Dry[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Dry[["edgelist1"]][["adja"]]
)

Genus_spring_Zone_Dry_Top100=Genus_spring_Zone_Dry_Top100[order(Genus_spring_Zone_Dry_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_spring_Zone_Dry_Top100,"Associations/Genus_spring_Zone_Dry_Top100.csv",row.names = FALSE)

############################
















#Intermediate Zone******************************************************************************
Intermediate_data <- readRDS("Zone/Intermediate Zone.rds")

top_Intermediate <- prune_taxa(names(sort(taxa_sums(Intermediate_data),TRUE)[1:50]), Intermediate_data)
#plot_heatmap(top_Intermediate)
Intermediate_data

net_single_Intermediate <- netConstruct(Intermediate_data, verbose = 3,
                               filtTax = "highestFreq",
                               filtTaxPar = list(highestFreq = 100),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "none", normMethod = "none",
                               measure = "spring",
                               measurePar = list(nlambda = 20, 
                                                 rep.num = 20,
                                                 thresh = 0.1, 
                                                 subsample.ratio = 0.8, # 10*sqrt(n)/n for n > 144
                                                 lambda.min.ratio = 0.01,
                                                 lambdaseq = "data-specific",
                                                 ncores=1),
                               sparsMethod = "none",  dissFunc = "signed", 
                               seed = 123456)

saveRDS(net_single_Intermediate, "Zone/Intermediate Zone/Intermediate_spring network.rds")

props_single_Intermediate <- netAnalyze(net_single_Intermediate, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Intermediate, "Zone/Intermediate Zone/Intermediate_spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Intermediate)
net.summary

p_Intermediate <- plot(props_single_Intermediate,
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
              title1 = "Intermediate Zone Network on Genus level with spring Method.", 
              showTitle = TRUE,
              cexTitle = 1.5)

saveRDS(p_Intermediate, "Zone/Intermediate Zone/Intermediate_spring plot.rds")

from <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$from]
to <- p_Intermediate$labels$labels1[p_Intermediate$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Intermediate$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Intermediate$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Intermediate$q1$Edgelist$from)))
edges$from <- p_Intermediate$q1$Edgelist$from
edges$to <- p_Intermediate$q1$Edgelist$to
edges$weight <- p_Intermediate$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Zone/Intermediate Zone/Intermediate edge data_spring.csv", row.names = FALSE)

hubs <- props_single_Intermediate$hubs$hubs1
write(hubs, "Zone/Intermediate Zone/Intermediate Zone Hubs_spring.txt")

node.lables <- p_Intermediate$labels$labels1
clust <- props_single_Intermediate$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Intermediate$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Zone/Intermediate Zone/Intermediate node data_spring.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Zone/Intermediate Zone/Intermediate hub data_spring.csv", row.names = FALSE)

##########################################################################################################
#For Associations

Genus_spring_Zone_Intermediate_Top100=data.frame(Node1=net_single_Intermediate[["edgelist1"]][["v1"]],
                                        Node2=net_single_Intermediate[["edgelist1"]][["v2"]],
                                        Association=net_single_Intermediate[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Intermediate[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Intermediate[["edgelist1"]][["adja"]]
)

Genus_spring_Zone_Intermediate_Top100=Genus_spring_Zone_Intermediate_Top100[order(Genus_spring_Zone_Intermediate_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_spring_Zone_Intermediate_Top100,"Associations/Genus_spring_Zone_Intermediate_Top100.csv",row.names = FALSE)

############################





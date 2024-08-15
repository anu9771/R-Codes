# Loading required packages
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(BiocManager)
library(microbiome)
library(NetCoMi)
library(WGCNA)

# Setting the working directory
setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology")


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


#aggregating taxa in the genus level
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

#Save as .rds file
saveRDS(filt_complete_agg, "Phase/All data.rds")

#Sparcc(Top100) Network************************************************************************
All_data <- readRDS("Phase/All data.rds")

top_All <- prune_taxa(names(sort(taxa_sums(All_data),TRUE)[1:50]), All_data)
#plot_heatmap(top_All)
All_data

#SparCC
net_single_All <- netConstruct(All_data,
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

saveRDS(net_single_All, "Phase/All Data/All_Sparcc network.rds")

props_single_All <- netAnalyze(net_single_All, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_All, "Phase/All Data/All_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_All)
net.summary

p_All <- plot(props_single_All,
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
                title1 = "Total Network on genus level with SparCC Method.", 
                showTitle = TRUE,
                cexTitle = 1.5)

saveRDS(p_All, "Phase/All Data/All_Sparcc plot.rds")

from <- p_All$labels$labels1[p_All$q1$Edgelist$from]
to <- p_All$labels$labels1[p_All$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_All$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_All$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_All$q1$Edgelist$from)))
edges$from <- p_All$q1$Edgelist$from
edges$to <- p_All$q1$Edgelist$to
edges$weight <- p_All$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/All Data/All edge data.csv", row.names = FALSE)

hubs <- props_single_All$hubs$hubs1
write(hubs, "Phase/All Data/All Phase Hubs.txt")

node.lables <- p_All$labels$labels1
clust <- props_single_All$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_All$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Phase/All Data/All node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/All Data/All hub data.csv", row.names = FALSE)

#Spring(Top100) Network**********************************************************************
All_data <- readRDS("Phase/All data.rds")

top_All <- prune_taxa(names(sort(taxa_sums(All_data),TRUE)[1:50]), All_data)
#plot_heatmap(top_All)
All_data

#SpRING
net_single_All_Spring <- netConstruct(All_data, verbose = 3,
                                        filtTax = "highestFreq",
                                        filtTaxPar = list(highestFreq = 100),
                                        #filtTax = "none",
                                        #filtTaxPar = "totalReads",
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

saveRDS(net_single_All_Spring, "Phase/All Data/All_Spring network.rds")


# ?netAnalyze
props_single_All_Spring <- netAnalyze(net_single_All_Spring, 
                                        centrLCC = TRUE,
                                        clustMethod = "cluster_fast_greedy",
                                        hubPar = "eigenvector",
                                        weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_All_Spring, "Phase/All Data/All_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_All_Spring)
net.summary

# ?plot.microNetProps
p_All_Spring <- plot(props_single_All_Spring,
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
                       title1 = "Total Network on genus level with Spring Method.", 
                       showTitle = TRUE,
                       cexTitle = 1.5)

saveRDS(p_All_Spring, "Phase/All Data/All_Spring plot.rds")

from <- p_All_Spring$labels$labels1[p_All_Spring$q1$Edgelist$from]
to <- p_All_Spring$labels$labels1[p_All_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_All_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_All_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_All_Spring$q1$Edgelist$from)))
edges$from <- p_All_Spring$q1$Edgelist$from
edges$to <- p_All_Spring$q1$Edgelist$to
edges$weight <- p_All_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/All Data/All_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_All_Spring$hubs$hubs1
write(hubs, "Phase/All Data/All_Spring Hubs.txt")

node.lables <- p_All_Spring$labels$labels1
clust <- props_single_All_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_All_Spring$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Phase/All Data/All_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/All Data/All_Spring hub data.csv", row.names = FALSE)



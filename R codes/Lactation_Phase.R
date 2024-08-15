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


# select data of Different Lactation Phases.
Early_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Early")
Late_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Late")
Mid_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Mid")
Dry_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Dry")

# Saving the phyloseq objects
saveRDS(complete_physeq, "Phase/Complete data.rds")
saveRDS(complete_agg, "Phase/Aggregated data.rds")
saveRDS(filt_complete_agg, "Phase/Filtered Aggregated data.rds")
saveRDS(Early_physeq, "Phase/Early Phase.rds")
saveRDS(Late_physeq, "Phase/Late Phase.rds")
saveRDS(Mid_physeq, "Phase/Mid Phase.rds")
saveRDS(Dry_physeq, "Phase/Dry Phase.rds")

#Early Lactation************************************************************************
Early_data <- readRDS("Phase/Early Phase.rds")

top_Early <- prune_taxa(names(sort(taxa_sums(Early_data),TRUE)[1:50]), Early_data)
plot_heatmap(top_Early)
Early_data

# ?netConstruct

#SparCC
net_single <- netConstruct(Early_data,
                           verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "sparcc",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 20221031)

saveRDS(net_single, "Phase/Early Phase/Early_Sparcc network.rds")


# CCLasso
net_single <- netConstruct(Early_data, verbose = 3,
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 200),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "pseudo", normMethod = "none",
                           measure = "cclasso",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 20221031)

saveRDS(net_single, "Phase/Early Phase/Early_Sparcc network.rds")


props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Phase/Early Phase/Early_Sparcc analysis.rds")


p <- plot(props_single,
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
          title1 = "Early Phase Network on genus level with CCLasso Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Early Phase/Early_Sparcc plot.rds")



#SpRING
net_single <- netConstruct(Early_data, verbose = 3,
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

saveRDS(net_single, "Phase/Early Phase/Early_Spring network.rds")


# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single, "Phase/Early Phase/Early_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)
net.summary
tax_table(Early_data)

# ?plot.microNetProps
p <- plot(props_single,
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
          title1 = "Early Phase Network on genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Early Phase/Early_Spring plot.rds")


# SpiecEasi
net_single <- netConstruct(Early_data, verbose = 3,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           #filtTax = "none",
                           #filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spieceasi",
                           measurePar = list(method = "glasso", 
                                             nlambda = 20,
                                             lambda.min.ratio = 0.01,
                                             pulsar.params = list(rep.num = 20,
                                                                 thresh = 0.1, 
                                                                 subsample.ratio = 0.8)),
                           sparsMethod = "none", dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Early Phase/Early_SpiecEasi network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single, "Phase/Early Phase/Early_SpiecEasi analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)


# Identify strongest positive and negative correlations





# ?plot.microNetProps
p <- plot(props_single,
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
          title1 = "Early Phase Network on genus level with SPIEC-EASI Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Early Phase/Early_SpiecEasi plot.rds")

from <- p$labels$labels1[p$q1$Edgelist$from]
to <- p$labels$labels1[p$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p$q1$Edgelist$from)))
edges$from <- p$q1$Edgelist$from
edges$to <- p$q1$Edgelist$to
edges$weight <- p$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/Early Phase/Early edge data.csv", row.names = FALSE)

hubs <- props_single$hubs$hubs1
write(hubs, "Phase/Early Phase/Early Phase Hubs.txt")

node.lables <- p$labels$labels1
clust <- props_single$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Phase/Early Phase/Early node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/Early Phase/Early hub data.csv", row.names = FALSE)

#Mid Lactation**************************************************************************
Mid_data <- readRDS("Phase/Mid Phase.rds")
Mid_data
top_Mid <- prune_taxa(names(sort(taxa_sums(Mid_data),TRUE)[1:50]), Mid_data)
top_Mid
#plot_heatmap(top_Mid)

# ?netConstruct

#SparCC
net_single <- netConstruct(Mid_data,
                           verbose = 3,
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "sparcc",
                           sparsMethod = "threshold", thresh = 0.4,
                           dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Mid Phase/Mid network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Phase/Mid Phase/Mid analysis.rds")

p <- plot(props_single,
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
          title1 = "Mid Phase Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Mid Phase/Mid plot.rds")

#Still has to do !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#SpRING
net_single <- netConstruct(Mid_data, verbose = 3,
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spring",
                           measurePar = list(nlambda = 20, 
                                             rep.num = 20,
                                             thresh = 0.1, 
                                             subsample.ratio = 0.7, # 10*sqrt(n)/n for n > 144
                                             lambda.min.ratio = 0.01,
                                             lambdaseq = "data-specific",
                                             ncores=1),
                           sparsMethod = "none",  dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Mid Phase/Mid network.rds")


# SpiecEasi
net_single <- netConstruct(Mid_data, verbose = 3,
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spieceasi",
                           measurePar = list(method = "glasso", 
                                             nlambda = 50,
                                             lambda.min.ratio = 0.01,
                                             pulsar.params = list(rep.num = 50,
                                                                  thresh = 0.1, 
                                                                  subsample.ratio = 0.7)),
                           sparsMethod = "none", dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Mid Phase/Mid network.rds")



# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single, "Phase/Mid Phase/Mid analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)

# ?plot.microNetProps
p <- plot(props_single,
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
          title1 = "Mid Phase Network on genus level with SPIEC-EASI Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Mid Phase/Mid plot.rds")


#Late Phase***********************************************************************
Late_data <- readRDS("Phase/Late Phase.rds")
Late_data
top_Late <- prune_taxa(names(sort(taxa_sums(Late_data),TRUE)[1:50]), Late_data)
plot_heatmap(top_Late)
Late_data

# ?netConstruct

#SPARCC
net_single <- netConstruct(Late_data,
                           verbose = 3,
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "sparcc",
                           sparsMethod = "threshold", thresh = 0.4,
                           dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Late Phase/Late network.rds")


props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Phase/Late Phase/Late analysis.rds")


p <- plot(props_single,
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
          title1 = "Late Phase Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Late Phase/Late plot.rds")


#Still has to do !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#SpRING
net_single <- netConstruct(Late_data, verbose = 3,
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spring",
                           measurePar = list(nlambda = 20, 
                                             rep.num = 20,
                                             thresh = 0.1, 
                                             subsample.ratio = 0.7, # 10*sqrt(n)/n for n > 144
                                             lambda.min.ratio = 0.01,
                                             lambdaseq = "data-specific",
                                             ncores=1),
                           sparsMethod = "none",  dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Late Phase/Late network.rds")


#Still has to do !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# SpiecEasi
net_single <- netConstruct(Late_data, verbose = 3,
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spieceasi",
                           measurePar = list(method = "glasso", 
                                             nlambda = 50,
                                             lambda.min.ratio = 0.01,
                                             pulsar.params = list(rep.num = 50,
                                                                  thresh = 0.1, 
                                                                  subsample.ratio = 0.7)),
                           sparsMethod = "none", dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Late Phase/Late network.rds")


# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single, "Phase/Late Phase/Late analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)

# ?plot.microNetProps
p <- plot(props_single,
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
          title1 = "Late Phase Network on genus level with SPIEC-EASI Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Late Phase/Late plot.rds")



#Dry Phase*************************************************************************
Dry_data <- readRDS("Phase/Dry Phase.rds")
Dry_data
top_Dry <- prune_taxa(names(sort(taxa_sums(Dry_data),TRUE)[1:50]), Dry_data)
plot_heatmap(top_Dry)

# ?netConstruct
net_single <- netConstruct(Dry_data,
                           verbose = 3,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "sparcc",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Dry Phase/Dry network.rds")

#SpRING
net_single <- netConstruct(Dry_data, verbose = 3,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spring",
                           measurePar = list(nlambda = 20, rep.num = 10,
                                             ncores=1),
                           sparsMethod = "none",  dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Dry Phase/Late network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "degree",
                           weightDeg = FALSE,
                           normDeg = TRUE,
                           lnormFit = TRUE,
                           verbose = 3)

saveRDS(props_single, "Phase/Dry Phase/Dry analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)

# ?plot.microNetProps
p <- plot(props_single,
          shortenLabels = "none",
          # labelLength = 16,
          # charToRm = "g__",
          labelScale = FALSE,
          rmSingles = "all",
          nodeSize = "degree",
          nodeColor = "cluster",
          hubBorderCol = "blue",
          cexNodes = 1,
          cexLabels = 1,
          edgeWidth = 1,
          highlightHubs = TRUE,
          cexHubs = 1.5,
          # cexHubLabels = 2,
          title1 = "Dry Phase Network on genus level with SparCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/Dry Phase/Dry plot.rds")

#All Phase************************************************************************
# select data of Different Lactation Phases.
filt_complete_agg
Phase_physeq <- subset_samples(filt_complete_agg, LactationPhase=="Early" | LactationPhase=="Mid" | LactationPhase =="Late" | LactationPhase=="Dry")
Phase_physeq

saveRDS(Phase_physeq, "Phase/All Phase.rds")

Phase_data <- readRDS("Phase/All Phase.rds")

top_Phase <- prune_taxa(names(sort(taxa_sums(Phase_data),TRUE)[1:50]), Phase_data)
plot_heatmap(top_Phase)

# ?netConstruct

#t-test
net_single <- netConstruct(Phase_data,
                           measure = "sparcc",
                           measurePar = list(iter=20),
                           sparsMethod = "t-test",
                           alpha = 0.05,
                           adjust = "adaptBH",
                           trueNullMethod = "convest",
                           verbose = 3,
                           seed = 123456)

saveRDS(net_single, "Phase/All Phase/All network.rds")

#Thereshold
net_single <- netConstruct(Phase_data,
                           verbose = 3,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "sparcc",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/All Phase/All network.rds")

#SpRING
net_single <- netConstruct(Phase_data, verbose = 3,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spring",
                           measurePar = list(nlambda = 20, rep.num = 10,
                                             ncores=1),
                           sparsMethod = "none",  dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Phase/Dry Phase/Late network.rds")


# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "degree",
                           weightDeg = FALSE,
                           normDeg = TRUE,
                           lnormFit = TRUE,
                           verbose = 3)

saveRDS(props_single, "Phase/All Phase/All analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)
net.summary

# ?plot.microNetProps
p <- plot(props_single,
          shortenLabels = "none",
          # labelLength = 16,
          # charToRm = "g__",
          labelScale = FALSE,
          rmSingles = "all",
          nodeSize = "degree",
          nodeColor = "cluster",
          hubBorderCol = "blue",
          cexNodes = 1,
          cexLabels = 0.5,
          edgeWidth = 1,
          highlightHubs = TRUE,
          cexHubs = 1.5,
          # cexHubLabels = 2,
          title1 = "All Phase Network on genus level with SparCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Phase/All Phase/All plot.rds")

from <- p$labels$labels1[p$q1$Edgelist$from]
to <- p$labels$labels1[p$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p$q1$Edgelist$from)))
edges$from <- p$q1$Edgelist$from
edges$to <- p$q1$Edgelist$to
edges$weight <- p$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/All Phase/All edge data.csv", row.names = FALSE)

hubs <- props_single$hubs$hubs1
write(hubs, "Phase/All Phase/All Phase Hubs.txt")

node.lables <- p$labels$labels1
clust <- props_single$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Phase/All Phase/All node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/All Phase/All hub data.csv", row.names = FALSE)



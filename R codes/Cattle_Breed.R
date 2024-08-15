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
physeq

#Reading the meta data
metx=read.delim("datafiles/metadata_table.txt",sep=",",header=T,row.names=sample_names(physeq))
head(metx)

#Meta data as sample data
metx=metx[,-1]
metx=sample_data(metx)
metx

#Merging the sample
complete_physeq=merge_phyloseq(physeq,metx)
complete_physeq

#**************************************************************************************************************************

#Visualization
sample_data(complete_physeq)
otu_table(complete_physeq)
tax_table(complete_physeq)

#aggregating taxa in the genus level
complete_agg <- tax_glom(complete_physeq, 'Genus')
complete_agg
tax_table(complete_agg)
complete_agg <- aggregate_taxa(complete_agg, 'Genus')
complete_agg


# keep only taxa that were observed at least twice
minTotRelAbun <- 5e-5
sums <- taxa_sums(complete_agg)
keepTaxa <- taxa_names(complete_agg)[which((sums / sum(sums)) > minTotRelAbun)]
filt_complete_agg <- prune_taxa(keepTaxa, complete_agg)

# select data of Different Cattle Breeds.
FJC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Friesian-Jersey_cross")
FC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Friesian_cross")
SC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Sahiwal_cross")
JSC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Jersey-Sahiwal_cross")
BC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Local_crossbreds(Batu_cross)")
FSC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Friesian-Sahiwal_cross")
AFS_physeq <- subset_samples(filt_complete_agg, CattleBreed =="AFS")
JC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Jersey_cross")
S_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Sahiwal")
A_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Ayrshire")
J_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Jersey")

J_physeq

# Saving the phyloseq objects
saveRDS(FJC_physeq, "Breed/FJC Breed.rds")
saveRDS(FC_physeq, "Breed/FC Breed.rds")
saveRDS(SC_physeq, "Breed/SC Breed.rds")
saveRDS(JSC_physeq, "Breed/JSC Breed.rds")
saveRDS(BC_physeq, "Breed/BC Breed.rds")
saveRDS(FSC_physeq, "Breed/FSC Breed.rds")
saveRDS(AFS_physeq, "Breed/AFS Breed.rds")
saveRDS(JC_physeq, "Breed/JC Breed.rds")
saveRDS(S_physeq, "Breed/S Breed.rds")
saveRDS(A_physeq, "Breed/A Breed.rds")
saveRDS(J_physeq, "Breed/J Breed.rds")

#FJC Breed************************************************************************
FJC_data <- readRDS("Breed/FJC Breed.rds")

top_FJC <- prune_taxa(names(sort(taxa_sums(FJC_data),TRUE)[1:50]), FJC_data)
plot_heatmap(top_FJC)

# ?netConstruct

#SparCC
net_single <- netConstruct(FJC_data,
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

saveRDS(net_single, "Breed/FJC Breed/FJC network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/FJC Breed/FJC analysis.rds")


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
          title1 = "FJC Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/FJC Breed/FJC plot.rds")






#SpRING
net_single <- netConstruct(FJC_data, verbose = 3,
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

saveRDS(net_single, "Breed/FJC Breed/FJC network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/FJC Breed/FJC analysis.rds")

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
          title1 = "FJC Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/FJC Breed/FJC plot.rds")


# SpiecEasi
net_single <- netConstruct(FJC_data, verbose = 3,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spieceasi",
                           measurePar = list(method = "glasso", 
                                             nlambda = 50,
                                             lambda.min.ratio = 0.01,
                                             pulsar.params = list(rep.num = 50,
                                                                  thresh = 0.1, 
                                                                  subsample.ratio = 0.8)),
                           sparsMethod = "none", dissFunc = "signed", 
                           seed = 123456)

saveRDS(net_single, "Breed/FJC Breed/FJC network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single, "Breed/FJC Breed/FJC analysis.rds")

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
          title1 = "FJC Breed Network on genus level with SPIEC-EASI Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/FJC Breed/FJC plot.rds")




#FC Breed************************************************************************
#4 Samples
#Gave Threshold Output
#No SPRING output

FC_data <- readRDS("Breed/FC Breed.rds")

top_FC <- prune_taxa(names(sort(taxa_sums(FC_data),TRUE)[1:50]), FC_data)
plot_heatmap(top_FC)

# ?netConstruct

#SparCC
net_single <- netConstruct(FC_data,
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

saveRDS(net_single, "Breed/FC Breed/FC network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/FC Breed/FC analysis.rds")


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
          title1 = "FC Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/FC Breed/FC plot.rds")



#SpRING
net_single <- netConstruct(FC_data, verbose = 3,
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

saveRDS(net_single, "Breed/FC Breed/FC network.rds")
#No Spring Output

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/FC Breed/FC analysis.rds")

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
          title1 = "FC Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/FC Breed/FC plot.rds")

#SC Breed************************************************************************
#2 Samples
#Abnormal Threshold Output
#Error in SPRING Output


SC_data <- readRDS("Breed/SC Breed.rds")

top_SC <- prune_taxa(names(sort(taxa_sums(SC_data),TRUE)[1:50]), SC_data)
plot_heatmap(top_SC)
SC_data

# ?netConstruct
net_single <- netConstruct(SC_data,
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

saveRDS(net_single, "Breed/SC Breed/SC network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/SC Breed/SC analysis.rds")


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
          title1 = "SC Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/SC Breed/SC plot.rds")




#SpRING
net_single <- netConstruct(SC_data, verbose = 3,
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

saveRDS(net_single, "Breed/SC Breed/SC network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/SC Breed/SC analysis.rds")

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
          title1 = "SC Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/SC Breed/SC plot.rds")

#JSC Breed************************************************************************
# 15 +4 Samples
# Threshold Output given
# SPRING Output given


JSC_data <- readRDS("Breed/JSC Breed.rds")

top_JSC <- prune_taxa(names(sort(taxa_sums(JSC_data),TRUE)[1:50]), JSC_data)
plot_heatmap(top_JSC)

# ?netConstruct

#SparCC
net_single <- netConstruct(JSC_data,
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

saveRDS(net_single, "Breed/JSC Breed/JSC network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/JSC Breed/JSC analysis.rds")


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
          title1 = "JSC Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/JSC Breed/JSC plot.rds")




#SpRING
net_single <- netConstruct(JSC_data, verbose = 3,
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

saveRDS(net_single, "Breed/JSC Breed/JSC network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/JSC Breed/JSC analysis.rds")

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
          title1 = "JSC Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/JSC Breed/JSC plot.rds")


#BC Breed************************************************************************
# 4 Samples
# Threshold Output given
# SPRING Output given


BC_data <- readRDS("Breed/BC Breed.rds")

top_BC <- prune_taxa(names(sort(taxa_sums(BC_data),TRUE)[1:50]), BC_data)
plot_heatmap(top_BC)

# ?netConstruct

#Sparcc
net_single <- netConstruct(BC_data,
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

saveRDS(net_single, "Breed/BC Breed/BC network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/BC Breed/BC analysis.rds")


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
          title1 = "BC Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/BC Breed/BC plot.rds")




#SpRING
net_single <- netConstruct(BC_data, verbose = 3,
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

saveRDS(net_single, "Breed/BC Breed/BC network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/BC Breed/BC analysis.rds")

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
          title1 = "BC Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/BC Breed/BC plot.rds")



#FSC Breed************************************************************************
# 2 Samples
# Threshold Output Abnormal
# SPRING Output 


FSC_data <- readRDS("Breed/FSC Breed.rds")

top_FSC <- prune_taxa(names(sort(taxa_sums(FSC_data),TRUE)[1:50]), FSC_data)
plot_heatmap(top_FSC)

# ?netConstruct

#SparCC
net_single <- netConstruct(FSC_data,
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

saveRDS(net_single, "Breed/FSC Breed/FSC network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/FSC Breed/FSC analysis.rds")


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
          title1 = "FSC Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/FSC Breed/FSC plot.rds")




#SpRING
net_single <- netConstruct(FSC_data, verbose = 3,
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

saveRDS(net_single, "Breed/FSC Breed/FSC network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/FSC Breed/FSC analysis.rds")

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
          title1 = "FSC Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/FSC Breed/FSC plot.rds")


#AFS Breed************************************************************************
# 2 Samples
# Threshold Output abnormal
# No SPRING Output 

AFS_data <- readRDS("Breed/AFS Breed.rds")

top_AFS <- prune_taxa(names(sort(taxa_sums(AFS_data),TRUE)[1:50]), AFS_data)
plot_heatmap(top_AFS)

# ?netConstruct
net_single <- netConstruct(AFS_data,
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

saveRDS(net_single, "Breed/AFS Breed/AFS network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/AFS Breed/AFS analysis.rds")


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
          title1 = "AFS Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/AFS Breed/AFS plot.rds")




#SpRING
net_single <- netConstruct(AFS_data, verbose = 3,
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

saveRDS(net_single, "Breed/AFS Breed/AFS network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/AFS Breed/AFS analysis.rds")

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
          title1 = "AFS Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/AFS Breed/AFS plot.rds")


#JC Breed************************************************************************
# 6 Samples
# Threshold Output given
# SPRING Output: Network has no edges 


JC_data <- readRDS("Breed/JC Breed.rds")

top_JC <- prune_taxa(names(sort(taxa_sums(JC_data),TRUE)[1:50]), JC_data)
plot_heatmap(top_JC)

# ?netConstruct

#SparCC
net_single <- netConstruct(JC_data,
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

saveRDS(net_single, "Breed/JC Breed/JC network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/JC Breed/JC analysis.rds")


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
          title1 = "JC Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/JC Breed/JC plot.rds")



#SpRING
net_single <- netConstruct(JC_data, verbose = 3,
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


saveRDS(net_single, "Breed/JC Breed/JC network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/JC Breed/JC analysis.rds")

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
          title1 = "JC Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/JC Breed/JC plot.rds")


# S Breed************************************************************************
# 2 Samples
# Threshold Output abnormal
# Error in SPRING Output 


S_data <- readRDS("Breed/S Breed.rds")

top_S <- prune_taxa(names(sort(taxa_sums(S_data),TRUE)[1:50]), S_data)
plot_heatmap(top_S)

# ?netConstruct
net_single <- netConstruct(S_data,
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

saveRDS(net_single, "Breed/S Breed/S network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/S Breed/S analysis.rds")


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
          title1 = "S Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/S Breed/S plot.rds")



#SpRING
net_single <- netConstruct(S_data, verbose = 3,
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

saveRDS(net_single, "Breed/S Breed/S network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/S Breed/S analysis.rds")

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
          title1 = "S Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/S Breed/S plot.rds")

# A Breed************************************************************************
# 0 Samples after filtering
# Threshold Output is null 
# No SPRING Output 


A_data <- readRDS("Breed/A Breed.rds")

top_A <- prune_taxa(names(sort(taxa_sums(A_data),TRUE)[1:50]), A_data)
plot_heatmap(top_A)

# ?netConstruct
net_single <- netConstruct(A_data,
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

saveRDS(net_single, "Breed/A Breed/A network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/A Breed/A analysis.rds")


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
          title1 = "A Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/A Breed/A plot.rds")



#SpRING
net_single <- netConstruct(A_data, verbose = 3,
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

saveRDS(net_single, "Breed/A Breed/A network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/A Breed/A analysis.rds")

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
          title1 = "A Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/A Breed/A plot.rds")


# J Breed************************************************************************
# 5 Samples 
# Threshold Output is given 
# No SPRING Output 


J_data <- readRDS("Breed/J Breed.rds")

top_J <- prune_taxa(names(sort(taxa_sums(J_data),TRUE)[1:50]), J_data)
plot_heatmap(top_J)

# ?netConstruct
net_single <- netConstruct(J_data,
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

saveRDS(net_single, "Breed/J Breed/J network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/J Breed/J analysis.rds")


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
          title1 = "J Breed Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/J Breed/J plot.rds")




#SpRING
net_single <- netConstruct(J_data, verbose = 3,
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

saveRDS(net_single, "Breed/J Breed/J network.rds")

# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/J Breed/J analysis.rds")

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
          title1 = "J Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


saveRDS(p, "Breed/J Breed/J plot.rds")


#All Breeds************************************************************************
# select data of Different Breeds.
# select data of Different Cattle Breeds.

filt_complete_agg
Breed_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Friesian-Jersey_cross" | CattleBreed =="Friesian_cross" |
                                 CattleBreed =="Sahiwal_cross" | CattleBreed =="Jersey-Sahiwal_cross" | 
                                 CattleBreed =="Local_crossbreds(Batu_cross)" | CattleBreed =="Friesian-Sahiwal_cross" | 
                                 CattleBreed =="AFS" | CattleBreed =="Jersey_cross" | CattleBreed =="Sahiwal" | 
                                 CattleBreed =="Ayrshire" | CattleBreed =="Jersey")
Breed_physeq

saveRDS(Breed_physeq, "Breed/All Breed.rds")

Breed_data <- readRDS("Breed/All Breed.rds")

top_Breed <- prune_taxa(names(sort(taxa_sums(Breed_data),TRUE)[1:50]), Breed_data)
plot_heatmap(top_Breed)

# ?netConstruct

#t-test
net_single <- netConstruct(Breed_data,
                           measure = "sparcc",
                           measurePar = list(iter=20),
                           sparsMethod = "t-test",
                           alpha = 0.05,
                           adjust = "adaptBH",
                           trueNullMethod = "convest",
                           verbose = 3,
                           seed = 123456)

saveRDS(net_single, "Breed/All Breed/All network.rds")

#Thereshold
net_single <- netConstruct(Breed_data,
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

saveRDS(net_single, "Breed/All Breed/All network.rds")

props_single <- netAnalyze(net_single, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single, "Breed/All Breed/All analysis.rds")


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
          title1 = "All Network on genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/All Breed/All plot.rds")




#SpRING
net_single <- netConstruct(Breed_data, verbose = 3,
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

saveRDS(net_single, "Breed/All Breed/All network.rds")


# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)


saveRDS(props_single, "Breed/All Breed/All analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)
net.summary

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
          title1 = "All Breed Network on genus level with SPRING Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Breed/All Breed/All plot.rds")









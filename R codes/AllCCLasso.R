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

# select data of different Cattle Breeds
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


# Saving the phyloseq objects for lactation phases
saveRDS(complete_physeq, "Phase/Complete data.rds")
saveRDS(complete_agg, "Phase/Aggregated data.rds")
saveRDS(filt_complete_agg, "Phase/Filtered Aggregated data.rds")
saveRDS(Early_physeq, "Phase/Early Phase.rds")
saveRDS(Late_physeq, "Phase/Late Phase.rds")
saveRDS(Mid_physeq, "Phase/Mid Phase.rds")
saveRDS(Dry_physeq, "Phase/Dry Phase.rds")

# Saving the phyloseq objects for carrle breeds
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


#Early Lactation************************************************************************
Early_data <- readRDS("Phase/Early Phase.rds")

top_Early <- prune_taxa(names(sort(taxa_sums(Early_data),TRUE)[1:50]), Early_data)
#plot_heatmap(top_Early)
Early_data

# CCLasso
net_single <- netConstruct(Early_data, verbose = 3,
                            # filtTax = "highestFreq",
                            # filtTaxPar = list(highestFreq = 100),
                            filtTax = "none",
                            filtTaxPar = "totalReads",
                            filtSamp = "totalReads",
                            filtSampPar = list(totalReads = 1000),
                            zeroMethod = "pseudo", normMethod = "none",
                            measure = "cclasso",
                            sparsMethod = "threshold", thresh = 0.2,
                            dissFunc = "signed", 
                            seed = 123456)

props_single <- netAnalyze(net_single, clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", hubQuant = 0.95)

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


#Mid Lactation**************************************************************************
Mid_data <- readRDS("Phase/Mid Phase.rds")
Mid_data
top_Mid <- prune_taxa(names(sort(taxa_sums(Mid_data),TRUE)[1:50]), Mid_data)
top_Mid
#plot_heatmap(top_Mid)

# CCLasso
net_single <- netConstruct(Mid_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "pseudo", normMethod = "none",
                           measure = "cclasso",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 123456)

props_single <- netAnalyze(net_single, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95)

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
          nodeSize = "eigenvector",
          nodeColor = "cluster",
          hubBorderCol = "blue",
          cexNodes = 1,
          cexLabels = 0.5,
          edgeWidth = 1,
          highlightHubs = TRUE,
          cexHubs = 1.5,
          # cexHubLabels = 2,
          title1 = "Mid Phase Network on genus level with CCLasso Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

#Late Lactation**************************************************************************
Late_data <- readRDS("Phase/Late Phase.rds")
Late_data
top_Late <- prune_taxa(names(sort(taxa_sums(Late_data),TRUE)[1:50]), Late_data)
top_Late
#plot_heatmap(top_Mid)

# CCLasso
net_single <- netConstruct(Late_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "pseudo", normMethod = "none",
                           measure = "cclasso",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 123456)

props_single <- netAnalyze(net_single, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95)

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
          nodeSize = "eigenvector",
          nodeColor = "cluster",
          hubBorderCol = "blue",
          cexNodes = 1,
          cexLabels = 0.5,
          edgeWidth = 1,
          highlightHubs = TRUE,
          cexHubs = 1.5,
          # cexHubLabels = 2,
          title1 = "Late Phase Network on genus level with CCLasso Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


#FJC Breed************************************************************************
FJC_data <- readRDS("Breed/FJC Breed.rds")

top_FJC <- prune_taxa(names(sort(taxa_sums(FJC_data),TRUE)[1:50]), FJC_data)
#plot_heatmap(top_FJC)

# CCLasso
net_single <- netConstruct(FJC_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "pseudo", normMethod = "none",
                           measure = "cclasso",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 123456)

props_single <- netAnalyze(net_single, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95)

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
          nodeSize = "eigenvector",
          nodeColor = "cluster",
          hubBorderCol = "blue",
          cexNodes = 1,
          cexLabels = 0.5,
          edgeWidth = 1,
          highlightHubs = TRUE,
          cexHubs = 1.5,
          # cexHubLabels = 2,
          title1 = "FJC Breed Network on genus level with CCLasso Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

#FC Breed************************************************************************
FC_data <- readRDS("Breed/FC Breed.rds")

top_FC <- prune_taxa(names(sort(taxa_sums(FC_data),TRUE)[1:50]), FC_data)
#plot_heatmap(top_FJC)

# CCLasso
net_single <- netConstruct(FC_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
                           filtTax = "none",
                           filtTaxPar = "totalReads",
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "pseudo", normMethod = "none",
                           measure = "cclasso",
                           sparsMethod = "threshold", thresh = 0.2,
                           dissFunc = "signed", 
                           seed = 123456)

props_single <- netAnalyze(net_single, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95)

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
          nodeSize = "eigenvector",
          nodeColor = "cluster",
          hubBorderCol = "blue",
          cexNodes = 1,
          cexLabels = 0.5,
          edgeWidth = 1,
          highlightHubs = TRUE,
          cexHubs = 1.5,
          # cexHubLabels = 2,
          title1 = "FC Breed Network on genus level with CCLasso Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)


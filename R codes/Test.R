library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
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

#install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)

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

Early_data <- readRDS("Phase/Early Phase.rds")

top_Early <- prune_taxa(names(sort(taxa_sums(Early_data),TRUE)[1:50]), Early_data)
plot_heatmap(top_Early)
Early_data

?spiec.easi


spiec.out=spiec.easi(Early_data, method="mb",
                     pulsar.params = list(rep.num = 50,
                                          thresh = 0.1, 
                                          subsample.ratio = 0.8),
                     icov.select.params=list(rep.num=20))

spiec.graph=adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(Early_data)))
spiec.graph

plot_network(spiec.graph, Early_data, type='taxa',label="Genus")



library("SPRING")
?SPRING

spring.out=SPRING(Early_data,
                  quantitative = FALSE,
                  method = "mb",
                  ambda.min.ratio = 0.01,
                  nlambda = 20,
                  lambdaseq = exp(seq(log(0.6), log(0.6 * lambda.min.ratio), length.out = nlambda)),
                  seed = 123456,ncores = 1,
                  thresh = 0.1,
                  subsample.ratio = 0.8,
                  rep.num = 20,
                  Rtol = 1e-06,
                  verbose = TRUE,
                  verboseR = FALSE,
                  Rmethod = "original"
  
)

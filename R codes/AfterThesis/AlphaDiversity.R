
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

#**************************************************************************************************************************

# #aggregating taxa in the Genus level
# complete_agg <- tax_glom(complete_physeq, 'Genus')
# #complete_agg
# tax_table(complete_agg)
# complete_agg <- aggregate_taxa(complete_agg, 'Genus')
# complete_agg


# Without aggregate ....................................................
minTotRelAbun <- 5e-5
sums <- taxa_sums(complete_physeq)
#sums
keepTaxa <- taxa_names(complete_physeq)[which((sums / sum(sums)) > minTotRelAbun)]
#keepTaxa
filt_complete_phy <- prune_taxa(keepTaxa, complete_physeq)
filt_complete_phy

plot_richness(filt_complete_phy, x="LactationPhase", measures=c("Observed","Chao1", "Shannon","Simpson"))
plot_richness(filt_complete_phy, x="LactationPhase", measures=c("Observed","Chao1", "Shannon","Simpson"))

# # With aggregate ....................................................
# minTotRelAbun <- 5e-5
# sums <- taxa_sums(complete_agg)
# #sums
# keepTaxa <- taxa_names(complete_agg)[which((sums / sum(sums)) > minTotRelAbun)]
# #keepTaxa
# filt_complete_agg <- prune_taxa(keepTaxa, complete_agg)
# filt_complete_agg
# 
# 
# plot_richness(filt_complete_agg, x="LactationPhase", measures=c("Observed","Chao1", "Shannon","Simpson"))
# 


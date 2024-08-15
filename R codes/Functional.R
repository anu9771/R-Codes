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


constructNetwork=function(method,par,subpar){
  
  Data <- readRDS(paste(par,"/",subpar," ",par,".rds",sep=""))
  Data
  
  if(method=="sparcc"){
    net_single <- netConstruct(Data,
                               verbose = 3,
                               # filtTax = "highestFreq",
                               # filtTaxPar = list(highestFreq = 100),
                               filtTax = "none",
                               filtTaxPar = "totalReads",
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 1000),
                               zeroMethod = "none", normMethod = "none",
                               measure = "sparcc",
                               sparsMethod = "threshold", thresh = 0.4,
                               dissFunc = "signed", 
                               seed = 123456)
  }
  
  else if(method=="spring"){
    net_single <- netConstruct(Data, verbose = 3,
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
  }
  
  
  saveRDS(net_single,paste(par,"/",subpar," ",par,"/",subpar,"_",method," network.rds",sep=""))
  
  props_single <- netAnalyze(net_single, 
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector", 
                             hubQuant = 0.95)
  
  saveRDS(props_single,paste(par,"/",subpar," ",par,"/",subpar,"_",method," analysis.rds",sep=""))
  
  #?summary.microNetProps
  net.summary <- summary(props_single)
  net.summary
  
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
            title1 = paste(subpar," ",par," Network on genus level with SparCC Method.",sep=""), 
            showTitle = TRUE,
            cexTitle = 1.5)
  
  saveRDS(p, paste(par,"/",subpar," ",par,"/",subpar,"_",method," plot.rds",sep=""))
  
}

####################################################################################################

constructNetwork("sparcc","breed","J")

















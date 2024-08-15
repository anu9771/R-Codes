# devtools::install_github("stefpeschel/NetCoMi", 
#                          dependencies = c("Depends", "Imports", "LinkingTo"),
#                          repos = c("https://cloud.r-project.org/",
#                                    BiocManager::repositories()))



# devtools::install_github("GraceYoon/SPRING",force = T)
# devtools::install_github("GraceYoon/SPRING")


# 1. Australian Friesian Sahiwal (AFS)
# 2. Ayrshire (AYS)
# 3. Friesian (FR)
# 4. Friesian-cross (FRC)
# 5. Friesian-Jersey_Cross (FRJRC)
# 6. Friesian-Sahiwal-Cross (FRSHC)
# 7. Jersey (JR)
# 8. Jersey-cross (JRC)
# 9. Jersey-Sahiwal-cross (JRSHC)
# 10. Local-Crossbreds: Batu-Cross (BC)
# 11. Sahiwal (SH)
# 12. Sahiwal-cross (SHC)
# 
# Categorize as below
# 1. Jersey (JR) and Jersey Cross (JRC)
# - Jersey-Sahiwal-Cross (JRSHC)
# 2. Friesian (FR) and Friesian-cross (FRC)
# - Friesian-Sahiwal-Cross (FRSHC)
# - Friesian-Jersey_Cross (FRJRC
#                          - Australian Friesian Sahiwal (AFS)
#                          3. Ayrshire (AYS)
#                          4. Sahiwal (SH) and Sahiwal-Cross (SHC)
#                          5. Local-crossbreds: Batu-Cross (BC)



# Loading required packages
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(BiocManager)
library(microbiome)
library(NetCoMi)
library(WGCNA)



# Setting the working directory
setwd("C:/Users/audar/Desktop/Anushka/Research/Methodology")


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


# keep only taxa that are ....................................................
# minTotRelAbun <- 5e-5
# sums <- taxa_sums(complete_agg)
# keepTaxa <- taxa_names(complete_agg)[which((sums / sum(sums)) > minTotRelAbun)]
# filt_complete_agg <- prune_taxa(keepTaxa, complete_agg)
# filt_complete_agg

# 
# # select data of Different Lactation Phases.
# Early_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Early")
# Late_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Late")
# Mid_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Mid")
# Dry_physeq <- subset_samples(filt_complete_agg, LactationPhase =="Dry")
# 
# # select data of different Cattle Breeds
# FJC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Friesian-Jersey_cross")
# FC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Friesian_cross")
# SC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Sahiwal_cross")
# JSC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Jersey-Sahiwal_cross")
# BC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Local_crossbreds(Batu_cross)")
# FSC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Friesian-Sahiwal_cross")
# AFS_physeq <- subset_samples(filt_complete_agg, CattleBreed =="AFS")
# JC_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Jersey_cross")
# S_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Sahiwal")
# A_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Ayrshire")
# J_physeq <- subset_samples(filt_complete_agg, CattleBreed =="Jersey")


#Select Updated Cattle Breeds
JRJRC_physeq=subset_samples(filt_complete_agg, CattleBreed =="Jersey" | CattleBreed=="Jersey_cross" | CattleBreed=="Jersey-Sahiwal_cross")
FRFRC_physeq=subset_samples(filt_complete_agg, CattleBreed =="Friesian_cross" | CattleBreed=="Friesian-Sahiwal_cross" | CattleBreed=="Friesian-Jersey_cross" | CattleBreed=="AFS" | CattleBreed=="Friesian")
Asian_physeq=subset_samples(filt_complete_agg, CattleBreed =="Sahiwal" | CattleBreed=="Sahiwal_cross" | CattleBreed =="Local_crossbreds(Batu_cross)")


# # Saving the phyloseq objects for lactation phases
# saveRDS(complete_physeq, "Phase/Complete data.rds")
# saveRDS(complete_agg, "Phase/Aggregated data.rds")
# saveRDS(filt_complete_agg, "Phase/Filtered Aggregated data.rds")
# saveRDS(Early_physeq, "Phase/Early Phase.rds")
# saveRDS(Late_physeq, "Phase/Late Phase.rds")
# saveRDS(Mid_physeq, "Phase/Mid Phase.rds")
# saveRDS(Dry_physeq, "Phase/Dry Phase.rds")
# 
# # Saving the phyloseq objects for cattle breeds
# saveRDS(FJC_physeq, "Breed/FJC Breed.rds")
# saveRDS(FC_physeq, "Breed/FC Breed.rds")
# saveRDS(SC_physeq, "Breed/SC Breed.rds")
# saveRDS(JSC_physeq, "Breed/JSC Breed.rds")
# saveRDS(BC_physeq, "Breed/BC Breed.rds")
# saveRDS(FSC_physeq, "Breed/FSC Breed.rds")
# saveRDS(AFS_physeq, "Breed/AFS Breed.rds")
# saveRDS(JC_physeq, "Breed/JC Breed.rds")
# saveRDS(S_physeq, "Breed/S Breed.rds")
# saveRDS(A_physeq, "Breed/A Breed.rds")
# saveRDS(J_physeq, "Breed/J Breed.rds")

# Saving the phyloseq objects for updated cattle  breeds
saveRDS(JRJRC_physeq, "Breed/JRJRC.rds")
saveRDS(FRFRC_physeq, "Breed/FRFRC.rds")
saveRDS(Asian_physeq, "Breed/Asian.rds")
# saveRDS(LCBC_physeq, "Breed/LCBC.rds")


comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)



#Sparcc#######################################################################################################

#JRJRC********************************************************************************************************
#30 Samples
JRJRC_data <- readRDS("Breed/JRJRC.rds")

top_JRJRC <- prune_taxa(names(sort(taxa_sums(JRJRC_data),TRUE)[1:50]), JRJRC_data)
#plot_heatmap(top_JRJRC)
JRJRC_data

#SparCC
net_single_JRJRC <- netConstruct(JRJRC_data,
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

saveRDS(net_single_JRJRC, "Breed/JRJRC Breed/JRJRC_Sparcc network.rds")

props_single_JRJRC <- netAnalyze(net_single_JRJRC, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)


saveRDS(props_single_JRJRC, "Breed/JRJRC Breed/JRJRC_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_JRJRC)
net.summary

p_JRJRC <- plot(props_single_JRJRC,
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
                title1 = "JRJRC Breed Network on Genus level with SparCC Method.", 
                showTitle = TRUE,
                cexTitle = 1.5)

saveRDS(p_JRJRC, "Breed/JRJRC Breed/JRJRC_Sparcc plot.rds")

from <- p_JRJRC$labels$labels1[p_JRJRC$q1$Edgelist$from]
to <- p_JRJRC$labels$labels1[p_JRJRC$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_JRJRC$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_JRJRC$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_JRJRC$q1$Edgelist$from)))
edges$from <- p_JRJRC$q1$Edgelist$from
edges$to <- p_JRJRC$q1$Edgelist$to
edges$weight <- p_JRJRC$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/JRJRC Breed/JRJRC edge data.csv", row.names = FALSE)

hubs <- props_single_JRJRC$hubs$hubs1
write(hubs, "Breed/JRJRC Breed/JRJRC Breed Hubs.txt")

node.lables <- p_JRJRC$labels$labels1
clust <- props_single_JRJRC$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_JRJRC$centralities$degree1[nodes$lable]
evs=props_single_JRJRC$centralities$eigenv1[nodes$lable]
betweennesses=props_single_JRJRC$centralities$between1[nodes$lable]
closenesses=props_single_JRJRC$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/JRJRC Breed/JRJRC node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/JRJRC Breed/JRJRC hub data.csv", row.names = FALSE)

#For Associations---------------------------------------------------------------------------------------------


Genus_Sparcc_Breed_JRJRC_Top100=data.frame(Node1=net_single_JRJRC[["edgelist1"]][["v1"]],
                                         Node2=net_single_JRJRC[["edgelist1"]][["v2"]],
                                         Association=net_single_JRJRC[["edgelist1"]][["asso"]],
                                         Dissimilarity=net_single_JRJRC[["edgelist1"]][["diss"]],
                                         Adjecency=net_single_JRJRC[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Breed_JRJRC_Top100=Genus_Sparcc_Breed_JRJRC_Top100[order(Genus_Sparcc_Breed_JRJRC_Top100$Association,
                                                                  decreasing = T),]

write.csv(Genus_Sparcc_Breed_JRJRC_Top100,"Associations/Genus_Sparcc_Breed_JRJRC_Top100.csv",row.names = FALSE)

#Functions----------------------------------------------------------------------------------------------------


funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/JRJRC Breed/JRJRC node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/JRJRCCluster',paste("JRJRC_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/JRJRC_Cluster1.csv",header=T)
# tax[2]
# 
# 
# a=list()  
# 
# for(i in 1:ncol(tax)) {      
#   if(i==2){
#     a[[i]] <- tax[ , i]
#   }    
# }
# 
# names(a)=colnames(tax)  
# print(a)


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)
comp_tax_csv

# for(i in a$x){
#   print(subset(comp_tax_csv, comp_tax_csv$genus==i))
#   path=normalizePath(file.path('./New folder/Otu',paste("JRJRC_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }
# 
# 
# 
# 
# 
# otu=read.csv("New folder/Otu/JRJRC_otuid_Papillibacter.csv",header=T)
# otu[2]
# 
# a=list()  
# 
# for(i in 1:ncol(otu)) {      
#   if(i==2){
#     a[[i]] <- otu[ , i]
#   }    
# }
# 
# print(a)
# 
# func_df=data.frame(funcFile)
# func_df[1]
# 
# 
# for(i in a[[2]]){
#   print(subset(func_df, func_df[1]==i))
#   # path=normalizePath(file.path('./New folder',paste("JRJRC_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 
# 
# 
# 
# 
# 
# 
# 


getwd()
fileList=list.files(path="New folder\\JRJRCCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/JRJRCCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/JRJRCOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\JRJRCOtu",pattern = ".txt")
fileList


func_df=data.frame(funcFile)
func_df[1]

# for(i in fileList){
#   cat(i,"*****************************************************************")
#   path=normalizePath(file.path('./New folder/Otu',paste(i, sep='')))
#   file=(read.table(path))
#   print(file[2])
# }




func_df=data.frame(funcFile)


fileList=list.files(path="New folder\\JRJRCOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/JRJRCOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/JRJRCOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:7){
  path1=normalizePath(file.path('./New folder/JRJRCOtu/Res',paste("JRJRC_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/JRJRCOtu',paste("JRJRC_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("JRJRC_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/JRJRC_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/JRJRC_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)


#FRFRC********************************************************************************************************
#26 Samples
FRFRC_data <- readRDS("Breed/FRFRC.rds")

top_FRFRC <- prune_taxa(names(sort(taxa_sums(FRFRC_data),TRUE)[1:50]), FRFRC_data)
#plot_heatmap(top_FRFRC)
FRFRC_data

#SparCC
net_single_FRFRC <- netConstruct(FRFRC_data,
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

saveRDS(net_single_FRFRC, "Breed/FRFRC Breed/FRFRC_Sparcc network.rds")

props_single_FRFRC <- netAnalyze(net_single_FRFRC, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_FRFRC, "Breed/FRFRC Breed/FRFRC_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_FRFRC)
net.summary

p_FRFRC <- plot(props_single_FRFRC,
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
                title1 = "FRFRC Breed Network on Genus level with SparCC Method.", 
                showTitle = TRUE,
                cexTitle = 1.5)

saveRDS(p_FRFRC, "Breed/FRFRC Breed/FRFRC_Sparcc plot.rds")

from <- p_FRFRC$labels$labels1[p_FRFRC$q1$Edgelist$from]
to <- p_FRFRC$labels$labels1[p_FRFRC$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_FRFRC$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_FRFRC$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_FRFRC$q1$Edgelist$from)))
edges$from <- p_FRFRC$q1$Edgelist$from
edges$to <- p_FRFRC$q1$Edgelist$to
edges$weight <- p_FRFRC$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/FRFRC Breed/FRFRC edge data.csv", row.names = FALSE)

hubs <- props_single_FRFRC$hubs$hubs1
write(hubs, "Breed/FRFRC Breed/FRFRC Breed Hubs.txt")

node.lables <- p_FRFRC$labels$labels1
clust <- props_single_FRFRC$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_FRFRC$centralities$degree1[nodes$lable]
evs=props_single_FRFRC$centralities$eigenv1[nodes$lable]
betweennesses=props_single_FRFRC$centralities$between1[nodes$lable]
closenesses=props_single_FRFRC$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/FRFRC Breed/FRFRC node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/FRFRC Breed/FRFRC hub data.csv", row.names = FALSE)


#For Associations---------------------------------------------------------------------------------------------


Genus_Sparcc_Breed_FRFRC_Top100=data.frame(Node1=net_single_FRFRC[["edgelist1"]][["v1"]],
                                           Node2=net_single_FRFRC[["edgelist1"]][["v2"]],
                                           Association=net_single_FRFRC[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_FRFRC[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_FRFRC[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Breed_FRFRC_Top100=Genus_Sparcc_Breed_FRFRC_Top100[order(Genus_Sparcc_Breed_FRFRC_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Sparcc_Breed_FRFRC_Top100,"Associations/Genus_Sparcc_Breed_FRFRC_Top100.csv",row.names = FALSE)


#Functions----------------------------------------------------------------------------------------------------


funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/FRFRC Breed/FRFRC node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/FRFRCCluster',paste("FRFRC_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/FRFRC_Cluster1.csv",header=T)
# tax[2]
# 
# 
# a=list()  
# 
# for(i in 1:ncol(tax)) {      
#   if(i==2){
#     a[[i]] <- tax[ , i]
#   }    
# }
# 
# names(a)=colnames(tax)  
# print(a)


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)
comp_tax_csv

# for(i in a$x){
#   print(subset(comp_tax_csv, comp_tax_csv$genus==i))
#   path=normalizePath(file.path('./New folder/Otu',paste("FRFRC_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }
# 
# 
# 
# 
# 
# otu=read.csv("New folder/Otu/FRFRC_otuid_Papillibacter.csv",header=T)
# otu[2]
# 
# a=list()  
# 
# for(i in 1:ncol(otu)) {      
#   if(i==2){
#     a[[i]] <- otu[ , i]
#   }    
# }
# 
# print(a)
# 
# func_df=data.frame(funcFile)
# func_df[1]
# 
# 
# for(i in a[[2]]){
#   print(subset(func_df, func_df[1]==i))
#   # path=normalizePath(file.path('./New folder',paste("FRFRC_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 
# 
# 
# 
# 
# 
# 
# 


getwd()
fileList=list.files(path="New folder\\FRFRCCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/FRFRCCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/FRFRCOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\FRFRCOtu",pattern = ".txt")
fileList


func_df=data.frame(funcFile)
func_df[1]

# for(i in fileList){
#   cat(i,"*****************************************************************")
#   path=normalizePath(file.path('./New folder/Otu',paste(i, sep='')))
#   file=(read.table(path))
#   print(file[2])
# }




func_df=data.frame(funcFile)


fileList=list.files(path="New folder\\FRFRCOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/FRFRCOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/FRFRCOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:6){
  path1=normalizePath(file.path('./New folder/FRFRCOtu/Res',paste("FRFRC_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/FRFRCOtu',paste("FRFRC_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("FRFRC_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/FRFRC_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/FRFRC_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)




#Asian********************************************************************************************************
#8 Samples
Asian_data <- readRDS("Breed/Asian.rds")

top_Asian <- prune_taxa(names(sort(taxa_sums(Asian_data),TRUE)[1:50]), Asian_data)
#plot_heatmap(top_Asian)
Asian_data

#SparCC
net_single_Asian <- netConstruct(Asian_data,
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

saveRDS(net_single_Asian, "Breed/Asian Breed/Asian_Sparcc network.rds")

props_single_Asian <- netAnalyze(net_single_Asian,
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector",
                                 hubQuant = 0.95)

saveRDS(props_single_Asian, "Breed/Asian Breed/Asian_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Asian)
net.summary

p_Asian <- plot(props_single_Asian,
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
                title1 = "Asian Breed Network on Genus level with SparCC Method.",
                showTitle = TRUE,
                cexTitle = 1.5)

saveRDS(p_Asian, "Breed/Asian Breed/Asian_Sparcc plot.rds")

from <- p_Asian$labels$labels1[p_Asian$q1$Edgelist$from]
to <- p_Asian$labels$labels1[p_Asian$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Asian$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1
  }
  else if (net_single_Asian$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1
  }
}

edges <- data.frame(row.names = c(1:length(p_Asian$q1$Edgelist$from)))
edges$from <- p_Asian$q1$Edgelist$from
edges$to <- p_Asian$q1$Edgelist$to
edges$weight <- p_Asian$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Asian Breed/Asian edge data.csv", row.names = FALSE)

hubs <- props_single_Asian$hubs$hubs1
write(hubs, "Breed/Asian Breed/Asian Breed Hubs.txt")

node.lables <- p_Asian$labels$labels1
clust <- props_single_Asian$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Asian$centralities$degree1[nodes$lable]
evs=props_single_Asian$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Asian$centralities$between1[nodes$lable]
closenesses=props_single_Asian$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/Asian Breed/Asian node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Asian Breed/Asian hub data.csv", row.names = FALSE)

#For Associations---------------------------------------------------------------------------------------------


Genus_Sparcc_Breed_Asian_Top100=data.frame(Node1=net_single_Asian[["edgelist1"]][["v1"]],
                                           Node2=net_single_Asian[["edgelist1"]][["v2"]],
                                           Association=net_single_Asian[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_Asian[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_Asian[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Breed_Asian_Top100=Genus_Sparcc_Breed_Asian_Top100[order(Genus_Sparcc_Breed_Asian_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Sparcc_Breed_Asian_Top100,"Associations/Genus_Sparcc_Breed_Asian_Top100.csv",row.names = FALSE)


#Functions----------------------------------------------------------------------------------------------------


funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/Asian Breed/Asian node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/AsianCluster',paste("Asian_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/Asian_Cluster1.csv",header=T)
# tax[2]
# 
# 
# a=list()  
# 
# for(i in 1:ncol(tax)) {      
#   if(i==2){
#     a[[i]] <- tax[ , i]
#   }    
# }
# 
# names(a)=colnames(tax)  
# print(a)


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)
comp_tax_csv

# for(i in a$x){
#   print(subset(comp_tax_csv, comp_tax_csv$genus==i))
#   path=normalizePath(file.path('./New folder/Otu',paste("Asian_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }
# 
# 
# 
# 
# 
# otu=read.csv("New folder/Otu/Asian_otuid_Papillibacter.csv",header=T)
# otu[2]
# 
# a=list()  
# 
# for(i in 1:ncol(otu)) {      
#   if(i==2){
#     a[[i]] <- otu[ , i]
#   }    
# }
# 
# print(a)
# 
# func_df=data.frame(funcFile)
# func_df[1]
# 
# 
# for(i in a[[2]]){
#   print(subset(func_df, func_df[1]==i))
#   # path=normalizePath(file.path('./New folder',paste("Asian_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 
# 
# 
# 
# 
# 
# 
# 


getwd()
fileList=list.files(path="New folder\\AsianCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/AsianCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/AsianOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\AsianOtu",pattern = ".txt")
fileList


func_df=data.frame(funcFile)
func_df[1]

# for(i in fileList){
#   cat(i,"*****************************************************************")
#   path=normalizePath(file.path('./New folder/Otu',paste(i, sep='')))
#   file=(read.table(path))
#   print(file[2])
# }




func_df=data.frame(funcFile)


fileList=list.files(path="New folder\\AsianOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/AsianOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/AsianOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:4){
  path1=normalizePath(file.path('./New folder/AsianOtu/Res',paste("Asian_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/AsianOtu',paste("Asian_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Asian_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/Asian_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/Asian_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)




# #LCBC********************************************************************************************************
# #4 Samples
# LCBC_data <- readRDS("Breed/LCBC.rds")
# 
# top_LCBC <- prune_taxa(names(sort(taxa_sums(LCBC_data),TRUE)[1:50]), LCBC_data)
# #plot_heatmap(top_LCBC)
# LCBC_data
# 
# #SparCC
# net_single_LCBC <- netConstruct(LCBC_data,
#                                  verbose = 3,
#                                  filtTax = "highestFreq",
#                                  filtTaxPar = list(highestFreq = 100),
#                                  # filtTax = "none",
#                                  # filtTaxPar = "totalReads",
#                                  filtSamp = "totalReads",
#                                  filtSampPar = list(totalReads = 1000),
#                                  zeroMethod = "none", normMethod = "none",
#                                  measure = "sparcc",
#                                  sparsMethod = "threshold", thresh = 0.4,
#                                  dissFunc = "signed", 
#                                  seed = 123456)
# 
# saveRDS(net_single_LCBC, "Breed/LCBC Breed/LCBC_Sparcc network.rds")
# 
# props_single_LCBC <- netAnalyze(net_single_LCBC, 
#                                  clustMethod = "cluster_fast_greedy",
#                                  hubPar = "eigenvector", 
#                                  hubQuant = 0.95)
# 
# saveRDS(props_single_LCBC, "Breed/LCBC Breed/LCBC_Sparcc analysis.rds")
# 
# #?summary.microNetProps
# net.summary <- summary(props_single_LCBC)
# net.summary
# 
# p_LCBC <- plot(props_single_LCBC,
#                 shortenLabels = "none",
#                 # labelLength = 16,
#                 # charToRm = "g__",
#                 labelScale = FALSE,
#                 rmSingles = "all",
#                 nodeSize = "eigenvector",
#                 nodeColor = "cluster",
#                 hubBorderCol = "blue",
#                 cexNodes = 1,
#                 cexLabels = 0.5,
#                 edgeWidth = 1,
#                 highlightHubs = TRUE,
#                 cexHubs = 1.5,
#                 # cexHubLabels = 2,
#                 title1 = "LCBC Breed Network on Genus level with SparCC Method.", 
#                 showTitle = TRUE,
#                 cexTitle = 1.5)
# 
# saveRDS(p_LCBC, "Breed/LCBC Breed/LCBC_Sparcc plot.rds")
# 
# from <- p_LCBC$labels$labels1[p_LCBC$q1$Edgelist$from]
# to <- p_LCBC$labels$labels1[p_LCBC$q1$Edgelist$to]
# direction <- vector(mode = "integer")
# 
# for (x in 1:length(from)){
#   if (net_single_LCBC$assoMat1[from[x],to[x]] < 0){
#     direction[x] <- -1 
#   }
#   else if (net_single_LCBC$assoMat1[from[x],to[x]] > 0){
#     direction[x] <- 1 
#   }
# } 
# 
# edges <- data.frame(row.names = c(1:length(p_LCBC$q1$Edgelist$from)))
# edges$from <- p_LCBC$q1$Edgelist$from
# edges$to <- p_LCBC$q1$Edgelist$to
# edges$weight <- p_LCBC$q1$Edgelist$weight
# edges$association <- direction
# 
# write.csv(edges, file = "Breed/LCBC Breed/LCBC edge data.csv", row.names = FALSE)
# 
# hubs <- props_single_LCBC$hubs$hubs1
# write(hubs, "Breed/LCBC Breed/LCBC Breed Hubs.txt")
# 
# node.lables <- p_LCBC$labels$labels1
# clust <- props_single_LCBC$clustering$clust1[node.lables]
# node.hubs <- integer(length((node.lables)))
# node.hubs[match(hubs,node.lables)] <- 1
# 
# 
# nodes <- data.frame(row.names = 1:length(node.lables))
# nodes$index <- 1:length(node.lables)
# nodes$lable <- node.lables
# nodes$cluster <- clust
# nodes$hubs <- node.hubs
# degrees <- props_single_LCBC$centralities$degree1[nodes$lable]
# nodes$degree <- degrees
# 
# write.csv(nodes, file = "Breed/LCBC Breed/LCBC node data.csv", row.names = FALSE)
# 
# hubdata <- subset(nodes, hubs==1, select = -hubs)
# write.csv(hubdata, file = "Breed/LCBC Breed/LCBC hub data.csv", row.names = FALSE)
# 
# 
# 

#Spring#############################################################################################################

#JRJRC**********************************************************************************************************
JRJRC_data <- readRDS("Breed/JRJRC.rds")

top_JRJRC <- prune_taxa(names(sort(taxa_sums(JRJRC_data),TRUE)[1:50]), JRJRC_data)
#plot_heatmap(top_JRJRC)

#SpRING
net_single_JRJRC_Spring <- netConstruct(JRJRC_data, verbose = 3,
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

saveRDS(net_single_JRJRC_Spring, "Breed/JRJRC Breed/JRJRC_Spring network.rds")


# ?netAnalyze
props_single_JRJRC_Spring <- netAnalyze(net_single_JRJRC_Spring, 
                                      centrLCC = TRUE,
                                      clustMethod = "cluster_fast_greedy",
                                      hubPar = "eigenvector",
                                      weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_JRJRC_Spring, "Breed/JRJRC Breed/JRJRC_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_JRJRC_Spring)
net.summary

# ?plot.microNetProps
p_JRJRC_Spring <- plot(props_single_JRJRC_Spring,
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
                     title1 = "JRJRC Breed Network on Genus level with Spring Method.", 
                     showTitle = TRUE,
                     cexTitle = 1.5)

saveRDS(p_JRJRC_Spring, "Breed/JRJRC Breed/JRJRC_Spring plot.rds")

from <- p_JRJRC_Spring$labels$labels1[p_JRJRC_Spring$q1$Edgelist$from]
to <- p_JRJRC_Spring$labels$labels1[p_JRJRC_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_JRJRC_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_JRJRC_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_JRJRC_Spring$q1$Edgelist$from)))
edges$from <- p_JRJRC_Spring$q1$Edgelist$from
edges$to <- p_JRJRC_Spring$q1$Edgelist$to
edges$weight <- p_JRJRC_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/JRJRC Breed/JRJRC_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_JRJRC_Spring$hubs$hubs1
write(hubs, "Breed/JRJRC Breed/JRJRC_Spring Hubs.txt")

node.lables <- p_JRJRC_Spring$labels$labels1
clust <- props_single_JRJRC_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_JRJRC_Spring$centralities$degree1[nodes$lable]
evs=props_single_JRJRC_Spring$centralities$eigenv1[nodes$lable]
betweennesses=props_single_JRJRC_Spring$centralities$between1[nodes$lable]
closenesses=props_single_JRJRC_Spring$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Breed/JRJRC Breed/JRJRC_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/JRJRC Breed/JRJRC_Spring hub data.csv", row.names = FALSE)


#Associations------------------------------------------------------------------------------------------

Genus_Spring_Breed_JRJRC_Top100=data.frame(Node1=net_single_JRJRC_Spring[["edgelist1"]][["v1"]],
                                      Node2=net_single_JRJRC_Spring[["edgelist1"]][["v2"]],
                                      Association=net_single_JRJRC_Spring[["edgelist1"]][["asso"]],
                                      Dissimilarity=net_single_JRJRC_Spring[["edgelist1"]][["diss"]],
                                      Adjecency=net_single_JRJRC_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Breed_JRJRC_Top100=Genus_Spring_Breed_JRJRC_Top100[order(Genus_Spring_Breed_JRJRC_Top100$Association,
                                                            decreasing = T),]

write.csv(Genus_Spring_Breed_JRJRC_Top100,"Associations/Genus_Spring_Breed_JRJRC_Top100.csv",row.names = FALSE)


#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/JRJRC Breed/JRJRC_Spring node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/Spring_JRJRCCluster',paste("Spring_JRJRC_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/JRJRC_Cluster1.csv",header=T)
# tax[2]
# 
# 
# a=list()  
# 
# for(i in 1:ncol(tax)) {      
#   if(i==2){
#     a[[i]] <- tax[ , i]
#   }    
# }
# 
# names(a)=colnames(tax)  
# print(a)


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)
comp_tax_csv

# for(i in a$x){
#   print(subset(comp_tax_csv, comp_tax_csv$genus==i))
#   path=normalizePath(file.path('./New folder/Otu',paste("JRJRC_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }

# otu=read.csv("New folder/Otu/JRJRC_otuid_Papillibacter.csv",header=T)
# otu[2]
# 
# a=list()  
# 
# for(i in 1:ncol(otu)) {      
#   if(i==2){
#     a[[i]] <- otu[ , i]
#   }    
# }
# 
# print(a)
# 
# func_df=data.frame(funcFile)
# func_df[1]
# 
# 
# for(i in a[[2]]){
#   print(subset(func_df, func_df[1]==i))
#   # path=normalizePath(file.path('./New folder',paste("JRJRC_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 



getwd()
fileList=list.files(path="New folder\\Spring_JRJRCCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/Spring_JRJRCCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_JRJRCOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\Spring_JRJRCOtu",pattern = ".txt")
fileList


func_df=data.frame(funcFile)
func_df[1]

# for(i in fileList){
#   cat(i,"*****************************************************************")
#   path=normalizePath(file.path('./New folder/Otu',paste(i, sep='')))
#   file=(read.table(path))
#   print(file[2])
# }

func_df=data.frame(funcFile)


fileList=list.files(path="New folder\\Spring_JRJRCOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/Spring_JRJRCOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_JRJRCOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:9){
  path1=normalizePath(file.path('./New folder/Spring_JRJRCOtu/Res',paste("Spring_JRJRC_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/Spring_JRJRCOtu',paste("Spring_JRJRC_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Spring_JRJRC_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/JRJRC_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/JRJRC_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)







#FRFRC**********************************************************************************************************
FRFRC_data <- readRDS("Breed/FRFRC.rds")

top_FRFRC <- prune_taxa(names(sort(taxa_sums(FRFRC_data),TRUE)[1:50]), FRFRC_data)
#plot_heatmap(top_FRFRC)

#SpRING
net_single_FRFRC_Spring <- netConstruct(FRFRC_data, verbose = 3,
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

saveRDS(net_single_FRFRC_Spring, "Breed/FRFRC Breed/FRFRC_Spring network.rds")


# ?netAnalyze
props_single_FRFRC_Spring <- netAnalyze(net_single_FRFRC_Spring, 
                                        centrLCC = TRUE,
                                        clustMethod = "cluster_fast_greedy",
                                        hubPar = "eigenvector",
                                        weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_FRFRC_Spring, "Breed/FRFRC Breed/FRFRC_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_FRFRC_Spring)
net.summary

# ?plot.microNetProps
p_FRFRC_Spring <- plot(props_single_FRFRC_Spring,
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
                       title1 = "FRFRC Breed Network on Genus level with Spring Method.", 
                       showTitle = TRUE,
                       cexTitle = 1.5)

saveRDS(p_FRFRC_Spring, "Breed/FRFRC Breed/FRFRC_Spring plot.rds")

from <- p_FRFRC_Spring$labels$labels1[p_FRFRC_Spring$q1$Edgelist$from]
to <- p_FRFRC_Spring$labels$labels1[p_FRFRC_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_FRFRC_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_FRFRC_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_FRFRC_Spring$q1$Edgelist$from)))
edges$from <- p_FRFRC_Spring$q1$Edgelist$from
edges$to <- p_FRFRC_Spring$q1$Edgelist$to
edges$weight <- p_FRFRC_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/FRFRC Breed/FRFRC_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_FRFRC_Spring$hubs$hubs1
write(hubs, "Breed/FRFRC Breed/FRFRC_Spring Hubs.txt")

node.lables <- p_FRFRC_Spring$labels$labels1
clust <- props_single_FRFRC_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_FRFRC_Spring$centralities$degree1[nodes$lable]
evs=props_single_FRFRC_Spring$centralities$eigenv1[nodes$lable]
betweennesses=props_single_FRFRC_Spring$centralities$between1[nodes$lable]
closenesses=props_single_FRFRC_Spring$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses

write.csv(nodes, file = "Breed/FRFRC Breed/FRFRC_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/FRFRC Breed/FRFRC_Spring hub data.csv", row.names = FALSE)


#Associations------------------------------------------------------------------------------------------

Genus_Spring_Breed_FRFRC_Top100=data.frame(Node1=net_single_FRFRC_Spring[["edgelist1"]][["v1"]],
                                           Node2=net_single_FRFRC_Spring[["edgelist1"]][["v2"]],
                                           Association=net_single_FRFRC_Spring[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_FRFRC_Spring[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_FRFRC_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Breed_FRFRC_Top100=Genus_Spring_Breed_FRFRC_Top100[order(Genus_Spring_Breed_FRFRC_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Spring_Breed_FRFRC_Top100,"Associations/Genus_Spring_Breed_FRFRC_Top100.csv",row.names = FALSE)


#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/FRFRC Breed/FRFRC_Spring node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/Spring_FRFRCCluster',paste("Spring_FRFRC_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/FRFRC_Cluster1.csv",header=T)
# tax[2]
# 
# 
# a=list()  
# 
# for(i in 1:ncol(tax)) {      
#   if(i==2){
#     a[[i]] <- tax[ , i]
#   }    
# }
# 
# names(a)=colnames(tax)  
# print(a)


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)
comp_tax_csv

# for(i in a$x){
#   print(subset(comp_tax_csv, comp_tax_csv$genus==i))
#   path=normalizePath(file.path('./New folder/Otu',paste("FRFRC_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }

# otu=read.csv("New folder/Otu/FRFRC_otuid_Papillibacter.csv",header=T)
# otu[2]
# 
# a=list()  
# 
# for(i in 1:ncol(otu)) {      
#   if(i==2){
#     a[[i]] <- otu[ , i]
#   }    
# }
# 
# print(a)
# 
# func_df=data.frame(funcFile)
# func_df[1]
# 
# 
# for(i in a[[2]]){
#   print(subset(func_df, func_df[1]==i))
#   # path=normalizePath(file.path('./New folder',paste("FRFRC_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 



getwd()
fileList=list.files(path="New folder\\Spring_FRFRCCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/Spring_FRFRCCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_FRFRCOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\Spring_FRFRCOtu",pattern = ".txt")
fileList


func_df=data.frame(funcFile)
func_df[1]

# for(i in fileList){
#   cat(i,"*****************************************************************")
#   path=normalizePath(file.path('./New folder/Otu',paste(i, sep='')))
#   file=(read.table(path))
#   print(file[2])
# }

func_df=data.frame(funcFile)


fileList=list.files(path="New folder\\Spring_FRFRCOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/Spring_FRFRCOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_FRFRCOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:7){
  path1=normalizePath(file.path('./New folder/Spring_FRFRCOtu/Res',paste("Spring_FRFRC_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/Spring_FRFRCOtu',paste("Spring_FRFRC_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Spring_FRFRC_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/FRFRC_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/FRFRC_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)









#Asian**********************************************************************************************************
Asian_data <- readRDS("Breed/Asian.rds")
#
top_Asian <- prune_taxa(names(sort(taxa_sums(Asian_data),TRUE)[1:50]), Asian_data)
#plot_heatmap(top_Asian)

#SpRING
net_single_Asian_Spring <- netConstruct(Asian_data, verbose = 3,
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
                                                           lambdaseq = exp(seq(log(0.6), log(0.9*0.1), length.out = 20)),
                                                           #lambdaseq="data-specific",
                                                           ncores=1),
                                         sparsMethod = "none",  dissFunc = "signed",
                                         seed = 123456)

saveRDS(net_single_Asian_Spring, "Breed/Asian Breed/Asian_Spring network.rds")
#
#
# ?netAnalyze
props_single_Asian_Spring <- netAnalyze(net_single_Asian_Spring,
                                        centrLCC = TRUE,
                                        clustMethod = "cluster_fast_greedy",
                                        hubPar = "eigenvector",
                                        weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_Asian_Spring, "Breed/Asian Breed/Asian_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Asian_Spring)
net.summary

# ?plot.microNetProps
p_Asian_Spring <- plot(props_single_Asian_Spring,
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
                       title1 = "Asian Breed Network on Genus level with Spring Method.",
                       showTitle = TRUE,
                       cexTitle = 1.5)

saveRDS(p_Asian_Spring, "Breed/Asian Breed/Asian_Spring plot.rds")


from <- p_Asian_Spring$labels$labels1[p_Asian_Spring$q1$Edgelist$from]
to <- p_Asian_Spring$labels$labels1[p_Asian_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Asian_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Asian_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Asian_Spring$q1$Edgelist$from)))
edges$from <- p_Asian_Spring$q1$Edgelist$from
edges$to <- p_Asian_Spring$q1$Edgelist$to
edges$weight <- p_Asian_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/Asian Breed/Asian_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_Asian_Spring$hubs$hubs1
write(hubs, "Breed/Asian Breed/Asian_Spring Hubs.txt")

node.lables <- p_Asian_Spring$labels$labels1
clust <- props_single_Asian_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Asian_Spring$centralities$degree1[nodes$lable]
evs=props_single_Asian_Spring$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Asian_Spring$centralities$between1[nodes$lable]
closenesses=props_single_Asian_Spring$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses

write.csv(nodes, file = "Breed/Asian Breed/Asian_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/Asian Breed/Asian_Spring hub data.csv", row.names = FALSE)


#Associations------------------------------------------------------------------------------------------

Genus_Spring_Breed_Asian_Top100=data.frame(Node1=net_single_Asian_Spring[["edgelist1"]][["v1"]],
                                           Node2=net_single_Asian_Spring[["edgelist1"]][["v2"]],
                                           Association=net_single_Asian_Spring[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_Asian_Spring[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_Asian_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Breed_Asian_Top100=Genus_Spring_Breed_Asian_Top100[order(Genus_Spring_Breed_Asian_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Spring_Breed_Asian_Top100,"Associations/Genus_Spring_Breed_Asian_Top100.csv",row.names = FALSE)


#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/Asian Breed/Asian_Spring node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/Spring_AsianCluster',paste("Spring_Asian_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/Asian_Cluster1.csv",header=T)
# tax[2]
# 
# 
# a=list()  
# 
# for(i in 1:ncol(tax)) {      
#   if(i==2){
#     a[[i]] <- tax[ , i]
#   }    
# }
# 
# names(a)=colnames(tax)  
# print(a)


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)
comp_tax_csv

# for(i in a$x){
#   print(subset(comp_tax_csv, comp_tax_csv$genus==i))
#   path=normalizePath(file.path('./New folder/Otu',paste("Asian_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }

# otu=read.csv("New folder/Otu/Asian_otuid_Papillibacter.csv",header=T)
# otu[2]
# 
# a=list()  
# 
# for(i in 1:ncol(otu)) {      
#   if(i==2){
#     a[[i]] <- otu[ , i]
#   }    
# }
# 
# print(a)
# 
# func_df=data.frame(funcFile)
# func_df[1]
# 
# 
# for(i in a[[2]]){
#   print(subset(func_df, func_df[1]==i))
#   # path=normalizePath(file.path('./New folder',paste("Asian_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 



getwd()
fileList=list.files(path="New folder\\Spring_AsianCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/Spring_AsianCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_AsianOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\Spring_AsianOtu",pattern = ".txt")
fileList


func_df=data.frame(funcFile)
func_df[1]

# for(i in fileList){
#   cat(i,"*****************************************************************")
#   path=normalizePath(file.path('./New folder/Otu',paste(i, sep='')))
#   file=(read.table(path))
#   print(file[2])
# }

func_df=data.frame(funcFile)


fileList=list.files(path="New folder\\Spring_AsianOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/Spring_AsianOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_AsianOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:7){
  path1=normalizePath(file.path('./New folder/Spring_AsianOtu/Res',paste("Spring_Asian_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/Spring_AsianOtu',paste("Spring_Asian_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Spring_Asian_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/Asian_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/Asian_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)







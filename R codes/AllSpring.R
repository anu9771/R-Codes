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
tax_table(filt_complete_agg)

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

comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)



#Early Lactation************************************************************************
Early_data <- readRDS("Phase/Early Phase.rds")

top_Early <- prune_taxa(names(sort(taxa_sums(Early_data),TRUE)[1:50]), Early_data)
#plot_heatmap(top_Early)
Early_data

#SpRING
net_single_Early_Spring <- netConstruct(Early_data, verbose = 3,
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

saveRDS(net_single_Early_Spring, "Phase/Early Phase/Early_Spring network.rds")


# ?netAnalyze
props_single_Early_Spring <- netAnalyze(net_single_Early_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_Early_Spring, "Phase/Early Phase/Early_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Early_Spring)
net.summary

# ?plot.microNetProps
p_Early_Spring <- plot(props_single_Early_Spring,
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
          title1 = "Early Phase Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_Early_Spring, "Phase/Early Phase/Early_Spring plot.rds")

from <- p_Early_Spring$labels$labels1[p_Early_Spring$q1$Edgelist$from]
to <- p_Early_Spring$labels$labels1[p_Early_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Early_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Early_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Early_Spring$q1$Edgelist$from)))
edges$from <- p_Early_Spring$q1$Edgelist$from
edges$to <- p_Early_Spring$q1$Edgelist$to
edges$weight <- p_Early_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/Early Phase/Early_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_Early_Spring$hubs$hubs1
write(hubs, "Phase/Early Phase/Early_Spring Phase Hubs.txt")

node.lables <- p_Early_Spring$labels$labels1
clust <- props_single_Early_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Early_Spring$centralities$degree1[nodes$lable]
evs=props_single_Early_Spring$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Early_Spring$centralities$between1[nodes$lable]
closenesses=props_single_Early_Spring$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses



write.csv(nodes, file = "Phase/Early Phase/Early_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/Early Phase/Early_Spring hub data.csv", row.names = FALSE)


#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Early Phase/Early_Spring node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/Spring_EarlyCluster',paste("Spring_Early_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/Early_Cluster1.csv",header=T)
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
#   path=normalizePath(file.path('./New folder/Otu',paste("Early_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }

# otu=read.csv("New folder/Otu/Early_otuid_Papillibacter.csv",header=T)
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
#   # path=normalizePath(file.path('./New folder',paste("Early_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 



getwd()
fileList=list.files(path="New folder\\Spring_EarlyCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/Spring_EarlyCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_EarlyOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\Spring_EarlyOtu",pattern = ".txt")
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


fileList=list.files(path="New folder\\Spring_EarlyOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/Spring_EarlyOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_EarlyOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:8){
  path1=normalizePath(file.path('./New folder/Spring_EarlyOtu/Res',paste("Spring_Early_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/Spring_EarlyOtu',paste("Spring_Early_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Spring_Early_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/Early_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/Early_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)




#Mid Lactation**************************************************************************
Mid_data <- readRDS("Phase/Mid Phase.rds")
Mid_data
top_Mid <- prune_taxa(names(sort(taxa_sums(Mid_data),TRUE)[1:50]), Mid_data)
top_Mid
#plot_heatmap(top_Mid)

#SpRING
net_single_Mid_Spring <- netConstruct(Mid_data, verbose = 3,
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

saveRDS(net_single_Mid_Spring, "Phase/Mid Phase/Mid_Spring network.rds")


# ?netAnalyze
props_single_Mid_Spring <- netAnalyze(net_single_Mid_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_Mid_Spring, "Phase/Mid Phase/Mid_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Mid_Spring)
net.summary

# ?plot.microNetProps
p_Mid_Spring <- plot(props_single_Mid_Spring,
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
          title1 = "Mid Phase Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_Mid_Spring, "Phase/Mid Phase/Mid_Spring plot.rds")

from <- p_Mid_Spring$labels$labels1[p_Mid_Spring$q1$Edgelist$from]
to <- p_Mid_Spring$labels$labels1[p_Mid_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Mid_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Mid_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Mid_Spring$q1$Edgelist$from)))
edges$from <- p_Mid_Spring$q1$Edgelist$from
edges$to <- p_Mid_Spring$q1$Edgelist$to
edges$weight <- p_Mid_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/Mid Phase/Mid_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_Mid_Spring$hubs$hubs1
write(hubs, "Phase/Mid Phase/Mid_Spring Phase Hubs.txt")

node.lables <- p_Mid_Spring$labels$labels1
clust <- props_single_Mid_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Mid_Spring$centralities$degree1[nodes$lable]
evs=props_single_Mid_Spring$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Mid_Spring$centralities$between1[nodes$lable]
closenesses=props_single_Mid_Spring$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses

write.csv(nodes, file = "Phase/Mid Phase/Mid_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/Mid Phase/Mid_Spring hub data.csv", row.names = FALSE)




#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Mid Phase/Mid_Spring node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/Spring_MidCluster',paste("Spring_Mid_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/Mid_Cluster1.csv",header=T)
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
#   path=normalizePath(file.path('./New folder/Otu',paste("Mid_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }

# otu=read.csv("New folder/Otu/Mid_otuid_Papillibacter.csv",header=T)
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
#   # path=normalizePath(file.path('./New folder',paste("Mid_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 



getwd()
fileList=list.files(path="New folder\\Spring_MidCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/Spring_MidCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_MidOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\Spring_MidOtu",pattern = ".txt")
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


fileList=list.files(path="New folder\\Spring_MidOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/Spring_MidOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_MidOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:15){
  path1=normalizePath(file.path('./New folder/Spring_MidOtu/Res',paste("Spring_Mid_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/Spring_MidOtu',paste("Spring_Mid_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Spring_Mid_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/Mid_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/Mid_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)




#Late Phase***********************************************************************
Late_data <- readRDS("Phase/Late Phase.rds")
Late_data
top_Late <- prune_taxa(names(sort(taxa_sums(Late_data),TRUE)[1:50]), Late_data)
#plot_heatmap(top_Late)
Late_data

#SpRING
net_single_Late_Spring <- netConstruct(Late_data, verbose = 3,
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

saveRDS(net_single_Late_Spring, "Phase/Late Phase/Late_Spring network.rds")


# ?netAnalyze
props_single_Late_Spring <- netAnalyze(net_single_Late_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_Late_Spring, "Phase/Late Phase/Late_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Late_Spring)
net.summary

# ?plot.microNetProps
p_Late_Spring <- plot(props_single_Late_Spring,
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
          title1 = "Late Phase Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_Late_Spring, "Phase/Late Phase/Late_Spring plot.rds")

from <- p_Late_Spring$labels$labels1[p_Late_Spring$q1$Edgelist$from]
to <- p_Late_Spring$labels$labels1[p_Late_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Late_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Late_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Late_Spring$q1$Edgelist$from)))
edges$from <- p_Late_Spring$q1$Edgelist$from
edges$to <- p_Late_Spring$q1$Edgelist$to
edges$weight <- p_Late_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/Late Phase/Late_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_Late_Spring$hubs$hubs1
write(hubs, "Phase/Late Phase/Late_Spring Phase Hubs.txt")

node.lables <- p_Late_Spring$labels$labels1
clust <- props_single_Late_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Late_Spring$centralities$degree1[nodes$lable]
evs=props_single_Late_Spring$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Late_Spring$centralities$between1[nodes$lable]
closenesses=props_single_Late_Spring$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses

write.csv(nodes, file = "Phase/Late Phase/Late_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/Late Phase/Late_Spring hub data.csv", row.names = FALSE)



#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Late Phase/Late_Spring node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/Spring_LateCluster',paste("Spring_Late_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}
# 
# tax=read.csv("New folder/Cluster/Late_Cluster1.csv",header=T)
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
#   path=normalizePath(file.path('./New folder/Otu',paste("Late_otuid_",i,'.csv', sep='')))
#   write.csv(subset(comp_tax_csv, comp_tax_csv$genus==i),path)
# }

# otu=read.csv("New folder/Otu/Late_otuid_Papillibacter.csv",header=T)
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
#   # path=normalizePath(file.path('./New folder',paste("Late_function_",i,'.csv', sep='')))
#   # write.csv(subset(func_df, func_df[1]==i),path)
# }
# 



getwd()
fileList=list.files(path="New folder\\Spring_LateCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/Spring_LateCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_LateOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\Spring_LateOtu",pattern = ".txt")
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


fileList=list.files(path="New folder\\Spring_LateOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/Spring_LateOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/Spring_LateOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:11){
  path1=normalizePath(file.path('./New folder/Spring_LateOtu/Res',paste("Spring_Late_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/Spring_LateOtu',paste("Spring_Late_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Spring_Late_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
  newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
  newDF=unique(newDF)
  write.table(newDF,path3)
  
}

# resulrDataFile=read.table("New folder/Otu/Res/Late_Cluster1.csv.txt_results.txt",col.names = c("1","2","3"))
# resultDataFrame=data.frame(resulrDataFile)
# 
# otuGeneraFile=read.table("New folder/Otu/Late_Cluster1.csv.txt",col.names = c("1","2","3"))
# otuGeneraDataFrame=data.frame(otuGeneraFile)
# 
# resultDataFrame$X2
# otuGeneraDataFrame$X2
# 
# newDF=merge(resultDataFrame,otuGeneraDataFrame,by="X2")
# newDF=unique(newDF)







#FJC Breed************************************************************************
FJC_data <- readRDS("Breed/FJC Breed.rds")

top_FJC <- prune_taxa(names(sort(taxa_sums(FJC_data),TRUE)[1:50]), FJC_data)
#plot_heatmap(top_FJC)

#SpRING
net_single_FJC_Spring <- netConstruct(FJC_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
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

saveRDS(net_single_FJC_Spring, "Breed/FJC Breed/FJC_Spring network.rds")


# ?netAnalyze
props_single_FJC_Spring <- netAnalyze(net_single_FJC_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_FJC_Spring, "Breed/FJC Breed/FJC_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_FJC_Spring)
net.summary

# ?plot.microNetProps
p_FJC_Spring <- plot(props_single_FJC_Spring,
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
          title1 = "FJC Breed Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_FJC_Spring, "Breed/FJC Breed/FJC_Spring plot.rds")

from <- p_FJC_Spring$labels$labels1[p_FJC_Spring$q1$Edgelist$from]
to <- p_FJC_Spring$labels$labels1[p_FJC_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_FJC_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_FJC_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_FJC_Spring$q1$Edgelist$from)))
edges$from <- p_FJC_Spring$q1$Edgelist$from
edges$to <- p_FJC_Spring$q1$Edgelist$to
edges$weight <- p_FJC_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/FJC Breed/FJC_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_FJC_Spring$hubs$hubs1
write(hubs, "Breed/FJC Breed/FJC_Spring Hubs.txt")

node.lables <- p_FJC_Spring$labels$labels1
clust <- props_single_FJC_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_FJC_Spring$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Breed/FJC Breed/FJC_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/FJC Breed/FJC_Spring hub data.csv", row.names = FALSE)


#FC Breed************************************************************************
#4 Samples
#Network has no edges


FC_data <- readRDS("Breed/FC Breed.rds")

top_FC <- prune_taxa(names(sort(taxa_sums(FC_data),TRUE)[1:50]), FC_data)
plot_heatmap(top_FC)

#SpRING
net_single_FC_Spring <- netConstruct(FC_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
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

saveRDS(net_single_FC_Spring, "Breed/FC Breed/FC_Spring network.rds")


# ?netAnalyze
props_single_FC_Spring <- netAnalyze(net_single_FC_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_FC_Spring, "Breed/FC Breed/FC_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_FC_Spring)
net.summary

# ?plot.microNetProps
p_FC_Spring <- plot(props_single_FC_Spring,
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
          title1 = "FC Breed Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_FC_Spring, "Breed/FC Breed/FC_Spring plot.rds")

from <- p_FC_Spring$labels$labels1[p_FC_Spring$q1$Edgelist$from]
to <- p_FC_Spring$labels$labels1[p_FC_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_FC_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_FC_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_FC_Spring$q1$Edgelist$from)))
edges$from <- p_FC_Spring$q1$Edgelist$from
edges$to <- p_FC_Spring$q1$Edgelist$to
edges$weight <- p_FC_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/FC Breed/FC_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_FC_Spring$hubs$hubs1
write(hubs, "Breed/FC Breed/FC_Spring Hubs.txt")

node.lables <- p_FC_Spring$labels$labels1
clust <- props_single_FC_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_FC_Spring$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Breed/FC Breed/FC_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/FC Breed/FC_Spring hub data.csv", row.names = FALSE)



#SC Breed************************************************************************
#2 Samples
#Abnormal Threshold Output



#JSC Breed************************************************************************
# 15 +4 Samples
# SPRING Output given


JSC_data <- readRDS("Breed/JSC Breed.rds")

top_JSC <- prune_taxa(names(sort(taxa_sums(JSC_data),TRUE)[1:50]), JSC_data)
#plot_heatmap(top_JSC)

#SpRING
net_single_JSC_Spring <- netConstruct(JSC_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
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

saveRDS(net_single_JSC_Spring, "Breed/JSC Breed/JSC_Spring network.rds")


# ?netAnalyze
props_single_JSC_Spring <- netAnalyze(net_single_JSC_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_JSC_Spring, "Breed/JSC Breed/JSC_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_JSC_Spring)
net.summary

# ?plot.microNetProps
p_JSC_Spring <- plot(props_single_JSC_Spring,
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
          title1 = "JSC Breed Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_JSC_Spring, "Breed/JSC Breed/JSC_Spring plot.rds")


from <- p_JSC_Spring$labels$labels1[p_JSC_Spring$q1$Edgelist$from]
to <- p_JSC_Spring$labels$labels1[p_JSC_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_JSC_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_JSC_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_JSC_Spring$q1$Edgelist$from)))
edges$from <- p_JSC_Spring$q1$Edgelist$from
edges$to <- p_JSC_Spring$q1$Edgelist$to
edges$weight <- p_JSC_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/JSC Breed/JSC_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_JSC_Spring$hubs$hubs1
write(hubs, "Breed/JSC Breed/JSC_Spring Hubs.txt")

node.lables <- p_JSC_Spring$labels$labels1
clust <- props_single_JSC_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_JSC_Spring$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Breed/JSC Breed/JSC_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/JSC Breed/JSC_Spring hub data.csv", row.names = FALSE)





#BC Breed************************************************************************
# 4 Samples
# Threshold Output given

BC_data <- readRDS("Breed/BC Breed.rds")

top_BC <- prune_taxa(names(sort(taxa_sums(BC_data),TRUE)[1:50]), BC_data)
#plot_heatmap(top_BC)

#SpRING
net_single_BC_Spring <- netConstruct(BC_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
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

saveRDS(net_single_BC_Spring, "Breed/BC Breed/BC_Spring network.rds")


# ?netAnalyze
props_single_BC_Spring <- netAnalyze(net_single_BC_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_BC_Spring, "Breed/BC Breed/BC_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_BC_Spring)
net.summary

# ?plot.microNetProps
p_BC_Spring <- plot(props_single_BC_Spring,
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
          title1 = "BC Breed Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_BC_Spring, "Breed/BC Breed/BC_Spring plot.rds")


from <- p_BC_Spring$labels$labels1[p_BC_Spring$q1$Edgelist$from]
to <- p_BC_Spring$labels$labels1[p_BC_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_BC_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_BC_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_BC_Spring$q1$Edgelist$from)))
edges$from <- p_BC_Spring$q1$Edgelist$from
edges$to <- p_BC_Spring$q1$Edgelist$to
edges$weight <- p_BC_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/BC Breed/BC_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_BC_Spring$hubs$hubs1
write(hubs, "Breed/BC Breed/BC_Spring Hubs.txt")

node.lables <- p_BC_Spring$labels$labels1
clust <- props_single_BC_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_BC_Spring$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Breed/BC Breed/BC_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/BC Breed/BC_Spring hub data.csv", row.names = FALSE)



#FSC Breed************************************************************************
# 2 Samples
# Threshold Output Abnormal


#AFS Breed************************************************************************
# 2 Samples
# Threshold Output abnormal


#JC Breed************************************************************************
# 6 Samples
# Network is empty

JC_data <- readRDS("Breed/JC Breed.rds")

top_JC <- prune_taxa(names(sort(taxa_sums(JC_data),TRUE)[1:50]), JC_data)
#plot_heatmap(top_JC)

#SpRING
net_single_JC_Spring <- netConstruct(JC_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
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

saveRDS(net_single_JC_Spring, "Breed/JC Breed/JC_Spring network.rds")


# ?netAnalyze
props_single_JC_Spring <- netAnalyze(net_single_JC_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_JC_Spring, "Breed/JC Breed/JC_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_JC_Spring)
net.summary

# ?plot.microNetProps
p_JC_Spring <- plot(props_single_JC_Spring,
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
          title1 = "JC Breed Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_JC_Spring, "Breed/JC Breed/JC_Spring plot.rds")


from <- p_JC_Spring$labels$labels1[p_JC_Spring$q1$Edgelist$from]
to <- p_JC_Spring$labels$labels1[p_JC_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_JC_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_JC_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_JC_Spring$q1$Edgelist$from)))
edges$from <- p_JC_Spring$q1$Edgelist$from
edges$to <- p_JC_Spring$q1$Edgelist$to
edges$weight <- p_JC_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/JC Breed/JC_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_JC_Spring$hubs$hubs1
write(hubs, "Breed/JC Breed/JC_Spring Hubs.txt")

node.lables <- p_JC_Spring$labels$labels1
clust <- props_single_JC_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_JC_Spring$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Breed/JC Breed/JC_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/JC Breed/JC_Spring hub data.csv", row.names = FALSE)



# S Breed************************************************************************
# 2 Samples
# Threshold Output abnormal


# A Breed************************************************************************
# 0 Samples after filtering
# Threshold Output is null 


# J Breed************************************************************************
# 5 Samples 
# Threshold Output is given 

J_data <- readRDS("Breed/J Breed.rds")

top_J <- prune_taxa(names(sort(taxa_sums(J_data),TRUE)[1:50]), J_data)
#plot_heatmap(top_J)

#SpRING
net_single_J_Spring <- netConstruct(J_data, verbose = 3,
                           # filtTax = "highestFreq",
                           # filtTaxPar = list(highestFreq = 100),
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

saveRDS(net_single_J_Spring, "Breed/J Breed/J_Spring network.rds")


# ?netAnalyze
props_single_J_Spring <- netAnalyze(net_single_J_Spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

saveRDS(props_single_J_Spring, "Breed/J Breed/J_Spring analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_J_Spring)
net.summary

# ?plot.microNetProps
p_J_Spring <- plot(props_single_J_Spring,
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
          title1 = "J Breed Network on Genus level with Spring Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_J_Spring, "Breed/J Breed/J_Spring plot.rds")


from <- p_J_Spring$labels$labels1[p_J_Spring$q1$Edgelist$from]
to <- p_J_Spring$labels$labels1[p_J_Spring$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_J_Spring$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_J_Spring$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_J_Spring$q1$Edgelist$from)))
edges$from <- p_J_Spring$q1$Edgelist$from
edges$to <- p_J_Spring$q1$Edgelist$to
edges$weight <- p_J_Spring$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Breed/J Breed/J_Spring edge data.csv", row.names = FALSE)

hubs <- props_single_J_Spring$hubs$hubs1
write(hubs, "Breed/J Breed/J_Spring Hubs.txt")

node.lables <- p_J_Spring$labels$labels1
clust <- props_single_J_Spring$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_J_Spring$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Breed/J Breed/J_Spring node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Breed/J Breed/J_Spring hub data.csv", row.names = FALSE)

###################################################################
#For Associations

Genus_Spring_Phase_Early_Top100=data.frame(Node1=net_single_Early_Spring[["edgelist1"]][["v1"]],
                                          Node2=net_single_Early_Spring[["edgelist1"]][["v2"]],
                                          Association=net_single_Early_Spring[["edgelist1"]][["asso"]],
                                          Dissimilarity=net_single_Early_Spring[["edgelist1"]][["diss"]],
                                          Adjecency=net_single_Early_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Phase_Early_Top100=Genus_Spring_Phase_Early_Top100[order(Genus_Spring_Phase_Early_Top100$Association,
                                                                    decreasing = T),]

write.csv(Genus_Spring_Phase_Early_Top100,"Associations/Genus_Spring_Phase_Early_Top100.csv",row.names = FALSE)

############################

Genus_Spring_Phase_Mid_Top100=data.frame(Node1=net_single_Mid_Spring[["edgelist1"]][["v1"]],
                                           Node2=net_single_Mid_Spring[["edgelist1"]][["v2"]],
                                           Association=net_single_Mid_Spring[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_Mid_Spring[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_Mid_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Phase_Mid_Top100=Genus_Spring_Phase_Mid_Top100[order(Genus_Spring_Phase_Mid_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Spring_Phase_Mid_Top100,"Associations/Genus_Spring_Phase_Mid_Top100.csv",row.names = FALSE)

############################

Genus_Spring_Phase_Late_Top100=data.frame(Node1=net_single_Late_Spring[["edgelist1"]][["v1"]],
                                           Node2=net_single_Late_Spring[["edgelist1"]][["v2"]],
                                           Association=net_single_Late_Spring[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_Late_Spring[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_Late_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Phase_Late_Top100=Genus_Spring_Phase_Late_Top100[order(Genus_Spring_Phase_Late_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Spring_Phase_Late_Top100,"Associations/Genus_Spring_Phase_Late_Top100.csv",row.names = FALSE)

############################

Genus_Spring_Breed_FJC_Top100=data.frame(Node1=net_single_FJC_Spring[["edgelist1"]][["v1"]],
                                           Node2=net_single_FJC_Spring[["edgelist1"]][["v2"]],
                                           Association=net_single_FJC_Spring[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_FJC_Spring[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_FJC_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Breed_FJC_Top100=Genus_Spring_Breed_FJC_Top100[order(Genus_Spring_Breed_FJC_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Spring_Breed_FJC_Top100,"Associations/Genus_Spring_Breed_FJC_Top100.csv",row.names = FALSE)

############################

Genus_Spring_Breed_JSC_Top100=data.frame(Node1=net_single_JSC_Spring[["edgelist1"]][["v1"]],
                                         Node2=net_single_JSC_Spring[["edgelist1"]][["v2"]],
                                         Association=net_single_JSC_Spring[["edgelist1"]][["asso"]],
                                         Dissimilarity=net_single_JSC_Spring[["edgelist1"]][["diss"]],
                                         Adjecency=net_single_JSC_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Breed_JSC_Top100=Genus_Spring_Breed_JSC_Top100[order(Genus_Spring_Breed_JSC_Top100$Association,
                                                                  decreasing = T),]

write.csv(Genus_Spring_Breed_JSC_Top100,"Associations/Genus_Spring_Breed_JSC_Top100.csv",row.names = FALSE)

############################

Genus_Spring_Breed_BC_Top100=data.frame(Node1=net_single_BC_Spring[["edgelist1"]][["v1"]],
                                         Node2=net_single_BC_Spring[["edgelist1"]][["v2"]],
                                         Association=net_single_BC_Spring[["edgelist1"]][["asso"]],
                                         Dissimilarity=net_single_BC_Spring[["edgelist1"]][["diss"]],
                                         Adjecency=net_single_BC_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Breed_BC_Top100=Genus_Spring_Breed_BC_Top100[order(Genus_Spring_Breed_BC_Top100$Association,
                                                                  decreasing = T),]

write.csv(Genus_Spring_Breed_BC_Top100,"Associations/Genus_Spring_Breed_BC_Top100.csv",row.names = FALSE)

############################

Genus_Spring_Breed_JC_Top100=data.frame(Node1=net_single_JC_Spring[["edgelist1"]][["v1"]],
                                         Node2=net_single_JC_Spring[["edgelist1"]][["v2"]],
                                         Association=net_single_JC_Spring[["edgelist1"]][["asso"]],
                                         Dissimilarity=net_single_JC_Spring[["edgelist1"]][["diss"]],
                                         Adjecency=net_single_JC_Spring[["edgelist1"]][["adja"]]
)

Genus_Spring_Breed_JC_Top100=Genus_Spring_Breed_JC_Top100[order(Genus_Spring_Breed_JC_Top100$Association,
                                                                  decreasing = T),]

write.csv(Genus_Spring_Breed_JC_Top100,"Associations/Genus_Spring_Breed_JC_Top100.csv",row.names = FALSE)



















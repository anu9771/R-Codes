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

# comp_tax=tax_table(complete_physeq)
# comp_tax_df=data.frame(rbind(comp_tax))
# write.csv(comp_tax_df,"Faprotax/tax/comp_tax.csv")

comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)




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

complete_agg <- aggregate_taxa(complete_agg, 'Genus')
complete_agg

agg_tax=tax_table(complete_agg)
agg_tax_df=data.frame(rbind(agg_tax))
write.csv(agg_tax_df,"Faprotax/tax/agg_tax.csv")

agg_otu=otu_table(complete_agg)
agg_otu_df=data.frame(rbind(agg_otu))
write.csv(agg_otu_df,"Faprotax/otu/agg_otu.csv")

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

#SparCC
net_single_Early <- netConstruct(Early_data,
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

saveRDS(net_single_Early, "Phase/Early Phase/Early_Sparcc network.rds")

props_single_Early <- netAnalyze(net_single_Early, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single_Early, "Phase/Early Phase/Early_Sparcc analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single_Early)
net.summary

p_Early <- plot(props_single_Early,
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
          title1 = "Early Phase Network on Genus level with SparCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_Early, "Phase/Early Phase/Early_Sparcc plot.rds")

from <- p_Early$labels$labels1[p_Early$q1$Edgelist$from]
to <- p_Early$labels$labels1[p_Early$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Early$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Early$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Early$q1$Edgelist$from)))
edges$from <- p_Early$q1$Edgelist$from
edges$to <- p_Early$q1$Edgelist$to
edges$weight <- p_Early$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/Early Phase/Early edge data.csv", row.names = FALSE)

hubs <- props_single_Early$hubs$hubs1
write(hubs, "Phase/Early Phase/Early Phase Hubs.txt")

node.lables <- p_Early$labels$labels1
clust <- props_single_Early$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Early$centralities$degree1[nodes$lable]
evs=props_single_Early$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Early$centralities$between1[nodes$lable]
closenesses=props_single_Early$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses


write.csv(nodes, file = "Phase/Early Phase/Early node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/Early Phase/Early hub data.csv", row.names = FALSE)

#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Early Phase/Early node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/EarlyCluster',paste("Early_Cluster",i,'.csv', sep='')))
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
# 
# 
# 
# 
# 
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
# 
# 
# 
# 
# 
# 
# 


getwd()
fileList=list.files(path="New folder\\EarlyCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/EarlyCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/EarlyOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\EarlyOtu",pattern = ".txt")
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


fileList=list.files(path="New folder\\EarlyOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/EarlyOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/EarlyOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:6){
  path1=normalizePath(file.path('./New folder/EarlyOtu/Res',paste("Early_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/EarlyOtu',paste("Early_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Early_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
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
comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)

Mid_data <- readRDS("Phase/Mid Phase.rds")
top_Mid <- prune_taxa(names(sort(taxa_sums(Mid_data),TRUE)[1:50]), Mid_data)
#plot_heatmap(top_Mid)

# ?netConstruct

#SparCC
net_single_Mid <- netConstruct(Mid_data,
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

saveRDS(net_single_Mid, "Phase/Mid Phase/Mid network.rds")

props_single_Mid <- netAnalyze(net_single_Mid, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single_Mid, "Phase/Mid Phase/Mid analysis.rds")

p_Mid <- plot(props_single_Mid,
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
          title1 = "Mid Phase Network on Genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_Mid, "Phase/Mid Phase/Mid plot.rds")

from <- p_Mid$labels$labels1[p_Mid$q1$Edgelist$from]
to <- p_Mid$labels$labels1[p_Mid$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Mid$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Mid$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Mid$q1$Edgelist$from)))
edges$from <- p_Mid$q1$Edgelist$from
edges$to <- p_Mid$q1$Edgelist$to
edges$weight <- p_Mid$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/Mid Phase/Mid edge data.csv", row.names = FALSE)

hubs <- props_single_Mid$hubs$hubs1
write(hubs, "Phase/Mid Phase/Mid Phase Hubs.txt")

node.lables <- p_Mid$labels$labels1
clust <- props_single_Mid$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Mid$centralities$degree1[nodes$lable]
evs=props_single_Mid$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Mid$centralities$between1[nodes$lable]
closenesses=props_single_Mid$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses

write.csv(nodes, file = "Phase/Mid Phase/Mid node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/Mid Phase/Mid hub data.csv", row.names = FALSE)


#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Mid Phase/Mid node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/MidCluster',paste("Mid_Cluster",i,'.csv', sep='')))
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
# 
# 
# 
# 
# 
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
# 
# 
# 
# 
# 
# 
# 


getwd()
fileList=list.files(path="New folder\\MidCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/MidCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/MidOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\MidOtu",pattern = ".txt")
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


fileList=list.files(path="New folder\\MidOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/MidOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/MidOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:6){
  path1=normalizePath(file.path('./New folder/MidOtu/Res',paste("Mid_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/MidOtu',paste("Mid_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Mid_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
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






#Late Phase***********************************************************************
comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)

Late_data <- readRDS("Phase/Late Phase.rds")
top_Late <- prune_taxa(names(sort(taxa_sums(Late_data),TRUE)[1:50]), Late_data)
#plot_heatmap(top_Late)

# ?netConstruct

#SPARCC
net_single_Late <- netConstruct(Late_data,
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

saveRDS(net_single_Late, "Phase/Late Phase/Late network.rds")


props_single_Late <- netAnalyze(net_single_Late, 
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", 
                           hubQuant = 0.95)

saveRDS(props_single_Late, "Phase/Late Phase/Late analysis.rds")


p_Late <- plot(props_single_Late,
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
          title1 = "Late Phase Network on Genus level with SPARCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p_Late, "Phase/Late Phase/Late plot.rds")

from <- p_Late$labels$labels1[p_Late$q1$Edgelist$from]
to <- p_Late$labels$labels1[p_Late$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
  if (net_single_Late$assoMat1[from[x],to[x]] < 0){
    direction[x] <- -1 
  }
  else if (net_single_Late$assoMat1[from[x],to[x]] > 0){
    direction[x] <- 1 
  }
} 

edges <- data.frame(row.names = c(1:length(p_Late$q1$Edgelist$from)))
edges$from <- p_Late$q1$Edgelist$from
edges$to <- p_Late$q1$Edgelist$to
edges$weight <- p_Late$q1$Edgelist$weight
edges$association <- direction

write.csv(edges, file = "Phase/Late Phase/Late edge data.csv", row.names = FALSE)

hubs <- props_single_Late$hubs$hubs1
write(hubs, "Phase/Late Phase/Late Phase Hubs.txt")

node.lables <- p_Late$labels$labels1
clust <- props_single_Late$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single_Late$centralities$degree1[nodes$lable]
evs=props_single_Late$centralities$eigenv1[nodes$lable]
betweennesses=props_single_Late$centralities$between1[nodes$lable]
closenesses=props_single_Late$centralities$close1[nodes$lable]
nodes$degree <- degrees
nodes$ev=evs
nodes$between=betweennesses
nodes$close=closenesses

write.csv(nodes, file = "Phase/Late Phase/Late node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Phase/Late Phase/Late hub data.csv", row.names = FALSE)


#Functions--------------------------------------------------------------------------------------

funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Late Phase/Late node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/LateCluster',paste("Late_Cluster",i,'.csv', sep='')))
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
# 
# 
# 
# 
# 
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
# 
# 
# 
# 
# 
# 
# 


getwd()
fileList=list.files(path="New folder\\LateCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/LateCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    #print(subset(comp_tax_csv, comp_tax_csv$genus==j))
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/LateOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\LateOtu",pattern = ".txt")
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


fileList=list.files(path="New folder\\LateOtu",pattern = ".txt")
fileList

for(file in fileList){
  cat(file,"*******************************************************************")
  path=normalizePath(file.path('./New folder/LateOtu',paste(file, sep='')))
  tab=read.table(path,header=F)
  tab=data.frame(otuid=tab[2],genus=tab[3])
  
  for(i in tab[1]$V2){
    print(subset(func_df, func_df[1]==i))
    x=subset(func_df, func_df[1]==i)
    print(x)
    path=normalizePath(file.path('./New folder/LateOtu/Res',paste(file,"_results.txt", sep='')))
    write.table(subset(func_df, func_df[1]==i),path,append=TRUE,col.names = F)
    
    
  }
}


for(i in 1:4){
  path1=normalizePath(file.path('./New folder/LateOtu/Res',paste("Late_Cluster",i,".csv.txt_results.txt", sep='')))
  path2=normalizePath(file.path('./New folder/LateOtu',paste("Late_Cluster",i,".csv.txt", sep='')))
  resulrDataFile=read.table(path1,col.names = c("1","2","3"))
  resultDataFrame=data.frame(resulrDataFile)
  
  otuGeneraFile=read.table(path2,col.names = c("1","2","3"))
  otuGeneraDataFrame=data.frame(otuGeneraFile)
  
  path3=normalizePath(file.path('./New folder/FinalRes',paste("Late_Cluster",i,".csv.txt_Finalresults.txt", sep='')))
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



##########################################################################################################
#For Associations


Genus_Sparcc_Phase_Early_Top100=data.frame(Node1=net_single_Early[["edgelist1"]][["v1"]],
                                    Node2=net_single_Early[["edgelist1"]][["v2"]],
                                    Association=net_single_Early[["edgelist1"]][["asso"]],
                                    Dissimilarity=net_single_Early[["edgelist1"]][["diss"]],
                                    Adjecency=net_single_Early[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Phase_Early_Top100=Genus_Sparcc_Phase_Early_Top100[order(Genus_Sparcc_Phase_Early_Top100$Association,
                                                        decreasing = T),]

write.csv(Genus_Sparcc_Phase_Early_Top100,"Associations/Genus_Sparcc_Phase_Early_Top100.csv",row.names = FALSE)

###########################
Genus_Sparcc_Phase_Mid_Top100=data.frame(Node1=net_single_Mid[["edgelist1"]][["v1"]],
                                        Node2=net_single_Mid[["edgelist1"]][["v2"]],
                                        Association=net_single_Mid[["edgelist1"]][["asso"]],
                                        Dissimilarity=net_single_Mid[["edgelist1"]][["diss"]],
                                        Adjecency=net_single_Mid[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Phase_Mid_Top100=Genus_Sparcc_Phase_Mid_Top100[order(Genus_Sparcc_Phase_Mid_Top100$Association,
                                                                decreasing = T),]

write.csv(Genus_Sparcc_Phase_Mid_Top100,"Associations/Genus_Sparcc_Phase_Mid_Top100.csv",row.names = FALSE)

###########################
Genus_Sparcc_Phase_Late_Top100=data.frame(Node1=net_single_Late[["edgelist1"]][["v1"]],
                                      Node2=net_single_Late[["edgelist1"]][["v2"]],
                                      Association=net_single_Late[["edgelist1"]][["asso"]],
                                      Dissimilarity=net_single_Late[["edgelist1"]][["diss"]],
                                      Adjecency=net_single_Late[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Phase_Late_Top100=Genus_Sparcc_Phase_Late_Top100[order(Genus_Sparcc_Phase_Late_Top100$Association,
                                                            decreasing = T),]

write.csv(Genus_Sparcc_Phase_Late_Top100,"Associations/Genus_Sparcc_Phase_Late_Top100.csv",row.names = FALSE)

###########################




#Author: M. Anushka Udara
#Date:   19/03/2023
#Title:  Microbial co-occurrence networks for different lactation phases
#Type:   R Script


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


#Creating the initial phyloseq object

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

#Reading taxa file: only important in functional analysis
comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)

#aggregating taxa in the Genus level
complete_agg <- tax_glom(complete_physeq, 'Genus')

complete_agg <- aggregate_taxa(complete_agg, 'Genus')

agg_tax=tax_table(complete_agg)
agg_tax_df=data.frame(rbind(agg_tax))
write.csv(agg_tax_df,"Faprotax/tax/agg_tax.csv")

agg_otu=otu_table(complete_agg)
agg_otu_df=data.frame(rbind(agg_otu))
write.csv(agg_otu_df,"Faprotax/otu/agg_otu.csv")

# Taxa filtering
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

# Saving the phyloseq objects for lactation phases
saveRDS(complete_physeq, "Phase/Complete data.rds")
saveRDS(complete_agg, "Phase/Aggregated data.rds")
saveRDS(filt_complete_agg, "Phase/Filtered Aggregated data.rds")
saveRDS(Early_physeq, "Phase/Early Phase.rds")
saveRDS(Late_physeq, "Phase/Late Phase.rds")
saveRDS(Mid_physeq, "Phase/Mid Phase.rds")
saveRDS(Dry_physeq, "Phase/Dry Phase.rds")


#Early Lactation Network Construction
Early_data <- readRDS("Phase/Early Phase.rds")

top_Early <- prune_taxa(names(sort(taxa_sums(Early_data),TRUE)[1:50]), Early_data)

# Setting up values for network constructing arguments. measure should be "sparcc"
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

#Analyzing the constructed network
props_single_Early <- netAnalyze(net_single_Early, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_Early, "Phase/Early Phase/Early_Sparcc analysis.rds")

#Summary of the analysis
net.summary <- summary(props_single_Early)
net.summary

#Visualizing the network
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

#Extracting edge data, hub data and node data
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

#Extracting the centrality values
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

#Functional analysis with FAPROTAX database
funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Early Phase/Early node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/EarlyCluster',paste("Early_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}

comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)


getwd()
fileList=list.files(path="New folder\\EarlyCluster",pattern = ".csv")
fileList

for(i in fileList){
  cat(i,"*****************************************************************")
  path=normalizePath(file.path('./New folder/EarlyCluster',paste(i, sep='')))
  file=(read.csv(path))
  print(file$x)
  
  for(j in file$x){
    x=subset(comp_tax_csv, comp_tax_csv$genus==j)
    print(x)
    path=normalizePath(file.path('./New folder/EarlyOtu',paste(i,".txt", sep='')))
    write.table(subset(comp_tax_csv, comp_tax_csv$genus==j),path,append=TRUE,col.names = F)
  }
  
}


fileList=list.files(path="New folder\\EarlyOtu",pattern = ".txt")
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


#Mid Lactation Network Construction
Mid_data <- readRDS("Phase/Mid Phase.rds")
top_Mid <- prune_taxa(names(sort(taxa_sums(Mid_data),TRUE)[1:50]), Mid_data)


#SparCC network construction for mid phase
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

#Mid lactation phase network analysis
props_single_Mid <- netAnalyze(net_single_Mid, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector", 
                               hubQuant = 0.95)

saveRDS(props_single_Mid, "Phase/Mid Phase/Mid analysis.rds")

#Network visualization
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

#Saving node data, edge data, and hub data.
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

#Saving the centrality values
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


#Functional analysis
funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Mid Phase/Mid node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/MidCluster',paste("Mid_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)

#Reading the individual genera of cluster files and checking the OTU similarity of samples and FAPROTAX
getwd()
fileList=list.files(path="New folder\\MidCluster",pattern = ".csv")

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

func_df=data.frame(funcFile)

fileList=list.files(path="New folder\\MidOtu",pattern = ".txt")

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


#Late Phase SparCC network construction
Late_data <- readRDS("Phase/Late Phase.rds")
top_Late <- prune_taxa(names(sort(taxa_sums(Late_data),TRUE)[1:50]), Late_data)

#Setting up values for network constructing arguments
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

#network analysis
props_single_Late <- netAnalyze(net_single_Late, 
                                clustMethod = "cluster_fast_greedy",
                                hubPar = "eigenvector", 
                                hubQuant = 0.95)

saveRDS(props_single_Late, "Phase/Late Phase/Late analysis.rds")

#Visualizing the network
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

#Saving node, edge and hub data
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

#saving the centrality values
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

#Functional analysis
funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Phase/Late Phase/Late node data.csv")
clustfile=data.frame(clustfile)


clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/LateCluster',paste("Late_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}


comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)

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

func_df=data.frame(funcFile)

fileList=list.files(path="New folder\\LateOtu",pattern = ".txt")

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




#########################################################################################################################

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

#Initial Phyloseq object
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


# Taxa filtering
minTotRelAbun <- 5e-5
sums <- taxa_sums(complete_agg)
#sums
keepTaxa <- taxa_names(complete_agg)[which((sums / sum(sums)) > minTotRelAbun)]
#keepTaxa
filt_complete_agg <- prune_taxa(keepTaxa, complete_agg)

#Select Updated Cattle Breeds
JRJRC_physeq=subset_samples(filt_complete_agg, CattleBreed =="Jersey" | CattleBreed=="Jersey_cross" | CattleBreed=="Jersey-Sahiwal_cross")
FRFRC_physeq=subset_samples(filt_complete_agg, CattleBreed =="Friesian_cross" | CattleBreed=="Friesian-Sahiwal_cross" | CattleBreed=="Friesian-Jersey_cross" | CattleBreed=="AFS" | CattleBreed=="Friesian")
Asian_physeq=subset_samples(filt_complete_agg, CattleBreed =="Sahiwal" | CattleBreed=="Sahiwal_cross" | CattleBreed =="Local_crossbreds(Batu_cross)")

# Saving the phyloseq objects for updated cattle  breeds
saveRDS(JRJRC_physeq, "Breed/JRJRC.rds")
saveRDS(FRFRC_physeq, "Breed/FRFRC.rds")
saveRDS(Asian_physeq, "Breed/Asian.rds")

comp_tax_csv=read.csv("Faprotax/tax/comp_tax.csv",header = T)
agg_tax_csv=read.csv("Faprotax/tax/agg_tax.csv",header = T)

#SparCC networks for Jersey and Jersey cross breeds
#30 Samples
JRJRC_data <- readRDS("Breed/JRJRC.rds")

top_JRJRC <- prune_taxa(names(sort(taxa_sums(JRJRC_data),TRUE)[1:50]), JRJRC_data)

#setting up sparcc arguments for network construction
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

#network analysis
props_single_JRJRC <- netAnalyze(net_single_JRJRC, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_JRJRC, "Breed/JRJRC Breed/JRJRC_Sparcc analysis.rds")

#summaryb of the analysis
net.summary <- summary(props_single_JRJRC)
net.summary

#network visualization
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

#saving edge data, node data, and hub data
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

#saving centrality values
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

#saving associations of top 100 nodes (both positive and negative)
Genus_Sparcc_Breed_JRJRC_Top100=data.frame(Node1=net_single_JRJRC[["edgelist1"]][["v1"]],
                                           Node2=net_single_JRJRC[["edgelist1"]][["v2"]],
                                           Association=net_single_JRJRC[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_JRJRC[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_JRJRC[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Breed_JRJRC_Top100=Genus_Sparcc_Breed_JRJRC_Top100[order(Genus_Sparcc_Breed_JRJRC_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Sparcc_Breed_JRJRC_Top100,"Associations/Genus_Sparcc_Breed_JRJRC_Top100.csv",row.names = FALSE)

#Functional analysis
funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/JRJRC Breed/JRJRC node data.csv")
clustfile=data.frame(clustfile)

#write node files for each cluster
clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/JRJRCCluster',paste("JRJRC_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}

#loading FAPROTAX results to a dataframe
comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)

#Reading each cluster file seperately
getwd()
fileList=list.files(path="New folder\\JRJRCCluster",pattern = ".csv")

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

#reading OTU files of each cluster
fileList=list.files(path="New folder\\JRJRCOtu",pattern = ".txt")

#getting FAPROTAX functional file as a dataframe
func_df=data.frame(funcFile)

#Looking for similar Genera/OTU in FAPROTAX file and our OTU files
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

#writing FAPROTAX functions for each cluster
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


#Frieasin and Frieasin ceoss breeds network
#26 Samples
FRFRC_data <- readRDS("Breed/FRFRC.rds")

top_FRFRC <- prune_taxa(names(sort(taxa_sums(FRFRC_data),TRUE)[1:50]), FRFRC_data)

#setting up arguments for sparcc frieasin and cross breeds network
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

#network analyzing
props_single_FRFRC <- netAnalyze(net_single_FRFRC, 
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector", 
                                 hubQuant = 0.95)

saveRDS(props_single_FRFRC, "Breed/FRFRC Breed/FRFRC_Sparcc analysis.rds")

#summary of the analysis
net.summary <- summary(props_single_FRFRC)
net.summary

#network visualization
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

#saving node data, edge data, and hub data
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

#saving centrality values
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


#Associations between top 100 nodes (both positive and negative)
Genus_Sparcc_Breed_FRFRC_Top100=data.frame(Node1=net_single_FRFRC[["edgelist1"]][["v1"]],
                                           Node2=net_single_FRFRC[["edgelist1"]][["v2"]],
                                           Association=net_single_FRFRC[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_FRFRC[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_FRFRC[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Breed_FRFRC_Top100=Genus_Sparcc_Breed_FRFRC_Top100[order(Genus_Sparcc_Breed_FRFRC_Top100$Association,
                                                                      decreasing = T),]

write.csv(Genus_Sparcc_Breed_FRFRC_Top100,"Associations/Genus_Sparcc_Breed_FRFRC_Top100.csv",row.names = FALSE)


#Functional analysis
funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/FRFRC Breed/FRFRC node data.csv")
clustfile=data.frame(clustfile)

#saving cluster files
clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/FRFRCCluster',paste("FRFRC_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}

#Creating dataframe for taxa data
comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)

#reading cluster file names
getwd()
fileList=list.files(path="New folder\\FRFRCCluster",pattern = ".csv")

#read cluster files
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


#reading OTU files
fileList=list.files(path="New folder\\FRFRCOtu",pattern = ".txt")

#reading FAPROTAX functional file
func_df=data.frame(funcFile)

#Comparison for genera and OTUs
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

#writing functions to file system
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


#Asian breeds network
#8 Samples
Asian_data <- readRDS("Breed/Asian.rds")
top_Asian <- prune_taxa(names(sort(taxa_sums(Asian_data),TRUE)[1:50]), Asian_data)


#Setting up arguments for network construction
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

#network analysis
props_single_Asian <- netAnalyze(net_single_Asian,
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = "eigenvector",
                                 hubQuant = 0.95)
saveRDS(props_single_Asian, "Breed/Asian Breed/Asian_Sparcc analysis.rds")


#summary of analysis
net.summary <- summary(props_single_Asian)
net.summary

#network visualizaton
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

#savinf edge data, node data, and hub data
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

#saving centrality values
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

#saving associations of top 100 nodes (both positive and negative)
Genus_Sparcc_Breed_Asian_Top100=data.frame(Node1=net_single_Asian[["edgelist1"]][["v1"]],
                                           Node2=net_single_Asian[["edgelist1"]][["v2"]],
                                           Association=net_single_Asian[["edgelist1"]][["asso"]],
                                           Dissimilarity=net_single_Asian[["edgelist1"]][["diss"]],
                                           Adjecency=net_single_Asian[["edgelist1"]][["adja"]]
)

Genus_Sparcc_Breed_Asian_Top100=Genus_Sparcc_Breed_Asian_Top100[order(Genus_Sparcc_Breed_Asian_Top100$Association,
                                                                      decreasing = T),]
write.csv(Genus_Sparcc_Breed_Asian_Top100,"Associations/Genus_Sparcc_Breed_Asian_Top100.csv",row.names = FALSE)


#Functional analysis
funcFile=read.csv("Faprotax/tax/1670601055.txt.faprotax.report.clean.otu_func.csv",header = F)
clustfile=read.csv("Breed/Asian Breed/Asian node data.csv")
clustfile=data.frame(clustfile)

#saving cluster files
clustfile=clustfile[order(clustfile$cluster,decreasing = F),]
for(i in clustfile$cluster){
  p=(clustfile$lable[clustfile$cluster==i])
  path=normalizePath(file.path('./New folder/AsianCluster',paste("Asian_Cluster",i,'.csv', sep='')))
  write.csv(p,path)
  
}

#reading taxa file
comp_tax_csv=data.frame(otuid=comp_tax_csv$otuid,genus=comp_tax_csv$Genus)

#listing cluster files
getwd()
fileList=list.files(path="New folder\\AsianCluster",pattern = ".csv")

#reading cluster files and taxa files to extract important data
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

#listing OTU files
fileList=list.files(path="New folder\\AsianOtu",pattern = ".txt")

#FAPROTAX functions into dataframe
func_df=data.frame(funcFile)
func_df[1]

#comparing OTU/genera with FAPROTAX
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

#saving up functions
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


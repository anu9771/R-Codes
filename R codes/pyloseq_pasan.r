### phyloseq on Saumya's data ###

## method via: https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

## setwd ##
# Linux
#setwd("/home/pbiggs/extraDrive2/data/metabarcoding/MGS00049_Saumya_Wickramasinghe_Delivery/MGS00049_1/results/")
# Windows on Dropbox
setwd("C:\\Users\\pbiggs\\Dropbox\\MasseyWork\\research_other\\nonMassey\\SaumyaW\\MGS00049MilkResults\\results\\")

## packages ##
library(tidyverse)
library("scales")
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(microbiome)

library(tidyverse)
library(phyloseq)
library(qiime2R)
library(BiocManager)
library(microbiome)
library(NetCoMi)
library(WGCNA)

setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology")

## import the data as per the method ##
# metadata
metadata<-read_tsv("datafiles/metadata_table.txt")
metadata
?read_tsv

# SVs
SVs<-read_qza("datafiles/table.qza")
#names(SVs)
#SVs$uuid
## [1] "3f97241d-4d46-4f6e-9b1d-8517a1235bb6"
#SVs$data[1:5,1:5] #show the first 5 samples and features

# taxonomy
taxonomy<-read_qza("datafiles/taxonomy.qza")
#taxonomy$uuid
##[1] "89a39f71-1344-448e-b779-d62512f81419"
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) 
#convert the table into a tabular split version
head(taxtable)
taxtable

# tree
tree<-read_qza("datafiles/rooted-tree.qza")
tree$uuid
##[1] "94a66a5d-5e49-4513-b8fb-3f6113bac6f2"
tree$data

# shannon
shannon<-read_qza("core-metrics-results-test/shannon_vector.qza")
shannon$uuid
##[1] "5b1c88b0-b409-473f-a6aa-c30f28a03d1f"
head(shannon$data)

# PCoA
pco<-read_qza("core-metrics-results-test/unweighted_unifrac_pcoa_results.qza")
pco$uuid
##[1] "073b4b56-eaa1-48d6-9a0b-f9879a9c922a"
head(pco$data$ProportionExplained) #this returns the variance explained
pco$data$Vectors[1:5, 1:3]


## make an object ##
phy<-qza_to_phyloseq(features="datafiles/table.qza", tree="datafiles/rooted-tree.qza", taxonomy="datafiles/taxonomy.qza")
#phy<-qza_to_phyloseq("table.qza", "rooted-tree.qza", "taxonomy.qza", "../setup/metafile82_R.txt")
phy
##phyloseq-class experiment-level object
##otu_table()   OTU Table:         [ 10918 taxa and 81 samples ]
##sample_data() Sample Data:       [ 81 samples by 16 sample variables ]
##tax_table()   Taxonomy Table:    [ 10918 taxa by 7 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 10918 tips and 10795 internal nodes ]

#Reading the meta data
metx=read.delim("datafiles/metadata_table.txt",sep=",",header=T,row.names=sample_names(phy))
head(metx)

#Meta data as sample data
metx=metx[,-1]
metx=sample_data(metx)
metx

#Merging the sample
phy=merge_phyloseq(phy,metx)
phy

##phyloseq-class experiment-level object
##otu_table()   OTU Table:         [ 10918 taxa and 81 samples ]
##sample_data() Sample Data:       [ 81 samples by 16 sample variables ]
##tax_table()   Taxonomy Table:    [ 10918 taxa by 7 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 10918 tips and 10795 internal nodes ]


head(taxa_names(phy))


## try alternative way to export and reimport data to sort out this 

write_phyloseq(phy, 'OTU', path = getwd())
write_phyloseq(phy, 'TAXONOMY', path = getwd())
write_phyloseq(phy, 'METADATA', path = getwd())

# Saving the phyloseq object as a RDATA file
save(phy, file = "datafiles/phy.RData")

# Opening the interactive Shiny program
shiny::runGitHub("shiny-phyloseq","joey711")

## hand over to phyloseq
library(phyloseq)
packageVersion("phyloseq")

# some basic playing #

sample_variables(phy)
head(taxa_names(phy))


top20 <- names(sort(taxa_sums(phy), decreasing=TRUE))[1:100]
sumtable=taxa_sums(phy)
head(sumtable)
head(top20)
ps.top20 <- transform_sample_counts(phy, function(OTU) OTU/sum(OTU))

# my modification to get percentages; use the above line for original implementation
#ps.top20 <- transform_sample_counts(phy, function(OTU) (OTU/sum(OTU))*100 )
# saving the otu table with percentages
write_phyloseq(ps.top20, 'OTU', path = getwd())

ps.top20 <- prune_taxa(top20, phy)
plot_bar(phy, x="AgriculturalZone", fill="Class") 
facet_wrap(~pH, scales="free_x")

plot_bar(ps.top20, x="CattleBreed", fill="Order") + facet_wrap(~AgriculturalZone, scales="free_x")

?facet_wrap

ord.nmds.bray <- ordinate(phy, method="NMDS", distance="bray")
plot_ordination(phy, ord.nmds.bray, color="pH", title="Bray NMDS")

barplot(sort(taxa_sums(phy), TRUE)[1:100]/nsamples(phy), las=2)

barplot(sort(taxa_sums(phy), TRUE)[1:100]/nsamples(phy), las=2)
plot_tree(ps.top20,shape="AgriculturalZone", color="CattleBreed", label.tips="Class")

barplot(sort(taxa_sums(phy), TRUE)[1:100]/nsamples(phy), las=2)
plot_tree(ps.top20,shape="AgriculturalZone", color="CattleBreed", label.tips="Order")

barplot(sort(taxa_sums(phy), TRUE)[1:100]/nsamples(phy), las=2)
plot_tree(ps.top20,shape="AgriculturalZone", color="CattleBreed", label.tips="Family")

barplot(sort(taxa_sums(phy), TRUE)[1:100]/nsamples(phy), las=2)
plot_tree(ps.top20,shape="AgriculturalZone", color="CattleBreed", label.tips="Genus")

barplot(sort(taxa_sums(phy), TRUE)[1:100]/nsamples(phy), las=2)
plot_tree(ps.top20, color="CattleBreed", label.tips="Species")
#############################
#
# now to try the microbiome package #
#
#############################

#library(BiocManager)
#BiocManager::install(version='devel')
#BiocManager::install("microbiome")


library(knitr)
library(dplyr)
library(ggpubr)

summarize_phyloseq(phy)


## global measures ##
# from https://microbiome.github.io/tutorials/Diversity.html

overall <- global(phy)

rich <- richness(phy)
kable(head(rich))

dominant <- dominance(phy, index = "all")
kable(head(dominant))

rare <- rarity(phy, index = "all")
kable(head(rare))

covMilk <- coverage(phy, threshold = 0.5)
kable(head(covMilk))

coreMilk <- core_abundance(phy, detection = .1/100, prevalence = 50/100)
kable(head(coreMilk))


## plot alpha diversity ##
# from: https://microbiome.github.io/tutorials/PlotDiversity.html

#ps1 <- prune_taxa(taxa_sums(phy) > 0, phy)
#tab <- alpha(ps1, index = "all")
#kable(head(tab))


## composition ##
# from: https://microbiome.github.io/tutorials/Composition.html

transform <- microbiome::transform

coreMilk <- core(phy, detection = .01/100, prevalence = 5/100)
coreMilk

p <- plot_composition(microbiome::transform(coreMilk, "compositional"),
                      plot.type = "heatmap",
                      sample.sort = "neatmap",
                      otu.sort = "neatmap")
print(p)
plot_taxa_prevalence(phy, "Phylum", detection = 10)


## core microbiota ##
# from: https://microbiome.github.io/tutorials/Core.html

phy.rel <- microbiome::transform(phy, "compositional")
head(prevalence(phy.rel, detection = 1/100, sort = TRUE))
head(prevalence(phy.rel, detection = 1/100, sort = TRUE, count = TRUE))

core.taxa.standard <- core_members(phy.rel, detection = 0, prevalence = 5/100)
phy.core <- core(phy.rel, detection = 0, prevalence = .05)
core.taxa <- taxa(phy.core)
core.abundance <- sample_sums(core(phy.rel, detection = 0.1/100, prevalence = 50/100))

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(phy.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

#### Pasan's codes ###########################################################################################

# plotting the abundance in a bar plot
plot_bar(phy, fill = "Family")

# plotting the heatmap
rank_names(phy)
top20 <- names(sort(taxa_sums(phy), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, phy)
plot_heatmap(phy)
             
plot_heatmap(ps.top20,taxa.label="Class", sample.label="AgriculturalZone", sample.order="AgriculturalZone")
head(sample_data(phy))



##################################################################################################################
# Tabulating the prevalance
prevelancedf = apply(X = otu_table(phy), MARGIN = 1,FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(phy),
                          tax_table(phy))
prevelancedf[1:10,]
write.csv(prevelancedf,'prevelancedf.csv')

prphyla=plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
write.csv(prphyla,'prevelancephylums.csv')

# filtering the taxa based on prevalance
phyla2Filter = c("p__Synergistetes", "p__OP8", "p__TM6",
                 "p__SR1","p__GN02","p__Fibrobacteres",
                 "p__Nitrospirae","p__WPS-2","p__Planctomycetes")
# Filter entries with unidentified Phylum.
phyfiltered = subset_taxa(phy, !Phylum %in% phyla2Filter)
phyfiltered

selectvs =c("CattleBreed","cleanliness","AgriculturalZone","LactationPhase","pooled","FarmSize","FarmCode"  )

# transforming to obtain relative abundances
trans <- transform_sample_counts(phyfiltered, function(x) x / sum(x))


###############################################################################################
# Calculating the alpha diversity

var1="CattleBreed"
p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
p <- p + geom_boxplot( alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
p <- p + geom_violin( alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
p <- p + geom_boxplot( notch=TRUE, alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
plot(p)

var1="LactationPhase"
p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
p <- p + geom_boxplot( alpha=0.2, aes(fill =LactationPhase, group=LactationPhase))
p <- p + geom_violin( alpha=0.2, aes(fill =LactationPhase, group=LactationPhase))
p <- p + geom_boxplot( notch=TRUE, alpha=0.2, aes(fill =LactationPhase, group=LactationPhase))
plot(p)

# a loop to go through all the variables
selectvs =c("CattleBreed","cleanliness","AgriculturalZone","LactationPhase","pooled","FarmSize","FarmCode","pH"  )
for (var1 in selectvs){
  #File name
  trfname= paste( "AD", var1,"notched_boxplots.pdf", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
  #p <- p + geom_boxplot( alpha=0.2) # converting string to variable name
  p <- p + geom_boxplot( notch= TRUE, alpha=0.2) # notched boxplot
  #p <- p + geom_violin( alpha=0.2) # violin plot
  
  print(plot(p))
  # Saving the pdf
  dev.off()
}

#generating richness values to a table
richtable =estimate_richness(phyfiltered)
# saving the table
write.csv(richtable,'ad_estimates.csv')


# list to combine wilcox test outputs
pwclist = list()

for (var1 in selectvs){
  pwc= pairwise.wilcox.test(richtable$Chao1, sample_data(phyfiltered)[[var1]],p.adj = "bonf")
  pwclist= c(pwclist,pwc)
}
# capturing the combined outputs and saving in a text file
chars= capture.output(print(pwclist))
writeLines(chars, con = file("_filtered_bonfadjustedpvalue.txt"))
# closing the opened file connection
close(file("_filtered_bonfadjustedpvalue.txt"))

# If you want to only save the p values as a table use the following code
#write.csv(pwc[["p.value"]],'pwc.csv')

###############################################################################################
###############################################################################################################
#### Beta diversity measures ####################################

# plotting the network ######
plot_net(phy, maxdist = 0.4)
ig <- make_network(phy, max.dist=0.8)
plot_network(ig, phy)
plot_network(ig, phy, color="FarmSize", shape="AgriculturalZone", line_weight=0.4)

igbray <- make_network(phyfiltered, dist.fun="bray", max.dist=0.9)
plot_network(igbray, phyfiltered, color="CattleBreed",  line_weight=0.4, label = "CattleBreed")
#plot_network(igbray, phy, color="SampleWell", line_weight=0.4)

for (i in selectvs){
  #File name
  trfname= paste(i,"network.pdf", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_network(igbray, phyfiltered, color=i, line_weight=0.4, label = "FarmCode"))
  # Saving the pdf
  dev.off()
}

# Ordination ######
# Calculate distances
DistBC = distance(trans, method = "bray")
DistUF = distance(trans, method = "wUniFrac")

DistBC = distance(phyfiltered, method = "bray")
DistUF = distance(phyfiltered, method = "wUniFrac")

# Calculating the ordinations
ordBC = ordinate(trans, method = "PCoA", distance = DistBC)
ordUF = ordinate(trans, method = "PCoA", distance = DistUF)

ordBC = ordinate(phyfiltered, method = "PCoA", distance = DistBC)
ordUF = ordinate(phyfiltered, method = "PCoA", distance = DistUF)

# Plotting scree plots
plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_scree(ordUF, "Scree Plot: Weighted UniFrac MDS")

# Plotting PCOAs
plot_ordination(trans, ordBC, color="FarmCode") +ggtitle("PCoA: Bray-Curtis")
plot_ordination(phyfiltered, ordUF, color="FarmCode", label = "FarmCode") +ggtitle("PCoA: Weigthed Unifrac")

# For bray-curtis
for (i in selectvs){
  #File name
  trfname= paste("Bray-curtis_ordinationplot",i,".pdf", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_ordination(trans, ordBC, color=i) +ggtitle("PCoA: Bray-Curtis"))
  # Saving the pdf
  dev.off()
}

# For unifrac
for (i in selectvs){
  #File name
  trfname= paste("Unifrac_ordinationplot",i,".pdf", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_ordination(trans, ordUF, color=i, label=i) +ggtitle("PCoA: Weigthed Unifrac"))
  # Saving the pdf
  dev.off()
}



###############################################################################################################
# one barplot to represent all the sample counts
#physamplesum = sample_sums(phy)
#physamplesum
#plot_bar(ps.top20, fill="Genus") 
# getting the percentages of total counts
#mphysum = mphysum/sum(mphysum)
# otable1 =read_csv('otu_table_percentages.csv', na = c("", "NA"))
# otable2 = as.matrix(otable1)

# getting taxa sums and converting to a matrix
physum = taxa_sums(phy)
physum = taxa_sums(phyfiltered)
mphysum = as.matrix(physum)
class(mphysum)
# adding the column names
colnames(mphysum) = c('X1')

# reading the taxonomy table
getwd()
ttable1 =read.csv('datafiles/taxonomy_table.csv', na = c("", "NA"))
class(ttable1)
ttable1=column_to_rownames(ttable1,'X1')
class(ttable1)
ttable2 = as.matrix(ttable1)
OTU1 = otu_table(mphysum, taxa_are_rows = TRUE)
TAX1 = tax_table(ttable2)

totsample1 = phyloseq(OTU1,TAX1)
totsample=merge_phyloseq( mphysum,tax_table(phy))

# generating the basic plot
top100 <- names(sort(taxa_sums(totsample1), decreasing=TRUE))[1:100]
ps.top100 <- prune_taxa(top100, totsample1)
plot_bar(ps.top100, fill="Order") 

# converting to RA values
transtot <- transform_sample_counts(totsample1, function(x) x / sum(x))
ps.top100 <- prune_taxa(top100, transtot)
plot_bar(ps.top100, fill="Order")

# defining a taxa level
t_level = "Genus"

#aggregating taxa
totGlommed = tax_glom(transtot, t_level, NArm=TRUE)

top100 <- names(sort(taxa_sums(totGlommed), decreasing=TRUE))[1:20]
ps.top100 <- prune_taxa(top100, totGlommed)
# filtering based on the percentage
cbrfr = filter_taxa(totGlommed, function(x)  max(x) > 0.005, TRUE)
# use the following plot margin for fat bars
plot_bar(cbrfr, fill=t_level)+theme(  plot.margin = margin(2, 9, 2, 9, "cm"))+ ylab("Relative abundance")
plot_bar(cbrfr, fill=t_level)+ylab("Relative abundance")
plot_bar(ps.top100, fill=t_level)+ylab("Relative abundance")

totmelted <- psmelt(totGlommed)

# concatenating string for csv file name
csvfname= paste(t_level, "_tot_RA.csv", sep="_")
write.csv(totmelted,csvfname)
############################################################################################
# Generating relative abundance bar plots


# aggregating taxa
aggrphy <- aggregate_taxa(trans, 'Phylum')
melted <- psmelt(toptax)
# Saving the otu table of the aggregated phyloseq object
#write_phyloseq(aggrphy, 'OTU', path = getwd())
write.csv(melted,'GEnus_percentage.csv')


# defining a taxa level
t_levellist = c("Species","Order","Phylum")

# The list of variables
selectvs =c("CattleBreed","cleanliness","AgriculturalZone","LactationPhase","pooled","FarmSize","FarmCode"  )

# # Another way to aggregate data
# phylumGlommed = tax_glom(trans, t_level, NArm=FALSE)
# 
# # getting the top abundant taxa
# topphy = names(sort(taxa_sums(phylumGlommed), TRUE)[1:30])
# toptax=prune_taxa(topphy,phylumGlommed)


# cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.01, TRUE)
plot_tree(toptax, color="CattleBreed", label.tips= "taxa_names", ladderize = TRUE, justify = "huha",nodelabf = nodeplotboot())

for(t_level in t_levellist){
  
  # Another way to aggregate data
  phylumGlommed = tax_glom(trans, t_level, NArm=TRUE)
  
  # getting the top abundant taxa
  topphy = names(sort(taxa_sums(phylumGlommed), TRUE)[1:30])
  toptax=prune_taxa(topphy,phylumGlommed)
  
  # Plotting the phylogenetic trees for all the variables using a loop
  for (i in selectvs){
    #File name
    trfname= paste( "Top30_justified",t_level,"tree", i,".pdf", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    #plotting
    print(plot_tree(toptax, color=i, label.tips= t_level,ladderize = "left", justify = "yes"))
    # Saving the pdf
    dev.off()
  }
  
  
}



# Saving the melted data file
phylummelted <- psmelt(phylumGlommed)

write.csv(phylummelted,'Order_glom_percentage.csv')


### Phylogenetic trees for different variable aggregates ###############
# different taxa ranks to iterate
taxal = c("Genus")
# Looping through the selected sample variables
for (i in selectvs){
  # merge based on variables
  mergedphy =merge_samples(phyfiltered,group = i)
  
  # transform to relative abundance
  cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
  # Aggregation based on given taxa level
  for(tl in t_levellist){
    phylumGlommed = tax_glom(cbreedphy, tl, NArm=FALSE)
    
    # filtering taxa based on percentage
    #cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.02, TRUE)
    
    topphy = names(sort(taxa_sums(phylumGlommed), TRUE)[1:30])
    toptax=prune_taxa(topphy,phylumGlommed)
    
    #File name
    trfname= paste( "top30",tl,"tree", i,".pdf", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    #plotting
    print(plot_tree(toptax, color="Sample", size="abundance",label.tips= tl,base.spacing=0.04, ladderize = TRUE))
    
    # Saving the pdf
    dev.off()
  }
}


# Assigning permenant set of colors for phylums; This can be used for any occassion that uses these data
library("RColorBrewer")
getPalette = colorRampPalette(brewer.pal(9, "Paired"))
familyList = unique(tax_table(trans)[,t_level])
familyPalette = getPalette(length(familyList))
names(familyPalette) = familyList


p = plot_bar(toptax,fill = t_level)+ facet_wrap(~cleanliness, scales="free_x") + ylab("Relative abundance")
p= p + scale_fill_manual(values= familyPalette)+theme(legend.position="bottom", legend.box.background = element_rect(colour = "black") )
plot(p)



###############################################################################################
#selective abundance
gp.ch = subset_taxa(phylumGlommed, Order == "o__Clostridiales")
mdf = psmelt(gp.ch)
nrow(mdf)
ncol(mdf)
colnames(mdf)
head(rownames(mdf))

p = ggplot(gp.ch, aes(x=cleanliness, y=Abundance, fill=Species))
#p = p + geom_bar(color="black", stat="identity", position="stack")
p = p + geom_bar(stat="identity", position="stack")
print(p)

write.csv(mdf,'Clostridiales.csv')

p = ggplot(mdf, aes(x=Sample, y=Abundance, fill=Order))
#p = p + geom_bar(color="black", stat="identity", position="stack")
p = p + geom_bar(stat="identity", position="stack")
print(p)

##################################################################################################################################
# Variable aggregation
# different taxa ranks to iterate
taxal = c("Genus")
# Looping through the selected sample variables
for (i in selectvs){
  # merge based on variables
  mergedphy =merge_samples(phyfiltered,group = i)
  
  # transform to relative abundance
  cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
  # Aggregation based on given taxa level
  for(tl in taxal){
    phylumGlommed = tax_glom(cbreedphy, tl, NArm=FALSE)
    
    # merging the OTU table with TAXA names and saving as a CSV
    # converting each table to dataframes
    ttdf=data.frame(as(tax_table(phylumGlommed), "matrix"))
    # have to transpose the OTU table as merging automatically transposes the matrix
    otdf=as.data.frame(t(otu_table(phylumGlommed)))
    # merging
    mtable= merge.data.frame(ttdf,otdf, by=0)
    # saving as a CSV
    csname= paste( tl, i,"relative_abandance.csv", sep="_")
    write.csv(mtable,csname)
   
    # filtering taxa based on percentage
    cbrfr = filter_taxa(phylumGlommed, function(x)  max(x) > 0.01, TRUE)
    
    #Plotting the stacked barplots
    #File name
    trfname= paste( tl, i,"barplot.pdf", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    # Plotting the stacked barplots
    p = plot_bar(cbrfr,fill = tl) + ylab("Relative abundance")
    p = p + geom_bar(color="black",stat="identity", position="stack")
    plot(p)
    # Saving the pdf
    dev.off()
    
    # Heatmaps
    #File name
    trfname= paste( tl, i,"heatmap.pdf", sep="_")
    # Command to save as a pdf
    pdf(file=trfname, paper="a4r",width = 0, height = 0)
    hm=plot_heatmap(cbrfr,method=NULL, taxa.label = tl,taxa.order=tl)+ scale_fill_gradient(name ="Relative abundance")
    plot(hm)
    # Saving the pdf
    dev.off()
  }
  
  
}


##################################################################################################################################

devtools::install_github("vmikk/metagMisc")
library(metagMisc)
phydf =phyloseq_to_df(phy, addtax = T, addtot = T, addmaxrank = F,sorting = "abundance")



##################################################################################################################################
# deseq differential abundance analysis
#First load DESeq2.

library("DESeq2")
packageVersion("DESeq2")

#converts your phyloseq-format microbiome data into a DESeqDataSet 
diagdds = phyloseq_to_deseq2(phyfiltered, ~ FarmCode)

# Workaround for "every gene contains zero" error
# defining a function to calculate geo mean
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}

geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

# Summarizing the results of the test
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

write.csv(sigtab,'deseq_significance.csv')
# nice ggplot2 summary of the results
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Species order
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))

ggplot(sigtab, aes(x=Species, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

resultsNames(diagdds)

devtools::install_github("umerijaz/microbiomeSeq")  # Install the package
library(microbiomeSeq)  #load the package

p <- boxplot_abundance(trans, x='CattleBreed', y='33cffd2c9eba73e4feee5abe2509b790')
plot(p)

####################################################################################
# Analyzing taxa subset

# extracting the selected phyla
sp2extract = c("s__agalactiae", "s__saprophyticus", "s__aerosaccus", "s__nasimurium")

spfiltered = subset_taxa(phy, Species %in% sp2extract)
#gp.ch = subset_taxa(phy, Species == "s__agalactiae")
for (i in selectvs){
  
  #File name
  plotname= paste( i,"pathogen_barplot.pdf", sep="_")
  # Command to save as a pdf
  pdf(file=plotname, paper="a4r",width = 0, height = 0)
  # Plotting the stacked barplots
  p= plot_bar(spfiltered, x=i, fill="Species")
  p =p+ geom_bar( stat="identity", position="stack")
  plot(p)
  # Saving the pdf
  dev.off()
}
p= plot_bar(spfiltered, x="CattleBreed", fill="Species")
p + geom_bar( stat="identity", position="stack")


mergedphy =merge_samples(phyfiltered,group = 'FarmCode')

# transform to relative abundance
cbreedphy <- transform_sample_counts(mergedphy, function(x) x / sum(x))
# Aggregation based on given taxa level

phylumGlommed = tax_glom(cbreedphy, 'Species', NArm=FALSE)


spfiltered = subset_taxa(phylumGlommed, Species %in% sp2extract)

# Plotting the stacked barplots
p = plot_bar(spfiltered, fill = 'Species') + ylab("Relative abundance")+ xlab('Farm Code')
p = p + geom_bar(color="black",stat="identity", position="stack")
plot(p)

#####################################################################################################
# Aggregated variable distributions

# Aggregating taxa based on taxonomic level
TaxGlommed = tax_glom(trans, "Species", NArm=FALSE)

# getting the top abundant taxa
topphy = names(sort(taxa_sums(TaxGlommed), TRUE)[1:100])
toptax=prune_taxa(topphy,TaxGlommed)

melted <- psmelt(toptax)
write.csv(melted,'Species_melted.csv')


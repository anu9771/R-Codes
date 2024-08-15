################################################################################
#Trying to choose only wanted variables
################################################################################

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
  features = "datafiles2/table.qza",
  taxonomy="datafiles2/taxonomy.qza",
  tree="datafiles2/rooted-tree.qza"
)

#Reading the meta data
metx=read.delim("datafiles2/Editedmetadata_table.csv",sep=",",header=T,row.names=sample_names(physeq))
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


ps_drop_incomplete <- function(ps, vars = NA, verbose = FALSE) {
  df <- sample_data(ps)
  if (identical(vars, NA)) vars <- phyloseq::sample_variables(ps)
  df_sub <- df[, vars, drop = FALSE]
  df_sub <- df_sub[stats::complete.cases(df_sub), , drop = FALSE]
  if (!isFALSE(verbose)) {
    incomplete <- nrow(df) - nrow(df_sub)
    if (incomplete > 0 || identical(verbose, "max")) {
      message("Dropping samples with missings: ", incomplete)
    }
  }
  if (identical(verbose, "max")) {
    for (v in vars) {
      n_missings <- sum(is.na(df[[v]]))
      if (n_missings > 0) message(v, " has NAs: ", n_missings)
    }
  }
  keepers <- rownames(df_sub)
  phyloseq::prune_samples(samples = keepers, x = ps)
}

complete_physeq=ps_drop_incomplete(complete_physeq)


# plotting the heatmap
rank_names(complete_physeq)
top20 <- names(sort(taxa_sums(complete_physeq), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, complete_physeq)
plot_heatmap(complete_physeq)

plot_heatmap(ps.top20,taxa.label="Genus", sample.label="LactationPhase", sample.order="LactationPhase")
head(sample_data(complete_physeq))


##################################################################################################################
# Tabulating the prevalance
prevelancedf = apply(X = otu_table(complete_physeq), MARGIN = 1,FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(complete_physeq),
                          tax_table(complete_physeq))
prevelancedf[1:10,]
write.csv(prevelancedf,'impData/prevelancedf.csv')

prphyla=plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
write.csv(prphyla,'impData/prevelancephylums.csv')

# filtering the taxa based on prevalance
phyla2Filter = c("Synergistetes", "OP8", "TM6",
                 "SR1","GN02","Fibrobacteres",
                 "Nitrospirae","WPS-2","Planctomycetes")
# Filter entries with unidentified Phylum.
phyfiltered = subset_taxa(complete_physeq, !Phylum %in% phyla2Filter)
phyfiltered


selectvs =c("CattleBreed")

# transforming to obtain relative abundances
trans <- transform_sample_counts(phyfiltered, function(x) x / sum(x))


###############################################################################################
# Calculating the alpha diversity
var1="CattleBreed"
p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
#p <- p + geom_boxplot( alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
p <- p + geom_violin( alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
#p <- p + geom_boxplot( notch=TRUE, alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
plot(p)

var1="LactationPhase"
p =plot_richness(phyfiltered, x=var1, measures=c("Shannon", "Simpson","Chao1","Observed"), color = var1)
#p <- p + geom_boxplot( alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
p <- p + geom_violin( alpha=0.2, aes(fill =LactationPhase, group=LactationPhase))
#p <- p + geom_boxplot( notch=TRUE, alpha=0.2, aes(fill =CattleBreed, group=CattleBreed))
plot(p)


# a loop to go through all the variables
selectvs =c("CattleBreed","LactationPhase")
for (var1 in selectvs){
  setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology/impData")
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

setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology")


#generating richness values to a table
richtable =estimate_richness(phyfiltered)
# saving the table
write.csv(richtable,'impData/ad_estimates.csv')


# list to combine wilcox test outputs
pwclist = list()

for (var1 in selectvs){
  pwc= pairwise.wilcox.test(richtable$Chao1, sample_data(phyfiltered)[[var1]],p.adj = "bonf")
  pwclist= c(pwclist,pwc)
}
# capturing the combined outputs and saving in a text file
chars= capture.output(print(pwclist))
writeLines(chars, con = file("impData/_filtered_bonfadjustedpvalue.txt"))
# closing the opened file connection
close(file("impData/_filtered_bonfadjustedpvalue.txt"))


###############################################################################################
###############################################################################################################
#### Beta diversity measures ####################################

# plotting the network ######
plot_net(complete_physeq, maxdist = 0.4)
ig <- make_network(complete_physeq, max.dist=0.8)
plot_network(ig, complete_physeq)
plot_network(ig, complete_physeq, color="LactationPhase", shape="LactationPhase", line_weight=0.4)

igbray <- make_network(phyfiltered, dist.fun="bray", max.dist=0.9)
plot_network(igbray, phyfiltered, color="CattleBreed",  line_weight=0.4, label = "CattleBreed")
plot_network(igbray, phyfiltered, color="LactationPhase",  line_weight=0.4, label = "LactationPhase")
#plot_network(igbray, phy, color="SampleWell", line_weight=0.4)

for (i in selectvs){
  setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology/impData")
  #File name
  trfname= paste("Edited",i,"network.pdf", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_network(igbray, phyfiltered, color=i, line_weight=0.4, label = i))
  # Saving the pdf
  dev.off()
}

setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology")

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
plot_ordination(trans, ordBC, color="CattleBreed") +ggtitle("PCoA: Bray-Curtis")
plot_ordination(phyfiltered, ordUF, color="CattleBreed", label = "CattleBreed") +ggtitle("PCoA: Weigthed Unifrac")

plot_ordination(trans, ordBC, color="LactationPhase") +ggtitle("PCoA: Bray-Curtis")
plot_ordination(phyfiltered, ordUF, color="LactationPhase", label = "LactationPhase") +ggtitle("PCoA: Weigthed Unifrac")

# For bray-curtis
for (i in selectvs){
  setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology/impData")
  #File name
  trfname= paste("Edited_Bray-curtis_ordinationplot",i,".pdf", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_ordination(trans, ordBC, color=i) +ggtitle("PCoA: Bray-Curtis"))
  # Saving the pdf
  dev.off()
}

setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology")

# For unifrac
for (i in selectvs){
  setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology/impData")
  #File name
  trfname= paste("Edited_Unifrac_ordinationplot",i,".pdf", sep="_")
  # Command to save as a pdf
  pdf(file=trfname, paper="a4r",width = 0, height = 0)
  #plotting
  print(plot_ordination(trans, ordUF, color=i, label=i) +ggtitle("PCoA: Weigthed Unifrac"))
  # Saving the pdf
  dev.off()
}

setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology")


######################################################################
#Rest of the graphs
######################################################################

# getting taxa sums and converting to a matrix
physum = taxa_sums(complete_physeq)
physum = taxa_sums(phyfiltered)
mphysum = as.matrix(physum)
class(mphysum)
# adding the column names
colnames(mphysum) = c('X1')


# reading the taxonomy table
getwd()
ttable1 =read.csv('EditedData_alpha/EditedTaxa.csv', na = c("", "NA"))
class(ttable1)
ttable1=column_to_rownames(ttable1,'OTU')
class(ttable1)
ttable2 = as.matrix(ttable1)
OTU1 = otu_table(mphysum, taxa_are_rows = TRUE)
TAX1 = tax_table(ttable2)

totsample1 = phyloseq(OTU1,TAX1)
totsample=merge_phyloseq( mphysum,tax_table(complete_physeq))


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
setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology/impData")
csvfname= paste(t_level, "_tot_RA.csv", sep="_")
write.csv(totmelted,csvfname)


setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology")

# Generating relative abundance bar plots



##################################################
######### Check This_Might Be Incorrect ##########
##################################################


phylumGlommed = tax_glom(trans, t_level, NArm=FALSE)

# getting the top abundant taxa
topphy = names(sort(taxa_sums(phylumGlommed), TRUE)[1:30])
toptax=prune_taxa(topphy,phylumGlommed)


# aggregating taxa
aggrphy <- aggregate_taxa(trans, 'Phylum')
melted <- psmelt(toptax)
# Saving the otu table of the aggregated phyloseq object
#write_phyloseq(aggrphy, 'OTU', path = getwd())
setwd("C:/Users/Anushka/Desktop/Anushka/Research/Methodology/impData")
write.csv(melted,'GEnus_percentage.csv')

# defining a taxa level
t_levellist = c("Species","Order","Phylum")

# The list of variables
selectvs =c("CattleBreed","LactationPhase")

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

###############################################################
################# An Error Occurred ###########################
###############################################################


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


#########################################################################

##################################################################################################################################
setwd("C:\\Users\\Anushka\\Desktop")
#Variable aggregation
# different taxa ranks to iterate
taxal = c("Species")

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
    mtable= merge.data.frame(ttdf,otdf, by=0,all.x = F)
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
    p = plot_bar(subset_taxa(cbrfr,!is.na(Species)),fill = tl) + ylab("Relative abundance")
    p = p + geom_bar(color="black",stat="identity", position="stack")
    p=p+scale_fill_manual(values=c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
                                   "#9467bd","#8c564b","#e377c2","#7f7f7f",
                                   "#bcbd22","#17becf","#aec7e8","#ffbb78",
                                   "#98df8a","#ff9896","#c5b0d5","#c49c94"))
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

#For Visualization Purposes
p = plot_bar(subset_taxa(cbrfr,!is.na(Species)),fill = tl) + ylab("Relative abundance")
p = p + geom_bar(color="black",stat="identity", position="stack")
p=p+scale_fill_manual(values=c("#1f77b4","#ff7f0e","#2ca02c","#d62728",
                               "#9467bd","#8c564b","#e377c2","#7f7f7f",
                               "#bcbd22","#17becf","#aec7e8","#ffbb78",
                               "#98df8a","#ff9896","#c5b0d5","#c49c94"))
plot(p)




##################################################################################################################################

###############################################################
################# An Error Occurred ###########################
###############################################################


devtools::install_github("vmikk/metagMisc")
library(metagMisc)


#sample_names(complete_physeq) <- paste0("z_", sample_names(complete_physeq))
phydf =phyloseq_to_df(complete_physeq, addtax = T, addtot = T, addmaxrank = F,sorting = "abundance")
#colnames(phydf) <- gsub(pattern = "^z_", replacement = "", x = colnames(phydf))



















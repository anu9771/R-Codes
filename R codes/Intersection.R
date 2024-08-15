# setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder\\EarlyCluster")
# 
# F_E1=read.csv("Early_Cluster1.csv",header = T)
# P_E1=F_E1$x
# P_E1
# 
# F_E2=read.csv("Early_Cluster2.csv",header = T)
# P_E2=F_E2$x
# P_E2
# 
# F_E3=read.csv("Early_Cluster3.csv",header = T)
# P_E3=F_E3$x
# P_E3
# 
# F_E4=read.csv("Early_Cluster4.csv",header = T)
# P_E4=F_E4$x
# P_E4
# 
# F_E5=read.csv("Early_Cluster5.csv",header = T)
# P_E5=F_E5$x
# P_E5
# 
# F_E6=read.csv("Early_Cluster6.csv",header = T)
# P_E6=F_E6$x
# P_E6
# 
# 
# setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder\\MidCluster")
# 
# F_M1=read.csv("Mid_Cluster1.csv",header = T)
# P_M1=F_M1$x
# P_M1
# 
# F_M2=read.csv("Mid_Cluster2.csv",header = T)
# P_M2=F_M2$x
# P_M2
# 
# F_M3=read.csv("Mid_Cluster3.csv",header = T)
# P_M3=F_M3$x
# P_M3
# 
# F_M4=read.csv("Mid_Cluster4.csv",header = T)
# P_M4=F_M4$x
# P_M4
# 
# F_M5=read.csv("Mid_Cluster5.csv",header = T)
# P_M5=F_M5$x
# P_M5
# 
# 
# setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder\\LateCluster")
# 
# F_L1=read.csv("Late_Cluster1.csv",header = T)
# P_L1=F_L1$x
# P_L1
# 
# F_L2=read.csv("Late_Cluster2.csv",header = T)
# P_L2=F_L2$x
# P_L2
# 
# F_L3=read.csv("Late_Cluster3.csv",header = T)
# P_L3=F_L3$x
# P_L3
# 
# F_L4=read.csv("Late_Cluster4.csv",header = T)
# P_L4=F_L4$x
# P_L4



# > intersect(P_M1,P_E1)
# character(0)
# > intersect(P_M2,P_E1)
# [1] "Atopococcus"        "Staphylococcus"     "Micrococcus"        "Helcobacillus"
# [5] "Atopostipes"        "Granulicatella"     "Ornithinimicrobium" "Dietzia"
# [9] "Papillibacter"      "Lactonifactor"
# > intersect(P_M3,P_E1)
# [1] "Bifidobacterium" "Sporobacter"     "Barnesiella"     "Oscillospira"
# > intersect(P_M4,P_E1)
# [1] "Succinispira"
# [2] "Prevotella"
# [3] "[Clostridium]"
# [4] "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Clostridium"
# [5] "CF231"
# [6] "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Clostridium"
# [7] "Roseburia"
# [8] "Alistipes"
# > intersect(P_M5,P_E1)
# [1] "Succinivibrio"      "Treponema"          "Anaerobiospirillum"
# >
# ######################################################################################################
# 
# > intersect(P_M1,P_E2)
# character(0)
# > intersect(P_M2,P_E2)
# [1] "Caldimonas"
# > intersect(P_M3,P_E2)
# [1] "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Bacteroides"
# > intersect(P_M4,P_E2)
# character(0)
# > intersect(P_M5,P_E2)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_M1,P_E3)
# character(0)
# > intersect(P_M2,P_E3)
# [1] "Corynebacterium" "Millisia"        "Trueperella"     "GW-34"           "ph2"
# [6] "Weissella"       "Kribbia"         "Isobaculum"      "Serinibacter"
# > intersect(P_M3,P_E3)
# [1] "Helcococcus"   "Petrimonas"    "Campylobacter"
# > intersect(P_M4,P_E3)
# [1] "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Clostridium"
# [2] "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Ruminococcus"
# > intersect(P_M5,P_E3)
# [1] "Pediococcus"
# >
# ######################################################################################################
# 
# > intersect(P_M1,P_E4)
# [1] "Streptococcus"
# [2] "Kocuria"
# [3] "Elizabethkingia"
# [4] "Lactococcus"
# [5] "Acinetobacter"
# [6] "Snodgrassella"
# [7] "Ochrobactrum"
# [8] "Alkanindiges"
# [9] "Rhizobium"
# [10] "Pseudoxanthomonas"
# [11] "Pedobacter"
# [12] "Nonlabens"
# [13] "Brevundimonas"
# [14] "Bacteria_Proteobacteria_Gammaproteobacteria_Pseudomonadales_Pseudomonadaceae_Pseudomonas"
# > intersect(P_M2,P_E4)
# [1] "Rothia"        "Enhydrobacter"
# > intersect(P_M3,P_E4)
# character(0)
# > intersect(P_M4,P_E4)
# [1] "Bacteria_Firmicutes_Clostridia_Clostridiales_Peptostreptococcaceae_Clostridium"
# > intersect(P_M5,P_E4)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_M1,P_E5)
# [1] "Sphingobacterium"
# > intersect(P_M2,P_E5)
# [1] "Arthrobacter" "Deinococcus"
# > intersect(P_M3,P_E5)
# [1] "Sphingobium"
# > intersect(P_M4,P_E5)
# character(0)
# > intersect(P_M5,P_E5)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_M1,P_E6)
# character(0)
# > intersect(P_M2,P_E6)
# character(0)
# > intersect(P_M3,P_E6)
# [1] "Novosphingobium" "Lacticigenium"
# > intersect(P_M4,P_E6)
# [1] "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides"
# > intersect(P_M5,P_E6)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_E1,P_L1)
# [1] "Bifidobacterium" "Succinivibrio" 
# >
# ######################################################################################################
# 
# > intersect(P_E1,P_L2)
# [1] "Staphylococcus" "Succinispira"   "Alistipes"     
# [4] "[Clostridium]"  "Treponema"      "CF231"         
# [7] "Sporobacter"    "Papillibacter"  "Barnesiella"   
# [10] "Oscillospira"
#>
# ######################################################################################################
# 
# > intersect(P_E1,P_L3)
# [1] "Atopococcus"                                                             
# [2] "Micrococcus"                                                             
# [3] "Granulicatella"                                                          
# [4] "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Clostridium"
# [5] "Roseburia"                                                               
# [6] "Atopostipes"                                                             
# [7] "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Clostridium"
# [8] "Prevotella"                                                              
# [9] "Lactonifactor"                                                           
# >
# ######################################################################################################
# 
# > intersect(P_E1,P_L4)
# [1] "Helcobacillus"      "Ornithinimicrobium" "Dietzia"           
# > 
# ######################################################################################################
# 
# > intersect(P_E2,P_L1)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_E2,P_L2)
# [1] "Caldimonas"                                                              
# [2] "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Bacteroides"
# > 
# ######################################################################################################
# 
# > intersect(P_E2,P_L3)
# character(0)
# > 
# ######################################################################################################
# 
# > intersect(P_E2,P_L4)
# character(0)
# > 
# ######################################################################################################
# 
# > intersect(P_E3,P_L1)
# character(0)
# > 
# ######################################################################################################
# 
# > intersect(P_E3,P_L2)
# [1] "Weissella"  "Isobaculum"
# >
# ######################################################################################################
# 
# > intersect(P_E3,P_L3)
# [1] "Corynebacterium"                                                          
# [2] "Millisia"                                                                 
# [3] "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Clostridium"  
# [4] "Trueperella"                                                              
# [5] "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Ruminococcus"
# [6] "ph2"                                                                      
# [7] "GW-34"                                                                    
# [8] "Kribbia"                                                                  
# [9] "Petrimonas"                                                               
# [10] "Serinibacter"                                                             
# >
# ######################################################################################################
# 
# > intersect(P_E3,P_L4)
# character(0)
# > 
# ######################################################################################################
# 
# > intersect(P_E4,P_L1)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_E4,P_L2)
# [1] "Rothia"       "Ochrobactrum" "Sphingomonas"
# > 
# ######################################################################################################
# 
# > intersect(P_E4,P_L3)
# [1] "Bacteria_Firmicutes_Clostridia_Clostridiales_Peptostreptococcaceae_Clostridium"          
# [2] "Bacteria_Proteobacteria_Gammaproteobacteria_Pseudomonadales_Pseudomonadaceae_Pseudomonas"
# > 
# ######################################################################################################
# 
# > intersect(P_E4,P_L4)
# [1] "Streptococcus"   "Elizabethkingia" "Enhydrobacter"  
# [4] "Kocuria"         "Acinetobacter"   "Alkanindiges"   
# [7] "Rhizobium"      
# > 
# ######################################################################################################
# 
# > intersect(P_E5,P_L1)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_E5,P_L2)
# [1] "Deinococcus"    "Sphingobium"    "Bergeyella"    
# [4] "Deinobacterium"
# > 
# ######################################################################################################
# 
# > intersect(P_E5,P_L3)
# character(0)
# > 
# ######################################################################################################
# 
# > intersect(P_E5,P_L4)
# [1] "Arthrobacter"
# > 
# ######################################################################################################
# 
# > intersect(P_E6,P_L1)
# [1] "Lacticigenium"                                                              
# [2] "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides"
# >
# ######################################################################################################
# 
# > intersect(P_E6,P_L2)
# character(0)
# >
# ######################################################################################################
# 
# > intersect(P_E6,P_L3)
# character(0)
# > 
# ######################################################################################################
# 
# > intersect(P_E6,P_L4)
# [1] "Novosphingobium"
# > 
# ######################################################################################################


# 
# 
# 
# intersect(P_E6,P_L4)
# 
# 
# library("gplots")
# #install.packages("ggVennDiagram")
# library(ggVennDiagram)
# 
# use list as input
# x <-list('Early6'=P_E6,
#          'Late4'=P_L4,"Mid3"=P_M3)

# # create Venn diagram with two sets
# ggVennDiagram(x,label = "count")




#For SparCC Phases 

setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder")

fileListEarly=list.files(path="EarlyCluster",pattern = ".csv")
fileListEarly

fileListMid=list.files(path="MidCluster",pattern = ".csv")
fileListMid

fileListLate=list.files(path="LateCluster",pattern = ".csv")
fileListLate

for(i in fileListEarly){
  path1=normalizePath(file.path('./EarlyCluster',paste(i, sep='')))
  file1=(read.csv(path1))
  Early=(file1$x)
  
  for(j in fileListMid){
    path2=normalizePath(file.path('./MidCluster',paste(j, sep='')))
    file2=(read.csv(path2))
    Mid=(file2$x)
    
    cat(i,"-",j)
    print(intersect(Early,Mid))
    print("*************************************************************************")
    
    pathTowards=normalizePath(file.path('./FinalRes/Intersections',paste(i,"___",j,"_Intersections.txt", sep='')))
    write.table(intersect(Early,Mid), 
                file = pathTowards, 
                col.names = FALSE, 
                row.names = FALSE, 
    )
    
  }
  
  
}


#For SparCC Breeds

setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder")

fileListFRFRC=list.files(path="FRFRCCluster",pattern = ".csv")
fileListFRFRC

fileListJRJRC=list.files(path="JRJRCCluster",pattern = ".csv")
fileListJRJRC

fileListAsian=list.files(path="AsianCluster",pattern = ".csv")
fileListAsian

for(i in fileListFRFRC){
  path1=normalizePath(file.path('./FRFRCCluster',paste(i, sep='')))
  file1=(read.csv(path1))
  FRFRC=(file1$x)
  
  for(j in fileListJRJRC){
    path2=normalizePath(file.path('./JRJRCCluster',paste(j, sep='')))
    file2=(read.csv(path2))
    JRJRC=(file2$x)
    
    cat(i,"-",j)
    print(intersect(FRFRC,JRJRC))
    print("*************************************************************************")
    
    pathTowards=normalizePath(file.path('./FinalRes/Intersections',paste(i,"___",j,"_Intersections.txt", sep='')))
    write.table(intersect(FRFRC,JRJRC), 
                file = pathTowards, 
                col.names = FALSE, 
                row.names = FALSE, 
    )
    
  }
  
  
}


#For SPRING Phases 

setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder")

fileListSpringEarly=list.files(path="Spring_EarlyCluster",pattern = ".csv")
fileListSpringEarly

fileListSpringMid=list.files(path="Spring_MidCluster",pattern = ".csv")
fileListSpringMid

fileListSpringLate=list.files(path="Spring_LateCluster",pattern = ".csv")
fileListSpringLate


for(i in fileListSpringEarly){
  path1=normalizePath(file.path('./Spring_EarlyCluster',paste(i, sep='')))
  file1=(read.csv(path1))
  Early=(file1$x)
  
  for(j in fileListSpringMid){
    path2=normalizePath(file.path('./Spring_MidCluster',paste(j, sep='')))
    file2=(read.csv(path2))
    Mid=(file2$x)
    
    cat(i,"-",j)
    print(intersect(Early,Mid))
    print("*************************************************************************")
    
    pathTowards=normalizePath(file.path('./FinalRes/Intersections',paste(i,"___",j,"_Intersections.txt", sep='')))
    write.table(intersect(Early,Mid), 
                file = pathTowards, 
                col.names = FALSE, 
                row.names = FALSE, 
    )
    
  }
  
  
}



#For SPRING Breeds 

setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder")

fileListSpringFRFRC=list.files(path="Spring_FRFRCCluster",pattern = ".csv")
fileListSpringFRFRC

fileListSpringJRJRC=list.files(path="Spring_JRJRCCluster",pattern = ".csv")
fileListSpringJRJRC

fileListSpringAsian=list.files(path="Spring_AsianCluster",pattern = ".csv")
fileListSpringAsian


for(i in fileListSpringJRJRC){
  path1=normalizePath(file.path('./Spring_JRJRCCluster',paste(i, sep='')))
  file1=(read.csv(path1))
  JRJRC=(file1$x)
  
  for(j in fileListSpringAsian){
    path2=normalizePath(file.path('./Spring_AsianCluster',paste(j, sep='')))
    file2=(read.csv(path2))
    Asian=(file2$x)
    
    cat(i,"-",j)
    print(intersect(JRJRC,Asian))
    print("*************************************************************************")
    
    pathTowards=normalizePath(file.path('./FinalRes/Intersections',paste(i,"___",j,"_Intersections.txt", sep='')))
    write.table(intersect(JRJRC,Asian), 
                file = pathTowards, 
                col.names = FALSE, 
                row.names = FALSE, 
    )
    
  }
  
  
}





#Functions to Intersections

setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder")

AllFunctionFiles=list.files(path="FinalRes",pattern = ".txt")
AllFunctionFiles

for(i in AllFunctionFiles){
  path1=normalizePath(file.path('./FinalRes',paste(i, sep='')))
  file1=(read.table(path1))
  df=data.frame(file1$X3.y,file1$X3.x)
  write.table(df,"FinalRes\\Intersections\\Intersection_Functions\\All_Functions.txt",append = TRUE,row.names = FALSE,col.names = FALSE)
}


AllFuncDF=read.table("FinalRes\\Intersections\\Intersection_Functions\\All_Functions.txt")
AllFuncDF=unique(AllFuncDF)
write.table(AllFuncDF,"FinalRes\\Intersections\\Intersection_Functions\\All_Functions_Unique.txt",append = FALSE,row.names = FALSE,col.names = c("Genus","Function"))


fTable=read.table("FinalRes\\Intersections\\Intersection_Functions\\All_Functions_Unique.txt")
#testTable=read.table("FinalRes\\Intersections\\Mid_Cluster2.csv___Late_Cluster4.csv_Intersections.txt")

intersectionFiles=list.files(path="FinalRes\\Intersections",pattern = ".txt")
intersectionFiles

for(i in intersectionFiles){
  df1=data.frame(Genus=c(),Function=c())
  path1=normalizePath(file.path('./FinalRes/Intersections',paste(i, sep='')))
  print(path1)
  
  tryCatch(
    expr={
      file=read.table(path1)
      #print(file$V1)
      for(a in file$V1){
        j=0
        while(j<259){
          j=j+1
          try(if(a==fTable$V1[j]){
            cat(a,":",fTable$V2[j],"\n")
            df2=data.frame(Genus=c(a),Function=c(fTable$V2[j]))
            df1=rbind(df1,df2)
            pathTowards=normalizePath(file.path('./FinalRes/Intersections/Intersection_Functions',paste(i,"_Intersection_Functions.txt", sep='')))
            write.table(df1,file = pathTowards,append = FALSE)
            
          })
          
        }
      }
    },
    error = function(e){         
      # print("There was an error message.")
    }
    
  )
  
  
  
  
  
}

#Iterate through all files
# Before that, add SPRING intersections
# for(i in testTable$V1){
#   j=0
#   while(j<259){
#     j=j+1
#     try(if(i==fTable$V1[j]){
#       cat(i,":",fTable$V2[j],"\n")
#     })
# 
#   }
# }



df1=data.frame(a=c(5,10),b=c(10,20))
df1

df2=data.frame(a=c(100,200),b=c(1000,2000))


df1=rbind(df1,df2)
df1


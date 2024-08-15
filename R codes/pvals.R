library(dplyr)

getwd()
setwd("C:\\Users\\Anushka\\PycharmProjects\\MilkMicrobiome")

fileList=list.files(path="pval",pattern = ".txt")
fileList

for(i in fileList){
  setwd("C:\\Users\\Anushka\\PycharmProjects\\MilkMicrobiome\\pval")
  DataFrame=read.table(i)
  DataFrame=arrange(DataFrame,V2)
  print(DataFrame)
  setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder\\FinalRes\\Pvals")
  write.table(DataFrame,i)
}


#Extra: Unnecessary

setwd("C:\\Users\\Anushka\\PycharmProjects\\MilkMicrobiome\\pval")
DataFrame=read.table("Mid_Cluster1.txt")
DataFrame=arrange(DataFrame,V2)
print(DataFrame)
write.csv(DataFrame,"C:\\Users\\Anushka\\Desktop\\MyFile.csv", row.names = FALSE)

#Extra: Unnecessary







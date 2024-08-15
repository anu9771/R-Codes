# install.packages("rrtable")
# install.packages("R2wd")


library(ggplot2)
library(rrtable)


setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology\\New folder")


funcFile=list.files(path="FinalRes",pattern = ".txt")
funcFile


for(i in funcFile){
  
  path=normalizePath(file.path('./FinalRes',paste(i, sep='')))
  File=(read.table(path,header = T))
  Table=as.data.frame(table(File$X3.x))
  print(i)
  print(Table)
  pathTowards=normalizePath(file.path('./FinalRes/Func_Frequency',paste(i,"_FuncFreq.txt", sep='')))
  write.table(Table, 
              file = pathTowards, 
              col.names = TRUE, 
              row.names = FALSE, 
              )
  
  
  try(
    
    plot2docx(ggplot(Table, aes(x = Var1, y = Freq)) +
                geom_col(fill = "#0099f9") +
                labs(
                  title = i
                )+
                theme_classic()+coord_flip(),title=i,append=TRUE)
  )
  
  
}




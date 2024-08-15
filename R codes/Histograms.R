# install.packages("rrtable")
# install.packages("R2wd")


library(ggplot2)
library(rrtable)


setwd("C:\\Users\\Anushka\\Desktop\\Anushka\\Research\\Methodology")


funcFile=list.files(path="./New folder/FinalRes",pattern = ".txt")
funcFile


for(i in funcFile){
  
  path=normalizePath(file.path('./New folder/FinalRes',paste(i, sep='')))
  File=(read.table(path,header = T))
  Table=as.data.frame(table(File$X3.x))
  print(i)
  print(Table)
  pathTowards=normalizePath(file.path('./New folder/FinalRes/Func_Frequency',paste(i,"_FuncFreq.txt", sep='')))
  write.table(Table, 
              file = pathTowards, 
              col.names = TRUE, 
              row.names = FALSE, 
  )
  
  try(
    
    plot2docx(ggplot(Table, aes(x = Var1, y = Freq)) +
                geom_col(fill = "#0099f9") +
                labs(
                  x="Functions",
                  y="Frequencies of Functions"
                )+
                theme_classic()+coord_flip(),title=i,append=TRUE)
  )
  
  
  
    

                
  
  
  
}



#########################################################################
# Install required package (if not already installed)
# install.packages("readxl")
library(readxl)

setwd("C:\\Users\\audar\\Desktop\\Anushka\\Research\\Methodology\\New folder\\FinalRes\\Func_Frequency\\All_Clusters_Togather")

# Load the ggplot2 library
library(ggplot2)

# Create a data frame with the provided data
phases=c("EarlyPhase.xlsx","MidPhase.xlsx","LatePhase.xlsx")

for(i in phases){
  data=read_excel(i)
  
  # Create a bar graph using ggplot2
  p=ggplot(data, aes(x = ID, y = Frequency, fill = Function)) +
    geom_bar(stat = "identity") +
    labs(title = "Frequency of Functions by ID", x = "ID", y = "Frequency") +
    scale_fill_manual(values=colorRampPalette(c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
                                                "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
                                                "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
                                                "#17becf", "#9edae5"))(32))
  p=p+theme(legend.text = element_text(size = 8),
            legend.title = element_text(size = 10),
            legend.key.size = unit(0.5, "cm"))
  
  plot(p)
  
  
  
}


  




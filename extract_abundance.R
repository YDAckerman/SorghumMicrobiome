## extract the abundance data

library(plyr)
lines <- readLines("~/Documents/PurdomGroup/data/abundance.txt")

data <- as.data.frame(llply(lines, function(line){
    splitLine <- strsplit(line, split = "\t")
    colname <- splitLine[[1]][1]
    data <- splitLine[[1]][2:3611]
    if(any(!is.na(as.numeric(data)))){
        data <- as.numeric(data)
    }
    col <- list(data)
    names(col) <- colname
    col
}), stringsAsFactors = FALSE)


## Looks like abundance data starts around 213, however it's in a
## format that is different from out sorghum-related data: the abundances
## are decimal numbers:

summary(rowSums(data[,213:3513]))
  ##  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  ## 600.0   784.7   791.7   788.4   796.5   800.0 
summary(colSums(data[,213:3513]))
  ## Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  ##  0.0      0.2      5.6    862.2     87.7 358800.0 

## I'll have to check out documentation?

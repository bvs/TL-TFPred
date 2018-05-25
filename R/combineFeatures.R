combineFeatures <- function(seq,csv) {
source("seqCTD.R")
ctd <- seqCTD(seq)
seqCSV <- read.csv("csv",header=FALSE)
seqCSV <- t(seqCSV[,-1])
mat <- rbind(seqCSV,ctd)
return(mat)
}

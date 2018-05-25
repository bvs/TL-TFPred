seqCTD <- function(seq) {

library(protr)

ctd <- c()
ctd1 <- c()
ctd2 <- c()
ctd3 <- c()

group1 = list(hydrophobicity  = c('R', 'K', 'E', 'D', 'Q', 'N'),
		normwaalsvolume = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
		polarizability  = c('G', 'A', 'S', 'D', 'T'),
		secondarystruct = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
		solventaccess   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))
     
group2 = list(hydrophobicity  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
		normwaalsvolume = c('N', 'V', 'E', 'Q', 'I', 'L'),
		polarizability  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
		secondarystruct = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
		solventaccess   = c('R', 'K', 'Q', 'E', 'N', 'D'))
     
group3 = list(hydrophobicity  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
		normwaalsvolume = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
		polarizability  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
		secondarystruct = c('G', 'N', 'P', 'S', 'D'),
		solventaccess   = c('M', 'S', 'P', 'T', 'H', 'Y'))

for(i in 1:length(seq)) {
	ctd1 <- cbind(ctd1, as.vector(extractCTDCClass(seq[[1]], aagroup1 = group1, aagroup2 = group2, aagroup3 = group3)))
}

for(i in 1:length(seq)) {
	ctd2 <- cbind(ctd2, as.vector(extractCTDTClass(seq[[1]], aagroup1 = group1, aagroup2 = group2, aagroup3 = group3)))
}

for(i in 1:length(seq)) {
	ctd3 <- cbind(ctd3, as.vector(extractCTDDClass(seq[[1]], aagroup1 = group1, aagroup2 = group2, aagroup3 = group3)))
}

ctd <- rbind(ctd1,ctd2,ctd3)

return(ctd)

}

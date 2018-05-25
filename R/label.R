### Function used to create binary labels ####
#vect1 <- c(4, 5, 10, 3, 1)
#rep(1:5, vect1)
biLabel <- function (Cs) {
	mat <- list()
	for(label in levels(Cs)) {
		#print(label)
		rows <- (Cs == label)
		Cn <- rep(0, length(Cs))
		Cn[rows] <- 1
		#print(Cn)
		mat[[label]] <- Cn
	}
	#riceLabel2 <- matrix(c(unlist(riceLabel)),ncol=56,nrow=2408)
	return(mat)
}
### End of biLabel ###

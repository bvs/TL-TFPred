### GTL: Graph Co-regularized Transfer Learning using Collective Matrix Factorization
### Usage: GTL(Xs=X_src,Xt=X_tar,Ys=Y_src,Yt=Y_tar,sigma=10,p=10,lambda=1.0,gama=1.0,iters=100)
GTL <- function(Xs,Xt,Ys,Yt,p,lambda,gama,sigma,iters) {
	print("Step: 1 Setting predefined variables");
	X <- as.matrix(cbind(as.matrix(Xs),as.matrix(Xt)))
	Y <- as.matrix(rbind(as.matrix(Ys),as.matrix(Yt)))
	m <- dim(X)[1]
	c <- length(unique(Y))
	ns <- dim(Xs)[2]
	nt <- dim(Xt)[2]
	YY <- c()
	for(i in sort(as.vector(array(unique(as.vector(Y)), c(1, length(unique(Y))))))) {
		YY <- cbind(YY, as.vector(Y) == i)
		YY[YY == FALSE] <- 0
	}
	Y <- as.matrix(apply(YY,1,which.max))
	Ys <- YY[1:ns,];
	Yt <- YY[(ns+1):dim(YY)[1],];
	print("Step: 2 Construction of Graph Laplacian for Xs and Xt (This will take a while!)");
	k <- p;
	Metric <- "Cosine";
	NeighborMode <- "KNN";
	WeightMode <- "Cosine";
	bNormalizeGraph <- 0;
	Wus <- affinity(as.matrix(Xs),k);
	Dus <- diagnalGTL(Wus,bNormalizeGraph);
	Wut <- affinity(as.matrix(Xt),k);
	Dut <- diagnalGTL(Wut,bNormalizeGraph);
	Wvt <- affinity(t(as.matrix(Xt)),k);
	Dvt <- diagnalGTL(Wvt,bNormalizeGraph);
	U <- t(matrix(runif(c*m),c));
	Vs <- 0.1+0.8*Ys;
	Vt <- matrix(runif(c*nt),nt);
	print("Step: 3 Starting Graph Co-Regularized Transfer Learning (GTL) Algorithm");
	Acc <- c()
	Obj <- c()
	prob1 <- c()
	prob2 <- c()
	cls2 <- c()
	cls <- list()
	Vt2 <- list()
	for(i in 1:iters) {
		U <- U * sqrt((Xs%*%Vs+Xt%*%Vt+lambda*Wus%*%U+lambda*Wut%*%U)/(U%*%(t(Vs)%*%Vs)+U%*%(t(Vt)%*%Vt)+lambda*Dus%*%U+lambda*Dut%*%U+.Machine$double.eps));

		Vs <- Vs * sqrt(Vs/(Vs%*%(t(Vs)%*%Vs)+.Machine$double.eps))			
		
		Vt <- Vt * sqrt((t(Xt)%*%U+gama*Wvt%*%Vt+sigma*Vt)/(Vt%*%(t(U)%*%U)+gama*Dvt%*%Vt+sigma*Vt%*%(t(Vt)%*%Vt)+.Machine$double.eps))	
		
		Vt2[[i]] <- Vt

		prob1 <- append(prob1, as.vector(Vt[,1]))
		prob2 <- append(prob2, as.vector(Vt[,2]))	
		Cls <- as.matrix(apply(Vt,1,which.max))
		cls2 <- rbind(cls2,Cls)
		Lbl <- as.matrix(apply(Yt,1,which.max))
		acc <- length(which(Cls == Lbl))/nt*100;
		Acc <- rbind(Acc,acc)
		if(i %% 10 ==0) {
			#print(paste("Iteration:",i, "Accuracy:",round(acc, digits=3)));
		}
	}
	cls[['Vt']] <- Vt2; cls[['Prob_No']] <- prob1; cls[['Prob_Yes']] <- prob2; cls[['Class Label']] <- cls2; cls[['Accuracy']] <- Acc;
	return(cls)
	print("Algorithm GTL terminated!!!");
}
### End of GTL ###
###
### affinity ###
affinity <- function(fea,k) {
#	k <- 10;
library(methods)
	Metric <- "Cosine";
	NeighborMode <- "KNN";
	WeightMode <- "Cosine";
	bNormalized <- 0;
	bSelfConnected <- 1;
	nSmp <- dim(fea)[1]
	maxM <- 62500000 # 500M
	BlockSize <- floor(maxM/(nSmp*3));
	if((NeighborMode == "KNN") && (k > 0)) {
		if(Metric == "Cosine") {
			if(!(bNormalized)) {
				fea <- as(fea, "matrix");
				nSmp <- dim(fea)[1];
				nFea <- dim(fea)[2];
				if(is.matrix(fea)) {
					fea2 <- t(fea);
					rm(fea);
					for(i in 1:nSmp) {
						fea2[,i] <- fea2[,i]/max(1e-10,sum(fea2[,i]^2)^0.5)
					}
					fea <- t(fea2);
					rm(fea2);
				} else { 
					feaNorm <- apply(fea^2,1,sum)^0.5;
					for(i in 1:nSmp) {
						fea[i,] <- fea[i,]/max(1e-12, feaNorm[i]);
					}
				}
			}
		}
		fea <- as(fea, "sparseMatrix");
		G <- array(0,c(nSmp*(k+1),3));
		for(i in 1:ceiling(nSmp/BlockSize)) {
			if(i == ceiling(nSmp/BlockSize)) {
				smpIdx <- ((i-1)*BlockSize+1):nSmp;
				fea <- as(fea,"sparseMatrix")
				dist <- fea[smpIdx,] %*% t(fea);
				fea <- as(fea,"matrix")
				dump <- t(apply(-dist,1,sort));
				idx <- t(apply(-dist,1,order));
				idx <- idx[,0:k+1];
				dump <- -dump[,0:k+1];
				G[((i-1)*BlockSize*(k+1)+1):(nSmp*(k+1)),1] <- kronecker(matrix(1,k+1,1),smpIdx);
				G[((i-1)*BlockSize*(k+1)+1):(nSmp*(k+1)),2] <- as.vector(idx);
				G[((i-1)*BlockSize*(k+1)+1):(nSmp*(k+1)),3] <- as.vector(dump);
			} else {
				smpIdx <- ((i-1)*BlockSize+1):(i*BlockSize);
				fea <- as(fea,"sparseMatrix")
				dist <- fea[smpIdx,] %*% t(fea);
				fea <- as(fea,"matrix")	
				dump <- t(apply(-dist,1,sort));
				idx <- t(apply(-dist,1,order));
				idx <- idx[,0:k+1];
				dump <- -dump[,0:k+1];
				G[((i-1)*BlockSize*(k+1)+1):((i*BlockSize)*(k+1)),1] <- kronecker(matrix(1,k+1,1),smpIdx);
				G[((i-1)*BlockSize*(k+1)+1):((i*BlockSize)*(k+1)),2] <- as.vector(idx);
				G[((i-1)*BlockSize*(k+1)+1):((i*BlockSize)*(k+1)),3] <- as.vector(dump);
			}
		}
		library("Matrix")
#		W <- sparseMatrix(i=G[,1],j=G[,2],x=G[,3])
		W <- sparseMatrix(i=G[,1],j=G[,2],x=G[,3], dims=c(nSmp,nSmp))
		if(!(bSelfConnected)) {
       			for(i in 1:dim(W)[1]) {
				W[i,i] = 0;
			}
		}
		W <- pmax.sparse(W,t(W));
		return(W);
	}
}
### End of affinity ###
###
### pmax.sparse ### 
pmax.sparse <- function(..., na.rm = FALSE) {
	#check that all matrices have conforming sizes
	num.rows <- unique(sapply(list(...), nrow))
	num.cols <- unique(sapply(list(...), ncol))
	stopifnot(length(num.rows) == 1)
	stopifnot(length(num.cols) == 1)
	cat.summary <- do.call(rbind, lapply(list(...), summary))
	out.summary <- aggregate(x ~ i + j, data = cat.summary, max, na.rm)
	sparseMatrix(i = out.summary$i, j = out.summary$j, x = out.summary$x, dims = c(num.rows, num.cols))
}
### End of pmax.sparse ###
### diagnalGTL ###
diagnalGTL <- function(W, bNormalizeGraph) {
	if(bNormalizeGraph) {
		D <- as(1/sqrt(rowSums(W)), "spareMatrix");
		D[is.infinite(D)] <- 0;
		D <- diag(c(as.matrix(D)));
		W <- D %*% W %*% D;
		D <- as(D, "matrix");
		D[D>0] <- 1;
		D <- as(D,"sparseMatrix");
	} else {
		D <- as(diag(c(as.matrix(rowSums(W)))), "sparseMatrix");
		}
	return(D)
}
### End of diagnalGTL ###

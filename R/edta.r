edta <- function(TP,FN,TN,FP){

	Se <- (TP+.5)/(TP+FN+1)
	Fp <- 1 - (TN+.5)/(TN+FP+1)

	Y1 <- log(Se/(1-Se))
	Y2 <- log(Fp/(1-Fp))

	n1 <- TP + FN + 1
	n0 <- TN + FP + 1

	N <- length(TP)   # number of trial
	p <- 2  # dimension of the multivariate meta-analysis

	V1 <- (n1*Se*(1-Se))^-1
	V2 <- (n0*Fp*(1-Fp))^-1

	y <- cbind(Y1,Y2)
	S <- cbind(V1,rep(0,N),V2)

	R1 <- list(y=y,S=S,Se=Se,Fp=Fp)
	
	return(R1)
	
}

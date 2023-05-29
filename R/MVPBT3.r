# Parametric resampling function from the estimated distribution (by constrained ML)

PBS3 <- function(y,S,b0,V0){

	N <- dim(y)[1]
	p <- dim(y)[2]

	y.pb <- matrix(numeric(N*p),N)

	for(i in 1:N){

		yi <- y[i,]
		Si <- matrix(c(S[i,1],S[i,2],S[i,2],S[i,3]),p)

		Vi <- Si + V0

		Xi <- diag( sqrt(diag(Si) + diag(V0))^-1 )
		Psii <- Xi %*% Vi %*% t(Xi)

		mui <- Xi %*% b0

		y.pb[i,] <- MASS::ginv(Xi) %*% MASS::mvrnorm(1, mui, Psii)

	}
		
	return(y.pb)
	
}

# Parametric bootstrap test for publication bias of DTA meta-analysis

MVPBT3 <- function(y,S,B=2000){

	V0 <- mvmeta::mvmeta(y,S)$Psi

	Q0 <- MVPBT2(y,S)

	T.b <- numeric(B)
	
		for(b in 1:B){
	
			y.pb <- PBS3(y,S,Q0$b0,V0)
			mm.pb <- mvmeta::mvmeta(y.pb,S)
			Q.b <- MVPBT2(y.pb,S)
			T.b[b] <- Q.b$T
			
			print1 <- paste0("The ",b,"th bootstrap (/",B,") is completed.")
			if(b%%100==0) print(print1)
		
		}


	
	QT <- function(x,x0){

		x1 <- sort(c(x,x0))
		w1 <- which(x1==as.numeric(x0))
		qt <- 1 - w1/(length(x)+1)
		return(qt)
 
	}

	P <- QT(T.b,Q0$T)	# p-value for the bootstrap test

	R1 <- list(T.b=T.b,T=Q0$T,P=P)
	
	return(R1)

}

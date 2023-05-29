# Noma's publication bias test for DTA meta-analysis

MVPBT2 <- function(y,S){

	V0 <- mvmeta::mvmeta(y,S)$Psi

	N <- dim(y)[1]
	p <- dim(V0)[1]

	Q1 <- matrix(numeric(p*p),p)
	Q2 <- numeric(p)

	for(i in 1:N){

		yi <- y[i,]
		Si <- matrix(c(S[i,1],S[i,2],S[i,2],S[i,3]),p)

		Vi <- Si + V0

		Xi <- diag( sqrt(diag(Si) + diag(V0))^-1 )
		Psii <- Xi %*% Vi %*% t(Xi)
		
		zi <- Xi %*% yi

		Q1 <- Q1 + t(Xi) %*% MASS::ginv(Psii) %*% Xi
		Q2 <- Q2 + t(Xi) %*% MASS::ginv(Psii) %*% zi
	
	}

	b0 <- as.numeric( MASS::ginv(Q1) %*% Q2 )		# constrained MLE of b under H_0

	Q3 <- Q4 <- numeric(p)
	Q5 <- Q6 <- Q7 <- matrix(numeric(p*p),p)

	for(i in 1:N){

		yi <- y[i,]
		Si <- matrix(c(S[i,1],S[i,2],S[i,2],S[i,3]),p)

		Vi <- Si + V0

		Xi <- diag( sqrt(diag(Si) + diag(V0))^-1 )
		Psii <- Xi %*% Vi %*% t(Xi)
		
		zi <- Xi %*% yi

		Q3 <- Q3 + MASS::ginv(Psii) %*% (zi - Xi %*% b0)
		Q4 <- Q4 + t(Xi) %*% MASS::ginv(Psii) %*% (zi - Xi %*% b0)

		Q5 <- Q5 + MASS::ginv(Psii)
		Q6 <- Q6 + t(Xi) %*% MASS::ginv(Psii)
		Q7 <- Q7 + t(Xi) %*% MASS::ginv(Psii) %*% Xi
	
	}

	#U <- c(Q3,Q4)	# score vector
	#J <- cbind(rbind(Q5,Q6),rbind(Q6,Q7))	# information matrix

	U <- Q3		# score vector for a
	J <- (Q5 - t(Q6)%*%MASS::ginv(Q7)%*%Q6)
	
	T <- t(U) %*% MASS::ginv(J) %*% U		# score statistic
	P <- 1 - pchisq(T,df=p)			# P-value

	R1 <- list(T=T,P=P,b0=b0)
	
	return(R1)

} 



sdta <- function(Se,Fp,Secl=NULL,Secu=NULL,Fpcl=NULL,Fpcu=NULL){

	Y1 <- log(Se/(1-Se))
	Y2 <- log(Fp/(1-Fp))

	N <- length(Y1)   # number of trial

	if(is.null(Secl)==FALSE){
		cl1 <- log(Secl/(1-Secl))
		se1 <- (Y1 - cl1)/qnorm(.975)
	}
		
	if(is.null(Fpcl)==FALSE){
		cl2 <- log(Fpcl/(1-Fpcl))
		se2 <- (Y2 - cl2)/qnorm(.975)
	}

	if(is.null(Secu)==FALSE){
		cu1 <- log(Secu/(1-Secu))
		se1 <- (cu1 - Y1)/qnorm(.975)
	}

	if(is.null(Fpcu)==FALSE){
		cu2 <- log(Fpcu/(1-Fpcu))
		se2 <- (cu2 - Y2)/qnorm(.975)
	}

	V1 <- se1*se1
	V2 <- se2*se2

	y <- cbind(Y1,Y2)
	S <- cbind(V1,rep(0,N),V2)

	R1 <- list(y=y,S=S,Se=Se,Fp=Fp)
	
	return(R1)
	
}

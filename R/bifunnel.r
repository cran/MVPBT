bifunnel <- function(y,S){

	oldpar <- par(mfrow=c(1,1))

	par(mfrow=c(1,2))

	res1 <- rma(y[,1], S[,1])
	funnel(res1,main="(a) Funnel plot for logit(Se)")

	res2 <- rma(y[,2], S[,3])
	funnel(res2,main="(b) Funnel plot for logit(FPR)")

	par(oldpar)    # Reset the graphic parameter

} 



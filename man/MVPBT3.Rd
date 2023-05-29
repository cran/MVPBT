\name{MVPBT3}
\alias{MVPBT3}
\title{Generalized Egger test to detect publication bias in bivariate meta-analysis for diagnostic accuracy test (MSSET3)}
\usage{
MVPBT3(y,S,B=2000)
}
\arguments{
\item{y}{Summary outcome statistics}

\item{S}{Covariance estimates of \code{y}}

\item{B}{Number of bootstrap resampling (default: 2000)}
}
\value{
\itemize{
\item \code{T.b}: Bootstrap samples of the efficient score statistic.
\item \code{T}: The efficient score statistic.
\item \code{P}: P-value of the publication bias test (bootstrap test).
}
}
\description{
Generalized Egger test to detect publication bias in bivariate meta-analysis for diagnostic accuracy test (called MSSET3 in Noma (2020)). This test adequately consider the uncertainties of heterogeneity variance-covariance parameters by bootstrapping.
}
\references{
Noma, H. (2020). 
Discussion of “Testing small study effects in multivariate meta-analysis” by Chuan Hong, Georgia Salanti, Sally Morton, Richard Riley, Haitao Chu, Stephen E. Kimmel, and Yong Chen. 
\emph{Biometrics} \strong{76}: 1255-1259.
\doi{10.1111/biom.13343}

Noma, H. (2022). 
MVPBT: R package for publication bias tests in meta-analysis of diagnostic accuracy studies. 
arXiv:2209.07270.
\doi{10.48550/arXiv.2209.07270}
}
\examples{

require(metafor)
require(mada)

data(cervical)

LAG <- cervical[cervical$method==2,]

fit1 <- reitsma(LAG)
summary(fit1)     # results of the bivariate meta-analysis

###

attach(LAG)

dta1 <- edta(TP,FN,TN,FP)

oldpar <- par(mfrow=c(1,1))
par(mfrow=c(1,3))

plot(fit1, predict=TRUE, cex=1.5, pch=19, sroclty=1, sroclwd=1.5, lty=2, 
 main="(a) SROC plot", xlim=c(0,1), ylim=c(0,1))
points(dta1$Fp,dta1$Se,pch=20,col="blue")
#legend(0.4,0.1,legend=c("95% confidence region","95% prediction region"),lty=c(2,3))

###

attach(dta1)

res1 <- rma(y[,1], S[,1])
funnel(res1,main="(b) Funnel plot for logit(Se)")
regtest(res1, model="lm")	# univariate Egger's test

res2 <- rma(y[,2], S[,3])
funnel(res2,main="(c) Funnel plot for logit(FPR)")
regtest(res2, model="lm")	# univariate Egger's test

###

MVPBT3(y,S,B=20)   # Generalized Egger test (MSSET3)
# This is an example command for illustration. B should be >= 1000.

par(oldpar)    # Reset the graphic parameter
}

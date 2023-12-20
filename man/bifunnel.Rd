\name{bifunnel}
\alias{bifunnel}
\title{Funnel plots for the bivariate outcomes}
\usage{
bifunnel(y,S)
}
\arguments{
\item{y}{Summary outcome statistics}

\item{S}{Covariance estimates of \code{y}}
}
\value{
Funnel plots for the logit-transformed sensitivities and false positive rates are presented.
}
\description{
Funnel plots for the bivariate outcomes of diagnostic meta-analysis are created.
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

###

attach(dta1)

bifunnel(y,S)
}

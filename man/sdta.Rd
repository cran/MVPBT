\name{sdta}
\alias{sdta}
\title{Transforming diagnostic measures to summary statistics for meta-analysis of diagnostic studies}
\usage{
sdta(Se,Fp,Secl,Secu,Fpcl,Fpcu)
}
\arguments{
\item{Se}{A vector of the sensitivity estimates}

\item{Fp}{A vector of the false positive rate estimates}

\item{Secl}{A vector of the lower confidence limits of sensitivities}

\item{Secu}{A vector of the upper confidence limits of sensitivities}

\item{Fpcl}{A vector of the lower confidence limits of false positive rates}

\item{Fpcu}{A vector of the upper confidence limits of false positive rates}
}
\value{
Summary statistics for meta-analysis are generated.
\itemize{
\item \code{y}: Logit-transformed sensitivities and false positive rates.
\item \code{S}: Within-study variances and covariances.
\item \code{Se}: Sensitivities.
\item \code{Fp}: False positive rates.
}
}
\description{
Transforming diagnostic measures to summary statistics for meta-analysis of diagnostic studies.
}
\examples{
library("mada")

MRI <- cervical[cervical$method==3,]

MRIa <- MRI[,5:8]
MRIad <- madad(MRIa)

sdta(Se=MRIad$sens$sens,Fp=MRIad$fpr$fpr,
 Secl=MRIad$sens$sens.ci[,1],Fpcl=MRIad$fpr$fpr.ci[,1])
 
sdta(Se=MRIad$sens$sens,Fp=MRIad$fpr$fpr,
 Secu=MRIad$sens$sens.ci[,2],Fpcu=MRIad$fpr$fpr.ci[,2])
}

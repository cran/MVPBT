\name{edta}
\alias{edta}
\title{Transforming contingency table data to summary statistics in diagnostic studies}
\usage{
edta(TP,FN,TN,FP)
}
\arguments{
\item{TP}{A vector of the number of true positives (TP)}

\item{FP}{A vector of the number of false positives (FP)}

\item{FN}{A vector of the number of false negatives (FN)}

\item{TN}{A vector of the number of true negatives (TN)}
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
Transforming contingency table data to summary statistics in diagnostic studies.
}
\examples{
data(cervical)
LAG <- cervical[cervical$method==2,]

attach(LAG)

dta1 <- edta(TP,FN,TN,FP)
}

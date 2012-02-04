\name{DDD-package}
\alias{DDD-package}
\alias{DDD}
\docType{package}
\title{
Diversity dependent diversification}
\description{
This package is a likelihood-based statistical package to estimate parameters
under diversity-dependent diversification models.}
\details{
\tabular{ll}{
Package: \tab DDD\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2012-02-03\cr
License: \tab GPL 2.0 ?\cr
}
dd_loglik computes the likelihood of any of 4 different models of diversity-dependent
diversification, given a set of parameters and a data set of phylogenetic branching times.

dd_ML finds the parameters that maximizes the likelihood computed by dd_loglik.

dd_SR_loglik is the same as dd_loglik except that it allows a shift in the parameters at
time tshift.

dd_SR_ML finds the parameters that maximizes the likelihood computed by dd_SR_loglik
}
\author{
Rampal S. Etienne & Bart Haegeman
Maintainer: Rampal S. Etienne (r.s.etienne@rug.nl)
}
\references{Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439}
\keyword{ models }
\seealso{
 \code{ laser }
}
\examples{
\donttest{
cat("This will estimate the parameters for a set of branching times brts with the default settings:")
brts = 1:20
dd_ML(brts)

cat("This will estimate the parameters for a set of branching times brts with the default settings:")
brts = 1:20
dd_SR_ML(brts)
}
}

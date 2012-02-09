\name{DDD-package}
\alias{DDD-package}
\alias{DDD}
\docType{package}
\title{Diversity dependent diversification}
\description{
This package is a likelihood-based statistical package to estimate parameters
under diversity-dependent diversification models.\cr \cr
   New in version 0.3:\cr
   - Models with shifts in the parameters at a time tshift.\cr
   New in version 0.4:\cr
   - Bug fix of prefactor in likelihood - bug did not affect parameter estimates.\cr
   - There is an option to condition or not on non-extinction of the phylogeny.\cr
   - The functions dd_loglik and dd_SR_loglik now accept untransformed parameters;
   the transformation is done in the auxiliary functions.\cr
   - Improved help pages.\cr
   - Auxiliary functions no longer in help pages
}

\details{
\tabular{ll}{
Package: \tab DDD\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2012-02-09\cr
License: \tab GPL 2.0\cr
}
dd_loglik computes the likelihood of any of 4 different models of diversity-dependent
diversification, given a set of parameters and a data set of phylogenetic branching times.

dd_ML finds the parameters that maximizes the likelihood computed by dd_loglik.

dd_SR_loglik is the same as dd_loglik except that it allows a shift in the parameters at
time tshift.

dd_SR_ML finds the parameters that maximizes the likelihood computed by dd_SR_loglik.
}
\author{
Rampal S. Etienne & Bart Haegeman
Maintainer: Rampal S. Etienne (r.s.etienne@rug.nl)
}
\references{Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439}
\keyword{ models }
\seealso{
\code{laser}
}
\examples{
\donttest{
cat("Estimate the parameters of the diversity-dependent model")
cat("for a set of branching times:")
brts = 1:20
dd_ML(brts)

cat("Estimating the parameters of the diversity-dependent model with shiftting rates")
cat("for a set of branching times brts with the default settings:")
brts = 1:20
dd_SR_ML(brts)
}
}

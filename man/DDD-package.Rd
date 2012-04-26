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
   - Auxiliary functions no longer in help pages.\cr
   New in version 1.0:\cr
   - Another bug fix of prefactor in likelihood - bug did not affect parameter estimates.\cr
   - Models with key innovations defined as decoupling of diversity-dependent dynamics.\cr
   - Option to compute the likelihood of a set of branching times or of the phylogeny;
   these differ by a constant factor.
   New in version 1.1\cr
   - A model with diversity-dependence in speciation and extinction rate, currently only
   available for the likelihood without shifts or decoupling.
}

\details{
\tabular{ll}{
Package: \tab DDD\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-04-10\cr
License: \tab GPL 2.0\cr
}
dd_loglik computes the likelihood of any of 4 different models of diversity-dependent
diversification, given a set of parameters and a data set of phylogenetic branching times.

dd_ML finds the parameters that maximizes the likelihood computed by dd_loglik.

dd_SR_loglik is the same as dd_loglik except that it allows a shift in the parameters at
time tshift.

dd_SR_ML finds the parameters that maximizes the likelihood computed by dd_SR_loglik.

dd_KI_loglik is the same as dd_loglik except that it allows decoupling of a subclade's dynamics
from the main clade's dynamics potentially accompanied by a shift in parameters at time t_d.

dd_KI_ML finds the parameters that maximizes the likelihood computed by dd_KI_loglik.

}
\author{
Rampal S. Etienne & Bart Haegeman
Maintainer: Rampal S. Etienne (r.s.etienne@rug.nl)
}
\references{Etienne et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439\cr
Etienne & Haegeman 2012, Am. Nat. In press.}
\keyword{ models }
\seealso{
\code{laser}
}
\examples{
\donttest{
cat("Estimate the parameters of the diversity-dependent model")
cat("for a set of branching times with the default settings:")
brts = 1:20
dd_ML(brts)

cat("Estimating the parameters of the diversity-dependent model with shiftting rates")
cat("for a set of branching times brts with the default settings:")
brts = 1:20
dd_SR_ML(brts)

cat("Estimating the parameters of the diversity-dependent model with decoupling")
cat("for a set of branching times with the default settings:")
brtsM = 1:20
brtsS = seq(0.2,1.5,0.3)
tsplit = 5
dd_ML(brtsM,brtsS,tsplit)

}
}

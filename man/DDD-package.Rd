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
   these differ by a constant factor.\cr
   New in version 1.1\cr
   - A model with diversity-dependence in speciation and extinction rate, currently only
   available for the likelihood without shifts or decoupling.\cr
   New in version 1.2\cr
   - POssibility of automatically exploring multiple initial conditons for the SR models\cr
   - Suppressing unnecessary output\cr
   - Flushing computed likelihood values immediately (even with buffer on)\cr
   - Bug fix in default value of the initial time of decoupling\cr
   New in version 1.3\cr
   - Fix of default length of idparsopt\cr
   - Fix of conversion problems for small branching times\cr
   - Default initial parameter values in dd_ML based on birth-death model\cr
   New in version 1.4\cr
   - Bug fix for when there is a large number of missing species\cr
   New in version 1.5\cr
   - Bug fix for returned output of functions when there are fixed parameters.
   Textual output of DDD was correct.\cr
   - Return to Nelder-Mead optimization for dd_ML (was set to "nlm" in versions 1.3 & 1.4)
   - New format for output off all ML-functions (dd_ML, dd_SR_ML, dd_KI_ML)
}

\details{
\tabular{ll}{
Package: \tab DDD\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2012-10-02\cr
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
Rampal S. Etienne & Bart Haegeman\cr
Maintainer: Rampal S. Etienne (r.s.etienne@rug.nl)
}
\references{Etienne et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439\cr
Etienne & Haegeman 2012, Am. Nat. In press, doi: 10.1086/667574}
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

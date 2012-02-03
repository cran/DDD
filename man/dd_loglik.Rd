\name{dd_loglik}
\alias{dd_loglik}

\title{ loglikelihood for diversity-dependent diversification models
}
\description{
  This function computes loglikelihood of a diversity-dependent diversification model
  for a given set of branching times and parameter values
}
\usage{
dd_loglik(trpars1, pars2, brts, missnumspec)
}
\arguments{
  \item{trpars1}{ transformed vector of parameters, after transformation pars/(1+pars). See dd_ML for am explanation of the parameters. }
  \item{pars2}{
   pars2[1] sets the maximum number of species for which a probability must be computed, must be larger than 1 + length(brts) 
   pars2[2] sets the model of diversity-dependence:
   pars2[2] == 1 : linear dependence in speciation rate
   pars2[2] == 2 : exponential dependence in speciation rate
   pars2[2] == 3 : linear dependence in extinction rate
   pars2[2] == 4 : exponential dependence in extinction rate
}
  \item{brts}{ a set of branching times of a phylogeny, all positive, from present to past }
  \item{missnumspec}{ The number of species that are in the clade but missing in the phylogeny } 
}
\value{the loglikelihood}
\references{Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439}
\author{Rampal S. Etienne & Bart Haegeman}

\seealso{
 \code{\link{dd_ML}}
 \code{\link{dd_SR_loglik}}
}

\examples{ dd_loglik(trpars = c(0.2,0.1,0.98), pars2 = c(100,1), brts = 1:10, missnumspec = 0) }
\keyword{ models }


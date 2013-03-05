\name{dd_loglik}
\alias{dd_loglik}
\title{Loglikelihood for diversity-dependent diversification models }
\description{
  This function computes loglikelihood of a diversity-dependent diversification model
  for a given set of branching times and parameter values.
}
\usage{dd_loglik(pars1, pars2, brts, missnumspec)}
\arguments{
  \item{pars1}{Vector of parameters: \cr \cr
    \code{pars1[1]} corresponds to lambda (speciation rate) before the shift \cr
    \code{pars1[2]} corresponds to mu (extinction rate) before the shift \cr
    \code{pars1[3]} corresponds to K (clade-level carrying capacity) before the shift
}
  \item{pars2}{Vector of model settings: \cr \cr
    \code{pars2[1]} sets the maximum number of species for which a probability must be computed.
    This must be larger than 1 + length(brts). \cr \cr
    \code{pars2[2]} sets the model of diversity-dependence: \cr
    - \code{pars2[2] == 1} linear dependence in speciation rate \cr
    - \code{pars2[2] == 2} exponential dependence in speciation rate \cr
    - \code{pars2[2] == 3} linear dependence in extinction rate \cr
    - \code{pars2[2] == 4} exponential dependence in extinction rate \cr \cr
    \code{pars2[3]} sets the conditioning on non-extinction of the phylogeny (TRUE or FALSE) \cr \cr
    \code{pars2[4]} sets whether the likelihood is for the branching times (0) or the phylogeny (1) \cr \cr
    \code{pars2[5]} sets whether the parameters and likelihood should be shown on screen (1) or not (0)\cr \cr
    \code{pars2[6]} sets whether the first data point is stem age (1) or crown age (2)
}
  \item{brts}{A set of branching times of a phylogeny, all positive }
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny } 
}
\value{The loglikelihood}
\references{Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439}
\author{Rampal S. Etienne & Bart Haegeman}
\seealso{
 \code{\link{dd_ML}},
 \code{\link{dd_SR_loglik}},
 \code{\link{dd_KI_loglik}}
}
\examples{dd_loglik(pars1 = c(0.5,0.1,100), pars2 = c(100,1,TRUE,1), brts = 1:10, missnumspec = 0) }
\keyword{models}
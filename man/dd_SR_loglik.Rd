\name{dd_SR_loglik}
\alias{dd_SR_loglik}

\title{Loglikelihood for diversity-dependent diversification models with a shift in the parameters at time t = tshift }
\description{
  This function computes loglikelihood of a diversity-dependent diversification model
  for a given set of branching times and parameter values where the parameters are allowed to shift at time t = tshift
}
\usage{dd_SR_loglik(pars1, pars2, brts, missnumspec)}
\arguments{
  \item{pars1}{Vector of parameters: \cr \cr
   \code{pars1[1]} corresponds to lambda (speciation rate) before the shift \cr
   \code{pars1[2]} corresponds to mu (extinction rate) before the shift \cr
   \code{pars1[3]} corresponds to K (clade-level carrying capacity) before the shift \cr
   \code{pars1[4]} corresponds to lambda (speciation rate) after the shift \cr
   \code{pars1[5]} corresponds to mu (extinction rate) after the shift \cr
   \code{pars1[6]} corresponds to K (clade-level carrying capacity) after the shift \cr
   \code{pars1[7]} corresponds to tshift (the time of shift)
}
  \item{pars2}{Vector of model settings: \cr \cr
    \code{pars2[1]} sets the maximum number of species for which a probability must be computed.
    This must be larger than 1 + length(brts). \cr \cr
    \code{pars2[2]} sets the model of diversity-dependence: \cr
    - \code{pars2[2] == 1} linear dependence in speciation rate \cr
    - \code{pars2[2] == 2} exponential dependence in speciation rate \cr
    - \code{pars2[2] == 3} linear dependence in extinction rate \cr
    - \code{pars2[2] == 4} exponential dependence in extinction rate \cr \cr
    \code{pars2[3]} sets the conditioning on non-extinction of the phylogeny (TRUE or FALSE) \cr
    \code{pars2[4]} sets whether the likelihood is for the branching times (0) or the phylogeny (1)
}
  \item{brts}{A set of branching times of a phylogeny, all positive }
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny } 
}
\value{The loglikelihood}
\references{
- Etienne et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439
- Etienne & Haegeman 2012. Am. Nat. In press, doi: 10.1086/667574. \cr}
\author{Rampal S. Etienne & Bart Haegeman }
\seealso{
 \code{\link{dd_SR_ML}},
 \code{\link{dd_loglik}},
 \code{\link{dd_KI_loglik}}
}
\examples{dd_SR_loglik(pars1 = c(0.2,0.1,50,0.2,0.1,70,5), pars2 = c(100,1,TRUE,1), brts = 1:10, missnumspec = 0) }
\keyword{models}

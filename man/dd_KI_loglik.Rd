\name{dd_KI_loglik}
\alias{dd_KI_loglik}

\title{Loglikelihood for diversity-dependent diversification models with decoupling of a subclade from a main clade at time t = t_d }
\description{
  This function computes loglikelihood of a diversity-dependent diversification model
  for a given set of branching times and parameter values where the diversity-dependent dynamics
  of a subclade decouple from the dynamics of the main clade at time t_d, potentially
  accompanied by a shift in parameters.
}
\usage{dd_KI_loglik(pars1, pars2, brtsM, brtsS, missnumspec)}
\arguments{
  \item{pars1}{Vector of parameters: \cr \cr
   \code{pars1[1]} corresponds to lambda_M (speciation rate) of the main clade \cr
   \code{pars1[2]} corresponds to mu_M (extinction rate) of the main clade \cr
   \code{pars1[3]} corresponds to K_M (clade-level carrying capacity) of the main clade \cr
   \code{pars1[4]} corresponds to lambda_M (speciation rate) of the subclade \cr
   \code{pars1[5]} corresponds to mu_S (extinction rate) of the subclade \cr
   \code{pars1[6]} corresponds to K_S (clade-level carrying capacity) of the subclade \cr
   \code{pars1[7]} corresponds to t_d (the time of decoupling)
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
    \code{pars2[4]} sets the time of splitting of the branch that will decouple
}
  \item{brtsM}{A set of branching times of the main clade in the phylogeny, all positive }
  \item{brtsS}{A set of branching times of the subclade in the phylogeny, all positive }
  \item{missnumspec}{The number of species that are missing in the phylogeny } 
}
\value{The loglikelihood}
\references{
- Etienne et al. 2012. Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439
- Etienne & Haegeman 2012. Am. Nat. In press. \cr}
\author{Rampal S. Etienne & Bart Haegeman }
\seealso{
 \code{\link{dd_KI_ML}},
 \code{\link{dd_loglik}}
 \code{\link{dd_SR_loglik}}
}
\examples{
\donttest{
pars1 = c(0.25,0.12,25.51,1.0,0.16,8.61,9.8)
pars2 = c(200,1,TRUE,18.8)
missnumspec = 0
brtsM = c(25.2,24.6,24.0,22.5,21.7,20.4,19.9,19.7,18.8,17.1,15.8,11.8,9.7,8.9,5.7,5.2)
brtsS = c(9.6,8.6,7.4,4.9,2.5)
dd_KI_loglik(pars1,pars2,brtsM,brtsS,missnumspec)
}
}
\keyword{models}

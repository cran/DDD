\name{dd_SR_ML}
\alias{dd_SR_ML}
\title{ Maximization of the loglikelihood under a diversity-dependent diversification model with a shift in the parameters at time t = tshift
}
\description{
   This function computes the maximum likelihood estimates of the parameters of a
   diversity-dependent diversification model with shifting parameters at time t = tshift for a given set of phylogenetic branching times.
   It also outputs the corresponding loglikelihood that can be used in model comparisons.
}
\usage{
dd_SR_ML(brts, initparsopt = c(0.2, 0.1, 2 * (1 + length(brts)), 2 * (1 + length(brts)),
   max(brts)/2), parsfix = NULL, idparsopt = c(1:3, 6:7), idparsfix = NULL,
   idparsnoshift = (1:7)[c(-idparsopt, (-1)^(length(idparsfix) != 0) * idparsfix)],
   res = 10 * (1 + length(brts)), ddmodel = 1, missnumspec = 0, cond = TRUE)
}
\arguments{
  \item{brts}{A set of branching times of a phylogeny, all positive, from present to past }
  \item{initparsopt}{The initial values of the parameters that must be optimized }
  \item{parsfix}{The values of the parameters that should not be optimized }
  \item{idparsopt}{
   The ids of the parameters that must be optimized, e.g. 1:7 for all parameters.
   The ids are defined as follows: \cr
   id == 1 corresponds to lambda (speciation rate) before the shift \cr
   id == 2 corresponds to mu (extinction rate) before the shift \cr
   id == 3 corresponds to K (clade-level carrying capacity) before the shift \cr
   id == 4 corresponds to lambda (speciation rate) after the shift \cr
   id == 5 corresponds to mu (extinction rate) after the shift \cr
   id == 6 corresponds to K (clade-level carrying capacity) after the shift \cr
   id == 7 corresponds to tshift (the time of shift)
   }
  \item{idparsfix}{
The ids of the parameters that should not be optimized, e.g. c(1,3,4,6) if lambda and K
should not be optimized, but only mu. In that case idparsopt must be c(2,5,7). The default is to fix all parameters
not specified in idparsopt.
  }
  \item{idparsnoshift}{ The ids of the parameters that should not shift; 
This can only apply to ids 4, 5 and 6, e.g. idparsnoshift = c(4,5) means that lambda and mu have the same values before and after tshift }
  \item{res}{ sets the maximum number of species for which a probability must be computed, must be larger than 1 + length(brts) }
  \item{ddmodel}{ sets the model of diversity-dependence: \cr
   ddmodel == 1 : linear dependence in speciation rate \cr
   ddmodel == 2 : exponential dependence in speciation rate \cr
   ddmodel == 3 : linear dependence in extinction rate \cr
   ddmodel == 4 : exponential dependence in extinction rate
  }
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny } 
  \item{cond}{Conditioning on non-extinction of the phylogeny (TRUE or FALSE) }
}
\details{
The output is a dataframe where $par and $fvalues are the first elements.
There are more elements that are all from the optimization algorithm, optimx.
The computed loglikelihood contains the factor
N_q,m = (q-1)/(q + m - 1) (q + m + 1)!/((q + 1)! m!)
where q is the number of species in the phylogeny and m is the number of missing species,
as explained in the supplementary material to Etienne et al. 2011.
}
\value{
 \item{pars}{gives the list of parameter estimates}
 \item{fvalues}{gives the maximum loglikelihood}
}
\references{
- Etienne & Haegeman 2012. In review. \cr
- Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439 }
\author{Rampal S. Etienne & Bart Haegeman }
\note{The optimization may get trapped in local optima. Try different starting values to search for the global optimum.
}
\seealso{
 \code{\link{dd_SR_loglik}},
 \code{\link{optimx}}
}
\examples{
\donttest{
cat("Estimating parameters for a set of branching times brts with the default settings:")
brts = 1:20
dd_SR_ML(brts)
}
\dontshow{
dd_SR_ML(1:3,c(11.676626,0.001088,2.183581,3.221303,1.006468),idparsopt = c(1:3,6,7))
}
}
\keyword{ models }

\name{dd_SR_ML}
\alias{dd_SR_ML}
\title{ Maximization of the loglikelihood under a diversity-dependent diversification model with a shift in the parameters
}
\description{
   This function computes the maximum likelihood estimates of the parameters of a
   diversity-dependent diversification model with shifting parameters at time t = tshift for a given set of phylogenetic branching times.
   It also outputs the corresponding loglikelihood that can be used in model comparisons.
}
\usage{
dd_SR_ML(
   brts,
   initparsopt = c(0.5, 0.1, 2 * (1 + length(brts) + missnumspec),
2 * (1 + length(brts) + missnumspec), max(brts)/2),
   parsfix = NULL,
   idparsopt = c(1:3, 6:7),
   idparsfix = NULL,
   idparsnoshift = (1:7)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)],
   res = 10 * (1 + length(brts) + missnumspec),
   ddmodel = 1,
   missnumspec = 0,
   cond = 1,
   btorph = 1,
   soc = 2,
   allbp = FALSE,
   tol = c(1E-3,1E-4,1E-6),
   maxiter = 1000 * round((1.25)^length(idparsopt))
   )
}
\arguments{
  \item{brts}{A set of branching times of a phylogeny, all positive }
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
   ddmodel == 2.1 : variant of exponential dependence in speciation rate with offset at infinity\cr
   ddmodel == 2.2 : 1/n dependence in speciation rate\cr
   ddmodel == 3 : linear dependence in extinction rate \cr
   ddmodel == 4 : exponential dependence in extinction rate \cr
   ddmodel == 4.1 : variant of exponential dependence in extinction rate with offset at infinity
   ddmodel == 4.2 : 1/n dependence in extinction rate with offset at infinity
  }
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny } 
  \item{cond}{ Conditioning: \cr
    cond == 0 : no conditioning \cr
    cond == 1 : conditioning on non-extinction of the phylogeny \cr
    cond == 2 : conditioning on non-extinction of the phylogeny and on the total number of extant taxa (including missing species) \cr
    cond == 3 : conditioning on the total number of extant taxa (including missing species) \cr
(including missing species) \cr
   Note: cond == 3 assumes a uniform prior on stem age, as is the standard in constant-rate birth-death models, see
e.g. D. Aldous & L. Popovic 2004. Adv. Appl. Prob. 37: 1094-1115 and T. Stadler 2009. J. Theor. Biol. 261: 58-66.
This conditioning turns out, for the diversity-dependent model, to provide the least-biased parameter estimates when extinction is low, but may still contain considerable bias if extinction is high. The default value has therefore been changed from 1 to 3 starting with DDD version 2.3. 
}
  \item{btorph}{Sets whether the likelihood is for the branching times (0) or the phylogeny (1)}
  \item{soc}{Sets whether stem or crown age should be used (1 or 2) }
  \item{allbp}{Sets whether a search should be done with various initial conditions, with tshift at each of the branching points (TRUE/FALSE)}
  \item{tol}{Sets the tolerances in the optimization. Consists of: \cr
   reltolx = relative tolerance of parameter values in optimization \cr
   reltolf = relative tolerance of function value in optimization \cr
   abstolx = absolute tolerance of parameter values in optimization
  }
  \item{maxiter}{Sets the maximum number of iterations in the optimization}
}
\details{
The output is a dataframe containing estimated parameters and maximum loglikelihood.
The computed loglikelihood contains the factor
q! m!/(q + m)!
where q is the number of species in the phylogeny and m is the number of missing species,
as explained in the supplementary material to Etienne et al. 2012.
}
\value{
 \item{lambda_1}{ gives the maximum likelihood estimate of lambda before the shift}
 \item{mu_1}{ gives the maximum likelihood estimate of mu before the shift}
 \item{K_1}{ gives the maximum likelihood estimate of K before the shift}
 \item{lambda_2}{ gives the maximum likelihood estimate of lambda after the shift}
 \item{mu_2}{ gives the maximum likelihood estimate of mu after the shift}
 \item{K_2}{ gives the maximum likelihood estimate of K after the shift}
 \item{t_shift}{ gives the time of the shift}
 \item{loglik}{ gives the maximum loglikelihood}
 \item{df}{ gives the number of estimated parameters, i.e. degrees of feedom}
 \item{conv}{ gives a message on convergence of optimization; conv = 0 means convergence}}
\references{
- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
}
\author{Rampal S. Etienne & Bart Haegeman }
\note{The optimization may get trapped in local optima. Try different starting values to search for the global optimum.
}
\seealso{
 \code{\link{dd_SR_loglik}},
 \code{\link{dd_ML}},
 \code{\link{dd_KI_ML}},
}
\examples{
\donttest{
cat("Estimating parameters for a set of branching times brts with the default settings:")
brts = 1:20
dd_SR_ML(brts)
}
\dontshow{
dd_SR_ML(1:4,initparsopt = c(2.23,0.0002,1.58,4.41,3.01),idparsopt = c(1:3,6,7),tol = c(0.01,0.01,0.01))
}
}
\keyword{ models }

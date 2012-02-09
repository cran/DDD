\name{dd_ML}
\alias{dd_ML}
\title{Maximization of the loglikelihood under a diversity-dependent diversification model
}
\description{
   This function computes the maximum likelihood estimates of the parameters of a
   diversity-dependent diversification model for a given set of phylogenetic branching times.
   It also outputs the corresponding loglikelihood that can be used in model comparisons.
}
\usage{
dd_ML(brts, initparsopt = c(0.2, 0.1, 2 * length(brts)), idparsopt = 1:3,
   idparsfix = (1:3)[-idparsopt], parsfix = c(0.2, 0.1, 2 * length(brts))[-idparsopt],
   res = 10 * (1 + length(brts)), ddmodel = 1, missnumspec = 0, cond = TRUE)
}
\arguments{
  \item{brts}{A set of branching times of a phylogeny, all positive, from present to past }
  \item{initparsopt}{The initial values of the parameters that must be optimized}
  \item{idparsopt}{
   The ids of the parameters that must be optimized, e.g. 1:3 for all parameters.
   The ids are defined as follows: \cr
   id == 1 corresponds to lambda (speciation rate) \cr
   id == 2 corresponds to mu (extinction rate) \cr
   id == 3 corresponds to K (clade-level carrying capacity)
   }
  \item{idparsfix}{
The ids of the parameters that should not be optimized, e.g. c(1,3) if lambda and K
should not be optimized, but only mu. In that case idparsopt must be 2. The default is to fix all parameters
not specified in idparsopt. 
  }
  \item{parsfix}{ The values of the parameters that should not be optimized }
  \item{res}{ Sets the maximum number of species for which a probability must be computed, must be larger than 1 + length(brts) }
  \item{ddmodel}{ Sets the model of diversity-dependence: \cr
   ddmodel == 1 : linear dependence in speciation rate \cr
   ddmodel == 2 : exponential dependence in speciation rate \cr
   ddmodel == 3 : linear dependence in extinction rate \cr
   ddmodel == 4 : exponential dependence in extinction rate \cr
  }
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny } 
  \item{cond}{Conditioning on non-extinction of the phylogeny (TRUE or FALSE) }
}
\details{
The output is a dataframe where $par and $fvalues are the first elements.
There are more elements that are all from the optimization algorithm, optimx.
The computed loglikelihood contains the factor
N_q,m = (q - 1)/(q + m - 1) (q + m + 1)!/((q + 1)! m!)
where q is the number of species in the phylogeny and m is the number of missing species,
as explained in the supplementary material to Etienne et al. 2011.
}
\value{
 \item{pars}{ gives the list of parameter estimates}
 \item{fvalues}{ gives the maximum loglikelihood}
}
\references{ Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439}
\author{ Rampal S. Etienne & Bart Haegeman }

\seealso{
 \code{\link{dd_loglik}},
 \code{\link{optimx}}
}
\examples{
\donttest{
cat("Estimating parameters for a set of branching times brts with the default settings:")
brts = 1:20
dd_ML(brts)

cat("Estimating the intrinsic speciation rate lambda and the carrying capacity K")
cat("for a fixed extinction rate of 0.1:")
brts = 1:20
dd_ML(brts,initparsopt = c(0.2,100), idparsopt = c(1,3), parsfix = 0.1)

cat("This will estimate all parameters under the linear diversity-dependent")
cat("extinction model:")
brts = 1:20
dd_ML(brts, ddmodel = 3)

cat("This will estimate all parameters under the default model")
cat("when there are 2 missing species:")
brts = 1:20
dd_ML(brts, missnumspec = 2)
}
\dontshow{
dd_ML(1:3,c(1.151160,0.000001,3.997702))
}
}
\keyword{ models }

\name{dd_loglik_choosepar}
\alias{dd_loglik_choosepar}
\title{ Auxiliary function to dd_ML }
\description{
This function reorders the parameters so that dd_loglik can be called
with the correct parameters. It is called from dd_ML.
}
\usage{
dd_loglik_choosepar(trparsopt, trparsfix, idparsopt, idparsfix, pars2, brts, missnumspec)
}
\arguments{
  \item{trparsopt}{
These are the initial values of the transformed parameters that must be optimized. The transformation
is trpars = pars/(1+pars). See dd_ML.
}
  \item{trparsfix}{
These are the values of the transformed parameters that must not be optimized. The transformation
is trpars = pars/(1+pars). See dd_ML.
}
  \item{idparsopt}{
These are the ids of the parameters to be optimized. For example, to optimize all parameters, idparsopt = 1:3
}
  \item{idparsfix}{
These are the ids of the parameters to be fixed. For example, idparsfix = c(1,3). In that case idparsopt must be 2
}
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
\details{ This function should not be called directly from R. }
\value{ the loglikelihood }
\references{Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439}
\author{Rampal S. Etienne & Bart Haegeman}

\seealso{\code{\link{dd_loglik}}}

\examples{
dd_loglik_choosepar(c(0.15,0.10,0.99),NULL,1:3,NULL,c(100,1),1:10,0)
}
\keyword{ models }


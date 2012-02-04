\name{dd_SR_loglik_choosepar}
\alias{dd_SR_loglik_choosepar}
\title{
Auxiliary function to dd_SR_ML
}
\description{
This function reorders the parameters so that dd_loglik can be called
with the correct parameters. It is called from dd_SR_ML.
}
\usage{
dd_SR_loglik_choosepar(trparsopt, trparsfix, idparsopt, idparsfix, idparsnoshift, pars2, brts, missnumspec)
}
\arguments{
  \item{trparsopt}{
These are the initial values of the transformed parameters that must be optimized. The transformation
is trpars = pars/(1+pars). See dd_SR_ML.
}
  \item{trparsfix}{
These are the values of the transformed parameters that must not be optimized. The transformation
is trpars = pars/(1+pars). See dd_SR_ML.
}
  \item{idparsopt}{
These are the ids of the parameters to be optimized. For example, to optimize all parameters, idparsopt = 1:7
}
  \item{idparsfix}{
These are the ids of the parameters to be fixed. For example, idparsfix = c(1,3).
}
  \item{idparsnoshift}{
These are the ids of the parameters that are not to be fixed, for example, idparsnoshift = 4:5.
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
\references{ Etienne & Haegeman 2012. In review.}
\author{Rampal S. Etienne & Bart Haegeman}

\seealso{
\code{\link{dd_SR_loglik}}
\code{\link{dd_SR_ML}}

}

\examples{
dd_SR_loglik_choosepar(c(0.15,0.10,0.99,0.99,0.8),NULL,c(1:3,6:7),NULL,4:5,c(100,1),1:10,0)
}
\keyword{ models }



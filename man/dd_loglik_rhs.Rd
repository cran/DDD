\name{dd_loglik_rhs}
\alias{dd_loglik_rhs}
\title{ Auxiliary function to dd_loglik}
\description{
This function computes the right-hand side of the master equation for the
diversity-dependent model. It is called from dd_loglik and from dd_SR_loglik.
}
\usage{
dd_loglik_rhs(t, x, pars)
}
\arguments{
  \item{t}{ time }
  \item{x}{ the probability vector at time t }
  \item{pars}{ the parameters:
  pars[1] = lambda, intrinsic speciation rate
  pars[2] = mu, extinction rate
  pars[3] = K, clade-level carrying capacity
  pars[4] = kk, the number of lineages in the phylogeny at time t
  pars[5] = ddep, the divcersity-dependent model, see dd_loglik
}
}
\details{ This function should not be called directly from R. }
\value{ The function outputs the incremental probability vector needed in the master equation}
\references{ Etienne et al. 2011, Proc. Roy. Soc. B, doi: 10.1098/rspb.2011.1439}
\author{ Rampal S. Etienne & Bart Haegeman }

\seealso{
\code{\link{dd_loglik}}
\code{\link{dd_SR_loglik}}

}
\examples{ dd_loglik_rhs(-3,(1:10)/sum(1:10),c(0.2,0.1,100,3,1)) }
\keyword{ models }


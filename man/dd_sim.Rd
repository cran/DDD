\name{dd_sim}
\alias{dd_sim}
\title{ Function to simulate the diversity-dependent diversification process }
\description{ Simulating the diversity-dependent diversification process }
\usage{ dd_sim(pars, age, ddmodel = 1) }
\arguments{
  \item{pars}{Vector of parameters: \cr \cr
    \code{pars[1]} corresponds to lambda (speciation rate) \cr
    \code{pars[2]} corresponds to mu (extinction rate)  \cr
    \code{pars[3]} corresponds to K (clade-level carrying capacity) 
  }
  \item{age}{ Sets the crown age for the simulation }
  \item{ddmodel}{ Sets the model of diversity-dependence: \cr
   \code{ddmodel == 1} : linear dependence in speciation rate \cr
   \code{ddmodel == 2} : exponential dependence in speciation rate \cr
   \code{ddmodel == 3} : linear dependence in extinction rate \cr
   \code{ddmodel == 4} : exponential dependence in extinction rate \cr
   \code{ddmodel == 5} : linear dependence in speciation and extinction rate\cr
  }
}
\value{
  \item{ out }{ A list with the following three elements:
   The first element is the tree of extant species in phylo format \cr
   The second element is the tree of all species, including extinct species, in phylo format \cr
   The third element is a matrix of all species where \cr
   - the first column is the time at which a species is born \cr
   - the second column is the parent of the species \cr
   - the third column is the daughter species \cr
   - the fourth column is the time of extinction of the species \cr
   If the fourth element equals -1, then the species is still extant.
  }
}
\references{
- Etienne et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr
- Etienne & Haegeman 2012. Am. Nat. Am. Nat. 180: E75-E89, doi: 10.1086/667574}
\author{ Rampal S. Etienne }
\examples{ dd_sim(c(0.2,0.1,20),10) }
\keyword{models}

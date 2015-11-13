\name{DDD-package}
\alias{DDD-package}
\alias{DDD}
\docType{package}
\title{Diversity Dependent Diversification}
\description{
This package is a likelihood-based statistical package to estimate parameters
under diversity-dependent diversification models.\cr \cr
   New in version 0.3:\cr
   - Models with shifts in the parameters at a time tshift.\cr\cr
   New in version 0.4:\cr
   - Bug fix of prefactor in likelihood - bug did not affect parameter estimates.\cr
   - There is an option to condition or not on non-extinction of the phylogeny.\cr
   - The functions dd_loglik and dd_SR_loglik now accept untransformed parameters;
   the transformation is done in the auxiliary functions.\cr
   - Improved help pages.\cr
   - Auxiliary functions no longer in help pages.\cr\cr
   New in version 1.0:\cr
   - Another bug fix of prefactor in likelihood - bug did not affect parameter estimates.\cr
   - Models with key innovations defined as decoupling of diversity-dependent dynamics.\cr
   - Option to compute the likelihood of a set of branching times or of the phylogeny;
   these differ by a constant factor.\cr\cr
   New in version 1.1:\cr
   - A model with diversity-dependence in speciation and extinction rate, currently only
   available for the likelihood without shifts or decoupling.\cr\cr
   New in version 1.2:\cr
   - Possibility of automatically exploring multiple initial conditons for the SR models\cr
   - Suppressing unnecessary output\cr
   - Flushing computed likelihood values immediately (even with buffer on)\cr
   - Bug fix in default value of the initial time of decoupling\cr\cr
   New in version 1.3:\cr
   - Fix of default length of idparsopt\cr
   - Fix of conversion problems for small branching times\cr
   - Default initial parameter values in dd_ML based on birth-death model\cr\cr
   New in version 1.4:\cr
   - Bug fix for when there is a large number of missing species\cr\cr
   New in version 1.5:\cr
   - Bug fix for returned output of functions when there are fixed parameters.\cr
   Textual output of DDD was correct.\cr
   - Return to Nelder-Mead optimization for dd_ML (was set to "nlm" in versions 1.3 & 1.4)\cr
   - New format for output off all ML-functions (dd_ML, dd_SR_ML, dd_KI_ML)\cr\cr
   New in version 1.6 (not on CRAN)\cr
   - Bug fix in returned result when optimization does not converge.\cr
   - Works also when only branching times are entered in ascending order
   (previously error was reported)\cr\cr
   New in version 1.7\cr
   - Cleaning up code, fixing various small bugs\cr
   - New optimization algorithm because of occasional strange results with optim package\cr\cr
   New in version 1.8\cr
   - Bug fix for missing number of species in KI model\cr
   - Bug fix for optimization of one parameter (due to R's inconsistent handling of matrices)\cr
   - Model 5 (dd_ML and dd_loglik) now also allows values of -1 < r < 0;
   i.e. when extinction also decreases with diversity but less so than speciation\cr\cr
   New in version 1.9\cr
   - Bug fix in optimizer (did not stop when likelihood was positive)\cr
   - Bug fix for large numbers of missing species\cr
   - Conv = -1 is now always reported when optimization is not started or does not converge\cr
   - An error is reported when parameter values cause numerical problems (NA or NaN) in likelihood computation,
   e.g. when parameter values are very high; loglik is then set to -Inf\cr
   - For dd_KI_loglik and dd_KI_ML: if the number of missing species in main clade and subclade
   is known, this can now be specified\cr
   - Fixing K at Inf is now possible in SR and KI models\cr\cr
   New in version 1.10\cr
   - All non-convergences now get conv = -1 (this was conv = 1 for true nonconvergence and conv = -1 for optimization
   that did not start because of low likelihood values)\cr
   - Likelihoods can also be computed for a tree with a stem age rather than only a crown age\cr\cr
   New in version 1.11\cr
   - Various small tweaks to make the code run faster and some bug fixes\cr
   New in version 1.12\cr
   - Allows conditioning on survival of crown linages AND on the number of extant species by setting cond = 2\cr\cr
   New in version 1.13\cr
   - Fast computation of dd_loglik for the case K = Inf\cr\cr
   New in version 1.14\cr
   - Fast computation of dd_loglik for the case mu = 0\cr
   - Allows time-dependence instead of diversity-dependence; exponential decline in speciation rate and/or extinction\cr
   - Start of time-dependence can be specified to allow synchronisation across multiple clades\cr
   - Removed link to LASER\cr\cr
   New in version 2.0\cr
   - Includes function to simulate trees under diversity-dependence\cr\cr
   New in version 2.1\cr
   - Fixes numerical problems with the convolve-function (used to deal with missing species)\cr
   - Fixes a bug in dd_KI_loglik (ddep unknown) \cr
   - Changed third element of the output of dd_sim, such that all speciation and extinction times are in Mya\cr\cr
   New in version 2.2\cr
   - Allows conditioning on the number of taxa only for dd_loglik, asssuming a uniform prior on stem age; this is now the default\cr
   - An additional model of time-dependence, following deterministic logistic equation\cr\cr
   New in version 2.3\cr
   - Allows conditioning on the number of taxa only for dd_SR_loglik, assuming a uniform prior on stem age; this is now the default\cr
   - Fix of bug for stem age (soc = 1), introduced in v2.2.\cr
   - Allows Rabosky's exponential model, with la(n) = la_0 * n^(-x(la_0,mu,K)) as ddmodel = 2.1, mu(n) = mu_0 * n^x(la, mu_0, K) as ddmodel = 4.1\cr
   - Allows specific case of Rabosky's exponential model, with la(n) = la_0 * n^(-1), as ddmodel = 2.2, and mu(n) = mu_0 * n, as ddmodel = 4.2\cr\cr
   New in version 2.4\cr
   - Fix of bug in dd_KI_loglik ('m' was not found)\cr\cr
   New in version 2.5\cr
   - Fix of bug in bd_loglik where NA was output to screen for a nonexistent parameter under diversity-dependence without extinction\cr
   - Fix of bug in dd_KI_loglik when m is a two-element vector\cr
   - Fix of bug in dd_KI_loglik, introduced in v2.3\cr\cr
   New in version 2.6\cr
   - Allows linear dependence of speciation with parameter K' (= diversity where speciation = 0)\cr\cr
   New in version 2.7\cr
   - Fixes bugs in dd_KI_loglik, for cond = 0 introduced in v2.3, and for cond = 1 introduced in v2.0\cr
   - Solves (some) memory problems with dd_KI_ML and dd_KI_loglik when cond = 1\cr
   - Provides option to set loglikelihood to -Inf if optimization did not converge\cr\cr 
   New in version 3.0\cr
   - Includes a model of macroevolutionary succession in which there is no decoupling of diversity-dependence but parameters may change\cr
   - Computes LRb, a bootstrap version of the likelihood ratio\cr
   - Allows use of alternative ODE solvers\cr   
   - Allows use of subplex as optimization algorithm.\cr\cr
   Forthcoming:\cr
   - High-precision arithmetic in numerically unstable regions of bd_loglik\cr
}

\details{
\tabular{ll}{
Package: \tab DDD\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2015-11-13\cr
License: \tab GPL 3.0\cr
}
bd_loglik compues the likelihood under the constant-rate birth-death model of several specific 
time-dependent birth-death models.

bd_ML finds the parameters that maximizes the likelihood computed by bd_loglik.

dd_loglik computes the likelihood of any of 4 different models of diversity-dependent
diversification, given a set of parameters and a data set of phylogenetic branching times.

dd_ML finds the parameters that maximizes the likelihood computed by dd_loglik.

dd_SR_loglik is the same as dd_loglik except that it allows a shift in the parameters at
time tshift.

dd_SR_ML finds the parameters that maximizes the likelihood computed by dd_SR_loglik.

dd_KI_loglik is the same as dd_loglik except that it allows decoupling of a subclade's dynamics
from the main clade's dynamics potentially accompanied by a shift in parameters at time t_d.

dd_KI_ML finds the parameters that maximizes the likelihood computed by dd_KI_loglik.

dd_MS_loglik is the same as dd_KI_loglik except that the subclade's dynamics do not decouple
from the main clade's dynamics but they do get different parameters at time t_d.

dd_MS_ML finds the parameters that maximizes the likelihood computed by dd_MS_loglik.

dd_sim simulates the diversity-dependent process

dd_LR performs a bootstrap likelihood ratio test
}
\author{
Rampal S. Etienne & Bart Haegeman\cr
Maintainer: Rampal S. Etienne (r.s.etienne@rug.nl)
}
\references{
- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
}
\keyword{ models }
\seealso{
\code{TESS}
}

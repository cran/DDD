dd_SR_ML = function(brts,initparsopt=c(0.2,0.1,2*(1+length(brts)),2*(1+length(brts)),max(brts)/2),parsfix = NULL,idparsopt=c(1:3,6:7),idparsfix=NULL,idparsnoshift=(1:7)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)],res=10*(1 + length(brts)),ddmodel=1,missnumspec=0,cond = TRUE)
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - idparsopt contains the ids of the parameters to be optimized, e.g. to optimize la, mu, K, K2 and tshift idparsopt = c(1,2,3,6,7)
# - initparsopt contains the starting values of the parameters to be optimized
# - idparsfix contains the ids of the parameters that are fixed and must not be optimized
# - parsfix contains the values of the fixed parameters
# - idparsnoshift contains the ids of the parameters la2, mu2 and K2 that do not shift, i.e. that need to be set equal to la, mu and K
# - pars[1] = la = (initial) speciation rate before shift
# - pars[2] = mu = extinction rate before shift
# - pars[3] = K = carrying capacity before shift
# - pars[4] = la2 = (initial) speciation rate after shift
# - pars[5] = mu2 = extinction rate after shift
# - pars[6] = K2 = carrying capacity after shift
# - pars[7] = tshift = time of shift
# - res = resolution of the method; res should be larger than the total number of species
# - ddmodel = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
# - missnumspec = number of missing species    
# - cond = conditioning on non-extinction of the phylogeny

options(warn=-1)
if(is.numeric(brts) == FALSE) { cat("The branching times should be numeric") } else {
idparsnoshift = sort(idparsnoshift)
idpars = sort(c(idparsopt,idparsfix,idparsnoshift))
if(sum(idpars == (1:7)) != 7) {cat("The parameters to be optimized, fixed and not shifted are incoherent.") } else {
namepars = c("la","mu","K","la2","mu2","K2","tshift")
if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
cat("You are optimizing",optstr,"\n")
if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
cat("You are fixing",fixstr,"\n")
if(length(namepars[idparsnoshift]) == 0) { noshiftstr = "nothing" } else { noshiftstr = namepars[idparsnoshift] }
cat("You are not shifting",noshiftstr,"\n")
trparsopt = initparsopt/(1 + initparsopt)
trparsfix = parsfix/(1 + parsfix)
cat("Optimizing the likelihood - this may take a while.","\n")
out = optimx(trparsopt,dd_SR_loglik_choosepar,hess=NULL,method = "Nelder-Mead",control = list(maximize = TRUE,abstol = 1E-6,reltol = 1E-6,trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,brts = brts,pars2 = c(res,ddmodel,cond),missnumspec = missnumspec)
MLtrpars = unlist(out$par)
MLpars = MLtrpars/(1-MLtrpars)
out$par = list(MLpars)
MLpars1 = rep(0,7)
MLpars1[idparsopt] = MLpars
if(length(idparsfix) != 0) {MLpars1[idparsfix] = parsfix }
if(length(idparsnoshift) != 0) { MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3] }
s1 = sprintf('Maximum likelihood parameter estimates: %f %f %f %f %f %f %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7])
s2 = sprintf('Maxmimum loglikelihood: %f',out$fvalues)
cat("\n",s1,"\n",s2,"\n")
return(out)
}
}
}
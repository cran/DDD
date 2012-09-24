dd_KI_ML = function(brtsM,brtsS,tsplit,initparsopt=c(0.5,0.1,2*(1+length(brtsM)),2*(1+length(brtsS)),(tsplit + max(brtsS))/2),parsfix = NULL,idparsopt=c(1:3,6:7),idparsfix=NULL,idparsnoshift=(1:7)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)],res=10*(1 + length(c(brtsM,brtsS)) + missnumspec),ddmodel=1,missnumspec=0,cond = TRUE)
{
# brtsM, brtsS = branching times of main clade and subclade (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - tsplit = the branching time where the subclade branches off from the main clade
# - idparsopt contains the ids of the parameters to be optimized, e.g. to optimize la, mu, K, K2 and tshift idparsopt = c(1,2,3,6,7)
# - initparsopt contains the starting values of the parameters to be optimized
# - idparsfix contains the ids of the parameters that are fixed and must not be optimized
# - parsfix contains the values of the fixed parameters
# - idparsnoshift contains the ids of the parameters la2, mu2 and K2 that do not shift, i.e. that need to be set equal to la, mu and K
# - pars[1] = la_M = (initial) speciation rate in main clade
# - pars[2] = mu_M = extinction rate in main clade
# - pars[3] = K_M = carrying capacity in main clade
# - pars[4] = la_S = (initial) speciation rate in subclade
# - pars[5] = mu_S = extinction rate in subclade
# - pars[6] = K_S = carrying capacity in subclade
# - pars[7] = t_d = time of decoupling
# - res = resolution of the method; res should be larger than the total number of species
# - ddmodel = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
# - missnumspec = number of missing species    
# - cond = conditioning on non-extinction of the phylogeny

options(warn=-1)
if(is.numeric(brtsM) == FALSE || is.numeric(brtsS) == FALSE) { cat("The branching times should be numeric") } else {
idparsnoshift = sort(idparsnoshift)
idpars = sort(c(idparsopt,idparsfix,idparsnoshift))
if(sum(idpars == (1:7)) != 7) {cat("The parameters to be optimized, fixed and not shifted are incoherent.") } else {
namepars = c("la_M","mu_M","K_M","la_S","mu_S","K_S","t_d")
if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
cat("You are optimizing",optstr,"\n")
if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
cat("You are fixing",fixstr,"\n")
if(length(namepars[idparsnoshift]) == 0) { noshiftstr = "anything" } else { noshiftstr = namepars[idparsnoshift] }
cat("You are not shifting",noshiftstr,"\n")
trparsopt = initparsopt/(1 + initparsopt)
trparsfix = parsfix/(1 + parsfix)
cat("Optimizing the likelihood - this may take a while.","\n")
out = optimx2(trparsopt,dd_KI_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = 1E-6,reltol = 1E-6,trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,brtsM = brtsM, brtsS = brtsS, pars2 = c(res,ddmodel,cond,tsplit),missnumspec = missnumspec)
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
invisible(out)
}
}
}
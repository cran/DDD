dd_ML = function(brts,initparsopt=if(ddmodel < 5) {c(0.1+bd(as.numeric(brts))$r1/(1-bd(as.numeric(brts))$a),0.1,2*(length(brts) + missnumspec))} else {c(0.1+bd(as.numeric(brts))$r1/(1-bd(as.numeric(brts))$a),0.1,2*(length(brts) + missnumspec),0.01)},idparsopt = 1:length(initparsopt),idparsfix = (1:(3 + (ddmodel == 5)))[-idparsopt],parsfix = (ddmodel < 5) * c(0.2,0.1,2*(length(brts) + missnumspec))[-idparsopt] + (ddmodel == 5) * c(0.2,0.1,2*(length(brts) + missnumspec),0)[-idparsopt],res=10*(1+length(brts)+missnumspec),ddmodel=1,missnumspec=0,cond=TRUE, btorph = 1)
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - initpars[1] = la = (initial) speciation rate
# - initpars[2] = mu = extinction rate
# - initpars[3] = K = carrying capacity
# - res = resolution of the method; res should be larger than the total number of species
# - ddmodel = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
#  . ddmodel == 5 : linear dependence in speciation rate and in extinction rate
# - missnumspec = number of missing species    
# - cond = conditioning on non-extinction of the phylogeny
# - btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species

options(warn=-1)
if(is.numeric(brts) == FALSE) { cat("The branching times should be numeric") } else {
idpars = sort(c(idparsopt,idparsfix))
if(sum(idpars == (1:3)) != 3) {cat("The parameters to be optimized and fixed are incoherent.") } else {
namepars = c("lambda","mu","K")
if(ddmodel == 5) {namepars = namepars = c("lambda","mu","K","r")}
if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
cat("You are optimizing",optstr,"\n")
if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
cat("You are fixing",fixstr,"\n")
trparsopt = initparsopt/(1 + initparsopt)
trparsfix = parsfix/(1 + parsfix)
trparsfix[parsfix == Inf] = 1
cat("Optimizing the likelihood - this may take a while.","\n")
out = optimx2(trparsopt,dd_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = 1E-4,reltol = 1E-6,trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,brts = brts,pars2 = c(res,ddmodel,cond,btorph),missnumspec = missnumspec)
if(out$conv > 0) {cat("Optimization has not converged. Try again with different starting values.\n")} else {
MLtrpars = unlist(out$par)
MLpars = MLtrpars/(1-MLtrpars)
MLpars1 = rep(0,3)
if(ddmodel == 5) {MLpars1 = rep(0,4)}
if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
MLpars1[idparsopt] = MLpars
if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
out2 = data.frame(row.names = "results",lambda = MLpars1[1],mu = MLpars1[2],K = MLpars1[3], loglik = unlist(out$fvalues), df = length(initparsopt), conv = unlist(out$conv))
s1 = sprintf('Maximum likelihood parameter estimates: lambda: %f, mu: %f, K: %f',MLpars1[1],MLpars1[2],MLpars1[3])
if(ddmodel == 5)
{
   s1 = sprintf('%s, r: %f',s1,MLpars1[4])
   out2 = data.frame(row.names = "results",lambda = MLpars1[1],mu = MLpars1[2],K = MLpars1[3], r = MLpars1[4], loglik = unlist(out$fvalues), df = length(initparsopt), conv = unlist(out$conv))
}
s2 = sprintf('Maximum loglikelihood: %f',out$fvalues)
cat("\n",s1,"\n",s2,"\n")
}
invisible(out2)
}
}
}
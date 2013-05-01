bd_ML = function(brts, initparsopt = c(0.1,0.05), idparsopt = 1:length(initparsopt), idparsfix = (1:4)[-idparsopt], parsfix = c(0,0), missnumspec = 0, tdmodel = 0, cond = 1, btorph = 1, soc = 2, tol = c(1E-3, 1E-4, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)))
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - initpars[1] = la0 = (initial) speciation rate
# - initpars[2] = mu0 = (initial) extinction rate
# - initpars[3] = la1 = exponential decline parameter of speciation rate
# - initpars[4] = mu1 = exponential decline parameter of extinction rate
# - res = resolution of the method; res should be larger than the total number of species
# - missnumspec = number of missing species    
# - cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
#  . cond == 2 : conditioning on non-extinction of the phylogeny and on the total number of extant taxa (including missing species)
# - btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - tol = tolerance in optimization
#  . reltolx = relative tolerance of parameter values in optimization
#  . reltolf = relative tolerance of function value in optimization
#  . abstolx = absolute tolerance of parameter values in optimization
# - maxiter = the maximum number of iterations in the optimization

options(warn=-1)
brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
if(is.numeric(brts) == FALSE)
{
   cat("The branching times should be numeric.\n")
   out2 = data.frame(lambda0 = -1,mu0 = -1,lambda1 = -1, mu1 = -1, loglik = -1, df = -1, conv = -1)
} else {
idpars = sort(c(idparsopt,idparsfix))
if((sum(idpars == (1:4)) != 4) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
{
   cat("The parameters to be optimized and/or fixed are incoherent.\n")
   out2 = data.frame(lambda0 = -1,mu0 = -1,lambda1 = -1, mu1 = -1, loglik = -1, df = -1, conv = -1)
} else {
namepars = c("lambda0","mu0","lambda1","mu1")
if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
cat("You are optimizing",optstr,"\n")
if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
cat("You are fixing",fixstr,"\n")
cat("Optimizing the likelihood - this may take a while.","\n")
flush.console()
trparsopt = initparsopt/(1 + initparsopt)
trparsfix = parsfix/(1 + parsfix)
trparsfix[which(parsfix == Inf)] = 1
pars2 = c(tdmodel,cond,btorph,0,soc,tol,maxiter)
initloglik = bd_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,pars2 = pars2,brts = brts,missnumspec = missnumspec)
cat("The likelihood for the inital parameter values is",initloglik,"\n")
if(initloglik == -Inf)
{
   cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
   out2 = data.frame(lambda0 = -1,mu0 = -1,lambda1 = -1, mu1 = -1, loglik = -1, df = -1, conv = -1)
} else {
#code up to DDD v1.6: out = optimx2(trparsopt,bd_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = pars2[8],reltol = pars2[7],trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,brts = brts, pars2 = pars2,missnumspec = missnumspec)
out = bd_simplex(trparsopt,idparsopt,trparsfix,idparsfix,pars2,brts,missnumspec)
if(out$conv > 0)
{
   cat("Optimization has not converged. Try again with different initial values.\n")
   out2 = data.frame(lambda0 = -1,mu0 = -1,lambda1 = -1, mu1 = -1, loglik = -1, df = -1, conv = -1)
} else {
MLtrpars = as.numeric(unlist(out$par))
MLpars = MLtrpars/(1-MLtrpars)
MLpars1 = rep(0,4)
MLpars1[idparsopt] = MLpars
if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
ML = as.numeric(unlist(out$fvalues))
out2 = data.frame(lambda0 = MLpars1[1], mu0 = MLpars1[2], lambda1 = MLpars1[3], mu1 = MLpars1[4], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
s1 = sprintf('Maximum likelihood parameter estimates: lambda0: %f, mu0: %f, lambda1: %f, mu1: %f: ',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4])
s2 = sprintf('Maximum loglikelihood: %f',ML)
cat("\n",s1,"\n",s2,"\n")
}
}
}
}
invisible(out2)
}
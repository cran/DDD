dd_SR_ML = function(brts,initparsopt=c(0.5,0.1,2*(1+length(brts)+missnumspec),2*(1+length(brts)+missnumspec),max(brts)/2),parsfix = NULL,idparsopt=c(1:3,6:7),idparsfix=NULL,idparsnoshift=(1:7)[c(-idparsopt,(-1)^(length(idparsfix) != 0) * idparsfix)],res=10*(1 + length(brts) + missnumspec),ddmodel=1,missnumspec=0,cond = TRUE,btorph = 1,allbp = FALSE)
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
# - btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - allbp = optimize likelihood for fixed tshift at all bp (TRUE) or by letting tshift vary freely (FALSE)

brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
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
if(length(namepars[idparsnoshift]) == 0) { noshiftstr = "anything" } else { noshiftstr = namepars[idparsnoshift] }
cat("You are not shifting",noshiftstr,"\n")
trparsopt = initparsopt/(1 + initparsopt)
trparsfix = parsfix/(1 + parsfix)
cat("Optimizing the likelihood - this may take a while.","\n")
out = optimx2(trparsopt,dd_SR_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = 1E-6,reltol = 1E-6,trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,brts = brts,pars2 = c(res,ddmodel,cond,btorph),missnumspec = missnumspec)
MLtrpars = unlist(out$par)
MLpars = MLtrpars/(1-MLtrpars)
out$par = list(MLpars)
MLpars1 = rep(0,7)
MLpars1[idparsopt] = MLpars
ML = out$fvalues
if(sum(idparsfix == 7) == 0 && allbp == TRUE)
{ 
   idparsopt1 = idparsopt[1:(length(idparsopt) - 1)]
   idparsfix1 = c(idparsfix,7)
   for(bp in 2:length(brts))
   {
      for(ba in seq(-1,1,2))
      {
         initparsopt1 = initparsopt[1:length(idparsopt1)]
         parsfix1 = c(idparsfix,brts[bp] + ba * 1E-8)
         trparsopt1 = initparsopt1/(1 + initparsopt1)
         trparsfix1 = parsfix1/(1 + parsfix1)
         out = optimx2(trparsopt1,dd_SR_loglik_choosepar,hess=NULL,method = "Nelder-Mead",hessian = FALSE,control = list(maximize = TRUE,abstol = 1E-4,reltol = 1E-6,trace = 0,starttests = FALSE,kkt = FALSE),trparsfix = trparsfix1,idparsopt = idparsopt1,idparsfix = idparsfix1,idparsnoshift = idparsnoshift,brts = brts,pars2 = c(res,ddmodel,cond,btorph),missnumspec = missnumspec)
         if(as.numeric(out$fvalues) > ML)
         {
            MLtrpars = unlist(out$par)
            MLpars = MLtrpars/(1-MLtrpars)
            out$par = list(MLpars)
            ML = out$fvalues
         }
      }
   }
}
if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
if(MLpars1[6] > 10^7){MLpars1[6] = Inf}
if(length(idparsfix) != 0) {MLpars1[idparsfix] = parsfix }
if(length(idparsnoshift) != 0) { MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 3] }
s1 = sprintf('Maximum likelihood parameter estimates: %f %f %f %f %f %f %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7])
s2 = sprintf('Maxmimum loglikelihood: %f',ML)
cat("\n",s1,"\n",s2,"\n")
out$par = list(MLpars1)
out$fvalues = list(ML)
out2 = data.frame(row.names = "results",lambda_1 = MLpars1[1],mu_1 = MLpars1[2],K_1 = MLpars1[3],lambda_2 = MLpars1[4],mu_2 = MLpars1[5],K_2 = MLpars1[6],t_shift = MLpars1[7],loglik = unlist(out$fvalues),df = length(initparsopt),conv = unlist(out$conv))
}
invisible(out2)
}
}
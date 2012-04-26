dd_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec)
{
trpars1 = rep(0,3)
if(pars2[2] == 5) {trpars1 = rep(0,4)}
trpars1[idparsopt] = trparsopt
if(length(idparsfix) != 0) { trpars1[idparsfix] = trparsfix }
if(max(trpars1) > 1 || min(trpars1) < 0 || trpars1[1] <= trpars1[2])
{ loglik = -Inf } else {
pars1 = trpars1/(1 - trpars1)
loglik = dd_loglik(pars1,pars2,brts,missnumspec)}
return(loglik)
}
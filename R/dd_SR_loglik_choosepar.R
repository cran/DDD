dd_SR_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,idparsnoshift,pars2,brts,missnumspec)
{
trpars1 = rep(0,7)
trpars1[idparsopt] = trparsopt
if(length(idparsfix) != 0) { trpars1[idparsfix] = trparsfix }
if(length(idparsnoshift) != 0) { trpars1[idparsnoshift] = trpars1[idparsnoshift - 3] }
loglik = dd_SR_loglik(trpars1,pars2,brts,missnumspec)
return(loglik)
}
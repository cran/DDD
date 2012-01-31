dd_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec)
{
trpars1 = rep(0,3)
trpars1[idparsopt] = trparsopt
if(length(idparsfix) != 0) { trpars1[idparsfix] = trparsfix }
loglik = dd_loglik(trpars1,pars2,brts,missnumspec)
return(loglik)
}
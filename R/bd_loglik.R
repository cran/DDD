bd_loglik = function(pars1,pars2,brts,missnumspec)
# pars1 contains model parameters
# - pars1[1] = la0 = speciation rate
# - pars1[2] = mu0 = extinction rate
# - pars1[3] = la1 = parameter in exponential decay of speciation rate (default 0)
# - pars1[4] = mu2 = parameter in exponential decay of extinction rate (default 0)
# - pars1[5] = T0 = age at which lambda is lambda0 (default T0 = age of phylogeny)
# pars2 contains settings
# - pars2[1] = tdmodel = model of time-dependence
# . tdmodel == 0: no time-dependence
# . tdmodel == 1: exponential decline in speciation or extinction rate
# . tdmodel == 2: stepwise decline following diversity-dependence when extinction = 0.
# - pars2[2] = cond = conditioning on non-extinction of the clade (1) or not (0)
# - pars2[3] = btorph = likelihood for branching times (0) or phylogeny (1)
# - pars2[4] = printing of parameters and likelihood (1) or not (0)
# - pars2[5] = likelihood is for a tree with crown age (2) or stem age (1)
{
rhotaut = function(tau,t1,pars)
{
   t = t1
   la0 = pars[1]
   mu0 = pars[2]
   la1 = pars[3]
   mu1 = pars[4]  
   if(la1 > 0 & mu1 == 0) {rtt = mu0 * (tau - t) - la0/la1 * (exp(-la1*t) - exp(-la1*tau))}
   if(la1 == 0 & mu1 > 0) {rtt = mu0/mu1 * (exp(-mu1*t) - exp(-mu1*tau)) - la0 * (tau - t)}
   if(la1 > 0 & mu1 > 0) {rtt = mu0/mu1 * (exp(-mu1*t) - exp(-mu1*tau)) - la0/la1 * (exp(-la1*t) - exp(-la1*tau))}
   return(rtt)
}

PtTint = function(x,t1,pars)
{ 
   PtTint = pars[2] * exp(-pars[4] * x) * exp(rhotaut(x,t1,pars))
   return(PtTint)
}

tdmodel = pars2[1]
cond = pars2[2]
btorph = pars2[3]
m = missnumspec

brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
brts = c(brts[1],brts)
T = brts[1]
t = T - brts

S = length(brts)
N = S + m
PtT = rep(0,S)
ux = rep(0,S)

abstol = 1e-16
reltol = 1e-10 

if(length(pars1) > 2)
{
   la1 = pars1[3]
   K = pars1[3]
   mu1 = pars1[4]
   if(length(pars1) == 5) { T0 = pars1[5] } else { T0 = T }
} else {
   la1 = 0
   mu1 = 0
   pars1 = cbind(pars1,c(0,0))
   T0 = T
}
la0 = pars1[1] * exp(-la1 * (T0 - T))
mu0 = pars1[2] * exp(-mu1 * (T0 - T))

loglik = (btorph == 0) * lgamma(S)

if(tdmodel == 2)
{
   if(N > K + 1)
   {
      loglik = -Inf
      return(as.numeric(loglik))
   }
   PtT = rep(1,S)
   la = pmax(0,la0 * (1 - (2:S)/K))
   dx = c(abs(diff(brts)),brts[S])[2:S]
   ladx = la * dx
   for(i in 2:S) { ux[i] = 1 - exp(-sum(ladx[(i - 1):(S - 1)])) }
   ux[1] = ux[2]
   if(cond == 2)
   {
      eps = 1E-6
      K2 = K
      if(floor(K) == ceiling(K)) { K2 = K + eps } else {
      if(floor(K + eps) == ceiling(K)) { K2 = K - eps } else {
      if(ceiling(K - eps) == floor(K)) { K2 = K + eps }}}
      s = (1:N) * pmax(0,la0 * (1 - (1:N)/K2))
      logsdiff = rep(0,N)
      sgnsdiff = rep(0,N)
      for(i in 2:N)
      { 
         s2 = s
         s2[i] = s[i] + 1
         sgnsdiff[i] = prod(sign(s2[2:N] - s[i]))
         logsdiff[i] = sum(log(abs(s2[2:N] - s[i])))
      }
      logp = sum(log(s[2:(N-1)])) + log(sum(sgnsdiff[2:N] * exp(-s[2:N]*T - logsdiff[2:N])))
   } else { logp = 0 }
} else {
   if(la1 == 0 & mu1 == 0)
   {
      if(abs(la0 - mu0) < 1E-10)
      {
         lamu = max(la0,mu0)
         PtT = 1/(1 + lamu * (T - t))
         ux = lamu * (T - t)/(1 + lamu * (T - t))
      } else {
         PtT = (la0 - mu0)/(la0 - mu0 * exp(-(la0 - mu0) * (T - t)))
         ux = la0 * (1 - exp(-(la0 - mu0) * (T - t)))/(la0 - mu0 * exp(-(la0 - mu0) * (T - t)))
      }
   } else {
      for(i in 1:S)
      {
         PtT[i] = (1 + integrate(PtTint, lower = t[i], upper = T, t1 = t[i], pars = pars1)$value)^(-1)
         ux[i] = 1 - PtT[i] * exp(rhotaut(T,t[i],pars1))
      }
   }
   logp = log(S + m - 1) + 2 * log(PtT[1]) + 2 * log(1 - ux[1]) + (S + m - 2) * log(ux[1])
}

if(S > 2)
{
   if(tdmodel == 2)
   {
      loglik = loglik + sum(log(la0 * (1 - (2:(S - 1))/K)))
   } else {
      loglik = loglik + sum(log(la0 * exp(-la1 * t[3:S]))) 
   }
}
loglik = loglik + sum(log(PtT)) + sum(log(1 - ux)) - (cond == 1) * (2 * log(PtT[1])) - (cond == 2) * logp

if(m > 0)
{
   if(tdmodel == 2)
   {
      cat("Missing species in diversity-dependence models cannot be dealt with by bd_loglik.R\n" )
      return(-Inf)
   }
   x = rep(0,m + 1)
   x[1] = 1
   for(j in 1:S)
   {
       #x = convolve(x,rev((1:(m + 1)) * ux[j]^(0:m)),type = 'open')[1:(m + 1)]
       x = conv(x,(1:(m + 1)) * ux[j]^(0:m))[1:(m+1)]
   }
   loglik = loglik + lgamma(S + 1) + lgamma(m + 1) - lgamma(S + m + 1) + log(x[m + 1])
   #loglik = loglik - log(S + m + 1) + log(S + 1) + log(S + m - 1) - log(S - 1)
}
if(pars2[4] == 1)
{
    if(la1 == 0 & mu1 == 0)
    {
       s1 = sprintf('Parameters: %f %f',pars1[1],pars1[2])
    } else {
       s1 = sprintf('Parameters: %f %f %f %f',pars1[1],pars1[2],pars1[3],pars1[4])
    }
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    flush.console()
}

return(as.numeric(loglik))
}
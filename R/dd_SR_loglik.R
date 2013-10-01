dd_SR_loglik = function(pars1,pars2,brts,missnumspec)
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = la = (initial) speciation rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = carrying capacity
# - pars1[4] = la2 = (initial) speciation rate
# - pars1[5] = mu2 = extinction rate
# - pars1[6] = K2 = carrying capacity
# - pars1[7] = tshift = time of shift
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model, mode of diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate
#  . ddep == 2 : exponential dependence in speciation rate
#  . ddep == 2.1: variant with offset at infinity
#  . ddep == 2.2: 1/n dependence in speciation rate
#  . ddep == 3 : linear dependence in extinction rate
#  . ddep == 4 : exponential dependence in extinction rate
#  . ddep == 4.1: variant with offset at infinity
#  . ddep == 4.2: 1/n dependence in speciation rate
# - pars2[3] = cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
#  . cond == 2 : conditioning on non-extinction of the phylogeny and on the total number of extant taxa (including missing species)
#  . cond == 3 : conditioning on the total number of extant taxa (including missing species)
# - pars2[4] = btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - pars2[5] = parameters and likelihood should be printed (1) or not (0)
# - pars2[6] = likelihood is for a tree with crown age (2) or stem age (1)
# missnumspec = number of missing species    

if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
ddep = pars2[2]
abstol = 1e-16
reltol = 1e-10 
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0) { brts[length(brts) + 1] = 0 }
soc = pars2[6]
S = length(brts) + (soc - 2)
if(min(pars1) < 0 || -pars1[7] <= min(brts))
{
    loglik = -Inf
} else {
if(((pars1[2] == 0 || pars1[4] == 0) && (ddep == 2 | ddep == 2.1 | ddep == 2.2)) || ((pars1[1] == 0 || pars1[3] == 0) && (ddep == 4 | ddep == 4.1 | ddep == 4.2)) || pars1[1] <= pars1[2] || pars1[4] <= pars1[5])
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
    loglik = -Inf
} else {
    la = pars1[1]
    mu = pars1[2]
    K = pars1[3]
    la2 = pars1[4]
    mu2 = pars1[5]
    K2 = pars1[6]
    tshift = -pars1[7]
    if(sum(abs(brts - tshift) < 1E-14) == 1) { tshift = tshift - 1E-8 }
    kshift = 1 + max(which(brts < tshift))
    n0 = (ddep == 2 | ddep == 4)
    if(ddep == 1) 
    { 
       lx = min(max(1 + missnumspec,1 + ceiling(max(la/(la - mu) * K,la2/(la2 - mu2) * K2))),round(pars2[1]))
    } else {
       lx = round(pars2[1])
    }     
    cond = pars2[3]
    btorph = pars2[4]

    if(ddep == 1 && (ceiling(la/(la - mu) * K) < kshift || ceiling(la2/(la2 - mu2) * K2) < (S + missnumspec))) { loglik = -Inf } else
    {
       loglik = (btorph == 0) * lgamma(S)
       if(cond != 3)
       {
           probs = rep(0,lx)
           probs[1] = 1 # change if other species at stem/crown age   
      
           if(kshift > 2)
           {
              for(k in 2:(kshift-1))
              {
                 k1 = k + (soc - 2)
                 y = lsoda(probs,brts[(k-1):k],dd_loglik_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol)
                 probs = y[2,2:(lx+1)]
                 if(k<(S + 2 - soc))
                 {
                     probs = flavec(ddep,la,mu,K,0,lx,k1,n0) * probs # speciation event
                     if(sum(probs) <= 0)
                     { 
                        loglik = -Inf
                        break
                     } else {
                        loglik = loglik + log(sum(probs))
                     }
                     probs = probs/sum(probs)
                 }
              }
           }   
           k = kshift
           k1 = k + (soc - 2)
           y = lsoda(probs,c(brts[k-1],tshift),dd_loglik_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol)
           probs = y[2,2:(lx+1)]
           y = lsoda(probs,c(tshift,brts[k]),dd_loglik_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol)
           probs = y[2,2:(lx+1)] 
           if(k<(S + 2 - soc))
           {
               probs = flavec(ddep,la2,mu2,K2,0,lx,k1,n0) * probs # speciation event
               if(sum(probs) <= 0) { loglik = -Inf } else
               {
                  loglik = loglik + log(sum(probs))
               }
               probs = probs/sum(probs)
           }
           if((kshift + 1) <= (S + 2 - soc))
           {
              for(k in (kshift + 1):(S + 2 - soc))
              {
                 k1 = k + (soc - 2)
                 y = lsoda(probs,brts[(k-1):k],dd_loglik_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol)
                 probs = y[2,2:(lx+1)]
                 if(k<(S + 2 - soc))
                 {
                     probs = flavec(ddep,la2,mu2,K2,0,lx,k1,n0) * probs # speciation event
                     if(sum(probs) <= 0) { loglik = -Inf } else
                     {
                        loglik = loglik + log(sum(probs))
                     }
                     probs = probs/sum(probs)
                 }
              }
           }    
       } else {
           probs = rep(0,lx + 1)
           probs[1 + missnumspec] = 1
           if((kshift + 1) <= S + 2 - soc)  
           {      
              for(k in (S + 2 - soc):(kshift + 1))
              {
                 k1 = k + (soc - 2)
                 y = lsoda(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol)
                 probs = y[2,2:(lx+2)]
                 if(k > soc)
                 {
                     probs = c(flavec(ddep,la2,mu2,K2,0,lx,k1-1,n0),1) * probs # speciation event
                     if(sum(probs[1:lx]) <= 0)
                     {
                        loglik = -Inf
                        break
                     } else {
                        loglik = loglik + log(sum(probs[1:lx]))
                     }
                     probs[1:lx] = probs[1:lx]/sum(probs[1:lx])
                 }    
              }
           }
           k = kshift
           k1 = k + (soc - 2)
           y = lsoda(probs,-c(brts[k],tshift),dd_loglik_bw_rhs,c(pars1[4:6],k1,ddep),rtol = reltol,atol = abstol)
           probs = y[2,2:(lx+2)]
           y = lsoda(probs,-c(tshift,brts[k-1]),dd_loglik_bw_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol)
           probs = y[2,2:(lx+2)]
           if(k > soc)
           {
              probs = c(flavec(ddep,la,mu,K,0,lx,k1-1,n0),1) * probs # speciation event
              if(sum(probs[1:lx]) <= 0)
              {
                 loglik = -Inf
                 break
              } else {
                 loglik = loglik + log(sum(probs[1:lx]))
              }
              probs[1:lx] = probs[1:lx]/sum(probs[1:lx])
           }    
           for(k in (kshift-1):2)
           {
              k1 = k + (soc - 2)
              y = lsoda(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol)
              probs = y[2,2:(lx+2)]
              if(k > soc)
              {
                  probs = c(flavec(ddep,la,mu,K,0,lx,k1-1,n0),1) * probs # speciation event
                  if(sum(probs[1:lx]) <= 0)
                  {
                     loglik = -Inf
                     break
                  } else {
                     loglik = loglik + log(sum(probs[1:lx]))
                  }
                  probs[1:lx] = probs[1:lx]/sum(probs[1:lx])
              }    
           }
       }
 
       if(probs[1 + (cond != 3) * missnumspec] <= 0 || loglik == -Inf) { loglik = -Inf } else
       {
          loglik = loglik + (cond != 3 & soc == 2) * log(probs[1 + (cond != 3) * missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
                    
          logliknorm = 0
          if(cond == 1 | cond == 2)
          { 
             probs = rep(0,lx)
             probs[1] = 1 # change if other species at crown age
             k = soc
             y = lsoda(probs,c(brts[1],tshift),dd_loglik_rhs,c(pars1[1:3],k,ddep),rtol = reltol,atol = abstol);
             probs = y[2,2:(lx+1)]
             y = lsoda(probs,c(tshift,brts[length(brts)]),dd_loglik_rhs,c(pars1[4:6],k,ddep),rtol = reltol,atol = abstol);
             probs = y[2,2:(lx+1)]
             if(soc == 1) { aux = 1:lx }
             if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
             probsc = probs/aux
             if(cond == 1) { logliknorm = log(sum(probsc)) }
             if(cond == 2) { logliknorm = log(probsc[S + missnumspec - 1])}             
          }
          if(cond == 3)
          { 
             probsn = rep(0,lx + 1)
             probsn[S + missnumspec + 1] = 1
             T = max(1,1/abs(la - mu),1/abs(la2 - mu2)) * 100 * max(abs(brts)) # make this more efficient later
             y = lsoda(probsn,c(0,-tshift),dd_loglik_bw_rhs,c(pars1[4:6],0,ddep),rtol = reltol,atol = abstol)
             probsn = y[2,2:(lx+2)]
             y = lsoda(probsn,c(-tshift,T),dd_loglik_bw_rhs,c(pars1[1:3],0,ddep),rtol = reltol,atol = abstol)
             logliknorm = log(y[2,lx + 2])
             if(soc == 2)
             {
                probsn = rep(0,lx + 1)
                probsn[1:lx] = probs[1:lx]
                probsn = c(flavec(ddep,la,mu,K,0,lx,1,n0),1) * probsn # speciation event
                y = lsoda(probsn,c(max(abs(brts)),T),dd_loglik_bw_rhs,c(pars1[1:3],1,ddep),rtol = reltol,atol = abstol)
                logliknorm = logliknorm - log(y[2,lx + 2])
             }
          }
          loglik = loglik - logliknorm   
       }
    }
}}
if(pars2[5] == 1)
{
    s1 = sprintf('Parameters: %f %f %f %f %f %f %f, ',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6],pars1[7])
    s2 = sprintf('Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    flush.console()
}
return(as.numeric(loglik))
}
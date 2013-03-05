dd_loglik = function(pars1,pars2,brts,missnumspec)
{
# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = la = (initial) speciation rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = carrying capacity
# - pars1[4] = r = ratio of diversity-dependence in extinction rate over speciation rate
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
#  . ddep==1 : linear dependence in speciation rate
#  . ddep==2 : exponential dependence in speciation rate
#  . ddep==3 : linear dependence in extinction rate
#  . ddep==4 : exponential dependence in extinction rate
#  . ddep==5 : linear dependence in speciation and extinction rate
# - pars2[3] = cond = conditioning on non-extinction of the phylogeny
# - pars2[4] = btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - pars2[5] = parameters and likelihood should be printed (1) or not (0)
# - pars2[6] = likelihood is for a tree with crown age (2) or stem age (1)
# missnumspec = number of missing species    

if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
abstol = 1e-16
reltol = 1e-10 
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0) { brts[length(brts) + 1] = 0 }
soc = pars2[6]
S = length(brts) + (soc - 2)
ddep = pars2[2]
if(min(pars1) < 0 || pars1[1] <= pars1[2] || pars1[3] <= (S + missnumspec)) { loglik = -Inf } else
{
    la = pars1[1]
    mu = pars1[2]
    K = pars1[3]
    if(ddep == 5) {r = pars1[4]} else {r = 0}
    cond = pars2[3]
    btorph = pars2[4]


    if((ddep == 1 || ddep == 5) && ceiling(la/(la - mu) * (r + 1) * K) < (S + missnumspec)) { loglik = -Inf } else
    {
       if(ddep == 1 || ddep == 5) { lx = min(max(1 + missnumspec,1 + ceiling(la/(la - mu) * (r + 1) * K)),round(pars2[1])) } else { lx = round(pars2[1]) }
       probs = rep(0,lx)
       probs[1] = 1 # change if other species at stem/crown age  
       loglik = (btorph == 0) * lgamma(S)
       for(k in 2:(S + 2 - soc))
       {
          k1 = k + (soc - 2)
          y = lsoda(probs,brts[(k-1):k],dd_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol)
          probs = y[2,2:(lx+1)]
          if(k<(S + 2 - soc))
          {
              if(ddep == 1 || ddep == 5) { lavec = pmax(rep(0,lx),la - 1/(r + 1) * (la-mu)/K * ((0:(lx-1))+k1)) } 
              if(ddep == 2) { lavec = pmax(rep(0,lx),la * (((0:(lx-1))+k1) + 1)^(-log(la/mu)/log(K+1))) }
              if(ddep == 3 || ddep == 4) { lavec = la }    
              probs = lavec * probs # speciation event
              if(sum(probs) <= 0) { loglik = -Inf } else
              {
                 loglik = loglik + log(sum(probs))
              }
              probs = probs/sum(probs)
          }
       }    
       if(probs[1+missnumspec]<=0) { loglik = -Inf } else
       {        
          loglik = loglik + log(probs[1 + missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
  
          if(cond == TRUE)
          {
             probs = rep(0,lx)
             probs[1] = 1 # change if other species at stem or crown age
             k = soc
             t1 = brts[1] 
             t2 = brts[S + 2 - soc]
             y = lsoda(probs,c(t1,t2),dd_loglik_rhs,c(pars1,k,ddep),rtol = reltol,atol = abstol);
             probs = y[2,2:(lx+1)]
             if(soc == 1) { aux = 1:lx }
             if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
             logliknorm = log(sum(probs/aux))
          } else { logliknorm = 0 }
          loglik = loglik - logliknorm
       }
    }
}
if(pars2[5] == 1)
{
    s1 = sprintf('Parameters: %f %f %f',pars1[1],pars1[2],pars1[3])
    if(ddep == 5) {s1 = sprintf('%s %f',s1,pars1[4])}
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    flush.console()
}
return(as.numeric(loglik))
}
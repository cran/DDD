dd_KI_loglik = function(pars1,pars2,brtsM,brtsS,missnumspec)
{
# brtsM = branching times of main clade M (positive, from present to past)
# brtsS = branching times of subclade S (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = laM = (initial) speciation rate of main clade
# - pars1[2] = muM = extinction rate of main clade
# - pars1[3] = KM = carrying capacity of main clade
# - pars1[4] = laS = (initial) speciation rate of subclade
# - pars1[5] = muS = extinction rate of subclade
# - pars1[6] = KS = carrying capacity of subclade
# - pars1[7] = tinn = time of key innovation
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model, mode of diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate
#  . ddep == 2 : exponential dependence in speciation rate
#  . ddep == 3 : linear dependence in extinction rate
#  . ddep == 4 : exponential dependence in extinction rate
# - pars2[3] = cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
# - pars2[4] = tsplit = time of split of innovative branch
# missnumspec = number of missing species in main clade M and subclade S

if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
abstol = 1e-16
reltol = 1e-14
# order branching times
brts=-sort(abs(c(brtsM,brtsS)),decreasing = TRUE)
if(sum(brts == 0) == 0) { brts[length(brts) + 1] = 0 }
soc = pars2[6]
S = length(brts) + (soc - 2)
brtsM =-sort(abs(brtsM),decreasing = TRUE)
if(sum(brtsM == 0) == 0) { brtsM[length(brtsM) + 1] = 0 }
brtsS =-sort(abs(brtsS),decreasing = TRUE)
if(sum(brtsS == 0) == 0) { brtsS[length(brtsS) + 1] = 0 }

if(min(pars1) < 0 || -pars1[7] <= min(brtsM) || -pars1[7] >= min(brtsS))
{
    loglik = -Inf
} else {
if(((pars1[2] == 0 || pars1[4] == 0) && ddep == 2) || ((pars1[1] == 0 || pars1[3] == 0) && ddep == 4) || pars1[1] <= pars1[2] || pars1[4] <= pars1[5])
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
    loglik = -Inf
} else {
    laM = pars1[1]
    muM = pars1[2]
    KM = pars1[3]
    laS = pars1[4]
    muS = pars1[5]
    KS = pars1[6]
    tinn = -pars1[7]
    lmax = pars2[1]
    ddep = pars2[2]
    cond = pars2[3]
    tsplit = -pars2[4]
    m = missnumspec
    S1 = length(brtsM)-1 + (soc - 2)
    if(sum(brtsS == tinn) == 0) { brtsS = c(tinn,brtsS) }
    S2 = length(brtsS)-1
    if(ddep == 1 && (ceiling(laM/(laM - muM) * KM) < S1 || (ceiling(laS/(laS - muS) * KS) < S2 ))) { loglik = -Inf } else {

    # avoid coincidence of branching time and key innovation time
    if(sum(abs(brtsM - tinn) < 1E-14) == 1) { tinn = tinn - 1E-8 }

    # compute likelihood of clade M
    loglikM = 0
    if(ddep == 1) { lx = min(max(1 + m[1],1 + ceiling(laM/(laM - muM) * KM)),round(lmax)) } else { lx = round(lmax) }
    probs = rep(0,lx)
    probs[1] = 1 # change if other species at crown age

    ka = sum(brtsM<tinn);
    for(k in 2:(ka+1))
    {
       k1 = k + (soc - 2)
       t1 = brtsM[k-1]; t2 = min(c(tinn,brtsM[k]))
       y = lsoda(probs,c(t1,t2),dd_loglik_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol)
       probs = y[2,2:(lx+1)]
       if(t2<tinn)
       {
           if(ddep == 1) { lavec = pmax(rep(0,lx),laM - (laM-muM)/KM * ((0:(lx-1))+k1)) } 
           if(ddep == 2) { lavec = pmax(rep(0,lx),laM * (((0:(lx-1))+k1) + 1)^(-log(laM/muM)/log(KM+1))) }
           if(ddep == 3 || ddep == 4) { lavec = laM }
           probs = lavec * probs # speciation event
           if(sum(probs) <= 0)
           { 
              loglik = -Inf
              break
           } else {
              loglikM = loglikM + log(sum(probs))
           }
           probs = probs/sum(probs)
       }
    }
    for(k in (ka+1):(S1+1))
    {
       k1 = k + (soc - 2)
       t1 = max(tinn,brtsM[k-1]); t2 = brtsM[k];
       y = lsoda(probs,c(t1,t2),dd_loglik_rhs,c(pars1[1:3],k1-1,ddep),rtol = reltol,atol = abstol)
       probs = y[2,2:(lx+1)]
       if(k<(S1+1))
       {
           if(ddep == 1) { lavec = pmax(rep(0,lx),laM - (laM-muM)/KM * ((0:(lx-1))+k1-1)) } 
           if(ddep == 2) { lavec = pmax(rep(0,lx),laM * (((0:(lx-1))+k1	-1) + 1)^(-log(laM/muM)/log(KM+1))) }
           if(ddep == 3 || ddep == 4) { lavec = laM }    
           probs = lavec * probs # speciation event
           if(sum(probs) <= 0) { loglik = -Inf } else
           {
              loglikM = loglikM + log(sum(probs))
           }
           probs = probs/sum(probs)
       }
    }
    if(length(m) == 1)
    { 
       loglikM = loglikM + log(probs[1 + (0:m)])   
    } else {
       loglikM = loglikM + log(probs[1 + m[1]])   
    }
    # compute likelihood of clade S
    loglikS = 0
    if(ddep == 1) { lx = min(max(1 + m[length(m)],1 + ceiling(laS/(laS - muS) * KS)),round(lmax)) } else { lx = round(lmax) }
    probs = rep(0,lx)
    probs[1] = 1
    for(k in 1:S2)
    {
       t1 = brtsS[k]; t2 = brtsS[k+1]
       y = lsoda(probs,c(t1,t2),dd_loglik_rhs,c(pars1[4:6],k,ddep),rtol = reltol,atol = abstol)
       probs = y[2,2:(lx+1)]
       if(k<S2)
       {
           if(ddep == 1) { lavec = pmax(rep(0,lx),laS - (laS-muS)/KS * ((0:(lx-1))+k)) } 
           if(ddep == 2) { lavec = pmax(rep(0,lx),laS * (((0:(lx-1))+k) + 1)^(-log(laS/muS)/log(KS+1))) }
           if(ddep == 3 || ddep == 4) { lavec = laS }    
           probs = lavec * probs # speciation event
           if(sum(probs) <= 0) { loglik = -Inf } else
           {
              loglikS = loglikS + log(sum(probs))
           }
           probs = probs/sum(probs)
       }
    }
    if(length(m) == 1)
    {
       loglikS = loglikS + log(probs[1 + (0:m)])
    } else {
       loglikS = loglikS + log(probs[1 + m[2]])
    }
    # total likelihood = likelihood clade M x likelihood clade S
    if(length(m) == 1)
    {
       loglik = log(sum(exp(loglikM + loglikS[length(loglikS):1])))
    } else {
       loglik = loglikM + loglikS
    }
   
    if(cond == 0 || loglik == -Inf){logliknorm = 0} else {   
       # COMPUTE NORMALIZATION
       tcrown = brts[1]
       tpres = 0
       # compute survival probability of clade S
       lx = min(lmax,ceiling(laS/(laS-muS)*KS)+1)
       nx = -1:lx
       if(ddep == 1) 
       { 
           lavec = pmax(rep(0,lx + 2),laS - (laS-muS)/KS * nx)
           muvec = muS * rep(1,lx + 2)
       } 
       if(ddep == 2)
       { 
           lavec = pmax(rep(0,lx),laS * (nx + 1)^(-log(laS/muS)/log(KS+1)))
           muvec = muS * rep(1,lx + 2)
       }
       if(ddep == 3)
       {
           lavec = laS * rep(1,lx + 2)
           muvec = muS + (laS - muS)/KS * nx
       }    
       if(ddep == 4)
       {
           lavec = laS * rep(1,lx + 2)
           muvec = (nx + 1)^(log(laS/muS)/log(KS+1))
       }    
       m1 = lavec[1:lx] * nx[1:lx]
       m2 = muvec[3:(lx+2)] * nx[3:(lx+2)]
       m3 = (lavec[2:(lx+1)] + muvec[2:(lx+1)]) * nx[2:(lx+1)]
       probs = rep(0,lx) # probs[1] = extinction probability
       probs[2] = 1 # clade S starts with one species
       y = lsoda(probs,c(tinn,tpres),dd_logliknorm_rhs1,c(m1,m2,m3),rtol = reltol,atol = abstol)
       probs = y[2,2:(lx+1)]   
       PS = 1 - probs[1]
   
       # compute survival probability of clade M
       lx = min(lmax,ceiling(laM/(laM-muM)*KM)+1)
       nx1 = rep(-1:lx,lx + 2)
       dim(nx1) = c(lx + 2,lx + 2) # row index = number of species in first group 
       nx2 = t(nx1) # column index = number of species in second group
       nxt = nx1 + nx2
       if(ddep == 1) 
       { 
           lavec = pmax(matrix(0,lx + 2,lx + 2),laM - (laM-muM)/KM * nxt)
           muvec = muM * matrix(1,lx + 2,lx + 2)
       } 
       if(ddep == 2)
       { 
           lavec = pmax(matrix(0,lx + 2,lx + 2),laM * (nxt + 1)^(-log(laM/muM)/log(KM+1)))
           muvec = muM * matrix(1,lx + 2,lx + 2)
       }
       if(ddep == 3)
       {
           lavec = laM * matrix(1,lx + 2,lx + 2)
           muvec = muM + (laM - muM)/KM * nxt
       }    
       if(ddep == 4)
       {
           lavec = laM * matrix(1,lx + 2,lx + 2)
           muvec = (nxt + 1)^(log(laM/muM)/log(KM+1))
       }    
       m1 = lavec[1:lx,2:(lx+1)] * nx1[1:lx,2:(lx+1)]
       m2 = muvec[3:(lx+2),2:(lx+1)] * nx1[3:(lx+2),2:(lx+1)]
       ma = lavec[2:(lx+1),2:(lx+1)] + muvec[2:(lx+1),2:(lx+1)]
       m3 = ma * nx1[2:(lx+1),2:(lx+1)]
       m4 = lavec[2:(lx+1),1:lx] * nx2[2:(lx+1),1:lx]
       m5 = muvec[2:(lx+1),3:(lx+2)] * nx2[2:(lx+1),3:(lx+2)]
       m6 = ma * nx2[2:(lx+1),2:(lx+1)]
       probs = matrix(0,lx,lx)

       # probs[1,1] = probability of extinction of both lineages
       # sum(probs[1:lx,1]) = probability of extinction of second lineage
       probs[2,2] = 1 # clade M starts with two species
       # STEP 1: integrate from tcrown to tinn
       dim(probs) = c(lx*lx,1)
       y = lsoda(probs,c(tcrown,tinn),dd_logliknorm_rhs2,c(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol)
       probs = y[2,2:(lx * lx + 1)]
       dim(probs) = c(lx,lx)
       probs[1,1:lx] = 0
       probs[1:lx,1] = 0
       # STEP 2: transformation at tinn
       nx1a = nx1[2:(lx+1),2:(lx+1)]
       nx2a = nx2[2:(lx+1),2:(lx+1)]
       probs = probs * nx1a/(nx1a+nx2a)
       probs = t(c(t(probs[2:lx,1:lx]), rep(0,lx)))
       dim(probs) = c(lx,lx)
       # STEP 3: integrate from tinn to tpres
       dim(probs) = c(lx*lx,1);
       y = lsoda(probs,c(tinn,tpres),dd_logliknorm_rhs2,c(m1,m2,m3,m4,m5,m6),rtol = reltol,atol = abstol)
       probs = y[2,2:(lx * lx + 1)]
       dim(probs) = c(lx,lx)
       PM12 = sum(sum(probs[2:lx,2:lx]))
       PM2 = sum(probs[2:lx,1])
       logliknorm = log(2)+log(PM12+PS*PM2)
    }
    loglik = loglik - logliknorm - lgamma(S + m + 1) + lgamma(S + 1) + lgamma(m + 1)
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

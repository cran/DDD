td_loglik_rhs = function(t,x,pars)
{
lx = length(x)
lp = pars[length(pars)]
lrs = (lx - lp)/2
la = pars[1]
mu = pars[2]
K = pars[3]
ddep = pars[length(pars) - 1]
r = (ddep == 5) * pars[4]
n0 = (ddep == 2 | ddep == 4)

nn = 0:(lp + 1)
lnn = length(nn)

p = c(0,x[1:lp],0)
rho = x[(lp + 1):(lp + lrs)]
sig = x[(lp + lrs + 1):(lp + 2 * lrs)]

nntd = sum(nn * p)

if(ddep == 1)
{
    lavec = pmax(rep(0,lnn),la - (la-mu)/K * nn)
    latd = max(0,la - (la - mu)/K * nntd)
    muvec = mu * rep(1,lnn)
    mutd = mu
} else {
if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
{
    
    y = -(log(la/mu)/log(K+n0))^(ddep != 2.2)
    lavec = pmax(rep(0,lnn),la * (nn + n0)^y)
    latd = max(0,la * (nntd + n0)^y)
    muvec = mu * rep(1,lnn)
    mutd = mu
} else {
if(ddep == 3)
{
    lavec = la * rep(1,lnn)
    muvec = mu + (la - mu) * nn/K
    latd = la
    mutd = mu + (la - mu) * nntd/K
} else {
if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
{
    y = (log(la/mu)/log(K+n0))^(ddep != 4.2)
    lavec = la * rep(1,lnn)
    muvec = mu * (nn + n0)^y
    latd = la
    mutd = mu * (nntd + n0)^y
} else {
if(ddep == 5)
{ 
    lavec = pmax(rep(0,lnn),la - 1/(r+1)*(la-mu)/K * nn)
    muvec = muvec = mu + r/(r+1)*(la-mu)/K * nn
    latd = max(0,la - 1/(r+1)*(la-mu)/K * nntd)
    mutd = mu + r/(r+1)*(la-mu)/K * nntd
}}}}}

dp = lavec[(2:(lp+1))-1] * nn[(2:(lp + 1))-1] * p[(2:(lp + 1))-1] + muvec[(2:(lp+1))+1] * nn[(2:(lp+1))+1] * p[(2:(lp+1))+1] - (lavec[(2:(lp+1))] + muvec[(2:(lp+1))]) * nn[(2:(lp+1))] * p[(2:(lp+1))]

drho = rep(mutd - latd,lrs)

dsig = mutd * exp(rho)

return(list(c(dp,drho,dsig)))

}
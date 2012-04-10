dd_loglik_rhs = function(t,x,pars)
{
lx = length(x)
la = pars[1]
mu = pars[2]
K = pars[3]
if(length(pars) < 6)
{
   kk = pars[4]
   ddep = pars[5]
} else {
   r = pars[4]
   kk = pars[5]
   ddep = pars[6]
}

nn = -1:(lx+2*kk)
lnn = length(nn)

if(ddep == 1)
{
    lavec = pmax(rep(0,lnn),la - (la-mu)/K * nn)
    muvec = mu * rep(1,lnn)
}
if(ddep == 2)
{
    lavec = pmax(rep(0,lnn),la * (nn + 1)^(-log(la/mu)/log(K+1)))
    muvec = mu * rep(1,lnn)
}
if(ddep == 3)
{
    lavec = la * rep(1,lnn)
    muvec = mu + (la - mu) * nn/K
}
if(ddep == 4)
{
    lavec = la * rep(1,lnn)
    muvec = mu * (nn + 1)^(log(la/mu)/log(K+1))
}
if(ddep == 5)
{ 
    lavec = pmax(rep(0,lnn),la - 1/(r+1)*(la-mu)/K * nn)
    muvec = muvec = mu + r/(r+1)*(la-mu)/K * nn
}
xx = c(0,x,0)

dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]

return(list(dx))

}
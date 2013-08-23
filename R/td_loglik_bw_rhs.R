td_loglik_bw_rhs = function(t,x,pars)
{
lx = length(x)
la = pars[1]
mu = pars[2]
K = pars[3]
ddep = pars[length(pars)]
r = (ddep == 5) * pars[4]

nn = 0:(lx + 1)
lnn = length(nn)

p = c(0,x[1:lx],0)

if(ddep == 1)
{
    lavec = pmax(rep(0,lnn),la - (la-mu)/K * nn)
    muvec = mu * rep(1,lnn)
} else {
if(ddep == 2)
{
    lavec = pmax(rep(0,lnn),la * (nn + 1)^(-log(la/mu)/log(K+1)))
    muvec = mu * rep(1,lnn)
} else {
if(ddep == 3)
{
    lavec = la * rep(1,lnn)
    muvec = mu + (la - mu) * nn/K
} else {
if(ddep == 4)
{
    lavec = la * rep(1,lnn)
    muvec = mu * (nn + 1)^(log(la/mu)/log(K+1))
} else {
if(ddep == 5)
{ 
    lavec = pmax(rep(0,lnn),la - 1/(r+1)*(la-mu)/K * nn)
    muvec = muvec = mu + r/(r+1)*(la-mu)/K * nn
}}}}}

dx = lavec[(2:(lx+1))] * nn[(2:(lx + 1))] * p[(2:(lx + 1))+1] + muvec[(2:(lx+1))] * nn[(2:(lx+1))] * p[(2:(lx+1))-1] - (lavec[(2:(lx+1))-1] * nn[(2:(lx+1))-1] + muvec[(2:(lx+1))+1] * nn[(2:(lx+1))+1]) * p[(2:(lx+1))]
dG = p[2]

return(list(c(dx,dG)))

}
flavec = function(ddep,la,mu,K,r,lx,kk,n0)
{
   if(ddep == 1 | ddep == 5)
   { 
       lavec = pmax(rep(0,lx),la - 1/(r + 1) * (la-mu)/K * ((0:(lx-1))+kk))
   } 
   if(ddep == 2 | ddep == 2.1 | ddep == 2.2) 
   {
       x = -(log(la/mu)/log(K+n0))^(ddep != 2.2)
       lavec = pmax(rep(0,lx),la * (((0:(lx-1))+kk) + n0)^x)
   }
   if(ddep == 3 | ddep == 4 | ddep == 4.1 | ddep == 4.2)
   {
       lavec = la * rep(1,lx)
   }
   return(lavec)
}    

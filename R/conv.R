conv = function(x,y)
{
   lx = length(x)
   ly = length(y)
   lxy = length(x) + length(y)
   x = c(x,rep(0,lxy - lx))
   y = c(y,rep(0,lxy - ly))
   cvxy = rep(0,lxy)
   for(i in 2:lxy)
   {
      cvxy[i] = crossprod(x[(i-1):1],y[1:(i-1)])
   }
   return(cvxy[2:lxy])
}
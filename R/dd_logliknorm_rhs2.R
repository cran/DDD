dd_logliknorm_rhs2 = function(t,x,m)
{
   nx = sqrt(length(x))
   dim(x) = c(nx,nx)     
   m1 = m[[1]]
   m2 = m[[2]]
   m3 = m[[3]]
   m4 = m[[4]]
   m5 = m[[5]]
   m6 = m[[6]]
   ##m1 = m[1:(nx^2)]
   ##dim(m1) = c(nx,nx)
   ##m2 = m[(nx^2+1):(2*nx^2)]
   ##dim(m2) = c(nx,nx)
   ##m3 = m[(2*nx^2+1):(3*nx^2)]
   ##dim(m3) = c(nx,nx)
   ##m4 = m[(3*nx^2+1):(4*nx^2)]
   ##dim(m4) = c(nx,nx)
   ##m5 = m[(4*nx^2+1):(5*nx^2)]
   ##dim(m5) = c(nx,nx)
   ##m6 = m[(5*nx^2+1):(6*nx^2)]
   ##dim(m6) = c(nx,nx)

   xx = matrix(0,nx+2,nx+2)
   xx[2:(nx+1),2:(nx+1)] = x
   dx = m1 * xx[1:nx,2:(nx+1)] + m2 * xx[3:(nx+2),2:(nx+1)] - m3 * xx[2:(nx+1),2:(nx+1)] + m4 * xx[2:(nx+1),1:nx] + m5 * xx[2:(nx+1),3:(nx+2)] - m6 * xx[2:(nx+1),2:(nx+1)]
   dim(dx) = c(nx^2,1)
   return(list(dx))
}
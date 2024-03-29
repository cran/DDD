dd_KI_lamuN = function(ddmodel,pars,N)
{
    laM = pars[1]
    muM = pars[2]
    KM = pars[3]
    laS = pars[4]
    muS = pars[5]
    KS = pars[6]
    NM = N[1]
    NS = N[2]
    n0 = (ddmodel == 2 | ddmodel == 4)
    if(ddmodel == 1)
    {
        # linear dependence in speciation rate
        laMN = max(0,laM - (laM - muM) * NM/KM)
        muMN = muM
        laSN = max(0,laS - (laS - muS) * NS/KS)
        muSN = muS
    }
    if(ddmodel == 1.3)
    {
        # linear dependence in speciation rate
        laMN = max(0,laM * (1 - NM/KM))
        muMN = muM
        laSN = max(0,laS * (1 - NS/KS))
        muSN = muS
    }
    if(ddmodel == 2 | ddmodel == 2.1 | ddmodel == 2.2)
    {
        # exponential dependence in speciation rate
        al = (log(laM/muM)/log(KM+n0))^(ddmodel != 2.2)
        laMN = laM * (NM + n0)^(-al)
        muMN = muM
        al = (log(laS/muS)/log(KS+n0))^(ddmodel != 2.2)
        laSN = laS * (NS + n0)^(-al)
        muSN = muS
    }
    if(ddmodel == 2.3)
    {
        # exponential dependence in speciation rate
        al = KM
        laMN = laM * (NM + n0)^(-al)
        muMN = muM
        al = KS
        laSN = laS * (NS + n0)^(-al)
        muSN = muS
    }
    if(ddmodel == 3)
    {
        # linear dependence in extinction rate
        laMN = laM
        muMN = muM + (laM - muM) * NM/KM
        laSN = laS
        muSN = muS + (laS - muS) * NS/KS
    }
    if(ddmodel == 4 | ddmodel == 4.1 | ddmodel == 4.2)
    {
        # exponential dependence in extinction rate
        al = (log(laM/muM)/log(KM+n0))^(ddmodel != 4.2)
        laMN = laM
        muMN = muM * (NM + n0)^al
        al = (log(laS/muS)/log(KS+n0))^(ddmodel != 4.2)
        laSN = laS
        muSN = muS * (NS + n0)^al
    }
    return(c(laMN,muMN,laSN,muSN))
}



#' Function to simulate a key innovation in macro-evolution with the innovative
#' clade decoupling from the diversity-dependent diversification dynamics of
#' the main clade
#' 
#' Simulating a diversity-dependent diversification process where at a given
#' time a new clade emerges with different inherent speciation rate and
#' extinction rate and clade-level carrying capacity and with decoupled
#' dynamics
#' 
#' 
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to
#' lambda_M (speciation rate of the main clade) \cr \code{pars[2]} corresponds
#' to mu_M (extinction rate of the main clade) \cr \code{pars[3]} corresponds
#' to K_M (clade-level carrying capacity of the main clade) \code{pars[4]}
#' corresponds to lambda_S (speciation rate of the subclade) \cr \code{pars[5]}
#' corresponds to mu_S (extinction rate of the subclade) \cr \code{pars[5]}
#' corresponds to K_S (clade-level carrying capacity of the subclade) \cr
#' \code{pars[7]} tinn, the time the shift in rates occurs in the lineage
#' leading to the subclade
#' @param age Sets the crown age for the simulation
#' @param ddmodel Sets the model of diversity-dependence: \cr \code{ddmodel ==
#' 1} : linear dependence in speciation rate with parameter K (= diversity
#' where speciation = extinction)\cr \code{ddmodel == 1.3} : linear dependence
#' in speciation rate with parameter K' (= diversity where speciation = 0)\cr
#' \code{ddmodel == 2} : exponential dependence in speciation rate with
#' parameter K (= diversity where speciation = extinction)\cr \code{ddmodel ==
#' 2.1} : variant of exponential dependence in speciation rate with offset at
#' infinity\cr \code{ddmodel == 2.2} : 1/n dependence in speciation rate\cr
#' \code{ddmodel == 2.3} : exponential dependence in speciation rate with
#' parameter x (= exponent)\cr \code{ddmodel == 3} : linear dependence in
#' extinction rate \cr \code{ddmodel == 4} : exponential dependence in
#' extinction rate \cr \code{ddmodel == 4.1} : variant of exponential
#' dependence in extinction rate with offset at infinity \cr \code{ddmodel ==
#' 4.2} : 1/n dependence in extinction rate with offset at infinity
#' @return \item{ out }{ A list with the following elements: The first element
#' is the tree of extant species in phylo format \cr The second element is the
#' tree of all species, including extinct species, in phylo format \cr The
#' third element is a matrix of all species where \cr - the first column is the
#' time at which a species is born \cr - the second column is the label of the
#' parent of the species; positive and negative values only indicate whether
#' the species belongs to the left or right crown lineage \cr - the third
#' column is the label of the daughter species itself; positive and negative
#' values only indicate whether the species belongs to the left or right crown
#' lineage \cr - the fourth column is the time of extinction of the species \cr
#' If the fourth element equals -1, then the species is still extant.\cr - the
#' fifth column indicates whether the species belong to the main clade (0) or
#' the subclade (1)\cr The fourth element is the subclade tree of extant
#' species (without stem) \cr The fifth element is the subclade tree of all
#' species (without stem) \cr The sixth element is the same as the first,
#' except that it has attributed 0 for the main clade and 1 for the subclade\cr
#' The seventh element is the same as the Second, except that it has attributed
#' 0 for the main clade and 1 for the subclade\cr The sixth and seventh element
#' will be NULL if the subclade does not exist (because it went extinct).  }
#' @author Rampal S. Etienne
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#'  dd_KI_sim(c(0.2,0.1,20,0.1,0.05,30,4),10) 
#' @export dd_KI_sim
dd_KI_sim = function(pars,age,ddmodel = 1)
{
# Simulation of diversity-dependent process
#  . start from crown age
#  . no additional species at crown node
#  . no missing species in present
# pars = [laM muM K laS muS tinn]
# - pars1[1] = laM = (initial) speciation rate of main clade
# - pars1[2] = muM = extinction rate of main clade
# - pars1[3] = K = clade-level carrying capacity of main clade
# - pars1[4] = laS = (initial) speciation rate of subclade
# - pars1[5] = muS = extinction rate of subclade
# - pars1[6] = K = clade-level carrying capacity of subclade
# - pars1[7] = tinn = time of key innovation
# age = crown age
# ddmodel = mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate with parameter K
#  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 2.1: variant with offset at infinity
#  . ddmodel == 2.2: 1/n dependence in speciation rate
#  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
#  . ddmodel == 4.1: variant with offset at infinity
#  . ddmodel == 4.2: 1/n dependence in speciation rate
done = 0
if(pars[7] > age)
{
   stop('The key innovation time is before the crown age of the main clade.')
}
if((pars[1] < pars[2]) | (pars[4] < pars[5]))
{
   stop('lambda0 is smaller than mu for one or both clades')
}
if(min(pars) < 0)
{
   stop('One of the parameters is negative')
}
if(!(ddmodel %in% c(1,1.3,2,2.1,2.2,2.3,3,4,4.1,4.2)))
{
   stop('This diversity-dependence model does not exist or is not implemented')
}
while(done == 0)
{
    # number of species N at time t
    # i = index running through t and N
    t = rep(0,1)
    L = matrix(0,2,5)
    i = 1
    t[1] = 0
    NM = 2
    NS = 0
    # L = data structure for lineages,
    # . L[,1] = branching times
    # . L[,2] = index of parent species
    # . L[,3] = index of daughter species
    # . L[,4] = time of extinction
    # . L[,5] = main clade (0) or subclade (1)
    # j = index running through L
    L[1,1:5] = c(0,0,-1,-1,0)
    L[2,1:5] = c(0,-1,2,-1,0)
    linlistM = c(-1,2)
    linlistS = NULL
    newL = 2
    tinn = age - pars[7]
    ff = dd_KI_lamuN(ddmodel,pars,c(NM[i],NS[i]))
    laMN = ff[1]
    muMN = ff[2]
    laSN = ff[3]
    muSN = ff[4]
    denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
    t[i + 1] = t[i] + stats::rexp(1,denom)
    if(t[i + 1] > tinn & t[i] < tinn)
    {
         NM[i] = NM[i] - 1
         NS[i] = NS[i] + 1
         linlistS = sample2(linlistM,1)
         L[abs(linlistS),5] = 1
         linlistM = linlistM[-which(linlistM == linlistS)]
         ff = dd_KI_lamuN(ddmodel,pars,c(NM[i],NS[i]))
         laMN = ff[1]
         muMN = ff[2]
         laSN = ff[3]
         muSN = ff[4]
         denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
         t[i + 1] = tinn + stats::rexp(1,denom)
    }
    while(t[i + 1] <= age)
    {
        event = sample2(x = 1:4,size = 1,prob = c(laMN * NM[i], muMN * NM[i], laSN * NS[i], muSN * NS[i]))
        i = i + 1
        if(event == 1)
        {
            # speciation event in main clade
            ranL = sample2(linlistM,1)
            NM[i] = NM[i - 1] + 1
            NS[i] = NS[i - 1]
            newL = newL + 1
            L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1,0))
            linlistM = c(linlistM,sign(ranL) * newL)
        } else if(event == 3)
        {
            # speciation event in subclade
            ranL = sample2(linlistS,1)
            NM[i] = NM[i - 1]
            NS[i] = NS[i - 1] + 1
            newL = newL + 1
            L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1,1))
            linlistS = c(linlistS,sign(ranL) * newL)           
        } else if(event == 2)
        {
            # extinction event in main clade
            ranL = sample2(linlistM,1)
            NM[i] = NM[i - 1] - 1
            NS[i] = NS[i - 1]
            L[abs(ranL),4] = t[i]
            w = which(linlistM == ranL)
            linlistM = linlistM[-w]
            linlistM = sort(linlistM)
        } else if(event == 4)
        {
            # extinction event in subclade
            ranL = sample2(linlistS,1)
            NM[i] = NM[i - 1]
            NS[i] = NS[i - 1] - 1
            L[abs(ranL),4] = t[i]
            w = which(linlistS == ranL)
            linlistS = linlistS[-w]
            linlistS = sort(linlistS)        
        }
        if(sum(c(linlistM,linlistS) < 0) == 0 | sum(c(linlistM,linlistS) > 0) == 0)
        {
            t[i + 1] = Inf
        } else {
            ff = dd_KI_lamuN(ddmodel,pars,c(NM[i],NS[i]))
            laMN = ff[1]
            muMN = ff[2]
            laSN = ff[3]
            muSN = ff[4]
            denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
            t[i + 1] = t[i] + stats::rexp(1,denom)
            if(t[i + 1] > tinn & t[i] < tinn)
            {
               NM[i] = NM[i] - 1
               NS[i] = NS[i] + 1
               ff = dd_KI_lamuN(ddmodel,pars,c(NM[i],NS[i]))
               laMN = ff[1]
               muMN = ff[2]
               laSN = ff[3]
               muSN = ff[4]
               linlistS = sample2(linlistM,1)
               L[abs(linlistS),5] = 1
               linlistM = linlistM[-which(linlistM == linlistS)]
               denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
               t[i + 1] = tinn + stats::rexp(1,denom)
            }
        }
    }
    if(sum(c(linlistM,linlistS) < 0) == 0 | sum(c(linlistM,linlistS) > 0) == 0)
    {
       done = 0
    } else {
       done = 1
    }
}

L[,1] = age - c(L[,1])
notmin1 = which(L[,4] != -1)
L[notmin1,4] = age - c(L[notmin1,4])
L[which(L[,4] == age + 1),4] = -1
tes = L2phylo(L[,1:4],dropextinct = T)
tas = L2phylo(L[,1:4],dropextinct = F)
tesS = NULL
tes2 = NULL
tasS = NULL
tas2 = NULL
graphics::par(mfrow = c(1,2))
graphics::par(mai = c(0,0.5,0.5,0))
graphics::plot(tes)
graphics::title('Reconstructed tree', cex.main = 0.7)
graphics::plot(tas)
graphics::title('Complete tree', cex.main = 0.7)
cols = c("blue","red")
names(cols) = c(0,1)
if(length(linlistS) > 0)
{
   namesS = paste('t',abs(linlistS), sep = "")
   if(length(linlistS) == 1)
   {
      m = which(tes$tip.label == namesS)
      b2 = 0
   }
   else if(length(linlistS) > 1)
   {
      m = ape::getMRCA(phy = tes,tip = namesS)
      tesS = ape::extract.clade(phy = tes,node = m)
      b2 = age - ape::node.depth.edgelength(tes)[m]
   }  
   m0 = tes$edge[which(tes$edge[,2] == m),1]
   b1 = age - ape::node.depth.edgelength(tes)[m0]
   tes2 = phytools::paintSubTree(tes,node = m,state = "1",anc.state = "0",stem = (pars[7] - b2)/(b1 - b2))
   phytools::plotSimmap(tes2,cols,lwd = 3,pts = F, mar = c(0.1,0.1,1,0.1))
   graphics::title('Reconstructed tree', cex.main = 0.7)
}
allS = which(L[,5] == 1)
if(length(allS) > 0)
{
   namesS = paste('t',abs(allS), sep = "")
   if(length(allS) == 1)
   {
      m = which(tas$tip.label == namesS)
      b2 = 0
   }
   else if(length(allS) > 1)
   {
      m = ape::getMRCA(phy = tas,tip = namesS)
      tasS = ape::extract.clade(phy = tas,node = m)
      b2 = age - ape::node.depth.edgelength(tas)[m]
   }
   m0 = tas$edge[which(tas$edge[,2] == m),1]
   b1 = age - ape::node.depth.edgelength(tas)[m0]
   tas2 = phytools::paintSubTree(tas,node = m,state = "1",anc.state = "0", stem = (pars[7] - b2)/(b1 - b2))
   phytools::plotSimmap(tas2,cols,lwd = 3,pts = F, mar = c(0.1,0.1,1,0.1)) 
   graphics::title('Complete tree', cex.main = 0.7)
}
out = list(tes = tes,tas = tas,L = L,tesS = tesS,tasS = tasS,tes2 = tes2,tas2 = tas2)
return(out)

}

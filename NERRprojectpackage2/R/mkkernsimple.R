#' Internal function to create kernel
#'
#' This is an internal function that creates the growth kernel that gets used in EPR_NERR_simple.  Not meant to be run on its own.
#' @param x This is a matrix x that comes from kernmatsimple
#' @param y This is a matrix y that comes from kernmatsimple
#' @param F This is the fishing rate, set to 0 in EPR_NERR_simple
#' @param M This is the background non-fishing mortality rate
#' @param T This is the "step size", set to 1 in EPR_NERR_simple
#' @param fixparm This is Params that has Params$M, Params$Linf, Params$k (only one value for each) but there are other components added in EPR_NERR_simple before Params goes into this function
#' @keywords NERR kernel
#' @export
#' @examples
#' mkkernsimple(x,y,F,fixparm,T)

mkkernsimple <- function(x,y,F,M,fixparm,T){

#% adapted from Will's code adapted from Easterling's code, used for PISCO
#% rockfish data set

#% extract parameters
isjuv_mean <- 70
isjuv_sd <- 2

growthvar <- 0.1

#%DEFINE WHO CAN GET FISHED
#%isjuv = 1 - normcdf(x,fixparm(4),diff(x(1,1:2))/2);
isjuv <- 1 - pnorm(x,isjuv_mean,isjuv_sd);

#%define pr(reproductive) using a maturity ogive
#%not needed now because open population
#%ismat = 1 - normcdf(x,fixparm(5),diff(x(1,1:2))/2);

#%SURVIVAL PART OF KERNEL
#%add in stochasticity for process error, but keep constant across size
#%Assume no stochasticity in fishing
#%mortality.
#M <- fixparm$M #%normrnd(fixparm(3), varparm(1)); %CHECK VALUE FOR STD
dims <- size(x)
m <-matrix(1,dims[1],dims[2])*M + (1-isjuv)*F #%this is a matrix size x,
#%mortality for each size

p1 <- exp(-m*T) #% convert mortality rate to survivorship, iterate over time steps

#%GROWTH PART OF KERNEL
#%add variability in growth to k
Linf = fixparm$Linf
k = fixparm$k
#%k = max(realmin,normrnd(fixparm(2),varparm(2)));
#%growth
pmean1=Linf - (Linf - x)*exp(-k) #% (do not add in x0 for the one-step growth)
#%add variability around von Bertalanffy growth
psig1 = pmean1*growthvar #% make this a parameter; %.*varparm(3);

#%evaluate growth part of kernel
p2 = dnorm(y, pmean1, psig1)

p1 = pmax(p1,0)    #%to make sure no negatives
p2 = pmax(p2,0)    #%to make sure no negatives

kxy = p1*p2
#%add fecundity to kernel if a closed population

return(kxy)
}

#' Internal function to produce SSD and EPR
#'
#' This is the internal function that does all the work of producing the EPR and SSD for each site that is fed back to run_EPR_NERR_simple.
#' @param Params This function is written such that it only takes one set of M, Linf, and k at a time for a given site.  Params is still a list
#' with Params$M, Params$Linf, Params$k, but they only contain one value each.  This function also requires Params$veclength and Params$x,
#' which are created in run_EPR_NERR_simple.
#' @keywords NERR
#' @export
#' @examples
#' EPR_NERR_simple(Params)

EPR_NERR_simple <- function(Params){

  # Do EPR calculations using IPM for GTM NERR project
  # Started by JW White, April 2019
  # Edited LS Storch Nov 2019 for NERR impending deadline of Dec 2019

  #WARNING - there are a lot of hardcoded variables/parameters in here, taken
  #from oyster_PP_params or EPR_NERR (each parameter says where it came from)
  #WARNING - there are also a lot of hardcoded variables in mkkernLS, taken from
  #oyster_PP_params

  #the time step in this function is 6 months - this is enforced by importing
  #6-month growth and mortality parameters into the function (right?)

  #function [EPR, SSD] = EPR_NERR_LS(Params)
  #this code has been restructured to only take in ONE set of k, Linf, M variables at a time, i.e., one site at a time.

  library(pracma)
  #setwd("~/Box/Storch_oyster_postdoc/code/NERR_Hybrid_Code_LS")
  #source("makeSimpVec.R")
  source("kernmatsimple.R")
  source("mkkernsimple.R")
  #temporarily added this in just to be able to test the script and see if it works, post-matlab conversion
  #Params <- list();
  #Params$k <- 1.2/6
  #Params$Linf <- 47.6
  #Params$M <- 3.5/6

  #Params is a struct with cols ""SiteName"", ""M"", ""k"", ""Linf"".  These
  #parameters aren't calculated in MATLAB, they're imported from R so you
  #have to access a .mat file which has them saved in this format, or use
  #process_oyster_stats.m which takes the data and puts them in this format
  #and then calls EPR_NERR_LS.  this code has been restructured to only take
  #in ONE set of variables at a time, i.e., one site at a time.

  #PARAMETERS - a lot of these were taken from oyster_PP_params or EPR_NERR
  Params$T <- 1 #time step
  Params$F <- 0 #harvest rate
  #Params$veclength <- 200 #number of patches in the vector, made 200 so the size bins aren't too coarse grain - this is an input now
  # Params$x <- linspace(0,160,Params$veclength) # size bins vector, max value is 160 because that's double the biggest Linf - this is an input now
  Params$dx <- diff(Params$x[1:2])
  Params$LW_afdw <- 5.09e-5*Params$x^2.365 # Ash Free Dry Weight; units = g ; derived from Kimbro field samples in salinity zone 2 (taken from oyster_PP_params)
  #Params$Mat <- 40 #maturity parameter, taken from oyster_PP_params
  Params$Fec <- 19.86e6*Params$LW_afdw^(1.17) #fecundity parameter, taken from oyster_PP_params
  Params$Fec = Params$Fec * (Params$x>=35) # Impose size at maturity of 35 mm
  Params$Rmean <- 7 # mean spat size in week 1 - taken from EPR_NERR
  Params$Rstd <- 5 # std in spat size - taken from EPR_NERR
  Params$Rvec <- dnorm(Params$x,Params$Rmean,Params$Rstd) #taken from EPR_NERR, initializing recruitment vector
  #Params$DD <- 1e-3 # based on Puckett & Eggleston (2012), Fig 7. (taken from oyster_PP_params)"
  #Params$j_length <- 15 #juvenille mortality size limit, taken from oyster_PP_params"

  #S <- matrix(data = NA,nrow=1,ncol=Params$T) #S <- nan(1,Params$T) # reef structure, currently unused
  #N <- matrix(data = NA, nrow = Params$veclength, ncol = Params$T) #N <- nan(Params$veclength,Params$T) # oyster population size matrix (semi-annual time steps, based on incoming parameters M and k)
  #R <- matrix(data = NA, nrow = Params$veclength, ncol = Params$T) #R <- nan(Params$veclength,Params$T) # recruitment matrix
  #Harv <- matrix(data = NA, nrow = Params$veclength, ncol = Params$T) #Harv <- nan(Params$veclength,Params$T)# oysters harvested (semiannual), currently unused
  #Params$lambdaTAF <- 0.1/52 # Based on annual rate of 0.1 from Powell et al. (2012)
  #Params$density <- 0.849 # grams/cm^3. From aqua-calc.com

  # Initialize the size distributions
  #N[,1] <- 1
  #S[1] <- 1
  #R[,1] <- 1

  #create the kernels for fecundity, growth, and mortality using "fancy"
  #mkkernLS function, taken from ""fancy"" mkkern in oyster_IPM_statespace
  #folder (mkkernLS called from kernmatSimp)
  kmat <- kernmatsimple(Params$x,Params$F,Params$M,Params,Params$T)# get kernel
  kmatR <- kernmatsimple(Params$x,Params$F,Params$Mjuv,Params,Params$T)# get kernel



  Rdist <- Params$Rvec

  #% Ensure that Rdist integrates to 1:
  Rdist <- Rdist/(sum(Rdist*Params$dx));

  #% Find stable size distribution:
  TT <- 100 #how many times to iterate?
  N <- repmat(cbind(Rdist),1,TT);
  R1 <- N;
  R2 <- R1
  R3 <- R2
  A <- N

  for (t in 2:TT){

    # Is there a way to parse this out with a 3-month delay for the juvenile mortality?
    # That would be great
    R1[,t] = cbind(Rdist)
    R2[,t] = (kmatR%*%R1[,t-1])*Params$dx
    R3[,t] = (kmatR%*%R2[,t-1])*Params$dx

    A[,t] = (kmat%*%A[,t-1])*Params$dx + R3[,t-1]

    N[,t] =  A[,t] + R1[,t] +  R2[,t] +  R3[,t]
  }

  #% Stable size distribution
  SSD <- N[,TT];

  #% Calcluate EPR (relative, at this point) (this is just biomass at this
  #% point, fix to be eggs) (does not yet include possibility of sex change
  EPR <- sum(Params$Fec*SSD)

  result <- list(SSD,EPR) #R doesn't allow multiple outputs for a function WHATTTTTT?!?!?!
  return(result)

}







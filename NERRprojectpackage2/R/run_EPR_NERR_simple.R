#' NERR function to produce SSD and EPR
#'
#' Given a set of parameters per site, Mjuv (juvenile mortality), Madult (post-juvenile mortality), Linf (max size), and k(growth rate), this function will produce
#' the EPR and stable size distribution per site and graph the results with the 5th and 95th percentile uncertainties
#' as upper and lower bounds
#'
#' @param Params Params is a list with Params$SiteNames, Params$Mjuv, Params$Madult, Params$k, Params$Linf for all sites.  These can be uploaded by loading Params.Rdata
#' @param CovMats CovMats is a list with the covariance matrices for Mjuv, Madult, Linf, and k (in that order) for the sites.  These can be uploaded by loading CovMats.Rdata
#' @param n n is the number of random parameter draws to create uncertainty bounds
#' @keywords NERR
#' @export
#' @examples
#' run_EPR_NERR_simple(Params,CovMats)
#'
run_EPR_NERR_simple <- function(Params,CovMats,n=1000){

  #LSS December 2019
  #making another function that can deal with the struct of Params for the 7 sites since I restructured EPR_NERR_LS to only take in the parameters for one site at a time
  #here is also where the uncertainty analysis is done

  #the data for this has to be a struct (list) of params: Params$k, Params$Linf, Params$M , Params$sites and also the covariance matrices for each site, CovMats
  #NEW - Params$veclength is now also an input!!!!!
  library(MASS)
  library(ggplot2)
  library(egg)
  library(grDevices)

  # Check whether these 'source' commands are necessary - I think not????
  #load('NERRdata.Rdata')
  source("EPR_NERR_simple.R")
  source("kernmatsimple.R")
  source("mkkernsimple.R")
  # The following should be moved to an input option with default = 200
  Params$veclength <- 200 #moved this to the outside of the function because lower down I have to specify the length of the vectors storing the SSD outputs
  Params$x <- seq(from=0,to=160,length.out=Params$veclength) # space vector, max value is 160 because that's double the biggest Linf

  #everything will be run through EPR_NERR_LS, but we first need to generate n number of random combinations of variables using the covariance matrices to get an estimate of the uncertainty
  #n = 1000; #try 1000 random combos of parameters - make sure n is even

  #EDIT THIS TO CHUCK OUT NEGATIVE VALUES
  RandParams <- vector('list',length(Params$sites))
  for (i in 1:length(Params$sites)){ #need n number of random combos for each site
    invec <- c(Params$M[i],Params$Linf[i],Params$k[i]) #order of params in the covariance matrices is M, Linf, k, so mimic this in the order of params going into mvrnorm
    out <- mvrnorm(n, invec, CovMats[[i]][[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    #out <- mvrnorm(n, invec, matrix(0,3,3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE) #testing, delete later
    RandParams[[i]] <- out
  }

  # Note that at present we have good estimates of juvenile mortality (1st 3 months), but poor estimates of later mortality.
  # As a first approximately, older mortality rates appear to be ~10% of juvenile mortality rates.

  #Now that we have all of the parameters and random parameter combinations, have to multiply M and k by 30 because we have a daily rate and want a monthly
  #Column 1 of RandParams[[i]] is M and column 3 is k
  for (i in 1:length(Params$sites)){
    Params$M[i] <- Params$M[i]*30
    Params$k[i] <- Params$k[i]*30
    RandParams[[i]][,1] <- RandParams[[i]][,1]*30
    RandParams[[i]][,3] <- RandParams[[i]][,3]*30
  }

  #get the outputs for the 1000 random parameter combos for each site and find the upper and lower 5 and 95 percentile
  RandSSD <- vector('list',length(Params$sites)) #save all the random SSDs
  RandEPR <- vector('list',length(Params$sites)) #save all the random EPRs
  SSDmeans <- vector('list',length(Params$sites)) #save the means for all of the random SSDs, using for finding 5th and 95th percentile
  EPR_highlow <- vector('list',length(Params$sites)) #save 5th and 95th percentile EPR
  SSD_highlow <- vector('list',length(Params$sites)) #save 5th and 95th percentile SSD
  SSD_means_highlow <- vector('list',length(Params$sites)) #save 5th and 95th percentile SSD
  for (i in 1:length(Params$sites)){ #outer loop goes through the sites

    RandEPR[[i]] <- vector('numeric',n)
    SSDmeans[[i]] <- vector('numeric',n)
    EPR_index <- vector('numeric',2)
    RandSSD[[i]] <- matrix(0,Params$veclength,n) #column-wise is each SSD output for each of the random parameter combos
    SSD_highlow[[i]] <- matrix(0,Params$veclength,2) #save 5th and 95th percentile SSDs
    miniParams <- list()

    for (j in 1:n){ #inner loop goes through the n random parameter combos for each site

      miniParams$Mjuv <- RandParams[[i]][j,1]
      miniParams$M <- RandParams[[i]][j,1]*0.1 # notice reduction factor here
      miniParams$Linf <- RandParams[[i]][j,2]
      miniParams$k <- RandParams[[i]][j,3]
      miniParams$veclength <- Params$veclength
      miniParams$x <- Params$x
      out <- EPR_NERR_simple(miniParams)
      if (!is.na(out[[2]])){
        RandEPR[[i]][j] <- out[[2]]
      }
      if (!anyNA(out[[1]])){
        RandSSD[[i]][,j] <- cbind(out[[1]])
        SSDmeans[[i]][j] <- apply(cbind(RandSSD[[i]][,j]),2,mean)
      }

    }

    #95 percent confidence interval (WW redid this using quantile() to avoid sample size issue)
    #u95 <- n*0.95 #just make sure n is even or this wont be an integer
    #l95 <- n*0.05 #same as above
    #EPRsort <- sort(RandEPR[[i]])
    #EPR_highlow[[i]][1] <- EPRsort[u95]
    #EPR_highlow[[i]][2] <- EPRsort[l95]
    EPR_highlow[[i]][1] <- quantile(RandEPR[[i]],0.95,na.rm=TRUE)
    #EPR_highlow[[i]][2] <- EPRsort[l95]
    EPR_highlow[[i]][2] <- quantile(RandEPR[[i]],0.05,na.rm=TRUE)
    EPR_diff <- abs(RandEPR[[i]] - EPR_highlow[[i]][1])
    EPR_index[1] <- which(EPR_diff == min(EPR_diff))
    EPR_diff <- abs(RandEPR[[i]] - EPR_highlow[[i]][2])
    EPR_index[2] <- which(EPR_diff == min(EPR_diff))

    #SSD_means_sort <- sort(SSDmeans[[i]],index.return = TRUE) #the indices will be in SSD_means_sort$ix (notice the i in front)
    SSD_highlow[[i]][,1] = RandSSD[[i]][,EPR_index[1]] #grab the SSD from the 95th
    SSD_highlow[[i]][,2] = RandSSD[[i]][,EPR_index[2]] #grab the SSD from the 5th
    #SSD_means_highlow[[i]][1] = SSD_means_sort[u95];
    #SSD_means_highlow[[i]][2] = SSD_means_sort[l95];
    SSD_means_highlow[[i]][1] <- SSDmeans[[i]][EPR_index[1]]
    SSD_means_highlow[[i]][2] <- SSDmeans[[i]][EPR_index[2]]

  }

  #now that we have the 5th and 95th percentile, lets get the actual SSD and EPRs using the actual parameters so we can plot them with the 5th and 95th
  EPR <- vector('numeric',length(Params$sites))
  SSD <- vector('list',length(Params$sites))
  maxSSD_highlow <- vector(length=length(Params$sites))
  for (i in 1:length(Params$sites)){ #outer loop goes through the sites
    miniParams$Mjuv <- Params$M[i]
    miniParams$M <- Params$M[i]*0.1 # notice reduction factor here
    miniParams$Linf <- Params$Linf[i]
    miniParams$k <- Params$k[i]
    miniParams$veclength <- Params$veclength
    out <- EPR_NERR_simple(miniParams)
    EPR[i] <- out[[2]]
    SSD[[i]] <- out[[1]]

    # also find max value of SSD_highlow, for plotting
    maxSSD_highlow[i] = max(SSD_highlow[[i]][,1])
  }

  #Plotting - we will want to plot a given SSD[[i]] against SSD_highlow[[i]][,1] and SSD_highlow[[i]][,2] and compare EPR[i] against EPR_highlow[[i]][1] and EPR_highlow[[i]][2]
  # Loop over each site, create a ggobject, then plot them up in a grid
  GGp <- list() # pre-allocate
  Ymax <- max(maxSSD_highlow)*1.001
  for (i in 1:length(Params$sites)){

    # create a data frame for plotting SSD
    Data.sub <- data.frame(SSD=SSD[[i]],Length=Params$x,Hi=SSD_highlow[[i]][,1],Lo=SSD_highlow[[i]][,2])

    GGp[[i]] <- ggplot(data=Data.sub,aes(y=SSD,x=Length))+
      geom_line()+
      geom_line(aes(x=Length,y=Lo))+
      geom_line(aes(x=Length,y=Hi))+
      geom_ribbon(aes(x=Length,ymin=Lo,ymax=Hi),color='blue',fill='blue',alpha=0.5)+
      xlab('Length (mm)')+
      scale_x_continuous(limits=c(0,100))+
      scale_y_continuous(breaks=NULL,limits=c(0,Ymax))+
      ggtitle(Params$sites[i])+
      theme_bw()


  } # end loop over sites



  # Separately, plot EPR:
  EPR.sub <- data.frame(Site=Params$sites,EPR=EPR,EPRmin = NA,EPRmax=NA)
  for (i in 1:length(Params$sites)){
    EPR.sub$EPRmax[i] = EPR_highlow[[i]][1]
    EPR.sub$EPRmin[i] = EPR_highlow[[i]][2]
  }
  # rescale
  Max <- max(EPR.sub$EPRmax)
  EPR.sub$EPR = EPR.sub$EPR/Max
  EPR.sub$EPRmin = EPR.sub$EPRmin/Max
  EPR.sub$EPRmax = EPR.sub$EPRmax/Max


  EPRgg <- ggplot(data=EPR.sub,aes(x=Site,y=EPR))+
    geom_linerange(aes(ymin=EPRmin,ymax=EPRmax))+
    geom_point()+
    scale_x_discrete(limits=c('T','G','St. A','SR','B','M','P'))+
    xlab(NULL)+
    ylab('Relative eggs per oyster recruit')+
    theme_bw()

  # Display all results
  ggarrange(GGp[[1]],GGp[[2]],GGp[[3]],GGp[[4]],GGp[[5]],GGp[[6]],GGp[[7]],EPRgg, ncol = 2)


}



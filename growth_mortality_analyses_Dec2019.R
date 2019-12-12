# GROWTH ####
setwd("~/Box Sync/Storch_oyster_postdoc/St. Augustine Oysters- OSU")
library(ggpubr)
library(tidyverse)
library(readxl)
library(egg)
library(nlme)

#---------------------------------------------------------
# 1) Survival analysis. From 2019 outplant
# Load in data.
data = read_xlsx("2019_outplant(updated3dec).xlsx",sheet="oysters")

# clean up data formatting
data$reef.id = as.factor(data$reef.id) # which zone within GTM
data$treatment = as.factor(data$treatment) # cage, control, cagecontrol
data$meter = as.factor(data$meter) # within-site replication (at which transect meter was the experimental unit)

# define some organizational variables
zones = levels(data$reef.id)

# strip weird NA rows & bad values (zeros)
NArows = data$reef.id==""
data = data[!NArows,]
Zrows  =data$size == 0
data = data[!Zrows,]



# Compile data for analysis by GLM. Take count of number of remaining oysters 
# on each data, at each replicate experimental unit
data$alive <- data$status==1
data.agg <- aggregate(data$alive,by=list(data$reef.id,data$treatment,data$timedays,data$meter),FUN=sum)
colnames(data.agg) <- c('reef.id','treatment','timedays','meter','alive')

# pre-allocate variables for storing results
sitelist <- vector()
optionlist <- vector()
M <- vector()
Ms <- matrix(nrow=length(zones),ncol=2)
intercept <- vector()
Mod_mor <- list()
Mor_cov <- vector()
k <- 1 # this is a counter variable for use in a loop

treatments <- c('cage','control','cagecontrol') # for this analysis we do not consider the cagecon treatment
treatnames <- c('cage','con','cagecon') # short names

# Loop over Zones & Treatments to estimate survival curve
for (i in 1:length(zones)){ 
  for (j in 1:2){ #length(treatments)){ 
    
    data.sub <- subset(data.agg,reef.id == zones[i] & treatment == treatments[j])
    m = glm(alive ~ timedays, family=poisson(link='log'),data = data.sub)
    Mod_mor[[k]] <- m
    Mor_cov[k] <- vcov(m)[2,2]
    #mod_mortality[[i,j]] <- m
    sitelist[k] <- zones[i]
    optionlist[k] <- treatnames[j]
    M[k] <- coef(m)[2] # for export to csv
    Ms[i,j] <- coef(m)[2] # for use in plotting later
    intercept[k] <- coef(m)[1]
    k = k+1
    
  }
}

NERRmortalitycoeffs <- data.frame(sitelist,optionlist,M,intercept,Mor_cov)
write.csv(NERRmortalitycoeffs,'NERRmortalitycoeffs2019dec.csv')

# Plot results:
# Need to rescale data to proportion of initial density

Gp <- list()# variable for holding ggplot objects
Gp2 <- list()
Colors <- c('green','magenta','blue','purple','yellow','red','cyan')
Zone_names = c('Butler','Guana R.','Matanzas R.','Pellicer','Salt R.','St. Augustine','Tolomato R.')

for (i in 1:length(zones)){
  
  data.sub <- subset(data.agg,reef.id == zones[i] & (treatment=='control' | treatment==
                                                       'cage'))
  meters = levels(data.sub$meter)
  for (m in 1:length(meters)){
    # get number alive in that unit on day 0
    con0 = data.sub$alive[data.sub$meter==meters[m] & data.sub$timedays==0 & data.sub$treatment=='control']
    cage0 = data.sub$alive[data.sub$meter==meters[m] & data.sub$timedays==0 & data.sub$treatment=='cage']

    data.sub$alive[data.sub$meter==meters[m]& data.sub$treatment=='control']=data.sub$alive[data.sub$meter==meters[m]& data.sub$treatment=='control']/con0
    data.sub$alive[data.sub$meter==meters[m]& data.sub$treatment=='cage']=data.sub$alive[data.sub$meter==meters[m]& data.sub$treatment=='cage']/cage0
    
    } # end loop over meters
  
  # dummy data for curves
  timedays <- 0:max(data.sub$timedays)
  alive.con <- exp(Ms[i,2]*timedays)
  alive.cage <- exp(Ms[i,1]*timedays)
  dummy.data <- data.frame(timedays,alive.con,alive.cage)
  
  Gp[[i]] <- ggplot(data=data.sub,aes(x=timedays,y=alive))+
             geom_jitter(aes(shape=treatment),width=2,show.legend = FALSE,color=Colors[i])+
             scale_shape_manual(values = c(1,16))+
             geom_line(data=dummy.data,aes(x=timedays,y=alive.con),color='blue',lty=2)+
             geom_line(data=dummy.data,aes(x=timedays,y=alive.cage),color='blue',lty=1)+
             xlab('Time (d)')+
             ylab('Proportion surviving')+
             ggtitle(Zone_names[i])
  
  # create a separate plot with cages only
  data.sub2 = subset(data.sub,treatment='cage')
  
  Gp2[[i]] <- ggplot(data=data.sub2,aes(x=timedays,y=alive))+
    geom_jitter(shape=16,width=2,color=Colors[i])+
    geom_line(data=dummy.data,aes(x=timedays,y=alive.cage),color='blue',lty=1)+
    xlab('Time (d)')+
    ylab('Proportion surviving')+
    ggtitle(Zone_names[i])
  
} # end loop over zones

# Arrange zones in desired order:
# 7 (tolo), 2 (guana), 6 (staug), 5 (salt), 1 (butler), 3 (matanzas), 4 (pellicer)
quartz(width=6,height=6)
ggarrange(Gp[[7]],Gp[[2]],Gp[[6]],Gp[[5]],Gp[[1]],Gp[[3]],Gp[[4]], ncol = 3)

quartz(width=6,height=6)
ggarrange(Gp2[[7]],Gp2[[2]],Gp2[[6]],Gp2[[5]],Gp2[[1]],Gp2[[3]],Gp2[[4]], ncol = 3)

#---------------------------------
# Revisit this to split mortality into early vs. late
# 1.1) Do 'post-juvenile' mortality using later dates:
# Load in data.
data = read_xlsx("2019_outplant(updated3dec).xlsx",sheet="oysters")

# clean up data formatting
data$reef.id = as.factor(data$reef.id) # which zone within GTM
data$treatment = as.factor(data$treatment) # cage, control, cagecontrol
data$meter = as.factor(data$meter) # within-site replication (at which transect meter was the experimental unit)

# define some organizational variables
zones = levels(data$reef.id)

# strip weird NA rows & bad values (zeros)
NArows = data$reef.id==""
data = data[!NArows,]
Zrows  =data$size == 0
data = data[!Zrows,]

# Compile data for analysis by GLM. Take count of number of remaining oysters 
# on each data, at each replicate experimental unit
data$alive <- data$status==1
data.old <- subset(data,timedays>=75)
# Rescale time so that it runs from zero onwards at each site
for (i in 1:length(zones)){
 min.temp <- subset(data.old,reef.id=zones[i])
 min.date <- min(min.temp$timedays)
 data.old$timedays = data.old$timedays-min.date
}
data.agg <- aggregate(data.old$alive,by=list(data.old$reef.id,data.old$treatment,data.old$timedays,data.old$meter),FUN=sum)
colnames(data.agg) <- c('reef.id','treatment','timedays','meter','alive')

# Do GLMs:
# pre-allocate variables for storing results
sitelist <- vector()
optionlist <- vector()
M <- vector()
Ms <- matrix(nrow=length(zones),ncol=2)
intercept <- vector()
Mod_mor <- list()
Mor_cov <- vector()
k <- 1 # this is a counter variable for use in a loop

treatments <- c('cage','control','cagecontrol') # for this analysis we do not consider the cagecon treatment
treatnames <- c('cage','con','cagecon') # short names

# Loop over Zones & Treatments to estimate survival curve
for (i in 1:length(zones)){ 
  for (j in 1:2){ #length(treatments)){ 
    
    data.sub <- subset(data.agg,reef.id == zones[i] & treatment == treatments[j])
    m = glm(alive ~ timedays, family=poisson(link='log'),data = data.sub)
    Mod_mor[[k]] <- m
    Mor_cov[k] <- vcov(m)[2,2]
    #mod_mortality[[i,j]] <- m
    sitelist[k] <- zones[i]
    optionlist[k] <- treatnames[j]
    M[k] <- coef(m)[2] # for export to csv
    Ms[i,j] <- coef(m)[2] # for use in plotting later
    intercept[k] <- coef(m)[1]
    k = k+1
    
  }
}

NERRmortalitycoeffs_adult <- data.frame(sitelist,optionlist,M,intercept,Mor_cov)
write.csv(NERRmortalitycoeffs_adult,'NERRmortalitycoeffs2019dec_adult.csv')

# As of Dec 2019, this estimation is problematic -- too little variation in the data, too few samples after 75 d.
# For results with plausible values, adult mortality appears to be ~10% of juvenile mortalilty.

#---------------------------------



#---------------------------------------------------------
# 1) Growth analysis. From 2018 outplant
# Load in data.
data = read_xlsx("2018_outplant(updated_3dec).xlsx",sheet="cluster")

# clean up data formatting
data$site = as.factor(data$site) # specific study sites within each zone
data$zone = as.factor(data$zone) # zones within GTM
data$cl.no = as.factor(data$cl.no) # within-site replication (each oyster cluster)
data$time.days = data$`time(days)`

# define some organizational variables
zones = levels(data$zone)

# Fit nonlinear von Bertalanffy growth curves 
# allocate vectors for results
Linf <- rep(NA,length(zones)) 
k <- Linf; t0 <- Linf
Mod_growth <- list()
Growth_cov <- list()

# Need to provide initial guesses for the model parameters; these can vary by site
Linf.starts = rep(80,length(zones)); Linf.starts[4]=10 # try different starting values for different sites

# Dummy variables for curve-plotting purposes
age.dummy = 0:max(data$time.days)
L.dummy = matrix(0,nrow=length(age.dummy),ncol=length(zones))

# Loop over zones to get fits:
for (z in 1:length(zones)){ 
  if (z != 4){# pellicer zone doesn't have enough data to fit
    
    data.sub = subset(data,zone==zones[z])
    data.sub = data.sub[!is.na(data.sub$size),] # trim out NAs
    
    # option 1: nonlinear least squares.
    # The vB equation is L(t) = Linf * (1 - exp(-k*(t - t0)))
    m = nls(formula = data.sub$size ~ Linf*(1 - exp(-k*(data.sub$time.days-t0))), 
            start = list(Linf=Linf.starts[z],k=0.1,t0=0),data=data.sub,control=nls.control(minFactor=1/(2^15),maxiter = 200))
    
    # option 2: nonlinear mixed effects. Had issues with convergence at some sites (and saw little
    # difference in fit in sites that did converge), so did not use this method
    #f1 <- data.sub$size ~ Linf*(1-exp(-(k*(data.sub$time.days. - t0))))
    #m1 = nlme(f1,data=as.data.frame(data.sub),start=coef(m),fixed=k+Linf+t0~1,random=cl.no~1,na.action=na.omit)#,group=~reefid,na.action=na.omit)  
    
    Mod_growth[[z]] <- m
    Growth_cov[[z]] <- vcov(m)[1:2,1:2]
    Linf[z] = coef(m)[1]
    k[z] = coef(m)[2]
    t0[z] = coef(m)[3]}
  
  # Results for plotting:
  L.dummy[,z] = Linf[z]*(1-exp(-(k[z]*(age.dummy-t0[z]))))  }

#temporarily commented out by Laura
#L.dummy[110:405,5] = NA # blank out some of the Salt R harvest data past the point of all oysters dying

# Plot data &  growth curves:
Gp <- list()
Colors <- c('green','magenta','blue','purple','black','yellow','red','black','cyan')
Zone_names = c('Butler','Guana R.','Matanzas R.','Pellicer','Salt R. harvest','Salt R.','St. Augustine','Tolomato R. harvest','Tolomato R.')
for (z in 1:length(zones)){ 
  
  # subset data
  data.sub = subset(data,zone==zones[z])
  data.dummy = data.frame(x = age.dummy,y=L.dummy[,z])
  
  Gp[[z]]<-ggplot(data.sub,aes(x=time.days,y=size))+
    geom_jitter(color=Colors[z])+
    geom_line(data=data.dummy,aes(x=x,y=y),color='blue')+
    ylim(c(0,80))+
    xlab('Days')+
    ylab('Length (mm)')+
    ggtitle(Zone_names[z])
  theme_bw()}

# How to order sites:
# 9 (tolo.non), 8 (tolo.harvest),2 (guana), 7 (staug), 6 (salt.non), 5 (salt.harvest), 1 (butler), 3 (matanzas), 4 (pellicer)
#Order = c(9,8,2,7,6,5,1,3,4)
quartz(width=5,height=3)
ggarrange(Gp[[9]],Gp[[8]],Gp[[2]],Gp[[7]],Gp[[6]],Gp[[5]],Gp[[1]],Gp[[3]],Gp[[4]], ncol = 3)

# 7 (tolo), 2 (guana), 6 (staug), 5 (salt), 1 (butler), 3 (matanzas), 4 (pellicer)
Order = c(7,2,6,5,1,3,4)
quartz(width=6,height=6)
ggarrange(Gp[[7]],Gp[[2]],Gp[[6]],Gp[[5]],Gp[[1]],Gp[[3]],Gp[[4]], ncol = 3)


# Fit vB and get fitted line
Linf <- rep(NA,length(Zones)) # allocate vectors for results
k <- Linf; t0 <- Linf
# Put in some fake guesses for Pellicer?

#Linf.starts = rep(80,length(Zones)); Linf.starts[4]=10 # try different starting values for different sites
age.dummy = 0:105
L.dummy = matrix(0,nrow=length(age.dummy),ncol=length(Zones))
Mod_growth <- list()
for (z in 1:length(Zones)){ 
  #if (z != 4){# pellicer zone doesn't have enough data to fit
  data.sub = subset(data,reefid==Zones[z] & treatment == Treat[1]) #just trying to do cage for now since 2018 data file doesnt have treatment
  
#replace nls with nlme so you can have random effect of meter   
f1 <- data.sub$size ~ Linf*(1-exp(-(k*(data.sub$timedays - t0))))
m = nlme(f1,data=data.sub,start=list(k=0.01,Linf=100,t0=0),fixed=k+Linf+t0~1,random=meter~1,na.action=na.omit)#,group=~reefid,na.action=na.omit)
#m = nls(formula = data.sub$size ~ Linf*(1 - exp(-k*(data.sub$timedays -t0))), 
        #start = list(Linf=50,k=0.01,t0=0),data=data.sub,control=nls.control(minFactor=1/(2^15),maxiter = 500))
Mod_growth[[z]] <- m
Linf[z] = coef(m)[1]
k[z] = coef(m)[2]
t0[z] = coef(m)[3]#}
L.dummy[,z] = Linf[z]*(1-exp(-(k[z]*(age.dummy-t0[z]))))  }

#L.dummy[110:405,5] = NA # blank out some of the Salt R harvest data past the point of all oysters dying

# Plot vs. time to estimate growth, including growth curves:
Gp <- list()
Colors <- c('green','magenta','blue','purple','black','yellow','red','black','cyan')
Zone_names = c('Butler','Guana R.','Matanzas R.','Pellicer','Salt R.','St. Augustine','Tolomato R.')


for (z in 1:length(Zones)){ 
  data.sub = subset(data,reefid==Zones[z] & treatment == Treat[1])
  data.dummy = data.frame(x = age.dummy,y=L.dummy[,z])
Gp[[z]]<-ggplot(data.sub,aes(x=timedays,y=size))+
  geom_jitter(color=Colors[z])+
  geom_line(data=data.dummy,aes(x=x,y=y),color='blue')+
  ylim(c(0,70))+
  xlab('Days')+
  ylab('Length (mm)')+
  ggtitle(Zone_names[z])
  theme_bw()}

# How to order sites:
# 9 (tolo.non), 8 (tolo.harvest),2 (guana), 7 (staug), 6 (salt.non), 5 (salt.harvest), 1 (butler), 3 (matanzas), 4 (pellicer)
Order = c(7,2,6,5,1,3,4)
quartz(width=5,height=3)
ggarrange(Gp[[7]],Gp[[2]],Gp[[6]],Gp[[5]],Gp[[1]],Gp[[3]],Gp[[4]], ncol = 2)





# # SURVIVAL ####
# data = read.csv("outplant.data_updated2sept.oysters.csv")
# data$elapsed[data$elapsed<0]=0
# data$meter = factor(data$meter)
# 
# # remove anything with "2" in it from status
# data1 = data%>%
#   filter(status !=2) %>%
#   group_by(reef.id, elapsed, treatment, meter)%>%
#   summarize(
#     avg.surv = mean(status, na.rm = TRUE)
#   )
# 
# 
# # Plot vs. time to estimate growth:
# ggplot(data1,aes(x=elapsed,y=avg.surv))+
#   geom_jitter(aes(color=reef.id))+
#   theme_bw()
# 
# ZonesM = levels(data1$reef.id)
# # Fit exp(-M) and get fitted line
# M <- rep(0,length(ZonesM)) # allocate vector for results
# age.dummy = 0:105
# Surv.dummy = matrix(0,nrow=length(age.dummy),ncol=length(ZonesM))
# Mod_survival <- list()
# for (z in 1:length(ZonesM)){ 
#     data.sub = subset(data1,reef.id==ZonesM[z]&treatment=='control')
#     m = nls(formula = data.sub$avg.surv ~ exp(-M*elapsed), 
#             start = list(M=0.2),data=data.sub)
#     M[z] = coef(m)[1]
#     Mod_survival[[z]] <- m
#   Surv.dummy[,z] = exp(-M[z]*age.dummy)  }
# 
# # Fit exp(-M) and get fitted line (CAGE data)
# Mc <- rep(0,length(ZonesM)) # allocate vector for results
# Modc_survival <- list()
# age.dummy = 0:105
# Surv.dummy.c = matrix(0,nrow=length(age.dummy),ncol=length(ZonesM))
# for (z in 1:length(ZonesM)){ 
#   data.sub = subset(data1,reef.id==ZonesM[z]&treatment=='cage')
#   m = nls(formula = data.sub$avg.surv ~ exp(-M*elapsed), 
#           start = list(M=0.2),data=data.sub)
#   Modc_survival[[z]] <- m
#   Mc[z] = coef(m)[1]
#   Surv.dummy.c[,z] = exp(-Mc[z]*age.dummy)  }
# 
# # Plot vs. time to estimate growth, including growth curves (CONTROL ONLY):
# Gpm <- list()
# Colors <- c('green','magenta','blue','purple','yellow','red','cyan')
# Zone_names = c('Butler','Guana R.','Matanzas R.','Pellicer','Salt R.','St. Augustine','Tolomato R.')
# for (z in 1:length(ZonesM)){ 
#   data.sub = subset(data1,reef.id==ZonesM[z]&treatment=='control')
#   data.dummy = data.frame(x = age.dummy,y=Surv.dummy[,z])
#   Gpm[[z]]<-ggplot(data.sub,aes(x=elapsed,y=avg.surv))+
#     geom_jitter(color=Colors[z])+
#     geom_line(data=data.dummy,aes(x=x,y=y),color='blue')+
#     ylim(c(0,1))+
#     xlab('Days')+
#     ylab('Average survival')+
#     ggtitle(Zone_names[z])
#     theme_bw()}
# 
# # How to order sites:
# # 7 (tolo), 2 (guana), 6 (staug), 5 (salt), 1 (butler), 3 (matanzas), 4 (pellicer)
# Order = c(7,2,6,5,1,3,4)
# quartz(width=5,height=3)
# ggarrange(Gpm[[7]],Gpm[[2]],Gpm[[6]],Gpm[[5]],Gpm[[1]],Gpm[[3]],Gpm[[4]], ncol = 3)
# 
# # Plot vs. time to estimate growth, including growth curves (CAGE + CONTROL):
# Gpmc <- list()
# Colors <- c('green','magenta','blue','purple','yellow','red','cyan')
# Zone_names = c('Butler','Guana R.','Matanzas R.','Pellicer','Salt R.','St. Augustine','Tolomato R.')
# for (z in 1:length(ZonesM)){ 
#   data.sub1 = subset(data1,reef.id==ZonesM[z]&treatment=='control')
#   data.sub2 = subset(data1,reef.id==ZonesM[z]&treatment=='cage')
#   data.dummy = data.frame(x = age.dummy,y=Surv.dummy[,z])
#   data.dummyc = data.frame(x = age.dummy,y=Surv.dummy.c[,z])
#   Gpmc[[z]]<-ggplot(data.sub1,aes(x=elapsed,y=avg.surv))+
#     geom_jitter(color=Colors[z])+
#     geom_line(data=data.dummy,aes(x=x,y=y),color='blue')+
#     geom_jitter(data=data.sub2,aes(x=elapsed,y=avg.surv),color=Colors[z],shape=21)+
#     geom_line(data=data.dummyc,aes(x=x,y=y),color='blue',linetype=2)+
#     ylim(c(0,1))+
#     xlab('Days')+
#     ylab('Average survival')+
#     ggtitle(Zone_names[z])
#   theme_bw()}
# 
# # How to order sites:
# # 7 (tolo), 2 (guana), 6 (staug), 5 (salt), 1 (butler), 3 (matanzas), 4 (pellicer)
# Order = c(7,2,6,5,1,3,4)
# quartz(width=5,height=3)
# ggarrange(Gpmc[[7]],Gpmc[[2]],Gpmc[[6]],Gpmc[[5]],Gpmc[[1]],Gpmc[[3]],Gpmc[[4]], ncol = 3)
# 
# #Save all the coefficients and covariance matrices to transport over to Matlab
# Growth_Cell <- list()
# Coeff_Growth_Cell <- list()
# for (z in 1:length(Mod_growth)){
#     {if(z!=4) #the 4th is that site that starts with a p and there isn't a model fit 
#   Growth_Cell[[z]] <- vcov(Mod_growth[[z]])
#   Coeff_Growth_Cell[[z]] <- coef(Mod_growth[[z]])
#     }
# }
# 
# Survival_Cell <- list()
# CSurvival_Cell <- list()
# Coeff_Survival_Cell <- list()
# Coeff_CSurvival_Cell <- list()
# for (z in 1:length(Mod_survival)){
#     Survival_Cell[[z]] <- vcov(Mod_survival[[z]])
#     Coeff_Survival_Cell[[z]] <- coef(Mod_survival[[z]])
#     CSurvival_Cell[[z]] <- vcov(Modc_survival[[z]])
#     Coeff_CSurvival_Cell[[z]] <- coef(Modc_survival[[z]])
# }

  
  


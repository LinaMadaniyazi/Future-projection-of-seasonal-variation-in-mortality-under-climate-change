################################################################################
# R code for an example of the ALTERNATIVE approach in:
# "Projecting seasonality of mortality under climate change"

# The example datasets are adapted from the original data available from
# Epidemiology. 2019 May; 30(3): 321–329. doi: 10.1097/EDE.0000000000000982
################################################################################

################################################################################
# LOAD PACKAGES
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel);library(pbs)
library(mgcv) ; library(ggplot2) ; library(lubridate); library(mondate)


# A CRUDE FUNCTION TO PLOT A PREVIOUSLY FITTED CYCLIC SPLINE
# Detailed information is available in the tutorial paper below
#   Madaniyazi, Lina, et al. International Journal of Epidemiology (2022): dyac115.
#   DOI: https://doi.org/10.1093/ije/dyac115.
plotcyclic <- function(model,data,df=4,col="black") {
  # COEFFICIENTS FROM CYCLIC SPLINE
  coef1<- coef(model)[ 2:(df+1) ]
  # log RRS FOR EACH DAY OF YEAR (LOGRR, SE)
  logrr1 <- matrix(NA, 1, 366 )
  Trange <- 1:366
  bvar<- pbs(data$doy , df = df  ) # doy is the exposure indicator for seasonality
  l <- length(Trange)
  for(j in 1:l){ logrr1[1,j] <-  t(bvar[j,]) %*% as.numeric( coef1)  }
  # CENTERING LOGRR AT TROUGH
  maxcen<- apply(logrr1, 1, which.min )
  bvar_cen<- bvar[maxcen,]
  # CENTERED LOGRR
  l <- length(Trange)
  for(j in 1:l){ logrr1[1,j] <- as.numeric( t(bvar[j,]-bvar_cen) %*% as.numeric( coef1) ) }
  # plot
  plot(Trange, exp(logrr1), type='l', xlab="Day of year", ylab="RR (95% CI)",  ylim=c(0.9,1.4),
       lwd=2, xaxt='n',col=col )
  axis(1, at=seq(0,350,50))
  abline(h=1, col="black")
}


# FIND FEB 29 IN DATE
find_leap = function(x){
  day(x) == 29 & month(x) == 2
}
################################################################################


################################################################################
# LOAD DATA FOR OBSERVED SERIES
data<-read.csv("obs.csv",colClasses=c(date="Date"))
# DEFINE THE OUTCOME & DOW
data$death <- data$all
data <- data[!is.na(data$death),]
data$dow<-wday(data$date, week_start=1)

# LOAD PROJECTED TEMPERATURE SERIES UNDER RCP8.5
# The data has been calibrated by following the tutorial paper on projection below
# Epidemiology. 2019 May; 30(3): 321–329. doi: 10.1097/EDE.0000000000000982
tm<-read.csv("rcp8p5.csv",colClasses=c(date="Date"))[-1]
# WE USE MEAN TEMPERATURE ACROSS FIVE GCMS IN EXAMPLE FOR PREDICTION
tm$tmean<-rowMeans(tm[,c(2:6)])
tm<- tm[!find_leap(tm$date),]
################################################################################


################################################################################
# DEFINE THE MAIN PARAMETERS FOR THE ANALYSIS
################################################################################
# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun = "ns"
varper <- c(10,75,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 2
################################################################################


################################################################################
# FIRST-STAGE ANALYSIS: RUN THE MODEL IN EACH CITY, OBTAIN PROJECTED MORTALITY
################################################################################
# DEFINE TEMPERATURE CROSS-BASIS
argvar=list(fun=varfun,knots=quantile(data$tmean,varper/100,na.rm=T),Boundary.knots=range(data$tmean,na.rm=T))
data$cb<-crossbasis(data$tmean,lag=lag,argvar=argvar,
                    arglag=list(knots=logknots(lag,lagnk)))

# DEFINE SPLINE FOR LONG TERM TREND  (2DF/DECADE NOT PER YEAR)
data$trendspline <- ns(data$date,df=round(dfseas*length(unique(data$year))/10))

# RUN THE MODEL 
model <- glm(death~pbs(doy,df=4)+trendspline
             + cb,  data, family=quasipoisson,na.action="na.exclude")

##### CREATE NEWDATA FOR PREDICTION #####
# SET UP A DECADE OF 365-DAY YEARS 
datanew<-data[which(data$year=="2000"),]   # 365 ROWS
datanew<-datanew[rep(c(1:365),10),]        # 3560 ROWS TO MATCH WITH TEMPERATURE SERIES FOR EACH DECADE

# FIX LONG-TERM TREND TERM TO VALUE PREDICTED AT 2000-01-01  (FIX AT A SINGLE DATE, NOT A YEAR) 
trendspline_base<- data[which(data$date=="2000-01-01"),]$trendspline
datanew$trendspline<- matrix(trendspline_base, nrow=3650, ncol=length(trendspline_base), byrow=TRUE)    # THIS WILL NOW BE USED IN PREDICTIONS FROM THE ORIGINAL FITTED COEFFICIENTS

# ADD MODELLED TEMPERATURE SERIES TO NEWDATA (WE USE THE LAST PERIOD IN EXAMPLE) 
datanew$tmean<-tm[tm$date>"2089-12-31"&tm$date<"2100-01-01",]$tmean # SELECT TEMPERATURE SERIES BETWEEN 2090 AND 2099 AND ADD IT TO NEW DATA 
datanew$cb<-crossbasis(datanew$tmean,lag=lag,argvar=argvar, arglag=list(knots=logknots(lag,lagnk))) # EXACTLY THE SAME TEMP KNOTS AS BASELINE FIT
datanew_p10<-datanew

# PREDICT MORTALITY FOR EACH DAY-OF
datanew_p10$bm_p10<-predict(model,newdata = datanew_p10, type="response")
bm_rcp8p5<-datanew_p10[,c("doy","date","bm_p10")]

# PLOT PROJECTED DAILY MORTALITY (THE LAST PERIOD)
model2 <- glm(bm_p10~pbs(doy,df=4), data=bm_rcp8p5,family=quasipoisson,na.action="na.exclude")

plotcyclic(model2, data=bm_rcp8p5, col="red")
title("Projected seasonality under RCP8.5")
################################################################################

################################################################################
# SECOND-STAGE ANALYSIS: SEASONALITY ASSESSMENT
################################################################################
# WE FOLLOW THE CODES IN OUR PREVIOUS TUTORIAL PAPER TO ASSESS THE SEASONALITY 
# I.E. TO OBTAIN TIMINGS AND SIZE FROM MODEL 2 
# Madaniyazi L, Tobias A, Kim Y, Chung Y, Armstrong B, Hashizume M. 
# Assessing seasonality and the role of its potential drivers in environmental 
# epidemiology: a tutorial. Int J Epidemiol 2022; 51: 1677–86.

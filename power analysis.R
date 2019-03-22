library(lme4)
library(lmerTest)
library(plot3D)


# broom the model outputs & summary
# nSamps: number of samples per site
# nSites: number of sites (or groups)
# sd.obs: standard deviation in ecological metric (e.g., size)
# sd.group: standard deviation in how the mean of the metric varies among groups

sd.obs <- 4
sd.group <- 1

costSite <- 4 # hours of work/driving to add extra site
costSamp <- 0.5 # hours of work to add extra sample

statsig <- function(nSamps,nSites,effect,sd.obs,sd.group){
  
  # imagine you have only 1 measurement of x per stream/site
  x <- rnorm(nSites,mean=0,sd=1) # generate the predictor variable at the site
  x1 <- rep(x,each=nSamps) # repeat the x variable for each sample at the site
  group <- as.factor(rep(1:nSites,each=nSamps)) # track the site ID for each samples
  y <- .2 + effect*x + rnorm(length(x),mean=0,sd=sd.group) # create random ecological metric
  y1 <- rnorm(length(x1),mean=rep(y,each=nSamps),sd=sd.obs) # create random data for each group
  fit1 <- lm(y1~x1)
  fit2 <- suppressMessages(lmer(y1~x1+(1|group),REML=F))
  s1 <- summary(fit1)
  s2 <- summary(fit2) 
  
  pw1 <- effect <= (s1[[4]][2,1] + 1.96*s1[[4]][2,2]) & effect >= (s1[[4]][2,1]-1.96*s1[[4]][2,2]) & 0 <= (s1[[4]][2,1]+1.96*s1[[4]][2,2])*(s1[[4]][2,1]-1.96*s1[[4]][2,2])
  
  pw2 <- effect<=(s2$coefficients[2,1]+1.96*s2$coefficients[2,2])&
    effect>=(s2$coefficients[2,1]-1.96*s2$coefficients[2,2])&
    0<=(s2$coefficients[2,1]+1.96*s2$coefficients[2,2])*
    (s2$coefficients[2,1]-1.96*s2$coefficients[2,2])
  
  round(c(slope1=s1$coefficients[2,1],R1=s1$adj.r.squared,
          Type_I_error_1=s1[[4]][2,4]<=0.05,pw1=pw1,slope2=s2$coefficients[2,1],
          Type_I_error_2=s2$coefficients[2,5]<=0.05,pw2=pw2),3)
  }


effect_vec <- c(0,1)
nSites <- seq(from=5,to=50,by=5)
nSamps <- c(3,5,10,15,20,25,30,35,50,100)
results <- array(NA,dim=c(length(effect_vec),length(nSites),length(nSamps),7),dimnames=list("Effect size"=effect_vec,"N Sites"=nSites,"N Samps per Site"=nSamps,"Metric"=c("slope (lm)","R2 (lm)","Type I error (lm)","Power (lm)","Slope (lme4)","Type I error (lme4)","Power (lme4)")))

counter <- 1
progBar <- txtProgressBar(min = 0, max = length(effect_vec)*length(nSites)*length(nSamps),title="More power",style=3,initial=0)
ptm = Sys.time()
for(i in 1:length(effect_vec)){
  for(j in 1:length(nSites)){
    for(k in 1:length(nSamps)){
      out <- replicate(1000,statsig(nSamps=nSamps[k],nSites=nSites[j],effect=effect_vec[i],sd.obs=sd.obs,sd.group=sd.group),simplify=T)
      results[i,j,k,] <- rowMeans(out)
      setTxtProgressBar(progBar,counter)
      counter <- counter + 1
    }
  }
}
endtime <- Sys.time()-ptm
endtime

# analyze your cost-benefit ratios
power <- results[2,,,"Power (lme4)"]
power[power<0.8] <- NA
noType1 <- 1-results[1,,,"Type I error (lme4)"]
noType1[noType1<0.95] <- NA

cost <- sapply(costSamp*nSamps,function(x){x+costSite*nSites})
dimnames(cost) <- dimnames(power)

# find where you have minimal costs that maximizes your power and chance to avoid type 1 error
costBenefit <- which(abs(cost/power-min(cost/power,na.rm=T)+(cost/type1-min(cost/type1,na.rm=T)))==min(abs(cost/power-min(cost/power,na.rm=T)+(cost/type1-min(cost/type1,na.rm=T)))),arr.ind=T)
paste("Sites = ",nSites[costBenefit[1]]," & ",
        "Samples = ",nSamps[costBenefit[2]],sep="")

# what if there was no cost to more sites or more samples per site
noCost <- which(abs(power-max(power)+(type1-max(type1)))==min(abs(power-max(power)+(type1-max(type1)))),arr.ind=T)
paste("Sites = ",nSites[noCost[1]]," & ",
      "Samples = ",nSamps[noCost[2]],sep="")

# what if we didn't account for pseudoreplication?

# analyze your cost-benefit ratios
power <- results[2,,,"Power (lm)"]
noType1 <- 1-results[1,,,"Type I error (lm)"]

cost <- sapply(costSamp*nSamps,function(x){x+costSite*nSites})
dimnames(cost) <- dimnames(power)

# find where you have minimal costs that maximizes your power and chance to avoid type 1 error
costBenefit <- which(abs(cost/power-min(cost/power)+(cost/type1-min(cost/type1)))==min(abs(cost/power-min(cost/power)+(cost/type1-min(cost/type1)))),arr.ind=T)
paste("Sites = ",nSites[costBenefit[1]]," & ",
      "Samples = ",nSamps[costBenefit[2]],sep="")

# what if there was no cost to more sites or more samples per site
noCost <- which(abs(power-max(power)+(type1-max(type1)))==min(abs(power-max(power)+(type1-max(type1)))),arr.ind=T)
paste("Sites = ",nSites[noCost[1]]," & ",
      "Samples = ",nSamps[noCost[2]],sep="")

# make contour plot to highlight tradeoff

jpeg("cost-benefit.jpeg",units="in",res=800,height=7,width=7)

layout(1)
par(mar=c(5,4,1,1))
filled.contour(x=costSamp*nSamps, y=costSite*nSites, z=results[2,,,"Power (lme4)"],xlab='Time spent sampling each site',ylab='Time spent getting to sites',color.palette = colorRampPalette(c("red",'orange','dodgerblue',"yellow")))

dev.off()

jpeg("sample size tradeoff.jpeg",units="in",res=800,height=7,width=7)

layout(1)
par(mar=c(5,4,1,1))
filled.contour(x=nSamps, y=nSites, z=results[2,,,"Power (lme4)"],xlab='Number of samples',ylab='Number of sites',color.palette = colorRampPalette(c("red",'orange','dodgerblue',"yellow")))

dev.off()

jpeg("sample size cost-benefit.jpeg",units="in",res=800,height=7,width=7)

layout(1)
par(mar=c(5,4,1,1))
filled.contour(x=nSamps, y=nSites, z=cost/results[2,,,"Power (lme4)"],xlab='Number of samples',ylab='Number of sites',color.palette = colorRampPalette(c("red",'orange','dodgerblue',"yellow")))

dev.off()

jpeg("type 1 error.jpeg",units="in",res=800,height=7,width=7)

layout(1)
par(mar=c(5,4,1,1))
filled.contour(x=nSamps, y=nSites, z=noType1,xlab='Number of samples',ylab='Number of sites',color.palette = colorRampPalette(c("red",'orange','dodgerblue',"yellow")))
dev.off()

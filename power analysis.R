library(lme4)
library(lmerTest)
library(plot3D)
# ns: number of samples per site
# ng: number of sites (or groups)
# tsig: standard deviation in ecological metric (e.g., size)
# gsig: standard deviation in how the mean of the metric varies among groups

sd.obs <- 5
sd.group <- 1

statsig <- function(ns,ng,effect,tsig,gsig){
  x <- rnorm(ng,0,1) # generate the predictor variable
  x1 <- rep(x,each=ns) # generate the NS number of x variables
  group <- as.factor(rep(1:ng,each=ns)) # generate NG number of groups associated for each of NS samples
  y <- .2 + effect*x + rnorm(length(x),0,gsig) # create random ecological metric
  y1 <- rnorm(length(x1),rep(y,each=ns),tsig) # create random data for each group
  fit1 <- lm(y1~x1)
  fit2 <- suppressMessages(lmer(y1~x1+(1|group),REML=F))
  s1 <- summary(fit1)
  s2 <- summary(fit2) 
  
  pw1 <- effect<=(s1[[4]][2,1]+1.96*s1[[4]][2,2])&
    effect>=(s1[[4]][2,1]-1.96*s1[[4]][2,2])&
    0<=(s1[[4]][2,1]+1.96*s1[[4]][2,2])*(s1[[4]][2,1]-1.96*s1[[4]][2,2])
  
  pw2 <- effect<=(s2$coefficients[2,1]+1.96*s2$coefficients[2,2])&
    effect>=(s2$coefficients[2,1]-1.96*s2$coefficients[2,2])&
    0<=(s2$coefficients[2,1]+1.96*s2$coefficients[2,2])*
    (s2$coefficients[2,1]-1.96*s2$coefficients[2,2])
  
  round(c(slope1=s1$coefficients[2,1],R1=s1$adj.r.squared,
          Type_I_error_1=s1[[4]][2,4]<=0.05,pw1=pw1,slope2=s2$coefficients[2,1],
          Type_I_error_2=s2$coefficients[2,5]<=0.05,pw2=pw2),3)
  }


effect_vec <- seq(0,2,length=10)
nSites <- seq(from=5,to=50,by=5)
nSamps <- seq(from=5,to=50,by=5)
results <- array(NA,dim=c(length(effect_vec),length(nSites),length(nSamps),7),dimnames=list("Effect size"=effect_vec,"N Sites"=nSites,"N Samps per Site"=nSamps,"Metric"=c("slope (lm)","R2 (lm)","Type I error (lm)","Power (lm)","Slope (lme4)","Type I error (lme4)","Power (lme4)")))

counter <- 1
progBar <- txtProgressBar(min = 0, max = length(effect_vec)*length(nSites)*length(nSamps),title="More power",style=3,initial=0)
ptm = proc.time()
for(i in 1:length(effect_vec)){
  for(j in 1:length(nSites)){
    for(k in 1:length(nSamps)){
      out <- replicate(9,statsig(ns=nSamps[k],ng=nSites[j],effect=effect_vec[i],tsig=sd.obs,gsig=sd.group),simplify=T)
      results[i,j,k,] <- rowMeans(out)
      setTxtProgressBar(progBar,counter)
      counter <- counter + 1
    }
  }
}

results[1,,,"Type I error (lme4)"]

endtime <- proc.time()-ptm
endtime[3]/60
# let's plot an z ~ f(x,y), but have the perspective for the 3-D plot move.
plotdegrees <- seq(from=1,to=360,by=5)
for(i in plotdegrees){
  scatter3D(x=rep(nSites,length(nSamps),each=length(effect_vec)),y=rep(nSamps,each=length(nSites)*length(effect_vec)),z=rep(effect_vec,length(nSites)*length(nSamps)),colvar = as.vector(results[,,,"Power (lme4)"]),phi=25,theta=i,type="h",clab="Power",xlab="Number of sites",ylab="Samples per site",zlab="Effect size",bty="g",col=ramp.col(c("dodgerblue", "orange")),pch=19,lwd=1.5,ticktype= "detailed")
  Sys.sleep(0.25)
  }

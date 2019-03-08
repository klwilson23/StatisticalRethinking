
# Chapter 7
# Interactions

library(rethinking)
rstan_options(auto_write=TRUE)
options(mc.cores=2)

data(rugged)
d <- rugged

dd <- d[complete.cases(d$rgdppc_2000),]

dd$log_gdp <- log(dd$rgdppc_2000)

pairs(data = dd, ~ log_gdp + rugged + cont_africa + log(slave_exports) + dist_coast, lower.panel=panel.smooth)

X <- model.matrix(log_gdp ~ cont_africa + rugged, data=dd)

hist(dd$log_gdp,breaks=seq(-30,30,by=0.5),freq=FALSE,ylim=c(0,1),xlim=c(-20,20))
hist(dd$log_gdp[dd$cont_africa==1],breaks=seq(-30,30,by=0.5),freq=FALSE,add=TRUE,col="lightblue")
curve(dnorm(x,mean=mean(dd$log_gdp),sd=5),add=TRUE)

# cool code:
lookup(length,ReturnType = character())

stan_model <- "
data{
  int N;
  int K; // number of predictors in design matrix
  vector[N] log_gdp;
  matrix[N,K] design_mat;
}
parameters{
  vector[K] beta;
  real sigma;
}
transformed parameters{
  vector[N] mu;
  mu = design_mat*beta;
}
model{
  log_gdp ~ normal(mu, sigma);
  beta[1] ~ normal(mean(log_gdp),5);
  for(i in 2:K) {
    beta[i] ~ normal(0, 5);
  }
  sigma ~ cauchy(0, 10);
}
generated quantities{
  vector[N] log_lik;
  vector[N] log_gdp_new; // posterior predictive distributions
  for(i in 1:N){
    log_lik[i] = normal_lpdf(log_gdp[i] | mu[i], sigma);
    log_gdp_new[i] = normal_rng(mu[i],sigma);
  }
}
"
dat1 <- list(N = nrow(dd), K=ncol(X),log_gdp=dd$log_gdp,design_mat=X)
fit_dummy <- stan(model_code=stan_model,data=dat1,chains=2,iter=2000,cores=2,pars=c("beta","sigma","log_lik","log_gdp_new"))

summary(fit_dummy,pars=c("beta","sigma"),probs=c(0.1,0.9))$summary

# plots
post1 <- as.data.frame(fit_dummy)
plot(log_gdp~rugged, data=dd, main="Ln(GDP)~ Ruggedness, Africa dummy variable",pch=21,bg=ifelse(dd$cont_africa==1,"orange","grey"))
invisible(sapply(1:250,function(i){
  curve(post1$`beta[1]`[i] + post1$`beta[3]`[i]*x,add=TRUE,lwd=2,col=adjustcolor("grey",alpha=0.5))
}))
invisible(sapply(1:250,function(i){
  curve(post1$`beta[1]`[i] + post1$`beta[2]`[i]*1 + post1$`beta[3]`[i]*x,add=TRUE,lwd=2,col=adjustcolor("lightblue",alpha=0.5))
}))

curve(mean(post1$`beta[1]`) + mean(post1$`beta[2]`)*0 + mean(post1$`beta[3]`)*x,add=TRUE,lwd=2,col=adjustcolor("black",alpha=0.5))

curve(mean(post1$`beta[1]`) + mean(post1$`beta[2]`)*1 + mean(post1$`beta[3]`)*x,add=TRUE,lwd=2,col=adjustcolor("black",alpha=0.5))

# Example 2

# GDP v. ruggedness with an interaction effect of continent
X2 <- model.matrix(log_gdp ~ cont_africa*rugged,data=dd)
dat2 <- list(N=nrow(dd), K = ncol(X2), log_gdp=dd$log_gdp,design_mat=X2)
fit_int <- stan(model_code=stan_model,data=dat2,chains=2,iter=2000,cores=2,pars=c("beta","sigma","log_lik","log_gdp_new"))

summary(fit_int,pars=c("beta","sigma"),probs=c(0.1,0.9))$summary

rethinking::compare(fit_dummy,fit_int)

# plots
post2 <- as.data.frame(fit_int)
par(mfrow=c(1,2))
plot(log_gdp~rugged, data=dd[dd$cont_africa==0,], main="Non-Africa",pch=21)
invisible(sapply(1:250,function(i){
  curve(post2$`beta[1]`[i] + post2$`beta[2]`[i]*0 + post2$`beta[3]`[i]*x + post2$`beta[4]`[i]*x*0,add=TRUE,lwd=2,col=adjustcolor("grey",alpha=0.5))
}))
curve(mean(post2$`beta[1]`) + mean(post2$`beta[2]`)*0 + mean(post2$`beta[3]`)*x + mean(post2$`beta[4]`)*x*0,add=TRUE,lwd=2,col=adjustcolor("black",alpha=1))

plot(log_gdp~rugged, data=dd[dd$cont_africa==1,], main="Africa",pch=21)
invisible(sapply(1:250,function(i){
  curve(post2$`beta[1]`[i] + post2$`beta[2]`[i]*1 + post2$`beta[3]`[i]*x + post2$`beta[4]`[i]*x*1,add=TRUE,lwd=2,col=adjustcolor("lightblue",alpha=0.5))
}))

curve(mean(post2$`beta[1]`) + mean(post2$`beta[2]`)*1 + mean(post2$`beta[3]`)*x + mean(post2$`beta[4]`)*x*1,add=TRUE,lwd=2,col=adjustcolor("black",alpha=1))

# Visualize predicted distributions
ppd <- post2[ , grep("log_gdp_new",colnames(post2))]
mn_ppd <- apply(ppd, 2 ,mean)
ci_ppd <- apply(ppd, 2, HPDI, prob=0.89)

par(mfrow=c(1,1))
plot(dd$log_gdp, mn_ppd, pch=21, bg=ifelse(dat2$design_mat[,2]==1,"orange","lightblue"),ylim=range(ci_ppd),xlab="Observed ln(GDP)",ylab="Posterior predictive ln(GDP)")
segments(x0=dd$log_gdp, y0=ci_ppd[1,], y1=ci_ppd[2,],lwd=2)
abline(b=1,a=0,lwd=2,col="red")

boxplot(ppd, col=ifelse(dat2$design_mat[,2]==1,"orange","lightblue"),outline=FALSE)
log_gdp_mn <- aggregate(log_gdp~cont_africa,data=dd,FUN=mean)
abline(h=log_gdp_mn$log_gdp,col=c("lightblue","orange"),lwd=3)

# Example 3
# GDP ~ dist_coast + slave exports

dd$log_slave <- log(1+dd$slave_exports)

X3 <- model.matrix(log_gdp ~ log_slave*dist_coast,data=dd)
dat3 <- list(N=nrow(X3),K=ncol(X3),log_gdp=dd$log_gdp,design_mat=X3)

fit_cont <- stan(model_code = stan_model, data=dat3, chains=2, iter=2000,cores=2,pars=c("beta","sigma","log_lik","log_gdp_new"))

summary(fit_cont,pars=c("beta","sigma"),probs=c(0.1,0.9))$summary

rethinking::compare(fit_dummy,fit_int,fit_cont)

post3 <- as.data.frame(fit_cont)

# Visualize predicted distributions
ppd <- post3[ , grep("log_gdp_new",colnames(post3))]
mn_ppd <- apply(ppd, 2 ,mean)
ci_ppd <- apply(ppd, 2, HPDI, prob=0.89)

par(mfrow=c(1,1))
plot(dd$log_gdp, mn_ppd, pch=21, bg=ifelse(dat2$design_mat[,2]==1,"orange","lightblue"),ylim=range(ci_ppd),xlab="Observed ln(GDP)",ylab="Posterior predictive ln(GDP)")
segments(x0=dd$log_gdp, y0=ci_ppd[1,], y1=ci_ppd[2,],lwd=2)
abline(b=1,a=0,lwd=2,col="red")

boxplot(ppd, col=ifelse(dat2$design_mat[,2]==1,"orange","lightblue"),outline=FALSE)
log_gdp_mn <- aggregate(log_gdp~cont_africa,data=dd,FUN=mean)
abline(h=log_gdp_mn$log_gdp,col=c("lightblue","orange"),lwd=3)



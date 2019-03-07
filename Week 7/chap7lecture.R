#Chapter 7 
#Interactions

library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

data(rugged)
d <- rugged
#Eliminate NA values 
dd <- d[complete.cases(d$rgdppc_2000),]

pairs(data = dd , ~ rgdppc_2000 + rugged + cont_africa + slave_exports + dist_coast,
      lower.panel = panel.smooth)
pairs(data = dd, ~ log(rgdppc_2000) + rugged + cont_africa + log(slave_exports) + dist_coast,
      lower.panel = panel.smooth)

aggregate(slave_exports ~ cont_africa, data = d, FUN = mean)

#Use log(gdp)
dd$log_gdp <- log(dd$rgdppc_2000)

#Example 1 - Dummy Variable
#log_gdp vs ruggedness with cont_africa as dummy variable 

#Visualizing your prior
# dens(dd$log_gdp)
# curve(dnorm(x, mean = mean(dd$log_gdp), sd = 5), add = TRUE)

#Prior for beta[1] can just cover the first distribution
hist(dd$log_gdp, breaks = seq(0,30, by = 0.5), freq = FALSE, ylim = c(0, 1), xlim = c(-30,30))
curve(dnorm(x, mean = mean(dd$log_gdp), sd = 5), add = TRUE)
# #Plot hist for africa only
# #Prior needs to cover both/the difference the ditribution of africa only and overall dist.
# hist(dd$log_gdp[dd$cont_africa == 1], breaks = seq(-30,30, by = 0.5), freq = FALSE, add = TRUE, col = "light blue")
#Average difference between the two peaks
#Prior covers that value 
curve(dnorm(x, mean = 0, sd = 2), col = "orange", from = -5, to = 5, xlim = c(-10,10))
abline(v = mean(dd$log_gdp)-mean(dd$log_gdp[dd$cont_africa==1]))

X <- model.matrix(log_gdp ~ cont_africa + rugged, data = dd)

stan_model <- "
data{
  int N;
  int K; //number of predictors
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
  beta[1] ~ normal(mean(log_gdp), 5);
  for (i in 2:K){
    beta[i] ~ normal(0, 2);
  }
  sigma ~ cauchy(0,10);
}
generated quantities{
  vector[N] log_lik;
  vector[N] log_gdp_new;
  for (i in 1:N){
  log_lik[i] = normal_lpdf(log_gdp[i] | mu[i], sigma);
  log_gdp_new[i] = normal_rng(mu[i], sigma);
  }
}
"

dat1 <- list(N = nrow(dd), K = ncol(X), log_gdp = dd$log_gdp,
             design_mat = X)
?stan
fit_dummy <- stan(model_code = stan_model, data = dat1,
                  chains = 2, iter = 2000, cores = 1,
                  pars = c("beta", "sigma", "log_lik","log_gdp_new"))

# the better way to view our posterior samples
summary(fit_dummy,pars=c("beta","sigma"),probs=c(0.1,0.9))$summary


precis(fit_dummy, depth = 2)
#Beta 1 and 2 are the intercepts and beta 3 is the slope 
#Only the intercepts are allowed to change in this model
#Plot to visualize this
post1 <- as.data.frame(fit_dummy)
plot(log_gdp ~ rugged, data = dd, main = "Log-GDP ~ Ruggedness, Africa dummy variable")
invisible(sapply(1:250, function(i){
  curve(post1$`beta[1]`[i] + post1$`beta[3]`[i]*x, add = TRUE, lwd = 2, col = adjustcolor("grey", alpha = 0.1))
  curve(post1$`beta[1]`[i] + post1$`beta[2]`[i] + post1$`beta[3]`[i]*x, add = TRUE, lwd = 2, col = adjustcolor("light blue", alpha = 0.1))
}))
curve(mean(post1$`beta[1]`) + mean(post1$`beta[3]`)*x, add = TRUE, lwd = 2)
curve(mean(post1$`beta[1]`) + mean(post1$`beta[2]`) + mean(post1$`beta[3]`)*x, add = TRUE, lwd = 2, col = "dodger blue")

#Example 2
#GDP vs. ruggedness, with an interaction effect of being within/outside of Africa
X2 <- model.matrix(log_gdp ~ cont_africa*rugged, data = dd)

dat2 <- list(N = nrow(dd), K = ncol(X2), log_gdp = dd$log_gdp,
             design_mat = X2)

fit_int <- stan(model_code = stan_model, data = dat2,
                 chains = 2, iter = 2000, cores = 1,
                 pars = c("beta", "sigma", "log_lik","log_gdp_new"))

precis(fit_int, depth = 2)

post2 <- as.data.frame(fit_int)

#Plot to visualize the interactions 
par(mfrow = c(1,2))
#Non-African nations - Base case
plot(log_gdp ~ rugged, data = dd[dd$cont_africa == 0,], main = "Non-African Nations")
invisible(sapply(1:250, function(i){
  curve(post2$`beta[1]`[i] + post2$`beta[3]`[i]*x, add = TRUE, lwd = 2, col = adjustcolor("grey", alpha = 0.1))
}))
curve(mean(post2$`beta[1]`) + mean(post2$`beta[3]`)*x, add = TRUE, lwd = 2)
#African nations
plot(log_gdp ~ rugged, data = dd[dd$cont_africa == 1,], main = "African Nations", col = "light blue")
invisible(sapply(1:250, function(i){
  curve(post2$`beta[1]`[i] + post2$`beta[2]`[i] + post2$`beta[4]`[i]*x + post2$`beta[3]`[i]*x, add = TRUE, lwd = 2, col = adjustcolor("light blue", alpha = 0.1))
}))
curve(mean(post2$`beta[1]`) + mean(post2$`beta[2]`) + mean(post2$`beta[3]`)*x + mean(post2$`beta[4]`)*x, add = TRUE, lwd = 2, col = "dodger blue")


#Another way of visualizing
#Generate distributions of expected value for our observed values
# delete below
log_gdp_hat2 <- matrix(NA, nrow = nrow(post2), ncol= nrow(X2))
for (i in 1:nrow(post2)){
  estimates <- sapply(1:nrow(X2), FUN = function(x){as.numeric(post2[i, grep("beta", colnames(post2))]*X2[x,])})
  value <- colSums(estimates)
  log_gdp_hat2[i,] <- rnorm(length(value), mean = value, sd = post2[i, "sigma"])
}

mn_ppd <- apply(log_gdp_hat2,2,mean)
ci_ppd <- apply(log_gdp_hat2,2,HPDI,prob=0.89)
par(mfrow = c(1,1))
plot(dat2$log_gdp,mn_ppd,pch=21,bg="lightblue",ylim=range(ci_ppd),xlim=range(ci_ppd), main = "Africa*rugged")
segments(x0=dd$log_gdp,y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=2,col="orange")
# delete above

ppd <- post2[,grep("log_gdp_new",colnames(post2))]
mn_ppd <- apply(ppd,2,mean)
ci_ppd <- apply(ppd,2,HPDI,prob=0.89)
par(mfrow = c(1,1))
plot(dat2$log_gdp,mn_ppd,pch=21,bg=ifelse(dat2$design_mat[,2]==1,"orange","lightblue"),ylim=range(ci_ppd),xlim=range(ci_ppd), main = "Africa*rugged",xlab="Observed ln(GDP)",ylab="Posterior predictive ln(GDP)")
segments(x0=dd$log_gdp,y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=2,col="orange")

boxplot(ppd,col=ifelse(dat2$design_mat[,2]==1,"orange","lightblue"),outline=FALSE)
log_gdp_mn <- aggregate(log_gdp ~ cont_africa, data = dd, FUN = mean)
abline(h=log_gdp_mn$log_gdp,col=c("lightblue","orange"),lwd=3)

#Example 2
#GDP vs ruggedness and slave exports within Africa 
#Need to log slave exports (+1) because range is so wide 
dd$slave_logged <- log(dd$slave_exports + 1)
X3 <- model.matrix(log_gdp ~ slave_logged*dist_coast, data = dd)

dat3 <- list(N = nrow(dd), K = ncol(X3), log_gdp = dd$log_gdp,
              design_mat = X3)

fit_cont <- stan(model_code = stan_model, data = dat3,
                 chains = 2, iter = 2000, cores = 1,
                 pars = c("beta", "sigma", "log_lik","log_gdp_new"))
summary(fit_cont,pars=c("beta","sigma"),probs=c(0.1,0.9))$summary

precis(fit_cont, depth = 2)

post3 <- as.data.frame(fit_cont)
#Make triptych plots like in the textbook?

#Plot predicted vs observed
#Make an empty matrix to populate our expected values with
log_gdp_hat3 <- matrix(NA, nrow = nrow(post3), ncol= nrow(X3))

for (i in 1:nrow(post3)){
  estimates <- sapply(1:nrow(X3), FUN = function(x){as.numeric(post3[i, grep("beta", colnames(post3))]*X3[x,])})
  value <- colSums(estimates)
  log_gdp_hat3[i,] <- rnorm(length(value), mean = value, sd = post3[i, "sigma"])
}

mn_ppd <- apply(log_gdp_hat3,2,mean)
ci_ppd <- apply(log_gdp_hat3,2,HPDI,prob=0.89)
plot(dd$log_gdp,mn_ppd,pch=21,bg="lightblue",ylim=range(ci_ppd),xlim=range(ci_ppd), main = "Slave exports*Dist to coast")
segments(x0=dd$log_gdp,y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=2,col="orange")

ppd <- post3[,grep("log_gdp_new",colnames(post2))]
mn_ppd <- apply(ppd,2,mean)
ci_ppd <- apply(ppd,2,HPDI,prob=0.89)
par(mfrow = c(1,1))
plot(dat3$log_gdp,mn_ppd,pch=21,bg=ifelse(dat2$design_mat[,2]==1,"orange","lightblue"),ylim=range(ci_ppd),xlim=range(ci_ppd), main = "Africa*rugged",xlab="Observed ln(GDP)",ylab="Posterior predictive ln(GDP)")
segments(x0=dd$log_gdp,y0=ci_ppd[1,],y1=ci_ppd[2,],lwd=1)
abline(b=1,a=0,lwd=2,lty=2,col="orange")


#Compare our models using WAIC
rethinking::compare(fit_dummy, fit_int, fit_cont)


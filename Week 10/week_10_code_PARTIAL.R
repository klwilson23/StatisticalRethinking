# Statistical Rethinking Directed Readings
# Chapter 12 - Multi-level Models
# 2019-03-28

rm(list=ls())

library(rethinking)

rstan_options(auto_write = TRUE)
options(mc.cores=2)

# Why multi-level? We want to model the average of a group
# Multi-level model does two things;
# 1. What is the average across groups?
# 2. What is the group-level effect?

data("chimpanzees")
d <- chimpanzees
head(d)

#### Start with one intercept model and intercept for each chimp
# 1. Fit model with just one intercept ------------
model1 <- "
data {
  int N;
  int<lower=0, upper=1> pulled_left[N]; // pulled left
  int prosoc_left[N]; // prosocial / two plates of food option is on left side
  int condition[N]; // 1 means partner is present, 0 for control
}
parameters{
  real a; // intercept
  real bp; // beta for pro social
  real bpc; // beta for pro social left and condition interaction
}
transformed parameters{
  vector[N] p;
  for(i in 1:N){
    p[i] = a + bp * prosoc_left[i] + bpc * prosoc_left[i] * condition[i]; 
  }
}
model{
  a ~ normal(0,10);
  bp ~ normal(0,10);
  bpc ~ normal(0, 10);
  pulled_left ~ binomial_logit(1,p);
}
generated quantities{
  vector[N] log_lik;
  for(i in 1:N) {
  log_lik[i] = bernoulli_logit_lpmf(pulled_left[i] | p[i]);
  }
}
"
dat1 <- list(N=nrow(d), 
             pulled_left = d$pulled_left, 
             prosoc_left = d$prosoc_left,
             condition = d$condition)
fit1 <-stan(model_code = model1, data=dat1, iter=3000,chains=2,cores=2,
            pars=c("a", "bp", "bpc", "log_lik"))
summary(fit1,pars=c("a","bp","bpc"), probs=c(0.1, 0.9))$summary
#traceplot(fit1, pars=c("a", "bp", "bpc"))

# 2. Binomial model: individuals have specific intercepts ------------

model2 <- "
data {
  int N;
  int<lower=0, upper=1> pulled_left[N]; // pulled left
  int prosoc_left[N]; // prosocial / two plates of food option is on left side
  int condition[N]; // 1 means partner is present, 0 for control
  int<lower=1,upper=N> N_chimps;
  int<lower=1, upper=N_chimps> chimp[N]; // id of chimp / actor
}
parameters{
  real a[N_chimps];
  real bp;
  real bpc;
}
transformed parameters {
  vector[N] p;
  for(i in 1:N){
    p[i] = a[chimp[i]] + bp * prosoc_left[i] + bpc * condition[i] * prosoc_left[i];
  }
}
model{
  for(i in 1:N_chimps) {
    a[i] ~ normal(0, 10);
  }
  bp ~ normal(0,10);
  bpc ~ normal(0,10);
  pulled_left ~ binomial_logit(1,p);
}
generated quantities{
  vector[N] log_lik;
  for(i in 1:N){
    log_lik[i] = binomial_logit_lpmf(pulled_left[i] | 1, p[i]);
  }
}
"

dat2 <- list(N=nrow(d), 
             pulled_left = d$pulled_left, 
             prosoc_left = d$prosoc_left, 
             condition = d$condition, 
             N_chimps = length(unique(d$actor)), 
             chimp = d$actor)
fit2 <- stan(model_code = model2, data=dat2, iter=3000,chains=2,cores=2,control=list("adapt_delta"=0.81),
             pars=c("a", "bp", "bpc", "log_lik"))
summary(fit2, pars=c("a","bp","bpc"),probs=c(0.1, 0.9))$summary
#traceplot(fit2)

### START TYPING WORKSHOP HERE

# make a multi-level model here

model3 <- "
data{
  int N;
  int<lower=0,upper=1> pulled_left[N];
  int prosoc_left[N];
  int condition[N];
  int N_chimps;
  int<lower=1, upper=N_chimps> chimp[N];
}
parameters{
  real a_chimp[N_chimps];
  real bp; // effect of prosoc_left
  real bpc; // interaction effect
  real mu_chimp; // mean intercept across chimps
  real sigma_chimp; // standard deviation across chimps
}
transformed parameters {
  vector[N] p; // logit-scale
  for(i in 1:N) {
    p[i] = a_chimp[chimp[i]] + bp*prosoc_left[i] + bpc*condition[i]*prosoc_left[i];
  }
}

model {
  for(i in 1:N_chimps) {
    a_chimp[i] ~ normal(mu_chimp, sigma_chimp);
  }
  mu_chimp ~ normal(0,10);
  sigma_chimp ~ cauchy(0, 10);
  bp ~ normal(0, 10);
  bpc ~ normal(0, 10);
  pulled_left ~ binomial_logit(1, p); // likelihood function for data given model
}

generated quantities{
  vector[N] log_lik;
  for(i in 1:N) {
    log_lik[i] = binomial_logit_lpmf(pulled_left[i] | 1, p[i]);
  }
}

"

dat3 <- list(N=nrow(d), 
             pulled_left = d$pulled_left, 
             prosoc_left = d$prosoc_left, 
             condition = d$condition, 
             N_chimps = length(unique(d$actor)), 
             chimp = d$actor)
fit3 <- stan(model_code = model3, data=dat3, iter=3000,chains=2,cores=2,control=list("adapt_delta"=0.81),pars=c("mu_chimp","sigma_chimp","a_chimp", "bp", "bpc", "log_lik"))

summary(fit3, pars=c("mu_chimp","a_chimp","bp","bpc"),probs=c(0.1,0.9))$summary

# compare alphas between the models
# make a list of the posterior samples
post_list <- list(post1 <- as.data.frame(fit1),
                  post2 <- as.data.frame(fit2),
                  post3 <- as.data.frame(fit3))

# Make a function that does mean and quantiles
multi.fun <- function(x) {
  c(mean=mean(x), quantile= quantile(x, probs = c(0.1, 0.9)))
}

# Get mean and quantiles for each of the three models
ints <- as.list(rep(NA,3)) # make three object list, each object is blank
# fill the list with the mean and 80% intervals
for(i in 1:3) {
  ints[[i]] <- apply(X = as.matrix(post_list[[i]][ , grep("^a", names(post_list[[i]]))]), 2, multi.fun)
}

# Plot intercept estimates for the three models
# start with points for model 2
plot(ints[[2]][1,], col="black", pch=16,ylab="Posterior estimates of intercept", xlab="Chimp")
segments(x0= 1:7, y0=ints[[2]][3,], y1=ints[[2]][2, ], col="black")
points(x=1:7+0.2, y=ints[[3]][1,], col="dodger blue", pch=16) # add points for model 3
segments(x0= 1:7+ 0.2, y0=ints[[3]][3,], y1=ints[[3]][2, ], col="dodger blue")
abline(h= ints[[1]][1,], col="orange") # add line for model 1
abline(h= ints[[1]][2:3,], col=adjustcolor("orange", alpha=0.5))

# now add a block effect

model4 <- "
data{
  int N;
  int<lower=0,upper=1> pulled_left[N];
  int prosoc_left[N];
  int condition[N];
  int N_chimps;
  int<lower=1, upper=N_chimps> chimp[N]; // id of chimp
  int N_blocks;
  int<lower=1, upper=N_blocks> block_id[N]; // id of block
}
parameters{
  real a; // global intercept
  real a_block[N_blocks]; // block effect on intercepts
  real a_chimp[N_chimps]; // chimp effect on intercepts
  real bp; // effect of prosoc_left
  real bpc; // interaction effect
  real sigma_chimp; // standard deviation across chimps
  real sigma_block; // standard deviation for how mean varies across block
}
transformed parameters {
  vector[N] p; // logit-scale
  for(i in 1:N) {
    p[i] = a + a_chimp[chimp[i]] + a_block[block_id[i]] + bp*prosoc_left[i] + bpc*condition[i]*prosoc_left[i];
  }
}

model {
  for(i in 1:N_chimps) {
    a_chimp[i] ~ normal(0, sigma_chimp);
  }
  for(i in 1:N_blocks) {
    a_block[i] ~ normal(0, sigma_block);
  }
  a ~ normal(0,10);
  sigma_chimp ~ cauchy(0, 10);
  sigma_block ~ cauchy(0, 10);
  bp ~ normal(0, 10);
  bpc ~ normal(0, 10);
  pulled_left ~ binomial_logit(1, p); // likelihood function for data given model
}

generated quantities{
  vector[N] log_lik;
  vector[N] pulled_left_new; // generate new predictions
  for(i in 1:N) {
    log_lik[i] = binomial_logit_lpmf(pulled_left[i] | 1, p[i]);
    pulled_left_new[i] = binomial_rng(1, inv_logit(p[i]));
  }
}

"

dat4 <- list(N=nrow(d), 
             pulled_left = d$pulled_left, 
             prosoc_left = d$prosoc_left, 
             condition = d$condition, 
             N_chimps = length(unique(d$actor)), 
             chimp = d$actor,
             N_blocks = length(unique(d$block)),
             block_id = d$block)
fit4 <- stan(model_code = model4, data=dat4, iter=3000,chains=2,cores=2,control=list("adapt_delta"=0.81),pars=c("a","sigma_chimp","a_chimp", "a_block","sigma_block", "bp", "bpc", "log_lik","pulled_left_new"))

summary(fit4, pars=c("a","sigma_chimp","a_chimp", "a_block","sigma_block", "bp", "bpc"),probs=c(0.1,0.9))$summary


compare(fit1,fit2,fit3,fit4)

post4 <- as.data.frame(fit4)

# plot predicted v. observed
pred <- post4[, grep("pulled_left_new",colnames(post4))]
pred_ppd <- matrix(NA, nrow=nrow(pred), ncol=4)
colnames(pred_ppd) <- c("0/0", "0/1", "1/0", "1/1") # prosocial-left/condition

counter <- 0
for(i in 0:1) {
  for(j in 0:1) {
    counter = counter + 1
    pred_ppd[ , counter] = rowSums(pred[ , d$prosoc_left==i & d$condition==j])
  }
}

layout(1)
boxplot(pred_ppd)
d$pro_soc_cond <- paste(d$prosoc_left,d$condition,sep="/")
obs <- aggregate(pulled_left~pro_soc_cond, data=d, FUN=sum)
points(1:4, obs$pulled_left, pch=21, bg="dodgerblue", cex=2)

# plot predicted v. observed for individuals
pred_mean <- apply(pred, 2, mean)
pred_hpdi <- apply(pred, 2, HPDI, prob=0.89)
pulled_left_jitter <- d$pulled_left + rnorm(nrow(d), 0 , 0.05)
col_vec <- rep(c("red", "green", "orange", "blue","black","gray","purple"), each=nrow(d)/length(unique(d$actor)))

plot(y=pred_mean, x=pulled_left_jitter, ylim=c(0,1),pch=21, bg=col_vec)


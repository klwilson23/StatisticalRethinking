
library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores=2)


data("chimpanzees")
d <- chimpanzees

# binomial model with just an intercept

model10_1 <- "
data{
  int N;
  int<lower=0,upper=1> L[N]; // did they pull the left lever
}

parameters{
  real a;
}

transformed parameters{
  vector[N] p;
  for(i in 1:N){
    p[i] = inv_logit(a);
  }
}

model{
  a ~ normal(0,10);
  L ~ binomial(1,p);
}

generated quantities{
  vector[N] log_lik;
  for(i in 1:N) {
    log_lik[i] = bernoulli_lpmf(L[i] | p[i]);
  }
}

"


dat <- list(N=nrow(d), L = d$pulled_left)
fit10_1 <- stan(model_code = model10_1, data=dat, iter=2000,chains=2,cores=2)

summary(fit10_1, pars=c("a"),probs=c(0.1,0.9))$summary

model10_2 <- "
data{
  int N;
  int L[N]; // pulled left
  int P[N];
}
parameters{
  real a; // intercept
  real bp; // beta for pro social
}
transformed parameters{
  vector[N] p;
  for(i in 1:N){
    p[i] = a + bp * P[i]; 
  }
}
model{
  a ~ normal(0,10);
  bp ~ normal(0,10);
  L ~ binomial_logit(1,p);
}
generated quantities{
  vector[N] log_lik;
  for(i in 1:N) {
    log_lik[i] = bernoulli_logit_lpmf(L[i] | p[i]);
  }
}

"
dat2 <- list(N=nrow(d), L = d$pulled_left, P = d$prosoc_left)
fit10_2 <- stan(model_code = model10_2, data=dat2, iter=3000,chains=2,cores=2)

summary(fit10_2, pars=c("a","bp"),probs=c(0.1,0.9))$summary

# binomial model: individuals have specific intercepts

model10_4 <- "
data {
  int N;
  int<lower=0, upper=1> L[N]; // pulled left
  vector[N] P; // prosocial left
  vector[N] C; // condition
  int<lower=1,upper=N> MaxA;
  int<lower=1, upper=MaxA> A[N]; // unique actors
}
parameters{
  real a[MaxA];
  real bp;
  real bpc;
}
transformed parameters {
  vector[N] p;
  for(i in 1:N){
    p[i] = a[A[i]] + bp * P[i] + bpc * C[i] * P[i];
  }
}
model{
  for(i in 1:MaxA) {
    a[i] ~ normal(0, 10);
  }
  bp ~ normal(0,10);
  bpc ~ normal(0,10);
  L ~ binomial_logit(1,p);
}
generated quantities{
  vector[N] log_lik;
  for(i in 1:N){
    log_lik[i] = binomial_logit_lpmf(L[i] | 1, p[i]);
  }
}
"
dat3 <- list(N=nrow(d), L = d$pulled_left, P = d$prosoc_left, C = d$condition, MaxA = length(unique(d$actor)), A = d$actor)
fit10_4 <- stan(model_code = model10_4, data=dat3, iter=3000,chains=2,cores=2,control=list("adapt_delta"=0.81))

summary(fit10_4, pars=c("a","bp","bpc"),probs=c(0.1,0.9))$summary
#pairs(fit10_4,pars=c("a","bp","bpc"))
compare(fit10_1,fit10_2,fit10_4)

# aggregated binomial

dA <- aggregate(d$pulled_left,list(prosoc_left=d$prosoc_left,
                                   condition=d$condition,
                                   actor=d$actor),
                sum)
model10_5 <- "
data{
  int N;
  int<lower=0> x[N];
  int<lower=1> MaxX;
  vector[N] P;
  vector[N] C;
}
parameters{
  real a;
  real bp;
  real bpc;
}
transformed parameters{
  vector[N] p;
  for(i in 1:N){
    p[i] = a+bp*P[i] + bpc*C[i]*P[i];
  }
}
model {
  a ~ normal(0,10);
  bp ~ normal(0,10);
  bpc ~ normal(0,10);
  x ~ binomial_logit(MaxX,p);
}

generated quantities{
  vector[N] log_lik;
  for(i in 1:N){
    log_lik[i] = binomial_logit_lpmf(x[i] | MaxX, p[i]);
  }
}

"

dat4 <- list(N=nrow(dA), x = dA$x, P = dA$prosoc_left, C = dA$condition, MaxX = max(dA$x), A = dA$actor)
fit10_5 <- stan(model_code = model10_5, data=dat4, iter=3000,chains=2,cores=2,control=list("adapt_delta"=0.81))

summary(fit10_5, pars=c("a","bp","bpc"),probs=c(0.1,0.9))$summary
data(Kline)
d <- Kline
d$contact_high <- ifelse(d$contact=="high",1,0)


model10_6 <- "
data{
  int N;
  int<lower=0> P[N]; // population
  int<lower=0> T[N]; // total tools in population
  int<lower=0> C[N]; // contact
}
parameters{
  real a;
  real bp;
  real bc;
  real bpc;
}
transformed parameters{
  vector[N] lambda;
  for(i in 1:N){
    lambda[i] = a+bp*log(P[i]) + bpc*C[i]*log(P[i]) + bc*C[i];
  }
}
model {
  a ~ normal(0,100);
  bp ~ normal(0,1);
  bpc ~ normal(0,1);
  bc ~ normal(0,1);
  T ~ poisson_log(lambda);
}

generated quantities{
vector[N] log_lik;
  for(i in 1:N){
  log_lik[i] = poisson_log_lpmf(T[i] | lambda[i]);
  }
}

"

dat5 <- list(N=nrow(d),P=d$population,C=d$contact_high,T=d$total_tools)
fit10_6 <- stan(model_code = model10_6, data=dat5, iter=3000,chains=2,cores=2,control=list("adapt_delta"=0.81))

summary(fit10_6, pars=c("a","bp","bpc","bc"),probs=c(0.1,0.9))$summary


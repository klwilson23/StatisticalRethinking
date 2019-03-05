library(rethinking)
library(MASS)
library(loo)
options(mc.cores = 4)

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

data("Howell1")
d <- Howell1

# check for NAs
sum(is.na(d))

# basic model

stan_linear <- "
data {
  int<lower=0> N;
  vector[N] height;
  vector[N] weight_c;
}

parameters {
  real b; // intercept
  real m; // intercept
  real sigma; // intercept
}

model {
  // declare length of a numeric vector called mu
  vector[N] mu = m*weight_c + b;
  height ~ normal(mu, sigma);
  b ~ normal(mean(height), 20);
  m ~ normal(0, 10);
  sigma ~ cauchy(5, 100);
}

generated quantities {
  vector[N] log_lik; // need to code log-like into model
  for(i in 1:N) {
    log_lik[i] = normal_lpdf(height[i] | b + m*weight_c[i], sigma);
  }
}

"

dat_linear <- list(N = nrow(d),
                   height=d$height,
                   weight_c= d$weight-mean(d$weight))

fit_linear <- stan(model_code=stan_linear,
                   data=dat_linear,
                   chains=2,
                   iter=2000,
                   cores=3,
                   pars=c("b","m","sigma","log_lik"))
beepr::beep(sound=8)

# the log_lik we generate is a matrix of N number of datapoints for M number of posterior samples
# we can extract them to calculate using wAIC package via the loo package

ll_linear <- extract_log_lik(fit_linear)

kyle_wAIC <- sum(-2*apply(ll_linear,2,mean)) + sum(apply(ll_linear,2,var))

kyle_DIC <- sum(-2*apply(ll_linear,1,mean)) + -2*apply(ll_linear,1,mean)

# calculate WAIC
loo::waic(ll_linear)
rethinking::WAIC(fit_linear)

stan_poly <- "
data {
  int<lower=0> N;
  vector[N] height;
  vector[N] weight_c;
  vector[N] weight_sq; //squared centered weight
}

parameters {
  real b; // intercept
  real m; // slope
  real poly; // polynomial inflection point
  real sigma; // variance
}

model {
  // declare length of a numeric vector called mu
  vector[N] mu = m*weight_c + poly*weight_sq + b;
  height ~ normal(mu, sigma);
  b ~ normal(mean(height), 20);
  m ~ normal(0, 10);
  poly ~ normal(0,10);
  sigma ~ cauchy(5, 100);
}

generated quantities {
  vector[N] log_lik; // need to code log-like into model
  for(i in 1:N) {
  log_lik[i] = normal_lpdf(height[i] | b + m*weight_c[i] + poly*weight_sq[i], sigma);
}
}

"

dat_poly <- list(N = nrow(d),
                height=d$height,
                weight_c = d$weight-mean(d$weight),
                 weight_sq = (d$weight-mean(d$weight))^2)

fit_poly <- stan(model_code=stan_poly,
                   data=dat_poly,
                   chains=2,
                   iter=2000,
                   cores=3,
                   pars=c("b","m","poly","sigma","log_lik"))
beepr::beep(sound=8)

summary(fit_poly)

log_lik_poly <- extract_log_lik(fit_poly)

waic(log_lik_poly)
waic(ll_linear)

plot(height~weight_c,data=dat_poly)
precis(fit_poly)

# even more complex!!

stan_poly_cat <- "
data {
  int<lower=0> N;
  vector[N] height;
  vector[N] weight_c;
  vector[N] weight_sq; //squared centered weight
  int Nsex;
  int sex[N]; // read in whether sample is male or female
}

parameters {
  vector[Nsex] b; // intercept
  real m; // slope
  real poly; // polynomial inflection point
  real sigma; // variance
}

transformed parameters {
  vector[N] mu;
  for(i in 1:N) {
    mu[i] = m*weight_c[i] + poly*weight_sq[i] + b[sex[i]];
  }
}

model {
  // declare length of a numeric vector called mu
  //vector[N] mu;
  for(i in 1:N) {
    //mu[i] = m*weight_c[i] + poly*weight_sq[i] + b[sex[i]];
    height[i] ~ normal(mu[i], sigma);
  }
  b ~ normal(mean(height), 20);
  m ~ normal(0, 10);
  poly ~ normal(0,10);
  sigma ~ cauchy(5, 100);
}

generated quantities {
  vector[N] log_lik; // need to code log-like into model
  //vector[N] new_mu;
    for(i in 1:N) {
      //mu[i] = b[sex[i]] + m*weight_c[i] + poly*weight_sq[i];
      log_lik[i] = normal_lpdf(height[i] | mu[i], sigma);
  }
}

"

dat_poly_sex <- list(N = nrow(d),
                     Nsex = length(unique(d$male)),
                     height=d$height,
                     weight_c = d$weight-mean(d$weight),
                     weight_sq = (d$weight-mean(d$weight))^2,
                     sex = as.integer(as.factor(d$male)))

fit_poly_sex <- stan(model_code=stan_poly_cat,
                 data=dat_poly_sex,
                 chains=2,
                 iter=2000,
                 cores=3,
                 pars=c("b","m","poly","sigma","log_lik"))
beepr::beep(sound=8)

summary(fit_poly_sex)

precis(fit_poly_sex,2)
log_lik_poly_sex <- extract_log_lik(fit_poly_sex)

rethinking::compare(fit_linear,fit_poly,fit_poly_sex)

# let's look at LOO

loo_linear <- loo(fit_linear)
loo_poly <- loo(fit_poly)
loo_poly_sex <- loo(fit_poly_sex)
loo_poly_sex
waic(log_lik_poly_sex)


coefs <- coeftab(fit_linear,fit_poly,fit_poly_sex)

params <- c("b","b[1]","b[2]","m","poly","sigma")

par(mfrow=c(1,3))
coeftab_plot(coefs,pars=c("b","b[1]","b[2]"))
coeftab_plot(coefs,pars="m")
coeftab_plot(coefs,pars="poly")

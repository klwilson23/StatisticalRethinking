library(rethinking)
library(MASS)

options(mc.cores=2)
rstan_options(auto_write=TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

data(Howell1)
d <- Howell1
head(Howell1)
plot(height~weight,data=d)

d2 <- d[d$age >= 18, ]
plot(height~weight,data=d2)

layout(matrix(1:2,ncol=2))
dens(d2$height)
hist(d2$height)

layout(1)

# MAP - the posterior distribution via Gaussian approximiation
flist <- alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0, 50)
)

mFit <- map(flist, data=d2)
precis(mFit)

# change the starting values of MAP
# change the prior on mu & sigma

mFit2 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0, 50)
  ),
  data=d2, start = list("mu"=mean(d2$height),"sigma"=sd(d2$height))
)

precis(mFit2)

# extract the variance-covariance matrix
vcov(mFit2)
diag(vcov(mFit2))
mu_SD <- sqrt(diag(vcov(mFit2)))["mu"]

# convert covariance to CORRELATION: bounded between -1 and 1
cov2cor(vcov(mFit2))

# plot the effect of different priors

curve( dnorm(x, mean=178, sd=20), from= 100, to = 250, lwd=2,col="orange",ylim=c(0,4))
curve( dnorm(x, mean=178, sd=0.1), add=T, lwd=2, lty=2,col="dodgerblue")

# compare the strong prior to the posterior: dominated by the prior

curve( dnorm(x, mean=178, sd=0.1), from= 100, to = 250, lwd=2,col="orange",ylim=c(0,3))
curve( dnorm(x, mean=177.86, sd=0.1002), add=T, lwd=2, lty=2,col="dodgerblue")

# compare the vague prior to the posterior: dominated by the data
precis(mFit)
curve( dnorm(x, mean=178, sd=20), from= 100, to = 250, lwd=2,col="orange",ylim=c(0,1))
curve( dnorm(x, mean=154.61, sd=0.41), add=T, lwd=2, lty=2,col="dodgerblue")
hist(d2$height,freq=F,add=T)

# let's code this model in Stan!
# start with a data section
# move to an estimated parameters section
# then write the model section: likelihoods & priors

stanModel <- "
data {
  int<lower=0> N; // this is a comment
  vector[N] height;
}

parameters {
  real mu;
  real sigma;
}

model {
  // log-likelihoods or log-priors
  // the normal likelihood is: normal_lpdf(data/name | mean, standard deviation)
  target += normal_lpdf(height | mu, sigma);
  target += normal_lpdf(mu | 178, 20);
  target += uniform_lpdf(sigma | 0, 50);

}

"
dat <- list("N"=nrow(d2),
            "height"=d2$height)

fit <- stan(model_code=stanModel, data=dat,
            chains=2, iter=2000, cores=2,
            pars=c("mu","sigma"))
fit

# compare model fits with the values from the gaussian approximation



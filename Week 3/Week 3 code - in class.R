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
plot(fit)
traceplot(fit)

# store posterior samples
post <- as.data.frame(fit)
str(post)

# compare model fits with the values from the gaussian approximation

mean(post$mu)
sd(post$mu)
precis(mFit)

# this is dumb. we want to do actual regression
head(d2)
plot(height~weight,data=d2,xlab="Weight (kg)",ylab="Height (cm)",
     pch=21,bg=ifelse(male==1,"orange","dodgerblue"))

stanRegression <- "
data{
  int<lower=0> N;
  vector[N] height;
  vector[N] weight;
}
parameters{
  real b; // intercept
  real m; // slope
  real sigma; // standard deviation (variance = sigma^2)
}
model{
  vector[N] mu = weight*m + b; // define mu as y=mx+b
  target += normal_lpdf(height | mu, sigma); // likelihood based on data
  target += normal_lpdf(b | mean(height), sd(height));
  target += normal_lpdf(m | 0, 10); // set this prior on slope
  target += uniform_lpdf(sigma | 0, 100);
}
"

dat <- list("N"=nrow(d2),
            "height"=d2$height,
            "weight"=d2$weight)

fitRegression <- stan(model_code=stanRegression,
                      data=dat, chains=2, iter = 2000, cores = 2,
                      pars = c("b","m","sigma"))

fitRegression

post2 <- as.data.frame(fitRegression)
head(post2)
fitRegression@sim$samples[[1]][["b"]][1001]

colMeans(post2)
pairs(post2, lower.panel=panel.smooth)
cor(post2)

# 

stanCentered <- "
data{
  int<lower=0> N;
  vector[N] height;
  vector[N] weight_c;
}
parameters{
  real b; // intercept
  real m; // slope
  real sigma; // standard deviation (variance = sigma^2)
}
model{
  vector[N] mu = weight_c*m + b; // define mu as y=mx+b
  target += normal_lpdf(height | mu, sigma); // likelihood based on data

  target += normal_lpdf(b | mean(height), sd(height));
  target += normal_lpdf(m | 0, 10); // set this prior on slope
  target += uniform_lpdf(sigma | 0, 100);
}
"

datCentered <- list("N"=nrow(d2),
                    "height"=d2$height,
                    "weight_c"=d2$weight-mean(d2$weight))

fitCentered <- stan(model_code=stanCentered,
                    data=datCentered,chains=2,iter=2000,cores=2,
                    pars=c("b","m","sigma"))

postCentered <- as.data.frame(fitCentered)
colMeans(postCentered)
colMeans(post2)

pairs(postCentered,lower.panel=panel.smooth)

plot(height~weight_c,data=datCentered,pch=21,bg="dodgerblue")
invisible(sapply(1:250,function(i){
  curve(postCentered$b[i]+postCentered$m[i]*x,add=T,lwd=2,col=adjustcolor("grey",alpha=0.1))
}))

curve(mean(postCentered$b)+mean(postCentered$m)*x,add=T,lwd=5)

# create a posterior predictive distribution and plot the uncertainty in the model

weight.seq <- seq(-20,30,by=1)
post_pred <- matrix(NA,nrow=nrow(postCentered),ncol=length(weight.seq))

for(i in 1:nrow(post_pred))
{
  post_pred[i,] <- rnorm(length(weight.seq),
                         mean=postCentered$b[i] + postCentered$m[i]*weight.seq,sd=postCentered$sigma[i])
}

mu.mean <- apply(post_pred,2,mean)
mu.hdi <- apply(post_pred,2,HPDI,prob=0.89)

plot(height~weight_c,data=datCentered,pch=21,bg="dodgerblue")
polygon(c(weight.seq,rev(weight.seq)),c(mu.hdi[1,],rev(mu.hdi[2,])),col=adjustcolor("dodgerblue",0.7),border=F)
lines(weight.seq,mu.mean,col=adjustcolor("black",0.8),lwd=3)

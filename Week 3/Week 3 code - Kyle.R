library(rethinking)
library(MASS)

options(mc.cores = 2)

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

## load the dataset Howell1

data(Howell1)
d <- Howell1
# filter for adults: people > 18 years old
d2 <- d[ d$age >= 18 , ]

# look at the structure of the data
str( d2 )

layout(matrix(1:2,ncol=2))
dens(d2$height)
hist(d2$height)

layout(1)
## R code 4.25
flist <- alist(
  height ~ dnorm( mu , sigma ) ,
  mu ~ dnorm( 178 , 20 ) ,
  sigma ~ dunif( 0 , 50 )
)

## fit the model via Gaussian approximation
m4.1 <- map( flist , data=d2 )

## R code 4.27
precis( m4.1 )

## R code 4.28
start <- list(
  mu=mean(d2$height),
  sigma=sd(d2$height)
)

# fit the model again with a strong prior on mu: let's say from a previous study
m4.2 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 0.1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2,start=start)
precis( m4.2 )

# extract variance-covariance matrix
vcov( m4.1 )

## interact with variance-covariance matrix
diag( vcov( m4.1 ) )
# convert to correlation coefficients: bounded between -1 and 1
cov2cor( vcov( m4.1 ) )


## explore some priors of the data and how they influence the posterior
# compare vague prior to a strong prior
curve( dnorm( x , 178 , 20 ) , from=100 , to=250 ,lwd=2,col="orange")
curve( dnorm( x , 178 , 0.1 ) , from=100 , to=250 ,add=T,lwd=2,lty=2,col="dodgerblue")

# compare strong prior and posterior distribution - dominated by prior
curve( dnorm( x , 178 , 0.1 ) , from=100 , to=250 ,lwd=2,lty=2,col="dodgerblue")
curve( dnorm( x , coef(m4.2)["mu"] , sqrt(diag(vcov(m4.2))["mu"])),add=T,lwd=2,lty=1,col="purple")

# compare vague prior and posterior distribution - dominated by data
curve( dnorm( x , 178 , 20 ) , from=100 , to=250 ,lwd=2,lty=2,col="dodgerblue",ylim=c(0,0.3))
curve( dnorm( x , coef(m4.1)["mu"] , sqrt(diag(vcov(m4.1))["mu"])),add=T,lwd=2,lty=1,col="purple")

# let's code this model in Stan:
# start with the data section
# move to estimated parameters section
# then write the model section: likelihoods!
stanModel <- "
data {
  int<lower=0> N;
  vector[N] height;
}
parameters {
  real mu;
  real<lower=0,upper=50> sigma;
}
model {
  // comment log-likelihood for normal is: normal_lpdf(data | mean, standard deviation);
  target += normal_lpdf(height | mu, sigma);
  target += normal_lpdf(mu | 178, 20);
}

"

# make a list for the data: "N" and "height"

dat <- list(N = nrow(d2),
            height = d2$height)

fit <- stan(model_code=stanModel, data=dat,
            chains = 2, iter = 2000,cores=2,
            pars=c("mu","sigma"),
            init=function(){return(list("mu"=mean(dat$height),"sigma"=sd(dat$height)))})

print(fit, probs = c(0.10, 0.5, 0.9))
fit
plot(fit)
# store posterior samples
post <- as.data.frame(fit)
mean(post$mu);sd(post$mu)
precis(m4.1)

# okay, this is dumb. We want regression!!

plot(height~weight,data=d2,xlab="Weight (kg)",ylab="Height (cm)",pch=21,bg=ifelse(male==1,"dodgerblue","orange"))

stanRegression <- "
data {
  int<lower=0> N;
  vector[N] height;
  vector[N] weight;
}
parameters {
  real b; // intercept
  real m; // slope
  real<lower=0,upper=50> sigma;
}
model {
  // declare length of a numeric vector called mu
  vector[N] mu = m * weight + b;
  target += normal_lpdf(height | mu, sigma);
  target += normal_lpdf(b | mean(height), 20);
  target += normal_lpdf(m | 0, 10);
}

"

dat <- list(N = nrow(d2),
            height = d2$height,
            weight = d2$weight)

fit <- stan(model_code=stanRegression, data=dat,
            chains = 2, iter = 2000,cores=2,
            pars=c("b","m","sigma"),
            init=function(){return(list("b"=mean(dat$height),"m"=0,"sigma"=sd(dat$height)))})
fit
plot(fit)

post <- as.data.frame(fit)
colMeans(post)
pairs(post,lower.panel=panel.smooth)
cor(post)
# that's weird! the slope and intercept are correlated!

# to reduce the correlation, let's center the data on the average weight:
stanCentered <- "
data {
  int<lower=0> N;
  vector[N] height;
  vector[N] weight_c;
}
parameters {
  real b; // intercept
  real m; // slope
  real<lower=0,upper=50> sigma;
}
model {
  // declare length of a numeric vector called mu
  vector[N] mu = m * weight_c + b;
  target += normal_lpdf(height | mu, sigma);
  target += normal_lpdf(b | mean(height), 20);
  target += normal_lpdf(m | 0, 10);
}

"

datCentered <- list(N = nrow(d2),
                    height = d2$height,
                    weight_c = (d2$weight-mean(d2$weight)))

fitCentered <- stan(model_code=stanCentered, data=datCentered,
                    chains = 2, iter = 2000,cores=2,
                    pars=c("b","m","sigma"),
                    init=function(){return(list("b"=mean(datCentered$height),"m"=0,"sigma"=sd(datCentered$height)))})
fitCentered
plot(fitCentered)

postCentered <- as.data.frame(fitCentered)
colMeans(postCentered) # average height at the average weight
colMeans(post) # average height when weight is 0
pairs(postCentered,lower.panel=panel.smooth)
cor(postCentered)

# plot the average response: what are we missing?
plot(height~weight_c,data=datCentered,pch=21,bg="dodgerblue")
invisible(sapply(1:250,function(i){
  curve(postCentered$b[i]+postCentered$m[i]*x,add=T,lwd=2,col=adjustcolor("grey",alpha=0.1))
}))

curve(mean(postCentered$b)+mean(postCentered$m)*x,add=T,lwd=2)

plot(height~weight,data=dat,pch=21,bg="dodgerblue")
invisible(sapply(1:250,function(i){
  curve(post$b[i]+post$m[i]*x,add=T,lwd=2,col=adjustcolor("grey",alpha=0.1))
}))

NsubSamp <- 10
datv1 <- list(N = NsubSamp,
                    height = d2$height[1:NsubSamp],
                    weight_c = (d2$weight[1:NsubSamp]-mean(d2$weight[1:NsubSamp])))

fitv1 <- stan(model_code=stanCentered, data=datv1,
                    chains = 2, iter = 2000,cores=2,
                    pars=c("b","m","sigma"))
postv1 <- as.data.frame(fitv1)

layout(matrix(1:2,nrow=2))
plot(height~weight_c,data=datv1,pch=21,bg="dodgerblue")
invisible(sapply(1:250,function(i){
  curve(postv1$b[i]+postv1$m[i]*x,add=T,lwd=2,col=adjustcolor("grey",alpha=0.1))
}))

NsubSamp <- 25
datv2 <- list(N = NsubSamp,
              height = d2$height[1:NsubSamp],
              weight_c = (d2$weight[1:NsubSamp]-mean(d2$weight[1:NsubSamp])))

fitv2 <- stan(model_code=stanCentered, data=datv2,
              chains = 2, iter = 2000,cores=2,
              pars=c("b","m","sigma"))
postv2 <- as.data.frame(fitv2)

plot(height~weight_c,data=datv2,pch=21,bg="dodgerblue")
invisible(sapply(1:250,function(i){
  curve(postv2$b[i]+postv2$m[i]*x,add=T,lwd=2,col=adjustcolor("grey",alpha=0.1))
}))

# create a posterior predictive distribution and plot the uncertainty in the model
weight.seq <- seq(from=-20,to=30,by=1)
post_pred <- matrix(NA,nrow=nrow(postCentered),ncol=length(weight.seq))
for(i in 1:nrow(post_pred))
{
  post_pred[i,] <- rnorm(length(weight.seq),
                         mean=postCentered$b[i]+postCentered$m[i]*weight.seq,
                         sd=postCentered$sigma[i])
}

mu.mean <- apply(post_pred,2,mean)
mu.hdi <- apply(post_pred,2,HPDI,prob=0.89)

layout(1)
plot(height~weight_c,data=datCentered,pch=21,bg="dodgerblue")
polygon(c(weight.seq,rev(weight.seq)),c(mu.hdi[1,],rev(mu.hdi[2,])),col=adjustcolor("dodgerblue",0.2),border=F)
lines(weight.seq,mu.mean,col=adjustcolor("black",0.8),lwd=3)

# we've been looking at just adults for a few reasons: young people put on height more quickly than weight!
# but what about the entire dataset?

plot(height~weight,data=d)

# what about a polynomial regression?
# can use polynomial regression to understand height ~ weight for ALL ages of people

# to reduce the correlation, let's center the data on the average weight:
stanPoly <- "
data {
  int<lower=0> N;
  vector[N] height;
  vector[N] weight_c;
  vector[N] weight_sq;
}
parameters {
  real b; // intercept
  real m; // slope
  real poly; // polynomial inflection point
  real<lower=0,upper=50> sigma;
}
model {
  // declare length of a numeric vector called mu
  vector[N] mu = m * weight_c + poly * weight_sq + b;
  target += normal_lpdf(height | mu, sigma);
  target += normal_lpdf(b | mean(height), 20); // prior on intercept
  target += normal_lpdf(m | 0, 10); // prior on slope
  target += normal_lpdf(poly | 0, 10); // prior on polynomial term
}

"

datPoly <- list(N = nrow(d),
                height = d$height,
                weight_c = (d$weight-mean(d$weight)),
                weight_sq = (d$weight-mean(d$weight))^2)

fitPoly <- stan(model_code=stanPoly, data=datPoly,
                    chains = 2, iter = 2000,cores=2,
                    pars=c("b","m","poly","sigma"))
fitPoly

postPoly <- as.data.frame(fitPoly)
colMeans(postPoly) # average height at the average weight
pairs(postPoly,lower.panel=panel.smooth)
cor(postPoly)

# create a posterior predictive distribution and plot the uncertainty in the model
weight.seq <- seq(from=-40,to=40,by=1)
post_pred <- matrix(NA,nrow=nrow(postPoly),ncol=length(weight.seq))
for(i in 1:nrow(post_pred))
{
  post_pred[i,] <- rnorm(length(weight.seq),
                         mean=postPoly$b[i]+postPoly$m[i]*weight.seq + postPoly$poly[i]*weight.seq^2,
                         sd=postPoly$sigma[i])
}

mu.mean <- apply(post_pred,2,mean)
mu.hdi <- apply(post_pred,2,HPDI,prob=0.89)

layout(1)
plot(height~weight_c,data=datPoly,pch=21,bg="dodgerblue")
polygon(c(weight.seq,rev(weight.seq)),c(mu.hdi[1,],rev(mu.hdi[2,])),col=adjustcolor("dodgerblue",0.2),border=F)
lines(weight.seq,mu.mean,col=adjustcolor("black",0.8),lwd=3)


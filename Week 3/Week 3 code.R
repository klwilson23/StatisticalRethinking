library(rethinking)
library(MASS)

## R code 4.1
pos <- replicate( 1000 , sum( runif(16,-1,1) ) )

## R code 4.2
prod( 1 + runif(12,0,0.1) )

## R code 4.3
growth <- replicate( 10000 , prod( 1 + runif(12,0,0.1) ) )
dens( growth , norm.comp=TRUE )

## R code 4.4
big <- replicate( 10000 , prod( 1 + runif(12,0,0.5) ) )
small <- replicate( 10000 , prod( 1 + runif(12,0,0.01) ) )

## R code 4.5
log.big <- replicate( 10000 , log(prod(1 + runif(12,0,0.5))) )

## R code 4.6
w <- 6; n <- 9;
p_grid <- seq(from=0,to=1,length.out=100)
posterior <- dbinom(w,n,p_grid)*dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)

## R code 4.7

data(Howell1)
d <- Howell1

## R code 4.8
str( d )

## R code 4.9
d$height

## R code 4.10
d2 <- d[ d$age >= 18 , ]

## R code 4.11
curve( dnorm( x , 178 , 20 ) , from=100 , to=250 )

## R code 4.12
curve( dunif( x , 0 , 50 ) , from=-10 , to=60 )

## R code 4.13
sample_mu <- rnorm( 1e4 , 178 , 20 )
sample_sigma <- runif( 1e4 , 0 , 50 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h )

## R code 4.14
mu.list <- seq( from=140, to=160 , length.out=200 )
sigma.list <- seq( from=4 , to=9 , length.out=200 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )
post$LL <- sapply( 1:nrow(post) , function(i) sum( dnorm(
  d2$height ,
  mean=post$mu[i] ,
  sd=post$sigma[i] ,
  log=TRUE ) ) )
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
  dunif( post$sigma , 0 , 50 , TRUE )
post$prob <- exp( post$prod - max(post$prod) )

## R code 4.15
contour_xyz( post$mu , post$sigma , post$prob )

## R code 4.16
image_xyz( post$mu , post$sigma , post$prob )

## R code 4.17
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]

## R code 4.18
plot( sample.mu , sample.sigma , cex=0.5 , pch=16 , col=col.alpha(rangi2,0.1) )

## R code 4.19
dens( sample.mu )
dens( sample.sigma )

## R code 4.20
HPDI( sample.mu )
HPDI( sample.sigma )

## R code 4.21
d3 <- sample( d2$height , size=20 )

## R code 4.22
mu.list <- seq( from=150, to=170 , length.out=200 )
sigma.list <- seq( from=4 , to=20 , length.out=200 )
post2 <- expand.grid( mu=mu.list , sigma=sigma.list )
post2$LL <- sapply( 1:nrow(post2) , function(i)
  sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] ,
              log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
  dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
                        prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
      col=col.alpha(rangi2,0.1) ,
      xlab="mu" , ylab="sigma" , pch=16 )

## R code 4.23
dens( sample2.sigma , norm.comp=TRUE )

## R code 4.24
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

## R code 4.25
flist <- alist(
  height ~ dnorm( mu , sigma ) ,
  mu ~ dnorm( 178 , 20 ) ,
  sigma ~ dunif( 0 , 50 )
)

## R code 4.26
m4.1 <- map( flist , data=d2 )

## R code 4.27
precis( m4.1 )

## R code 4.28
start <- list(
  mu=mean(d2$height),
  sigma=sd(d2$height)
)

## R code 4.29
m4.2 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 0.1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2 )
precis( m4.2 )

## R code 4.30
vcov( m4.1 )

## R code 4.31
diag( vcov( m4.1 ) )
cov2cor( vcov( m4.1 ) )

## R code 4.32
post <- extract.samples( m4.1 , n=1e4 )
head(post)

## R code 4.33
precis(post)

## R code 4.34

post <- mvrnorm( n=1e4 , mu=coef(m4.1) , Sigma=vcov(m4.1) )

## R code 4.35
m4.1_logsigma <- map(
  alist(
    height ~ dnorm( mu , exp(log_sigma) ) ,
    mu ~ dnorm( 178 , 20 ) ,
    log_sigma ~ dnorm( 2 , 10 )
  ) , data=d2 )

## R code 4.36
post <- extract.samples( m4.1_logsigma )
sigma <- exp( post$log_sigma )

## R code 4.37
plot( d2$height ~ d2$weight )

## R code 4.38
# load data again, since it's a long way back

data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

# fit model
m4.3 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ dnorm( 156 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2 )

## R code 4.39
m4.3 <- map(
  alist(
    height ~ dnorm( a + b*weight , sigma ) ,
    a ~ dnorm( 178 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2 )

## R code 4.40
precis( m4.3 )

## R code 4.41
precis( m4.3 , corr=TRUE )

## R code 4.42
d2$weight.c <- d2$weight - mean(d2$weight)

## R code 4.43
m4.4 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*weight.c ,
    a ~ dnorm( 178 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2 )

## R code 4.44
precis( m4.4 , corr=TRUE )

## R code 4.45
plot( height ~ weight , data=d2 )
abline( a=coef(m4.3)["a"] , b=coef(m4.3)["b"] )

## R code 4.46
post <- extract.samples( m4.3 )

## R code 4.47
post[1:5,]

## R code 4.48
N <- 10
dN <- d2[ 1:N , ]
mN <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ dnorm( 178 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=dN )

## R code 4.49
# extract 20 samples from the posterior
post <- extract.samples( mN , n=20 )

# display raw data and sample size
plot( dN$weight , dN$height ,
      xlim=range(d2$weight) , ylim=range(d2$height) ,
      col=rangi2 , xlab="weight" , ylab="height" )
mtext(concat("N = ",N))

# plot the lines, with transparency
for ( i in 1:20 )
  abline( a=post$a[i] , b=post$b[i] , col=col.alpha("black",0.3) )

## R code 4.50
mu_at_50 <- post$a + post$b * 50

## R code 4.51
dens( mu_at_50 , col=rangi2 , lwd=2 , xlab="mu|weight=50" )

## R code 4.52
HPDI( mu_at_50 , prob=0.89 )

## R code 4.53
mu <- link( m4.3 )
str(mu)

## R code 4.54
# define sequence of weights to compute predictions for
# these values will be on the horizontal axis
weight.seq <- seq( from=25 , to=70 , by=1 )

# use link to compute mu
# for each sample from posterior
# and for each weight in weight.seq
mu <- link( m4.3 , data=data.frame(weight=weight.seq) )
str(mu)

## R code 4.55
# use type="n" to hide raw data
plot( height ~ weight , d2 , type="n" )

# loop over samples and plot each mu value
for ( i in 1:100 )
  points( weight.seq , mu[i,] , pch=16 , col=col.alpha(rangi2,0.1) )

## R code 4.56
# summarize the distribution of mu
mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

## R code 4.57
# plot raw data
# fading out points to make line and interval more visible
plot( height ~ weight , data=d2 , col=col.alpha(rangi2,0.5) )

# plot the MAP line, aka the mean mu for each weight
lines( weight.seq , mu.mean )

# plot a shaded region for 89% HPDI
shade( mu.HPDI , weight.seq )

## R code 4.58
post <- extract.samples(m4.3)
mu.link <- function(weight) post$a + post$b*weight
weight.seq <- seq( from=25 , to=70 , by=1 )
mu <- sapply( weight.seq , mu.link )
mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

## R code 4.59
sim.height <- sim( m4.3 , data=list(weight=weight.seq) )
str(sim.height)

## R code 4.60
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.61
# plot raw data
plot( height ~ weight , d2 , col=col.alpha(rangi2,0.5) )

# draw MAP line
lines( weight.seq , mu.mean )

# draw HPDI region for line
shade( mu.HPDI , weight.seq )

# draw PI region for simulated heights
shade( height.PI , weight.seq )

## R code 4.62
sim.height <- sim( m4.3 , data=list(weight=weight.seq) , n=1e4 )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.63
post <- extract.samples(m4.3)
weight.seq <- 25:70
sim.height <- sapply( weight.seq , function(weight)
  rnorm(
    n=nrow(post) ,
    mean=post$a + post$b*weight ,
    sd=post$sigma ) )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.64

data(Howell1)
d <- Howell1
str(d)

## R code 4.65
d$weight.s <- ( d$weight - mean(d$weight) )/sd(d$weight)

## R code 4.66
d$weight.s2 <- d$weight.s^2
m4.5 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b1*weight.s + b2*weight.s2 ,
    a ~ dnorm( 178 , 100 ) ,
    b1 ~ dnorm( 0 , 10 ) ,
    b2 ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d )

## R code 4.67
precis( m4.5 )

## R code 4.68
weight.seq <- seq( from=-2.2 , to=2 , length.out=30 )
pred_dat <- list( weight.s=weight.seq , weight.s2=weight.seq^2 )
mu <- link( m4.5 , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.89 )
sim.height <- sim( m4.5 , data=pred_dat )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

## R code 4.69
plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) )
lines( weight.seq , mu.mean )
shade( mu.PI , weight.seq )
shade( height.PI , weight.seq )

## R code 4.70
d$weight.s3 <- d$weight.s^3
m4.6 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b1*weight.s + b2*weight.s2 + b3*weight.s3 ,
    a ~ dnorm( 178 , 100 ) ,
    b1 ~ dnorm( 0 , 10 ) ,
    b2 ~ dnorm( 0 , 10 ) ,
    b3 ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d )

## R code 4.71
plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) , xaxt="n" )

## R code 4.72
at <- c(-2,-1,0,1,2)
labels <- at*sd(d$weight) + mean(d$weight)
axis( side=1 , at=at , labels=round(labels,1) )

## R code 4.73
plot( height ~ weight , data=Howell1 ,
      col=col.alpha(rangi2,0.4) )
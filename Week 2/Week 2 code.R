## R code 2.3

# define the size and values of the grid
Ngrid <- 20
p_grid <- seq( from=0 , to=1 , length.out=Ngrid )

# define prior
prior1 <- rep( 1 , Ngrid )

## R code 2.5
prior2 <- ifelse( p_grid < 0.5 , 0 , 1 )
prior3 <- exp( -5*abs( p_grid - 0.5 ) )

plot(prior1,ylim=range(prior1,prior2,prior3),pch=21,bg=rgb(1,0,0,0.7))
points(prior2,pch=21,bg=rgb(0,1,1,0.7))
points(prior3,pch=21,bg=rgb(0,0,1,0.7))

# compute likelihood at each value in grid
# binomial distribution: flip 9 coins, you observe 6 heads. What is the probability the coin is weighted at any particular value?
prior <- prior3
likelihood <- dbinom( 6 , size=9 , prob=p_grid )

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

## R code 2.4
plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

## R code 2.6

library(rethinking)
globe.qa <- map(
  alist(
    w ~ dbinom(9,p) ,  # binomial likelihood
    p ~ dunif(0,1)     # uniform prior
  ) ,
  data=list(w=6) ) # pass 'map()' the number of successes

# display summary of quadratic approximation
precis( globe.qa )

## R code 2.7
# analytical calculation
w <- 6 # number of successes
n <- 9 # number of trials
curve( dbeta( x , w+1 , n-w+1 ) , from=0 , to=1 )
# quadratic approximation
curve( dnorm( x , 0.67 , 0.16 ) , lty=2 , add=TRUE )


# let's see where this breaks down:
# let's approximate it via gaussian
globe.qa <- map(
  alist(
    w ~ dbinom(10,p) ,  # binomial likelihood
    p ~ dunif(0,1)     # uniform prior
  ) ,
  data=list(w=9), start=list(p=0.99)) # pass 'map()' the number of successes

# display summary of quadratic approximation
precis( globe.qa )
w <- 9 # number of successes
n <- 10 # number of trials
curve( dbeta( x , w+1 , n-w+1 ) , from=0 , to=1 )
# quadratic approximation
curve( dnorm( x , precis( globe.qa )@output$Mean , precis( globe.qa )@output$StdDev ) , lty=2 , add=TRUE )


## R code 3.1
# PrPV = probability of positive of vampire when its a vampire vampires
# PrPM = probability of positive of vampire when its a mortal person
# PrV = probability of vampire in population
# 1-PrV = probability of mortal in population
# PrP = probability of a positive for anyone

PrPV <- 0.95
PrPM <- 0.01
PrV <- 0.001
PrP <- PrPV*PrV + PrPM*(1-PrV) # whats the probability of a positive
( PrVP <- PrPV*PrV / PrP ) # whats the probability of a vampire with a  successful test

# instead say there are 100,000 people
# there are 100 vampires
# 95 vampires test positive
# 90,900 mortals, 999 will test positive for vampirism
# what's the probability a positive test is a vampire?
# Pr(vampire|positive)
95/(95+999)


## R code 3.2
Ngrid <- 1000
p_grid <- seq( from=0 , to=1 , length.out=Ngrid )
prior <- rep( 1 , Ngrid )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

## R code 3.3
NewSamp <- 1e5
samples <- sample( p_grid , prob=posterior , size=NewSamp , replace=TRUE )

## R code 3.4
plot( samples )

## R code 3.5
library(rethinking)
dens( samples )
hist(samples,freq=F,add=T,col=adjustcolor("dodgerblue", 0.4))

## R code 3.6
# add up posterior probability where p < 0.5
sum( posterior[ p_grid < 0.5 ] )

## R code 3.7
sum( samples < 0.5 ) / NewSamp

## R code 3.8
sum( samples > 0.5 & samples < 0.75 ) / NewSamp

# how would you calculate how much is greater than 0.75

dens(samples,lwd=2,xlab="proportion of water") # plot the density/probability of each value
shade(density(samples),lim=c(-Inf,0.5),col=adjustcolor("dodgerblue",0.5)) # ~ 17% probability here
shade(density(samples),lim=c(0.5,0.75),col=adjustcolor("orange",0.5)) # ~ 60% of the probability is here
shade(density(samples),lim=c(0.75,1),col=adjustcolor("grey50",0.5))

## R code 3.9
quantile( samples , 0.8 )

## R code 3.10
quantile( samples , c( 0.1 , 0.9 ) ) # called the 80% credible interval

## R code 3.11

# Now imagine we got 3 "waters" in a row - let's calculate the posterior

Ngrid <- 1000 # define how long our grid is
NewSamps <- 1e4
p_grid <- seq( from=0 , to=1 , length.out=Ngrid )
prior <- rep(1,Ngrid)
likelihood <- dbinom( 3 , size=3 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample( p_grid , size=NewSamps , replace=TRUE , prob=posterior )

dens(samples)

## R code 3.12
PI( samples , prob=0.5 )

## R code 3.13
HPDI( samples , prob=0.5 )

## R code 3.14
p_grid[ which.max(posterior) ]

## R code 3.15
chainmode( samples , adj=0.01 )

## R code 3.16
mean( samples )
median( samples )

## R code 3.17
sum( posterior*abs( 0.5 - p_grid ) )

## R code 3.18
loss <- sapply( p_grid , function(d) sum( posterior*abs( d - p_grid ) ) )

## R code 3.19
p_grid[ which.min(loss) ]

## R code 3.20
# simulate fake data which mirrors the "small world" assumptions

dbinom( 0:2 , size=2 , prob=0.7 )

## R code 3.21
rbinom( 1 , size=2 , prob=0.7 )

## R code 3.22
rbinom( 10 , size=2 , prob=0.7 )

## R code 3.23
Nsamps <- 1e5
dummy_w <- rbinom( Nsamps , size=2 , prob=0.7 )
table(dummy_w)/Nsamps

## R code 3.24
dummy_w <- rbinom( Nsamps , size=9 , prob=0.7 )
simplehist( dummy_w , xlab="dummy water count" ,xaxt="n")
axis(1,at=names(table(dummy_w)),labels=round(as.numeric(names(table(dummy_w)))/max(dummy_w),2))

## R code 3.25
# model checks and adequacy
Nsamps <- 1e4
w <- rbinom( Nsamps , size=9 , prob=0.6 )

## R code 3.26
w <- rbinom( Nsamps , size=9 , prob=samples )

hist(w[round(samples,1)==0.3],xlab="dummy water count",breaks=0:9,col=adjustcolor("orange",0.5),ylim=c(0,max(hist(w,plot=F)$count)),main="") # plot the density/probability of each value
hist(w[round(samples,1)==0.6],xlab="dummy water count",breaks=0:9,col=adjustcolor("dodgerblue",0.5),ylim=c(0,max(hist(w,plot=F)$count)),add=T) # plot the density/probability of each value
hist(w[round(samples,1)==0.9],xlab="dummy water count",breaks=0:9,col=adjustcolor("darkblue",0.5),ylim=c(0,max(hist(w,plot=F)$count)),add=T) # plot the density/probability of each value

hist(w,breaks=0:9,col="steelblue")
abline(v=6,lwd=3,col="brown")

# calculate number of switches and maximum run length
maxSwitch <- maxRun <- rep(NA,length(w))
for(i in 1:length(w))
{
  fakeDat <- rbinom(9,size=1,prob=samples[i])
  maxRun[i] <- max(rle(fakeDat)$lengths)
  maxSwitch[i] <- length(rle(fakeDat)$lengths)-1
}

simplehist(maxRun)
abline(v=max(rle(c(1,0,1,1,1,0,1,0,1))$lengths),lwd=4,col="orange")

simplehist(maxSwitch)
abline(v=length(rle(c(1,0,1,1,1,0,1,0,1))$lengths)-1,lwd=4,col="orange")

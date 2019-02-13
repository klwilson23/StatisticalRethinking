library(rethinking)
library(MASS)

options(mc.cores = 4)

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

## load the dataset milk

# does brain composition explain milk energy

data(milk)
milk

# use complete.cases()

milk_c <- milk[complete.cases(milk),]

# look at the data - check for correlation
pairs(data=milk_c, ~ kcal.per.g + neocortex.perc + mass,lower.panel=panel.smooth)
cor(milk_c$neocortex.perc, milk_c$mass)


# STAN code
stan_model <- "
  data{
    int N;
    vector[N] kcal_per_g;
    vector[N] neocortex_c; // centering the data
  }
  parameters {
    real a; // intercept
    real b_neo; // slope with neocortex
    real sigma; // variance with kcal_per_g
  }
  model{
    vector[N] mu = a + b_neo*neocortex_c; // mean response in kcal
    kcal_per_g ~ normal(mu, sigma); // likelihood
    a ~ normal(mean(kcal_per_g), sd(kcal_per_g)); // prior on alpha
    b_neo ~ normal(0, 1); // prior on slope
    sigma ~ cauchy(0, 1); // half-t distribution for variance
  }

"

dat_centered <- list("N"=nrow(milk_c),
                     "kcal_per_g"=milk_c$kcal.per.g,
                     "neocortex_c"=milk_c$neocortex.perc-mean(milk_c$neocortex.perc))
milk_fit <- stan(model_code=stan_model,
                 data=dat_centered,
                 chains=2,iter=2000,pars=c("a","b_neo","sigma"))

apply(summary(milk_fit)$summary,2,round,3)

# plot our model of kcal_per_g ~ neocortex_c
layout(matrix(1:4,nrow=2,ncol=2,byrow=T))

# get the posterior as a data.frame
post1 <- as.data.frame(milk_fit)
dim(post1)

# make the plot
plot(kcal_per_g~neocortex_c,data=dat_centered,pch=21,bg="dodgerblue",main="~ %Neocortex")
invisible(sapply(1:250, function(i){
  curve(post1[i,"a"]+post1[i,"b_neo"]*x,min=-15,max=15,add=T,lwd=2,col=adjustcolor("grey",alpha=0.2))}))
curve(mean(post1[,"a"])+mean(post1[,"b_neo"])*x,lwd=2,add=T,col="black")

# posterior predictive distribution

neocortex_seq <- seq(from=-20,to=20,by=1)

post_pred <- matrix(NA, nrow=nrow(post1),ncol=length(neocortex_seq))
for(i in 1:nrow(post1))
{
  post_pred[i,] <- rnorm(length(neocortex_seq),mean=post1[i,"a"]+post1[i,"b_neo"]*neocortex_seq,sd=post1[i,"sigma"])
}

mn_ppd <- colMeans(post_pred)
ci_ppd <- apply(post_pred,2,HPDI,0.89)

plot(kcal_per_g~neocortex_c,data=dat_centered,pch=21,bg="dodgerblue", main="Posterior predictive distribution",ylim=c(0,1))
lines(neocortex_seq,mn_ppd,lwd=2)
rethinking::shade(ci_ppd,neocortex_seq,0.89)

# multiple regression. MORE PREDICTORS!!!

# kcal_per_g ~ neocortex_c + logmass_c

stan_code2 <- "
  data{
  int N;
  vector [N] kcal_per_g;
  vector[N] neocortex_c;
  vector[N] logmass_c;
  }
  parameters{
  real a;
  real b_neo;
  real b_mass;
  real sigma;
  }
  model{
  vector[N] mu = a + b_neo*neocortex_c + b_mass*logmass_c;
  kcal_per_g ~ normal(mu,sigma);
  a ~ normal(0,100);
  b_neo ~ normal(0,100);
  b_mass ~ normal(0,100);
  sigma ~ cauchy(0,1);
  }

"

dat_Multi <- list("N"=nrow(milk_c),
                  "kcal_per_g"=milk_c$kcal.per.g,
                  "neocortex_c"=milk_c$neocortex.perc-mean(milk_c$neocortex),
                  "logmass_c"=log(milk_c$mass)-mean(log(milk_c$mass)))
fit <- stan(model_code=stan_code2,data=dat_Multi,iter=2000,chains=2,pars=c("a","b_neo","b_mass","sigma"))

apply(summary(fit)$summary,2,round,3)

post2 <- as.data.frame(fit)

post_pred <- matrix(NA, nrow=nrow(post2),ncol=nrow(milk_c))
for(i in 1:nrow(post2))
{
  post_pred[i,] <- rnorm(nrow(milk_c),mean=post2[i,"a"]+post2[i,"b_neo"]*dat_Multi$neocortex_c+post2[i,"b_mass"]*dat_Multi$logmass_c,sd=post2[i,"sigma"])
}

hpdi <- apply(post_pred,2,HPDI,0.89)
mn_ppd <- colMeans(post_pred)
plot(dat_Multi$kcal_per_g,mn_ppd,pch=21,bg="dodgerblue",ylim=range(hpdi))
segments(x0=dat_Multi$kcal_per_g,y0=hpdi[1,],y1=hpdi[2,])
abline(b=1,a=0)



post_pred2 <- matrix(NA, nrow=nrow(post2), ncol=length(neocortex_seq))
for(i in 1:nrow(post_pred)) {
  post_pred2[i,] <- rnorm(length(neocortex_seq),
                         mean=post2$a[i] + post2$b_neo[i]*neocortex_seq,
                         sd=post2$sigma[i])
}


mn_ppd <- colMeans(post_pred2)
ci_ppd <- apply(post_pred2,2,HPDI,0.89)

plot(kcal_per_g~neocortex_c,data=dat_centered,pch=21,bg="dodgerblue", main="Posterior predictive distribution",ylim=c(0,1))
lines(neocortex_seq,mn_ppd,lwd=2)
rethinking::shade(ci_ppd,neocortex_seq,0.89)


# categorical variables
str(milk)

# kcal.per.g ~ cladeA + cladeB + cladeC + cladeD

X <- model.matrix(kcal.per.g~clade,milk)

stanCat <- "
data{
  int N;
  int K; // number of clades
  vector[N] kcal_per_g;
  matrix[N,K] clades; // matrix of clades
}
parameters{
  vector[K] beta;
  real sigma;
}
model{
  kcal_per_g ~ normal(clades*beta,sigma);
  
  beta[1] ~ normal(mean(kcal_per_g),10); // prior for the intercept
  for(i in 2:K)
  {
    beta[i] ~ normal(0,1); // prior on intercepts
  }

  sigma ~ cauchy(0,1);
}

"

datCat <- list("N"=nrow(milk),
               "K"=length(unique(milk$clade)),
               "kcal_per_g"=milk$kcal.per.g,
               "clades"=X)
fitCat <- stan(model_code=stanCat,data=datCat,iter=2000,chains=2,pars=c("beta","sigma"))

fitCat

aggregate(kcal.per.g~clade,data=milk,mean)

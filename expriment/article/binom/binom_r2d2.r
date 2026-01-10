library(brms)
n.gibbs <- 500
burnin  <- 100
sim <- 30
n	<- 4000
p1 <- 100
p2 <- 100
p3 <- 600
p	<- p1+p2+p3

set.seed(1)

X.mat	<- array(rnorm(n*p,mean = 0, sd = sqrt(1/n)),dim=c(n,p))

beta.1	<- c(rep(-10,p1),rep(10,p2),rep(0,p3))
beta.2	<- rnorm(p,mean = 3, sd = 4)
beta.3	<- rbinom(p,1,0.5)*rnorm(p,mean = 7, sd = 1)

inv.log.odds <- function(Xbeta){
    return(1/(1+ exp(-Xbeta)))
}


# Scenario 1

set.seed(1)

y.1		  <- rbinom(n,1,inv.log.odds(X.mat %*% beta.1))


# Scenario 2

set.seed(1)

y.2		  <- rbinom(n,1,inv.log.odds(X.mat %*% beta.2))

set.seed(1)

y.3		  <- rbinom(n,1,inv.log.odds(X.mat %*% beta.3))

library(brms)

# give the columns of X.mat sensible names
p <- ncol(X.mat)
colnames(X.mat) <- paste0("x", 1:p)

dat1 <- data.frame(y = y.1, X.mat)
dat2 <- data.frame(y = y.2, X.mat)
dat3 <- data.frame(y = y.3, X.mat)
priors <- prior(
    R2D2(
        mean_R2  = 0.5,
        prec_R2  = 2,
        cons_D2  = 0.5,
        autoscale = TRUE
    ),
    class = "b"
)

fit1 <- brm(
    formula = y ~ -1 + .,
    data    = dat1,
    family  = bernoulli(link = "logit"),
    prior   = priors,
    chains  = 4,
    cores   = 4,
    iter    = 2000,
    seed    = 1
)
beta1 <- fixef(fit1)[,1]
fit2 <- brm(
    formula = y ~ -1 + .,
    data    = dat2,
    family  = bernoulli(link = "logit"),
    prior   = priors,
    chains  = 4,
    cores   = 4,
    iter    = 2000,
    seed    = 1
)
beta2 <- fixef(fit2)[,1]
fit3 <- brm(
    formula = y ~ -1 + .,
    data    = dat3,
    family  = bernoulli(link = "logit"),
    prior   = priors,
    chains  = 4,
    cores   = 4,
    iter    = 2000,
    seed    = 1
)

beta3 <- fixef(fit3)[,1]

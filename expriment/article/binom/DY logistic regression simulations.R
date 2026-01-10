

############################################################################

#	The script runs the logistic regression simulations of Sur and Candes

############################################################################


library(NPBayes)


#load("Logistic regression simulation.Rdata")


#  Generate X matrix and coefficient vectors


n	<- 4000
K	<- 800

set.seed(1)
X.mat	<- array(rnorm(n*K,mean = 0, sd = sqrt(1/n)),dim=c(n,K))

beta.1	<- c(rep(-10,100),rep(10,100),rep(0,600))
beta.2	<- rnorm(K,mean = 3, sd = 4)
beta.3	<- rbinom(K,1,0.5)*rnorm(K,mean = 7, sd = 1)



############################################################################

#	Simulation 1 -- logistic regression with the beta.1 coefficient

############################################################################


# 1.1 Initial simulation run to generate pi estimation plot and shrinkage plot -- takes about 12 minutes


set.seed(1)


sim1.mse.mat  <- array(dim = c(5,30))
dimnames(sim1.mse.mat) <- list(c("LSE","oracle","hBeta","lasso","ridge"),paste("Run",1:30))

set.seed(1)

for(i in 1:30)
{

    print(i)
    y.1		  <- rbinom(n,1,inv.log.odds(X.mat %*% beta.1))

    fit.1	    <- glm(y.1 ~ X.mat + 0,family = binomial)
    beta1.lse <- fit.1$coeff

    a.min <- min(c(-24,beta1.lse - 0.5))
    a.max <- max(c(24,beta1.lse + 0.5))
    hBeta.1 <- gibbs.logistic(500,y.1,X.mat,c(a.min,a.max),6,as.vector(beta1.lse),cpp=T)
    beta1.hBeta   <- apply(hBeta.1$beta.gibbs[101:500,],2,mean)

#    oracle.1 <- joint.EOD.gbbs.sampler(X.mat,sort(beta.1)[rank(beta1.lse)],y.1,50)
#    beta1.oracle  <- apply(oracle.1[1:50,],2,mean)

    fit.1.glmnet	<- cv.glmnet(X.mat,y.1,alpha = 1,intercept = FALSE, family = "binomial",type.measure = "class")
    beta1.lasso	<- as.numeric(coef(fit.1.glmnet, s = "lambda.min"))[2:801]

    fit.1.glmnet	<- cv.glmnet(X.mat,y.1,alpha = 0,intercept = FALSE, family = "binomial",type.measure = "class")
    beta1.ridge	<- as.numeric(coef(fit.1.glmnet, s = "lambda.min"))[2:801]

    sim1.mse.mat[1,i] <- mean((beta1.lse - beta.1)^2)
    sim1.mse.mat[2,i] <- mean((beta1.oracle - beta.1)^2)
    sim1.mse.mat[3,i] <- mean((beta1.hBeta - beta.1)^2)
    sim1.mse.mat[4,i] <- mean((beta1.lasso - beta.1)^2)
    sim1.mse.mat[5,i] <- mean((beta1.ridge - beta.1)^2)
}
date()


############################################################################

#	Simulation 2 -- logistic regression with the beta.2 coefficient

############################################################################

# 2.1 Initial simulation run to generate pi estimation plot and shrinkage plot -- takes about 12 minutes


##################################################################################

# 2.2  Perform 30 simulation runs and record MSE for the 5 estimates of beta

sim2.mse.mat  <- array(dim = c(5,30))
dimnames(sim2.mse.mat) <- list(c("LSE","oracle","hBeta","lasso","ridge"),paste("Run",1:30))

set.seed(100)

date()
for(i in 1:30)
{

  print(i)
  y.2		  <- rbinom(n,1,inv.log.odds(X.mat %*% beta.2))

  fit.2	    <- glm(y.2 ~ X.mat + 0,family = binomial)
  beta2.lse <- fit.2$coeff

  a.min <- min(c(-24,beta2.lse - 0.5))
  a.max <- max(c(24,beta2.lse + 0.5))
  hBeta.2 <- gibbs.logistic(500,y.2,X.mat,c(a.min,a.max),6,as.vector(beta2.lse),cpp=T)
  beta2.hBeta   <- apply(hBeta.2$beta.gibbs[101:500,],2,mean)

  oracle.2 <- joint.EOD.gbbs.sampler(X.mat,sort(beta.2)[rank(beta2.lse)],y.2,100)
  beta2.oracle <- apply(oracle.2[1:100,],2,mean)

  fit.2.glmnet	<- cv.glmnet(X.mat,y.2,alpha = 1,intercept = FALSE, family = "binomial",type.measure = "class")
  beta2.lasso	<- as.numeric(coef(fit.2.glmnet, s = "lambda.min"))[2:801]

  fit.2.glmnet	<- cv.glmnet(X.mat,y.2,alpha = 0,intercept = FALSE, family = "binomial",type.measure = "class")
  beta2.ridge	<- as.numeric(coef(fit.2.glmnet, s = "lambda.min"))[2:801]

  sim2.mse.mat[1,i] <- mean((beta2.lse - beta.2)^2)
  sim2.mse.mat[2,i] <- mean((beta2.oracle - beta.2)^2)
  sim2.mse.mat[3,i] <- mean((beta2.hBeta - beta.2)^2)
  sim2.mse.mat[4,i] <- mean((beta2.lasso - beta.2)^2)
  sim2.mse.mat[5,i] <- mean((beta2.ridge - beta.2)^2)
}


############################################################################

#	Simulation 3 -- logistic regression with the beta.3 coefficient

############################################################################


# 3.1 Initial simulation run to generate pi estimation plot and shrinkage plot -- takes about 12 minutes


set.seed(1)

date()
sim3.mse.mat  <- array(dim = c(5,30))
dimnames(sim3.mse.mat) <- list(c("LSE","oracle","hBeta","lasso","ridge"),paste("Run",1:30))

set.seed(100)

for(i in 1:30)
{

  print(i)
  y.3		  <- rbinom(n,1,inv.log.odds(X.mat %*% beta.3))

  fit.3	    <- glm(y.3 ~ X.mat + 0,family = binomial)
  beta3.lse <- fit.3$coeff

  a.min <- min(c(-24,beta3.lse - 0.5))
  a.max <- max(c(24,beta3.lse + 0.5))
  hBeta.3 <- gibbs.logistic(500,y.3,X.mat,c(a.min,a.max),6,as.vector(beta3.lse),cpp=T)
  beta3.hBeta   <- apply(hBeta.3$beta.gibbs[101:500,],2,mean)

#  oracle.3 <- joint.EOD.gbbs.sampler(X.mat,beta.3[order(beta3.lse)],y.3,100)
#  beta3.oracle  <- apply(oracle.3[51:100,],2,mean)

  fit.3.glmnet	<- cv.glmnet(X.mat,y.3,alpha = 1,intercept = FALSE, family = "binomial",type.measure = "class")
  beta3.lasso	<- as.numeric(coef(fit.3.glmnet, s = "lambda.min"))[2:801]

  fit.3.glmnet	<- cv.glmnet(X.mat,y.3,alpha = 0,intercept = FALSE, family = "binomial",type.measure = "class")
  beta3.ridge	<- as.numeric(coef(fit.3.glmnet, s = "lambda.min"))[2:801]

  sim3.mse.mat[1,i] <- mean((beta3.lse - beta.3)^2)
  sim3.mse.mat[2,i] <- mean((beta3.oracle - beta.3)^2)
  sim3.mse.mat[3,i] <- mean((beta3.hBeta - beta.3)^2)
  sim3.mse.mat[4,i] <- mean((beta3.lasso - beta.3)^2)
  sim3.mse.mat[5,i] <- mean((beta3.ridge - beta.3)^2)
}






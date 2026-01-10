##
# power simulation for mixed effect model
# D: 2022-02-15
# D: 2023-11-10 adding methods from https://arxiv.org/pdf/2208.10910.pdf
# D: 2025-03-04 oracle version
##
#rm(list=ls())
plot.fig = F
save.fig = F
save.data = T
graphics.off()
library(NPBayes)
source("qtl_util.R")
library(glmnet)
library(latex2exp)
library(zoo)
library(doParallel)
library(purrr)
library(PolyMixed)
library(ggplot2)
set.seed(12)
n.gibbs <-  20000 #10000 #samples for MCMC

if(plot.fig){
    q     <- 1  #we only generate a single sample to plot data.
    n_cores = 1
}else{
    q     <- 200
    n_cores = detectCores() - 6
}
distribution = c(2)  # 1 - gaussian random effect, 2- Laplace randeom effect, 3 -sparse
shift <- -2. #-0.02 # shift in the laplace distribution
nu <- 0.75 #assymetry in Laplace distribution

n_cores <- min(q, n_cores)
#thin =  how often to observe X (i.e. cm distnace)
#type =  1 - pure measurement error (\sigma^2I )
#        2 -  mixed effect     (\sigma^2I + XX^T \tau )


thin    = 1
active_propotion = 0.2 # if distribution 3- then how many non-zero signals
# parameters----------------------------------------------------------------

chr <- 10                 # chromosomes
m <-  150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
ns <-  400            # observations
d <- 1                   # distance between markers (cM)
sigma <- 0.1               # sd for error term
tau <- 0.005              # sd for polygenic effect
mu <- -rep(0.01, q)      # mean for polygenic effect




nc = NULL
qtl.pos <- ceiling(M * c(0.2,0.5,0.8))
mag.pos = 0.2
mag.neg = 0.2
qtl  <- c(mag.pos, qtl.pos[1], 1) # qtl (beta; position: marker, trait or traits)
qtl2 <- c(mag.pos, qtl.pos[2], 1)
qtl3 <- c(-mag.neg, qtl.pos[3], 1)
beta.true <- c(qtl[1], qtl2[1], qtl3[1])

# --------------------------------------------------------------------------

t.mix.forward.known <- matrix(0, q, M)
t.mix.forward <- matrix(0, q, M)
sigma.est.forward <- numeric(q)
tau.est.forward <- numeric(q)
crit.mix.forward <- numeric(q)
ignoreTau <- numeric(q) # LRT test


for(n in ns){
    string_out <- paste('\n\n',n)
    Xs <- Xs_new <- list()
    ys <- ys_new <- ys_new2 <- list()
    betas <- list()

    dupls <- list()
    SVDXs <- list()
    for(i in 1:q) {
        # Should we allow for duplitacte?
        X <- simulateDesignMatrix(chr, m, n, d)
        X_new <- simulateDesignMatrix(chr, m, n, d)

        Z <- rnorm(dim(X)[2])
        beta_          <- rep(0, dim(X)[2])
        if(distribution==1){
            beta_    = tau * Z
            beta_[qtl.pos] <- beta.true
        }else if(distribution==2){
            V <- rgamma(dim(X)[2],nu,nu)
            beta_ <- mu[1] + tau * (shift * (V-1) +  sqrt(V)*Z)
            beta_[qtl.pos] <- beta.true

        }else if(distribution == 3){
            index    = sample(1:dim(X)[2],ceiling(dim(X)[2]*active_propotion))
            V        = rep(0, dim(X)[2])
            V[index] = 1
            beta_    = V*mu[1]/active_propotion + tau/active_propotion * sqrt(V)*Z
            beta_[qtl.pos] <- beta.true
        }

        y <- X%*%(beta_) + sigma * rnorm(dim(X)[1])
        y_new <- X%*%(beta_) + sigma * rnorm(dim(X)[1])
        y_new_X_new <- X_new%*%(beta_) + sigma * rnorm(dim(X)[1])
        #y <- simulateTraits(X, q = 1, mu, tau, qtl = 0, qtl2 = 0, qtl3 = 0)

        thin_index   <- seq(1,dim(X)[2],by = thin)
        X            <- X[,thin_index]
        X_new        <- X_new[,thin_index]
        Xs[[i]]      <- X
        Xs_new[[i]]  <- X_new
        ys[[i]]      <- y
        ys_new[[i]]  <- y_new
        ys_new2[[i]] <- y_new_X_new
        betas[[i]]   <- beta_


    }

    cl <- makeCluster(n_cores, outfile="")
    registerDoParallel(cl)
    i0 = 1
    par.out <- foreach(i = i0:q)%dopar%{

        library(NPBayes)
        library(glmnet)
        #ridge regression
        cvfit <- cv.glmnet(Xs[[i]], ys[[i]], alpha=0)
        coef.ridge <- coef(cvfit, s = "lambda.1se")
        betas.ridge_unadj <- coef.ridge[-1]
        intercept.ridge_unadj <- coef.ridge[1]

        res <- gibbs.normal.permute.beta(n.gibbs,
                                        ys[[i]],
                                        Xs[[i]],
                                        betas[[i]],
                                        10)
        beta.oracle <-  base::colMeans(res$beta.gibbs[floor(n.gibbs/3):n.gibbs,])
        list(
            beta.oracle = beta.oracle
              )
    }
    list2env(purrr::transpose(par.out),globalenv())
    stopCluster(cl)
}


    if(plot.fig){
        window.size= 10
        res.beta.oracle <- plot.graph(beta.oracle[[q]], window.size, markers, beta_true = beta_,name="oracle")

        if(save.fig){
            ggsave('oracle.jpeg',res.beta.oracle$fig,width = 8, height = 6)
        }
    }


    ###
    # computing predictive power
    #
    ###
    MSE_XB_oracle<-  rep(0,q)
    MSE_beta_oracle  <- rep(0,q)
    beta.qtl.oracle <- matrix(0,nrow=q, ncol = length(qtl.pos))
    beta.smooth.qtl.oracle  <- matrix(0,nrow=q, ncol = length(qtl.pos))
    MSE_beta_smooth.oracle <- rep(0,q)

    for(i in 1:q){
        XB <- Xs[[i]]%*%betas[[i]]
        res.true <- smooth.beta(betas[[i]], betas[[i]], markers, window.size=10)

        MSE_XB_oracle[i] <- sqrt(mean((XB- Xs[[i]]%*%  beta.oracle[[i]] )^2))/sd(XB)
        MSE_beta_oracle[i] <- sqrt(mean((betas[[i]]  -beta.oracle[[i]] )^2))/sd(betas[[i]])

        res.oracle <- smooth.beta(beta.oracle[[i]], betas[[i]], markers, window.size=10)
        beta.smooth.qtl.oracle[i,] <- res.oracle$x[qtl.pos]
        MSE_beta_smooth.oracle[i] <- res.oracle$RMSE


    }

    Table <- cbind(c(mean(MSE_XB_oracle),mean(MSE_beta_oracle),mean(MSE_beta_smooth.oracle))
                   )


                   rownames(Table) <- c("RMSE of $ {\\bf X} \\beta $","RMSE of $ \\beta $","RMSE of $ \\tilde{\\beta} $")
    colnames(Table) <- c("oracle")

    library(xtable)
    print.xtable(xtable(Table, digits=3),type="latex",sanitize.text.function = function(x) x)


    sqrt(rowMeans(apply(beta.smooth.qtl.oracle, 1,function(x){(x-betas[[1]][qtl.pos])^2})))
    if(save.data){
    saveRDS(list(
        MSE_beta_oracle = MSE_beta_oracle,
        MSE_XB_oracle               = MSE_XB_oracle,
        MSE_beta_smooth.oracle                 = MSE_beta_smooth.oracle)
            ,"MSE_sim3.rds")}

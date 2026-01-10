##
# power simulation for mixed effect model
# D: 2022-02-15
# D: 2023-11-10 adding methods from https://arxiv.org/pdf/2208.10910.pdf
# D: 2024-01-24 updated format
##

#rm(list=ls())
plot.fig = T
save.fig = T
graphics.off()
source("qtl_util.R")
library(NPBayes)
library(zoo)
library(doParallel)
library(purrr)
library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(bigstep)
set.seed(12)

n.gibbs <-  20000 #10000 #samples for MCMC
if(plot.fig){
    q     <- 1  #we only generate a single sample to plot data.
    n_cores = 1
}else{
    q     <- 200
    n_cores = detectCores() - 2
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

        library(glmnet)
        library(SSLASSO)
        library(NPBayes)
        library(PolyMixed)
        cat('i=',i,'\n')
        set.seed(i)
        #unadjusted lasso
        cvfit <- cv.glmnet(Xs[[i]], ys[[i]])
        coef.lasso <- coef(cvfit, s = "lambda.1se")
        index.lasso <- which(coef.lasso!=0)
        index.lasso <- index.lasso[-1]
        lasso.beta <- rep(0,dim(X)[2])
        lasso.beta[index.lasso-1] <- coef.lasso[index.lasso]
        betas.lasso_unadj <- lasso.beta
        intercept.lasso_unadj <-coef.lasso[1]


        # Lasso
        cvfit <- cv.glmnet(cbind(rowSums(Xs[[i]]),Xs[[i]]), ys[[i]])
        coef.lasso <- coef(cvfit, s = "lambda.1se")
        index.lasso <- which(coef.lasso!=0)
        index.lasso <- index.lasso[-(1:2)]
        lasso.beta <- rep( coef.lasso[2],dim(X)[2])
        lasso.beta[index.lasso-2] <- lasso.beta[index.lasso-2]  + coef.lasso[index.lasso]
        betas.lasso <- lasso.beta
        intercept.lasso <-coef.lasso[1]

        #Spike and slab lasso
        sslasso.fit <- SSLASSO(cbind(rowSums(Xs[[i]]),Xs[[i]]), ys[[i]], variance = "unknown")
        index <- setdiff(sslasso.fit$model,1)
        sslasso.beta <- rep( sslasso.fit$beta[1, ncol(sslasso.fit$beta)], dim(X)[2])
        sslasso.beta[index-1] <- sslasso.beta[index-1] + sslasso.fit$beta[setdiff(sslasso.fit$model,1),ncol(sslasso.fit$beta)]
        betas.SSlasso <- sslasso.beta
        intercept.sslasso <- sslasso.fit$intercept[1,dim(sslasso.fit$intercept)[2]]

        sslasso.fit <- SSLASSO(Xs[[i]], ys[[i]], variance = "unknown")
        index <- sslasso.fit$model
        sslasso.beta <- rep( 0, dim(X)[2])
        sslasso.beta[index-1] <- sslasso.fit$beta[setdiff(sslasso.fit$model,1),ncol(sslasso.fit$beta)]
        betas.SSlasso_unadj <- sslasso.beta
        intercept.sslasso_unadj <- sslasso.fit$intercept[1,dim(sslasso.fit$intercept)[2]]


        #ridge regression
        cvfit <- cv.glmnet(cbind(rowSums(Xs[[i]]),Xs[[i]]), ys[[i]], alpha=0)
        coef.ridge <- coef(cvfit, s = "lambda.1se")
        betas.ridge <- coef.ridge[-(1:2)] + coef.ridge[2]
        intercept.ridge <- coef.ridge[1]


        cvfit <- cv.glmnet(Xs[[i]], ys[[i]], alpha=0)
        coef.ridge <- coef(cvfit, s = "lambda.1se")
        betas.ridge_unadj <- coef.ridge[-1]
        intercept.ridge_unadj <- coef.ridge[1]
        #hBayes
        res <- gibbs.normal(n.gibbs,
                            ys[[i]],
                            Xs[[i]],
                            c(-0.5,0.5),
                            8,
                            betas.ridge_unadj)
        #res <- gibbs.normal.fixed.sigma(n.gibbs,
        #                    ys[[i]],
        #                    Xs[[i]],
        #                    sigma,
        #                    c(-0.5,0.5),
        #                    8,
        #                    sslasso.fit$beta[,ncol(sslasso.fit$beta)])

        beta.npbayes <-  base::colMeans(res$beta.gibbs[floor(n.gibbs/3):n.gibbs,])



        #forward backward
        MixGeneObj <- SetupGeneMix('Y ~ 1',
                                   data = data.frame(Y=ys[[i]]),
                                   X=Xs[[i]],
                                   meanMuOff = F,
                                   tauOff = F)

        MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                                markers,
                                                dupl  = findDuplicate(Xs[[i]]))
        #compute posterior estimate of beta
        Q_2 = (t(Xs[[i]])%*%Xs[[i]])/MixGeneObj$sigma^2
        diag(Q_2) <- diag(Q_2)  + 1/MixGeneObj$tau^2
        Ebeta_2 = solve(Q_2, (t(Xs[[i]])%*%(ys[[i]] - MixGeneObj$beta[1] - Xs[[i]][,MixGeneObj$find,drop=F]%*%MixGeneObj$beta[-c(1:2)]))/MixGeneObj$sigma^2 + MixGeneObj$beta[2]/ MixGeneObj$tau^2 )
        Ebeta_2[MixGeneObj$find] =Ebeta_2[MixGeneObj$find] +MixGeneObj$beta[-c(1:2)]

        beta.forwardback <- Ebeta_2
        intercept.forwardback <-  MixGeneObj$beta[1]



        list( betas.lasso = betas.lasso,
              intercept.lasso =intercept.lasso,

              betas.lasso_unadj = betas.lasso_unadj,
              intercept.lasso_unadj = intercept.lasso_unadj,

              betas.SSlasso = betas.SSlasso,
              intercept.sslasso = intercept.sslasso,

              betas.SSlasso_unadj = betas.SSlasso_unadj,
              intercept.sslasso_unadj = intercept.sslasso_unadj,

              betas.ridge = betas.ridge,
              intercept.ridge = intercept.ridge,


              betas.ridge_unadj = betas.ridge_unadj,
              intercept.ridge_unadj = intercept.ridge_unadj,

              beta.npbayes = beta.npbayes,

              beta.forwardback=  beta.forwardback,
              intercept.forwardback = intercept.forwardback
              )
    }
    list2env(purrr::transpose(par.out),globalenv())
    stopCluster(cl)

}


if(plot.fig){
    window.size= 10
    #' @param x.sim simulations to get posterior conf int

    res.lasso       <- plot.graph(betas.lasso_unadj[[q]], window.size, markers, beta_true = beta_,name="Lasso")
    res.sslasso     <- plot.graph(betas.SSlasso_unadj[[q]], window.size, markers, beta_true = beta_,name="SSLasso")
    res.ridge       <- plot.graph(betas.ridge_unadj[[q]], window.size, markers, beta_true = beta_,name="Ridge")
    res.npbayes     <- plot.graph(beta.npbayes[[q]], window.size, markers, beta_true = beta_,name="HBayes")
    res.forwardback <- plot.graph(beta.forwardback[[q]], window.size, markers, beta_true = beta_,name="Mixed eff")

    if(save.fig){
        ggsave('lasso.jpeg',res.lasso$fig,width = 8, height = 6)
        ggsave('sslasso.jpeg',res.sslasso$fig,width = 8, height = 6)
        ggsave('ridge.jpeg',res.ridge$fig,width = 8, height = 6)
        ggsave('bayes.jpeg',res.npbayes$fig,width = 8, height = 6)
        ggsave('fb.jpeg',res.forwardback$fig,width = 8, height = 6)
    }

}else{


    ###
    # computing predictive power
    #
    ###

    smooth.beta <- function(x, beta_true, markers, window.size=10){
        x.window <- x
        beta.smooth <- beta_true
        marker.tot <- c(0, cumsum(markers))
        for(i in 2:length(marker.tot)){
            x.window[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(x[(1+marker.tot[i-1]):marker.tot[i]],
                                                                  window.size,
                                                                  na.pad=T,
                                                                  fill="extend")
            beta.smooth[(1+marker.tot[i-1]):marker.tot[i]] <- zoo::rollsum(beta_true[(1+marker.tot[i-1]):marker.tot[i]],
                                                                     window.size,
                                                                     na.pad=T,
                                                                     fill="extend")
        }
        return(list(x=x.window, RMSE = sqrt(mean( (x.window-beta.smooth)^2 ))/sd(beta.smooth)))
    }
    MSE_XB_lasso <- MSE_XB_SS <- MSE_XB_ridge <- MSE_XB_fb <- MSE_XB_npbayes <- rep(0,q)
    MSE_beta_lasso <- MSE_beta_SS <- MSE_beta_ridge <-MSE_beta_fb <- MSE_beta_npbayes <- rep(0,q)
    beta.qtl.ridge <-  beta.qtl.lasso <- beta.qtl.fb <- beta.qtl.bayes <- matrix(0,nrow=q, ncol = length(qtl.pos))
    beta.smooth.qtl.ridge <- beta.smooth.qtl.true <-  beta.smooth.qtl.lasso <- beta.smooth.qtl.fb <- beta.smooth.qtl.bayes <- matrix(0,nrow=q, ncol = length(qtl.pos))
    MSE_beta_smooth.ridge <- MSE_beta_smooth.lasso <- MSE_beta_smooth.fb <-MSE_beta_smooth.bayes <- rep(0,q)

    for(i in 1:q){
        XB <- Xs[[i]]%*%betas[[i]]
        MSE_XB_lasso[i] <- sqrt(mean((XB- Xs[[i]]%*%betas.lasso[[i]] - intercept.lasso[[i]])^2))/sd(XB)
        MSE_beta_lasso[i] <- sqrt(mean((betas[[i]]  -betas.lasso[[i]] )^2))/sd(betas[[i]])

        MSE_XB_SS[i] <- sqrt(mean((XB- Xs[[i]]%*%betas.SSlasso[[i]] - intercept.sslasso[[i]])^2))/sd(XB)
        MSE_beta_SS[i] <- sqrt(mean((betas[[i]]  -betas.SSlasso[[i]]  )^2))/sd(betas[[i]])


        MSE_XB_ridge[i] <- mean((XB- Xs[[i]]%*%betas.ridge_unadj[[i]] - intercept.ridge_unadj[[i]])^2)/sd(XB)
        MSE_beta_ridge[i] <- sqrt(mean((betas[[i]]  -betas.ridge_unadj[[i]] )^2))/sd(betas[[i]])

        MSE_XB_npbayes[i] <- sqrt(mean((XB- Xs[[i]]%*%  beta.npbayes[[i]] )^2))/sd(XB)
        MSE_beta_npbayes[i] <- sqrt(mean((betas[[i]]  -beta.npbayes[[i]] )^2))/sd(betas[[i]])


        MSE_XB_fb[i] <- sqrt(mean((XB - Xs[[i]]%*%beta.forwardback[[i]] - intercept.forwardback[[i]])^2))/sd(XB)
        MSE_beta_fb[i] <- sqrt(mean((betas[[i]]  -beta.forwardback[[i]] )^2))/sd(betas[[i]])


        beta.qtl.ridge[i, ] <- betas.ridge_unadj[[i]][qtl.pos]
        beta.qtl.lasso[i, ] <- betas.lasso[[i]][qtl.pos]
        beta.qtl.fb[i, ]    <- beta.forwardback[[i]][qtl.pos]
        beta.qtl.bayes[i, ] <- beta.npbayes[[i]][qtl.pos]

        res.ridge <- smooth.beta(betas.ridge_unadj[[i]], betas[[i]], markers, window.size=10)
        beta.smooth.qtl.ridge[i,] <- res.ridge$x[qtl.pos]
        MSE_beta_smooth.ridge[i] <- res.ridge$RMSE

        res.lasso <- smooth.beta(betas.lasso[[i]], betas[[i]], markers, window.size=10)
        beta.smooth.qtl.lasso[i,] <- res.lasso$x[qtl.pos]
        MSE_beta_smooth.lasso[i] <- res.lasso$RMSE

        res.fb <- smooth.beta(beta.forwardback[[i]], betas[[i]], markers, window.size=10)
        beta.smooth.qtl.fb[i,] <- res.fb$x[qtl.pos]
        MSE_beta_smooth.fb[i] <- res.fb$RMSE

        res.bayes <- smooth.beta(beta.npbayes[[i]], betas[[i]], markers, window.size=10)
        beta.smooth.qtl.bayes[i,] <- res.bayes$x[qtl.pos]
        MSE_beta_smooth.bayes[i] <- res.bayes$RMSE

        res.true <- smooth.beta(betas[[i]], betas[[i]], markers, window.size=10)
        beta.smooth.qtl.true[i, ] <- res.true$x[qtl.pos]

    }

    Table <- rbind(c(mean(MSE_XB_lasso),mean(MSE_XB_ridge),mean(MSE_XB_fb),mean(MSE_XB_npbayes)),
                   c(mean(MSE_beta_lasso),mean(MSE_beta_ridge),mean(MSE_beta_fb),mean(MSE_beta_npbayes)),
                   c(mean(MSE_beta_smooth.lasso),mean(MSE_beta_smooth.ridge),mean(MSE_beta_smooth.bayes),mean(MSE_beta_smooth.fb)))
    rownames(Table) <- c("RMSE of $ {\\bf X} \\beta $","RMSE of $ \\beta $","RMSE of $ \\tilde{\\beta} $")
    colnames(Table) <- c("lasso",  "ridge", "mix","hBeta")

    library(xtable)
    print.xtable(xtable(Table, digits=3),type="latex",sanitize.text.function = function(x) x)


    saveRDS(list(MSE_XB_lasso = MSE_XB_lasso,
                 MSE_XB_SS    = MSE_XB_SS,
                 MSE_XB_ridge = MSE_XB_ridge,
                 MSE_XB_npbayes = MSE_XB_npbayes,
                 MSE_XB_fb      = MSE_XB_fb,
                 MSE_beta_lasso = MSE_beta_lasso,
                 MSE_beta_SS    = MSE_beta_SS,
                 MSE_beta_ridge = MSE_beta_ridge,
                 MSE_beta_npbayes = MSE_beta_npbayes,
                 MSE_beta_fb    = MSE_beta_fb,
                 beta.qtl.ridge = beta.qtl.ridge,
                 beta.qtl.lasso = beta.qtl.lasso,
                 beta.qtl.fb  = beta.qtl.fb,
                 beta.qtl.bayes = beta.qtl.bayes,
                 beta.smooth.qtl.true = beta.smooth.qtl.true,
                 beta.smooth.qtl.bayes = beta.smooth.qtl.bayes,
                 beta.smooth.qtl.fb = beta.smooth.qtl.fb,
                 beta.smooth.qtl.lasso = beta.smooth.qtl.lasso,
                 beta.smooth.qtl.ridge = beta.smooth.qtl.ridge,
                 MSE_beta_smooth.ridge = MSE_beta_smooth.ridge,
                 MSE_beta_smooth.lasso = MSE_beta_smooth.lasso,
                 MSE_beta_smooth.fb = MSE_beta_smooth.fb,
                 MSE_beta_smooth.bayes = MSE_beta_smooth.bayes)
                ,"MSE_sim.rds")
}

#'
#' Haar measure on O(2) or Haar measure on the two last dimension
#' @param rho first parameter
#' @param b   second parameter
#' @param m   dimension of output matrix (if zero m=2)
#' @return sample matrix
R.rho.b	<- function(rho,b,m = NULL)
{
    if(is.null(m)) gamma.tmp <- cbind(c(cos(rho),-b * sin(rho)),c(sin(rho),b * cos(rho)))

    if(!is.null(m))
    {
        gamma.tmp <- diag(rep(1,m))
        gamma.tmp[(m-1):m,(m-1):m] <- cbind(c(cos(rho),-b * sin(rho)),c(sin(rho),b * cos(rho)))
    }
    return(gamma.tmp)
}


#' spherical standard normal N(vm, sm)
#' @param  p - dimension
#' @param  vm - center to sample around
#' @param  sm - standard dev
#' @return Z - return "normal random variable" of length 1
smp.v	<- function(p = NULL, vm = NULL, sm = 1)
{
    if(!is.null(vm))
        p = length(vm)
    rz	<- sm * rnorm(p,mean = 0, sd = 1)
    if(!is.null(vm))
        rz <- rz + vm
    return(rz / sqrt(sum(rz^2)))
}

#'
#' Building the Householder reflection
#' @param v.p
#' @param m
#'
house.reflctn	<- function(v.p,m = NULL)
{
    p	<- length(v.p)
    if(is.null(m))
    {
        e1	<- c(1,rep(0,p-1))
        v.m <- v.p
        m   <- p
    }
    if(!is.null(m))
    {
        e1	<- c(rep(0,m-p),1,rep(0,p-1))
        v.m <- c(rep(0,m-p),v.p)
    }
    den <- sqrt(sum((e1 - v.m)^2))
    if(den==0){
        x.m <- (e1 - v.m)
    }else{
        x.m	<- (e1 - v.m) / sqrt(sum((e1 - v.m)^2))
    }

    E = - 2* x.m %*% t(x.m)
    diag(E) <- diag(E) + 1
    return( E)
}


#'
#' multiply Gamma with the Householder reflection
#' @param Gamma
#' @param v.p
#'
house.reflctn.mult	<- function(Gamma,
                               v.p,
                               right = TRUE)
{
    p <- dim(Gamma)[1]
    m <- length(v.p)
    v.p[1]<- v.p[1]-1
    x.m	<- sqrt(2) * -v.p / sqrt(sum((v.p)^2))
    ind <- (1+p-m):p
    if(right){
        Gamma[ind,] <- Gamma[ind,] - x.m%*%(t(x.m)%*%(Gamma[ind,]))
    }else{
        Gamma[,ind] <- Gamma[,ind] - (Gamma[,ind]%*%x.m)%*%t(x.m)
    }
    return( Gamma)
}


acos.2pi <- function(cos.rho,sin.rho)
{
    if(0 <= sin.rho) return(acos(cos.rho))
    if(0 > sin.rho) return(2*pi - acos(cos.rho))
}

#' compute log likelihood of multivariate normal
#' @param R - cholesky of XtX
#' @param E - eigenvectors of Sigma
#' @param D - eigenvalues of Sigma
#' @param n - number of observations
oracle.gamma.R.loglik <- function(R,E, D, n)
{
    V <- R%*%E #p^2 (complexity)
    p <- dim(V)[1]
    return(-0.5*( sum(colSums(V^2)/D) + n* sum(log(D)) + n*p*log(2*pi)))
}

oracle.gamma.x.loglik <- function(x.tmp,gamma.tmp, D)
{
    x.trns	    <- x.tmp %*% gamma.tmp
    p <- dim(x.tmp)[2]
    n <- dim(x.tmp)[1]
    return(sum(dnorm(x.trns,mean = rep(0,p*n),sd = rep(sqrt(D),each = n),log = T)))
}

eval.orth.mat.2D	<- function(orth.mat)	c(acos.2pi(orth.mat[1,1],orth.mat[1,2]), orth.mat[2,2] / orth.mat[1,1])

vp.repr.orth.mat <- function(orth.mat)
{
    m 		<- dim(orth.mat)[1]
    v.p.mat	<- array(dim=c(m,m))

    if(m == 2) dimnames(v.p.mat) <- list(NULL,c("rho","b"))

    if(2 < m)
    {
        dimnames(v.p.mat) <- list(NULL,c(paste("v.",m:3,sep=""),"rho","b"))
        for(i in 1:(m-2))
        {
            v.p.mat[1:(m-i+1),i]	<- 	orth.mat[,1]
            orth.mat					<- (house.reflctn(v.p.mat[1:(m-i+1),i]) %*% orth.mat)[2:(m-i+1),2:(m-i+1)]
        }
    }
    if(is.na(sum(eval.orth.mat.2D(orth.mat)))==T){
        print('error')
    }

    v.p.mat[1,c(m-1,m)] <- 	 eval.orth.mat.2D(orth.mat)
    return(v.p.mat)
}

smp.rnd.orth.mat	<- function(n = 10)
{

    P.tmp	<- R.rho.b(runif(1,0,2*pi),1-2*rbinom(1,1,0.5))

    if(2 < n)
        for(i in 3:n)
        {
            v.i		<- smp.v(i)
            hr		<- house.reflctn(v.i)
            P.tmp	<- hr %*% cbind(c(1,rep(0,i-1)),rbind(0,P.tmp))
        }

    return(P.tmp)
}


#####################################################

# 1. Oracle Metropolis sampler assumes given D.oracle

#####################################################

#'
#' Sampling the eigenvectors
#'
#' @param E.init     - (p x p ) inital guess of eigenvectors
#' @param X          - (n x p) data assume centered
#' @param n.mtrp     - (int) number of metropolis steps
#' @param D          - (p x 1) eigenvalues
#' @param target.acc - (p x 1) target acceptance rate for the MH
#' @param MH.objs    - (list) adaptive object for MH
#' @param R          - (p xp) cholesky of XtX
oracle.metrop.sampler <- function(E.init, X, D, n.mtrp, target.acc = 0.23,
                                  MH.objs = NULL, R = NULL)
{
    n	<- dim(X)[1]
    p <- dim(X)[2]

    if(is.null(R)){
        XtX <- t(X)%*%X
        R <- chol(XtX)
    }



    v.p.ini <-  vp.repr.orth.mat(E.init)

    mtrp.loglik <- rep(NA,n.mtrp)
    mtrp.gamma <- array(dim=c(p,p,n.mtrp))
    acc.ind <- array(dim=c(n.mtrp-1,p))

    mtrp.gamma[,,1] <- E.init
    mtrp.loglik[1] <- oracle.gamma.R.loglik(R,mtrp.gamma[,,1], D, n)

    v.p.tmp <- v.p.ini

    loglik.old <- mtrp.loglik[1]


    if(is.null(MH.objs)){
        MH.objs <- list()
        for(i in 1:p){
            MH.objs[[i]] <- MH.init(0.3/(i+1), target.acc)
        }
    }

    for(i in 2:n.mtrp)
    {
        gamma.left <- diag(rep(1,p))
        gamma.right <- mtrp.gamma[,,i-1]

        # Run metropolis iterations for v.p vectors p to 3
        if(p>2){
            for(pp in p:3)
            {
                #h.v.p.m <- house.reflctn(v.p.tmp[1:pp,p-pp+1],m=p)
                #gamma.right_old <- h.v.p.m %*% gamma.right
                house_reflection_mult_inplace_cpp(gamma.right, v.p.tmp[1:pp,p-pp+1])
                v.prpsl <- smp.v(vm =v.p.tmp[1:pp,p-pp+1],
                                 sm = MH.objs[[pp]]$sigma)
                gamma.left.prop <- house_reflection_mult_cpp(gamma.left,v.prpsl, right =FALSE)
                gamma.prpsl <- gamma.left.prop %*% gamma.right
                prpsl.loglik <- oracle.gamma.R.loglik(R, gamma.prpsl, D, n)

                acc.ind[i-1,pp] <- runif(1) <= exp(prpsl.loglik - loglik.old)

                if(acc.ind[i-1,pp])
                {
                    v.p.tmp[1:pp,p-pp+1] <- v.prpsl
                    loglik.old <- prpsl.loglik
                    gamma.left <- gamma.left.prop
                }else{
                    #gamma.left <- gamma.left %*% h.v.p.m
                    house_reflection_mult_inplace_cpp(gamma.left, v.p.tmp[1:pp,p-pp+1], FALSE)
                }
                MH.objs[[pp]] <- MH.adaptive(acc.ind[i-1,pp], MH.objs[[pp]])

            }
            }

        # Sample metropolis sample for rho keeping b fixed as in MLE

        rho.prpsl <- rnorm(1,mean = v.p.tmp[1,p-1],sd = MH.objs[[2]]$sigma)
        R.rho.b.prop <- R.rho.b(rho.prpsl,v.p.tmp[1,p],p)
        gamma.prpsl <- gamma.left %*% R.rho.b.prop
        prpsl.loglik <- oracle.gamma.R.loglik(R,gamma.prpsl, D, n)

        acc.ind[i-1,2] <- runif(1) <= exp(prpsl.loglik - loglik.old)
        if(acc.ind[i-1,2])
        {
            v.p.tmp[1,p-1] <- rho.prpsl
            loglik.old <- prpsl.loglik
            mtrp.gamma[,,i] <- gamma.prpsl
        }
        if(!acc.ind[i-1,2])
            mtrp.gamma[,,i] <- gamma.left %*% R.rho.b(v.p.tmp[1,p-1],v.p.tmp[1,p],p)


        MH.objs[[2]] <- MH.adaptive(acc.ind[i-1,2], MH.objs[[2]])
        mtrp.loglik[i] <- loglik.old
    }

    return(list(Es = mtrp.gamma,
                loglik = mtrp.loglik,acc = acc.ind,
                MH.objs =MH.objs ))

}



#'
#' Sampling the eigenvectors
#'
#' @param E.init     - (p x p ) inital geuss of eigenvectors
#' @param X          - (n x p) data assume centered
#' @param n.mtrp     - (int) number of metropolis steps
#' @param D          - (p x 1) eigenvalues
#' @param target.acc - (p x 1) target acceptance rate for the MH
#' @param MH.objs    - (list) adaptive object for MH
oracle.metrop.sampler.old <- function(E.init, X, D, n.mtrp, target.acc = 0.23,
                                  MH.objs = NULL)
{
    n	<- dim(X)[1]
    p <- dim(X)[2]



    v.p.ini <-  vp.repr.orth.mat(E.init)

    mtrp.loglik <- rep(NA,n.mtrp)
    mtrp.gamma <- array(dim=c(p,p,n.mtrp))
    acc.ind <- array(dim=c(n.mtrp-1,p))

    mtrp.gamma[,,1] <- E.init
    mtrp.loglik[1] <- oracle.gamma.x.loglik(X,mtrp.gamma[,,1], D)

    v.p.tmp <- v.p.ini

    loglik.old <- mtrp.loglik[1]


    if(is.null(MH.objs)){
        MH.objs <- list()
        for(i in 1:p){
            MH.objs[[i]] <- MH.init(0.3/(i+1), target.acc)
        }
    }

    for(i in 2:n.mtrp)
    {
        gamma.left <- diag(rep(1,p))
        gamma.right <- mtrp.gamma[,,i-1]

        # Run metropolis iterations for v.p vectors p to 3

        for(pp in p:3)
        {
            h.v.p.m <- house.reflctn(v.p.tmp[1:pp,p-pp+1],m=p)
            gamma.right <- h.v.p.m %*% gamma.right
            v.prpsl <- smp.v(vm =v.p.tmp[1:pp,p-pp+1],
                             sm = MH.objs[[pp]]$sigma)
            h.v.prpsl <- house.reflctn(v.prpsl,m = p)
            gamma.left.prop <- gamma.left %*% h.v.prpsl
            gamma.prpsl <- gamma.left.prop %*% gamma.right
            prpsl.loglik <- oracle.gamma.x.loglik(X, gamma.prpsl, D)
            acc.ind[i-1,pp] <- runif(1) <= exp(prpsl.loglik - loglik.old)
            if(acc.ind[i-1,pp])
            {
                v.p.tmp[1:pp,p-pp+1] <- v.prpsl
                loglik.old <- prpsl.loglik
                gamma.left <- gamma.left.prop
            }else{
                gamma.left <- gamma.left %*% h.v.p.m
            }
            MH.objs[[pp]] <- MH.adaptive(acc.ind[i-1,pp], MH.objs[[pp]])

        }

        # Sample metropolis sample for rho keeping b fixed as in MLE

        rho.prpsl <- rnorm(1,mean = v.p.tmp[1,p-1],sd = MH.objs[[2]]$sigma)
        R.rho.b.prop <- R.rho.b(rho.prpsl,v.p.tmp[1,p],p)
        gamma.prpsl <- gamma.left %*% R.rho.b.prop
        prpsl.loglik <- oracle.gamma.x.loglik(X,gamma.prpsl, D)

        acc.ind[i-1,2] <- runif(1) <= exp(prpsl.loglik - loglik.old)

        if(acc.ind[i-1,2])
        {
            v.p.tmp[1,p-1] <- rho.prpsl
            loglik.old <- prpsl.loglik
            mtrp.gamma[,,i] <- gamma.prpsl
        }
        if(!acc.ind[i-1,2])
            mtrp.gamma[,,i] <- gamma.left %*% R.rho.b(v.p.tmp[1,p-1],v.p.tmp[1,p],p)


        MH.objs[[2]] <- MH.adaptive(acc.ind[i-1,2], MH.objs[[2]])
        mtrp.loglik[i] <- loglik.old
    }

    return(list(Es = mtrp.gamma,
                loglik = mtrp.loglik,acc = acc.ind,
                MH.objs =MH.objs ))

}

#'
#' Sampling the eigenvectors
#'
#' @param E.init - (p x p ) inital geuss of eigenvectors
#'@param X       - (n x p) data assume centered
#'@param n.mtrp  - (int) number of metropolis steps
#'@param D       - (p x 1) eigenvalues
oracle.metrop.sampler.old2 <- function(E.init, X, D, n.mtrp = 10)
{
    n	<- dim(X)[1]
    p <- dim(X)[2]



    v.p.ini <-  vp.repr.orth.mat(E.init)

    mtrp.loglik <- rep(NA,n.mtrp)
    mtrp.gamma <- array(dim=c(p,p,n.mtrp))
    acc.ind <- array(dim=c(n.mtrp-1,p))

    mtrp.gamma[,,1] <- E.init
    mtrp.loglik[1] <- oracle.gamma.x.loglik(X,mtrp.gamma[,,1], D)

    v.p.tmp <- v.p.ini

    loglik.old <- mtrp.loglik[1]

    for(i in 2:n.mtrp)
    {
        gamma.left <- diag(rep(1,p))
        gamma.right <- mtrp.gamma[,,i-1]

        # Run metropolis iterations for v.p vectors p to 3

        for(pp in p:3)
        {
            gamma.right <- house.reflctn(v.p.tmp[1:pp,p-pp+1],m=p) %*% gamma.right
            v.prpsl <- smp.v(vm =v.p.tmp[1:pp,p-pp+1],
                             sm = 0.2)
            gamma.prpsl <- gamma.left %*% house.reflctn(v.prpsl,m = p) %*% gamma.right
            prpsl.loglik <- oracle.gamma.x.loglik(X, gamma.prpsl, D)
            acc.ind[i-1,pp] <- runif(1) <= exp(prpsl.loglik - loglik.old)

            if(acc.ind[i-1,pp])
            {
                v.p.tmp[1:pp,p-pp+1] <- v.prpsl
                loglik.old <- prpsl.loglik
            }

            gamma.left <- gamma.left %*% house.reflctn(v.p.tmp[1:pp,p-pp+1],m=p)
        }

        # Sample metropolis sample for rho keeping b fixed as in MLE

        rho.prpsl <- rnorm(1,mean = v.p.tmp[1,p-1],sd = 0.2)
        gamma.prpsl <- gamma.left %*% R.rho.b(rho.prpsl,v.p.tmp[1,p],p)
        prpsl.loglik <- oracle.gamma.x.loglik(X,gamma.prpsl, D)

        acc.ind[i-1,2] <- runif(1) <= exp(prpsl.loglik - loglik.old)

        if(acc.ind[i-1,2])
        {
            v.p.tmp[1,p-1] <- rho.prpsl
            loglik.old <- prpsl.loglik
            mtrp.gamma[,,i] <- gamma.prpsl
        }
        if(!acc.ind[i-1,2])
            mtrp.gamma[,,i] <- gamma.left %*% R.rho.b(v.p.tmp[1,p-1],v.p.tmp[1,p],p)

        mtrp.loglik[i] <- loglik.old
    }

    return(list(Es = mtrp.gamma,loglik = mtrp.loglik,acc = acc.ind))

}



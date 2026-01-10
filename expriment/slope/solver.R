graphics.off()
iter= 500

#' Solver of for asymptotic distribution of generalized lasso
#' L(u) = \frac{1}{2} u^T C u - u^T w + f'(\beta_0;u)
#' where f' is directional derivative of f(\beta) = \lambda |A\beta|
#' Let A be m x p matrix, C pxp symmetric invertable, w a p x 1 vector
#' Let E_A be a diagonal m x m matrix with diagonal sign(A\beta) where
#' sign(0)=0
#' Then f'(\beta_0;u) = \lambda \sum_i (E_A A)_iu + \lambda|(I - abs(E_A)) Au|
Loss<-function(u, C, A, w, beta0, lambda ){
    E_1 <- diag(c(sign(A%*%beta0)))
    ind <- which(diag(E_1)==0)
    w_tilde <- w - lambda * colSums(E_1%*%A)

    F_0 <- A[ind,] #THIS is F in (6.3)
    return(0.5*t(u)%*%(C%*%u) - t(u)%*%w_tilde + lambda *sum(abs(F_0%*%u)))
}

#' admm p.34
soft.thersholding <- function(a, kappa){
    res <- rep(0,length(a))
    ind_pos <- a > kappa
    res[ind_pos] = a[ind_pos] - kappa
    ind_neg <- a < -kappa
    res[ind_neg] = a[ind_neg] + kappa
    return(res)
}
#run a sequence ofGLASS0
#' L(x) = \frac{1}{2} x^T C x - x^T w + f'(\beta_0;x)
#' @return sum(iter)+1 x p solution path of the algorithm
ADMM.GLASSO.outer <-function(C, A, w, beta0, lambda, rho = c(10,100,1000),x0 = NULL , iter=c(100,100,1000), proj.error=1e-3){

    x0 <- x0
    u0 <- NULL
    z0 <- NULL
    x <- NULL
    for( i in 1:length(iter) ){
        u.est <-ADMM.GLASSO(C, A, w, beta0, lambda,
                            x0 = x0,
                            z0 = z0,
                            u0 = u0,
                            rho = rho[i],
                            iter=iter[i])
        x0 <- u.est$x[iter[i],]
        u0 <- u.est$u[iter[i],]
        z0 <- u.est$z[iter[i],]
        if(is.null(x)){
            x <- u.est$x
        }else{
            x <- rbind(x,u.est$x)
        }
    }
    ind <- abs(A%*%x0)< proj.error
    if(sum(ind)>0){
        X <- A[ind,,drop=F]
        x_final <- c(x0 - t(X)%*%solve(X%*%t(X))%*%X%*%x0)
        x <- rbind(x,x_final)
    }
    return(x)

}

#' ADMM p.44 https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
#' L(x) = \frac{1}{2} x^T C x - x^T w + f'(\beta_0;x)
#' where f' is directional derivative of f(\beta) = \lambda |A\beta|
#' Let A be m x p matrix, C pxp symmetric invertable, w a p x 1 vector
#' Let E_A be a diagonal m x m matrix with diagonal sign(A\beta) where
#' sign(0)=0
#' Then f'(\beta_0;u) = \lambda \sum_i (E_A A)_iu + \lambda|(I - abs(E_A)) Au|
#' @return (x, u , z) by ADMM p.44
ADMM.GLASSO <- function(C, A, w, beta0, lambda, rho = 1.,x0 = NULL, u0=NULL, z0=NULL,  iter=100){

    if(is.null(x0))
        x0 <- rep(0, dim(C)[1])
    E_1 <- diag(c(sign(A%*%beta0)))
    ind <- which(diag(E_1)==0)
    w_tilde <- w - lambda * colSums(E_1%*%A)

    F_0 <- A[ind,,drop=F] #THIS is F in (6.3)

    #use A^TA = C (6.3)
    #    b = C^{-1/2}w_tilde
    # however only A^T b= w_tilde  is used in the algorithm
    Q <- solve(C + rho * t(F_0)%*%F_0)
    x_b <- Q%*%(w_tilde)
    QF <- rho * Q%*%t(F_0)
    x <-      matrix(0, nrow = iter,ncol=length(beta0))
    x[1,] <- x0
    z <- u <- matrix(0, nrow = iter,ncol=dim(F_0)[1])


    if(is.null(u0))
        u0 <- rep(0, dim(F_0)[1])
    if(is.null(z0))
        z0 <- rep(0, dim(F_0)[1])

    u[1,] <- u0
    z[1,] <- z0
    for(i in 2:iter){
        x[i,] <- x_b + QF%*%(z[i-1,] - u[i-1,])
        Fx <-F_0%*%x[i,]
        z[i,] <- soft.thersholding(Fx + u[i-1,], lambda/rho)
        u[i,] <-u[i-1,] + Fx - z[i,]
    }
    return(list(x=x, u=u, z = z))
}



set.seed(1)

beta0 <- c(0,0,0,1,1,1,2,2,2)
p <- length(beta0)
m <- p
sigma <- 1
lambda <- 4
A <- matrix(0, nrow= m, ncol = p)
w <- rnorm(p)
C <- matrix(0.8, p, p)
diag(C) <- 1
w <- sigma*t(chol(C))%*%rnorm(p)
for(i in 1:(m-1)){
    A[i,i] <- 1
    A[i,i+1] <- -1
}
A[m,1] <- 1
u <- beta0
x <- ADMM.GLASSO.outer(C, A, w, beta0, lambda, rho = c(10^2,10^3,10^4,10^5,10^6),x0 = NULL , iter=c(100,100,100,100,100))
Ls <- rep(0, dim(x)[1])
for(i in 1:dim(x)[1]){
    Ls[i] <- Loss(x[i,],C, A, w, beta0, lambda)
}
cat('final x = ',round(x[dim(x)[1],],2),'\n')
plot(Ls,xlab='iter',ylab='loss',type='l')



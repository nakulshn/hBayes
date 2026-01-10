


#' Gibbs sampler derived from Yekuteili code
#' then in the gibbs sampler we randomly permute the beta coeffients (one at the
#' time)
#' where \pi_q is a density estimated from data
#'
#' @paran n.gibbs - (int ) number of gibbs samples to be performed
#' @param Y       -  ( n x 1) observastions
#' @param X       -  ( n x p)  covariates vector
#' @param sigma   -  (double) standard devation
#' @param a.in    -  (2 x 1) lower and upper point of the density
#' @param L       -  (int)   2^L is number of bins in the tree
#' @param beta    -  (p x 1)  inital guess of beta
#' @param cpp     -  (bool) use C++
#'
#'@export
gibbs.normal.permute.beta.fixed.sigma          <- function(n.gibbs,
                                                               Y,
                                                               X,
                                                               sigma,
                                                               a.int,
                                                               L,
                                                               beta,
                                                               quite = T){



    p			<- dim(X)[2]
    n			<- dim(X)[1]
    tree.data <- build.tree(L)

    I = tree.data$I


    a.vec	<- seq(a.int[1],
                 a.int[2],
                 length = I + 1)

    #truncate beta
    beta[beta< a.int[1]] = a.int[1] + 1e-8
    beta[beta> a.int[2]] = a.int[2] - 1e-8


    pi.gibbs		<- array(dim = c(n.gibbs, I))
    beta.gibbs	    <- array(dim = c(n.gibbs, p))
    delta.gibbs    <-  array(dim = c(n.gibbs, p))



    beta.gibbs[1,]	<- beta
    beta.ind <- findInterval(beta, a.vec, left.open = T)

    pi.gibbs[1,]	<- log( (c(sum(beta.gibbs[1,] < a.vec[2]),diff(rank(c(a.vec[2:I],100,beta.gibbs[1,]))[1:I] -
                                                                    (1:I))) + 100/I) / 900)
    delta.gibbs[1, ] <- beta.ind
    Xbeta	<- X %*% beta.gibbs[1,]




    for (g in 2:n.gibbs)
    {
        if(quite==F){
            if (g %% 100 == 1)	print(paste(g))
        }


        ##
        # sample beta
        ##
        pi.gibbs_g       <-  pi.gibbs[g-1,]
        beta.gibbs.g     <- beta.gibbs[g-1,]
        res.beta <- permute_gibbs_beta_normal_cpp(Y,
                                           Xbeta,
                                           X,
                                         sigma,
                                          beta.gibbs.g,
                                         beta.ind)

        beta.ind <-         res.beta$beta.ind
        beta.gibbs[g,]   <- beta.gibbs.g

        delta.gibbs[g, ] <- beta.ind


        ##
        # sample pi
        ##
        nvec <- table(factor(beta.ind, levels=1:I))
        pi.gibbs[g,] <- polya.gibs(nvec, tree.data)
    }
    return(list(a.vec=  a.vec,
                pi.gibbs = pi.gibbs,
                beta.gibbs = beta.gibbs,
                delta.gibbs = delta.gibbs))

}

#' Gibbs Sampler with Beta Permutation and Sigma Sampling
#'
#' Performs Gibbs sampling for normal linear model with beta coefficient permutation
#' and variance parameter (sigma^2) sampling.
#'
#' @param n.gibbs Number of Gibbs samples (integer)
#' @param Y Response vector (n x 1)
#' @param X Design matrix (n x p)
#' @param a.int Interval boundaries for beta coefficients (vector length 2)
#' @param L Number of bins (2^L) for density estimation (integer)
#' @param beta Initial beta coefficients (vector length p)
#' @param cpp Use C++ implementation (logical, default = TRUE)
#' @param quiet Suppress iteration messages (logical, default = TRUE)
#'
#' @return List containing:
#' \itemize{
#'   \item a.vec - Bin boundaries
#'   \item pi.gibbs - Density estimates
#'   \item beta.gibbs - Sampled beta coefficients
#'   \item delta.gibbs - Bin indices
#'   \item sigma.gibbs - Sampled sigma values
#' }
#'
#' @export
gibbs.normal.permute.beta <- function(n.gibbs, Y, X, beta,gamma=1,quiet = TRUE) {

    # Initialize dimensions and data structures
    p <- ncol(X)
    n <- nrow(X)

    cumP <- compute_cumulative_probabilities(X, gamma)
    # Initialize storage arrays
    beta.gibbs <- array(dim = c(n.gibbs, p))
    sigma.gibbs <- numeric(n.gibbs)

    # Initialize beta and truncate values
    beta.gibbs[1,] <- beta

    # Initialize priors and first sample
    Xbeta <- X %*% beta
    sigma.gibbs[1] <- sqrt(1/rgamma(1, 0.5*n + 1, 0.5*sum((Y - Xbeta)^2)))
    betaind = 1:p
    # Main Gibbs loop
    for (g in 2:n.gibbs) {
        if (!quiet && g %% 100 == 1) message("Iteration ", g)

        # Sample sigma
        sigma <- sqrt(1/rgamma(1, 0.5*n + 1, 0.5*sum((Y - X %*% beta.gibbs[g-1,])^2)))
        sigma.gibbs[g] <- sigma

        # Permute beta coefficients
        beta.gibbs.g   <- beta.gibbs[g-1,]
        if(g%%2 == 0){
            res <- permute_gibbs_beta_normal_cpp(Y, Xbeta, X, sigma,
                                                     beta.gibbs.g,betaind)
        }else{
            res <- permute_gibbs_beta_normal_corr_cpp(Y, Xbeta, X, sigma,
                                                          beta.gibbs.g,betaind,cumP)
        }

        # Update beta values and indices
        beta.gibbs[g,] <- beta.gibbs.g

        # Update linear predictor
        #Xbeta <- X %*% beta.gibbs[g,]
    }

    return(list(
        beta.gibbs = beta.gibbs,
        sigma.gibbs = sigma.gibbs
    ))
}

compute_cumulative_probabilities <- function(X, gamma) {
    # Compute correlation matrix (columns are predictors)
    S <- cor(X)

    # Apply exponential transformation and zero out the diagonal
    P <- exp(gamma * S)
    diag(P) <- 0  # Ensure the diagonal is zero

    # Normalize each row so that probabilities sum to 1.
    # Using sweep to vectorize the row division.
    row_sums <- rowSums(P)
    P <- sweep(P, 1, row_sums, FUN = "/")
    # Replace any NaN (from division by zero) with zeros.
    P[is.nan(P)] <- 0

    # Compute cumulative sum for each row.
    cumP <- t(apply(P, 1, cumsum))

    return(cumP)
}

library(pracma)

# Gradient flow for EB regression

# compute the function [log N(0,tau^2)*g]'(phi)
dlog.conv <- function(centered, log.density, w, tau) {
  log.density.w = log.density + log(w)
  col.max = apply(log.density.w,2,max)
  density.w = t(exp(t(log.density.w)-col.max))
  num = colSums(-centered * density.w / tau^2)
  denom = colSums(density.w)
  return(num/denom)
}

# iteration on phi
langevin.step <- function(phi, centered, log.density, w, tau, XSinvX, XSinvy, eta.phi, Q, Q2) {
  if (is.null(Q)) {
    # no preconditioning
    phi = phi - eta.phi * (XSinvX %*% phi - XSinvy)
    phi = phi + eta.phi * dlog.conv(centered, log.density, w, tau)
    phi = phi + sqrt(2*eta.phi) * rnorm(length(phi))
  } else {
    # precondition by Q^2
    phi = phi - eta.phi * Q2 %*% (XSinvX %*% phi - XSinvy)
    phi = phi + eta.phi * Q2 %*% dlog.conv(centered, log.density, w, tau)
    phi = phi + sqrt(2*eta.phi) * Q %*% rnorm(length(phi))
  }
  return(c(phi))
}

# iteration on w
prior.step <- function(w, log.density, M, eta.w) {
  log.density.w = log.density + log(w)
  col.max = apply(log.density.w,2,max)
  density.w = exp(t(log.density.w)-col.max)
  update.w = colMeans(density.w/rowSums(density.w))
  # add spline penalty (if M != 0)
  grad.spline = as.numeric(M%*%w)
  update.w = update.w - w*grad.spline + w*sum(w*grad.spline)
  return(eta.w * update.w + (1-eta.w) * w)
}

# discrete difference approximation of second derivative
discrete.der = function(K,delta) {
  D = matrix(0,K-2,K)
  for (i in 1:(K-2)) {
    D[i,i] = 1/delta^2
    D[i,i+1] = -2/delta^2
    D[i,i+2] = 1/delta^2
  }
  return(D)
}

gradient.flow.EB.init <- function(X, y, sigma, grid, eta.w, eta.phi, max.iter, tausq.scale, w.init, phi.init, precondition, w.lower.bound) {
  n = nrow(X); p = ncol(X); K = length(grid)
  # Set tau
  X.svd = svd(X,nv=p)
  max.sing = max(X.svd$d)
  if (is.null(tausq.scale)) {
    tausq.scale = 0.5
  } else if (tausq.scale >= 1.0) {
    stop("Sigma = sigma^2 I - tau^2 XX' is not positive definite")
  }
  tau = sqrt(tausq.scale) * sigma / max.sing
  # Set eta.phi
  if (is.null(eta.phi)) {
    # Default: exp decay from 1 to 0.1
    eta.phi = logspace(0,-1,max.iter)
  }
  L = length(eta.phi)
  if (L < max.iter+1) {
    # Repeat last step size for all remaining iters
    eta.phi = c(eta.phi,rep(eta.phi[L],max.iter+1-L))
  } else {
    eta.phi = eta.phi[1:(max.iter+1)]
  }
  # Set eta.w
  if (is.null(eta.w)) {
    # Default: 0.01 * eta.phi
    eta.w = 0.01 * eta.phi
  } else {
    L = length(eta.w)
    if (L < max.iter+1) {
      # Repeat last step size for all remaining iters
      eta.w = c(eta.w,rep(eta.w[L],max.iter+1-L))
    } else {
      eta.w = eta.w[1:(max.iter+1)]
    }
  }
  # Precompute: X'Sigma^{-1}X
  #             X'Sigma^{-1}y
  #             Preconditioning matrices Q, Q^2
  D = X.svd$d^2/(rep(sigma^2,min(n,p))-tau^2*X.svd$d^2)
  if (p > n) { D = c(D,rep(0,p-n)) }
  XSinvX = X.svd$v%*%(D*t(X.svd$v))
  D2 = X.svd$d/(rep(sigma^2,min(n,p))-tau^2*X.svd$d^2)
  XSinvy = X.svd$v[,1:min(n,p)]%*%(D2*(t(X.svd$u)%*%y))
  if (precondition) {
    eta.phi.scale = 1
    Q = X.svd$v%*%(1/sqrt(D+1/tau^2)*t(X.svd$v))
    Q2 = X.svd$v%*%(1/(D+1/tau^2)*t(X.svd$v))
  } else {
    eta.phi.scale = 1/(max(D)+1/tau^2)
    Q = NULL
    Q2 = NULL
  }
  # Set w.init and phi.init
  if (is.null(w.init)) {
    w.init = rep(1/K,K)
  } else {
    w.init[which(w.init<w.lower.bound)] = w.lower.bound
    w.init = w.init/sum(w.init)
  }
  if (is.null(phi.init)) {
    phi.init = rep(0,p)
  }
  return(list(tau=tau, eta.w = eta.w, eta.phi=eta.phi,
              XSinvX=XSinvX, XSinvy=XSinvy,
	      eta.phi.scale=eta.phi.scale, Q=Q, Q2=Q2,
              w.init=w.init, phi.init=phi.init))
}



# ---- W1(CDF) between true empirical betas and current discrete w on grid ----
cdf_from_w_grid <- function(w, grid, x) {
  cw <- cumsum(w)
  idx <- findInterval(x, grid)   # 0..K
  out <- numeric(length(x))
  out[idx > 0] <- cw[idx[idx > 0]]
  out
}

w1_cdf_distance <- function(beta_true, w, grid, bounds = c(-3, 3), nbins = 10000) {
  xfine <- seq(bounds[1], bounds[2], length.out = nbins)
  dx <- xfine[2] - xfine[1]
  F_true <- ecdf(beta_true)(xfine)
  F_est  <- cdf_from_w_grid(w, grid, xfine)
  sum(abs(F_true - F_est)) * dx
}


# Main EBflow algorithm
#
# Inputs:
# X -- n x p design matrix 
# y -- response vector, length n
# sigma -- noise std dev
# grid -- prior support points, length K
# lambda -- spline smoothing penalty size
#
# Algorithm parameters:
# eta.w -- step size for prior update
#        Default: eta.phi * 0.01
# eta.phi -- step size for Langevin dynamics
#        Default: Exponential decay from 1.0 to 0.1
#        Will be scaled by 1/lambda_max without preconditioning,
#        or Q^{-1} with preconditioning
# tausq.scale -- reparametrize by phi = theta + N(0,tau^2)
#        tau^2 = tausq.scale * sigma^2 / ||X||_op^2
#        If NULL, will be initialized to default tausq.scale = 0.5
# w.init -- initial prior weights
#        Default: uniform over support points grid
# w.lower.bound -- lower bound for estimated prior weights
# phi.init -- initialization for Langevin dynamics
#        Default: all-0's vector
# max.iter -- total number of iterations to run
# precondition -- whether to precondition using
#        Q = X(sigma^2-tau^2 XX')^{-1}X+tau^{-2}
#
# Output/display options:
# save.iter -- save w and phi in increments of this many iterations
# verbose -- print progress
# w.true -- true w, used only to compute TV error when printing progress
# print.iter -- print progress in increments of this many iterations
#
# Returns:
# w -- final estimated prior weights, length K
#      (correpsonding to support points grid)
# hist -- saved prior weights w and samples phi every save.iter iterations
# sigma -- noise std dev (same as input)
# tau -- std dev in reparametrization by phi
# eta.w -- vector of prior update step sizes per iteration
# eta.phi -- vector of Langevin step sizes per iteration
gradient.flow.EB <- function(X, y, sigma, grid, lambda=1e-3,
        eta.w=NULL, eta.phi=NULL, tausq.scale=NULL,
        w.init=NULL, w.lower.bound=1e-5, phi.init=NULL,
        max.iter=10000, precondition=TRUE,
        verbose=TRUE, print.iter=100, beta.true = NULL, bounds = c(-3, 3)) {
  # Initialize default algorithm parameters
  K = length(grid)
  p = ncol(X)
  inits = gradient.flow.EB.init(X, y, sigma, grid, eta.w, eta.phi, max.iter, tausq.scale, w.init, phi.init, precondition, w.lower.bound)
  tau=inits$tau; eta.w=inits$eta.w; eta.phi=inits$eta.phi
  XSinvX=inits$XSinvX; XSinvy=inits$XSinvy
  eta.phi.scale=inits$eta.phi.scale; Q=inits$Q; Q2=inits$Q2
  w.init=inits$w.init; phi.init=inits$phi.init
  if (verbose) {
    print(sprintf("Gradient flow EB with sigma = %f, tau = %f", sigma, tau))
  }
  # Run gradient flow EB
  hist = list()
  w = w.init
  phi = phi.init
  # For smoothing spline
  delta = grid[2]-grid[1]
  D = discrete.der(K,delta)
  M = (lambda/delta)*t(D)%*%D
  # For first iteration
  centered = rep(1,K) %o% phi - grid %o% rep(1,p)
  log.density = dnorm(centered, sd=tau, log=TRUE)
  t.start = Sys.time()
  for (iter in 1:(max.iter+1)) {
    # Save history, print change in w
    t.curr = Sys.time()
    hist[[iter]] = list(w = w, phi = phi, time = t.curr - t.start)
    if ((iter %% print.iter == 1) | (print.iter == 1)) {
      if (verbose) {
        if (!is.null(beta.true)) {
          err <- w1_cdf_distance(beta.true, w, grid, bounds = bounds, nbins = 10000)
          print(sprintf("Iteration: %d, W1(CDF) from truth: %f", iter, err))
        } else {
          print(sprintf("Iteration: %d", iter))
        }
      }
    }
    # Precompute Gaussian density at coordinates of phi - grid
    if (max(abs(phi)) > 100 * max(abs(grid))) {
      print("WARNING: Diverging phi iterates, consider reducing eta.phi")
    }
    # Take one Langevin step in phi
    phi = langevin.step(phi, centered, log.density, w, tau, XSinvX, XSinvy,
		    eta.phi[iter]*eta.phi.scale, Q, Q2)
    # Update prior weights w
    centered = rep(1,K) %o% phi - grid %o% rep(1,p)
    log.density = dnorm(centered, sd=tau, log=TRUE)
    w = prior.step(w, log.density, M, eta.w[iter])
    w[which(w<w.lower.bound)] = w.lower.bound
    w = w/sum(w)
  }
  return(list(w = w, hist = hist, sigma = sigma, tau = tau, eta.w = eta.w, eta.phi = eta.phi))
}

# Compute posterior mean E[theta|X,y] from estimated prior
# and Langevin samples of phi
#
# Inputs:
# w -- estimated prior weights, length K
# phi -- Langevin samples, p x n.iters
# grid -- prior support points, length K
# tau -- std dev in reparametrization by phi
#
# Returns posterior mean of theta, length p
posterior.mean = function(w,phi,grid,tau) {
  p = nrow(phi); niters = ncol(phi); p.tot = niters*p; K = length(grid)
  dim(phi) = p.tot
  centered = rep(1,K) %o% phi - grid %o% rep(1,p.tot)
  log.density.w = dnorm(centered, sd=tau, log=TRUE) + log(w)
  col.max = apply(log.density.w,2,max)
  density.w = t(exp(t(log.density.w)-col.max))
  theta.mean = colMeans(grid*density.w)/colMeans(density.w)
  dim(theta.mean) = c(p,niters)
  return(rowMeans(theta.mean))
}

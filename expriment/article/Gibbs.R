library(CVXR)

stable.softmax = function(x) {
  exp_x = exp(x-max(x))
  return (exp_x / sum(exp_x))
}

### Minimize -E_qbar[log g] + lambda/2 int (log g)''
spline.optimize = function(qbar,M) {
  K = length(qbar)
  w = Variable(K)
  obj = -sum(qbar*log(w))+sum((M%*%w)^2)
  const = list(sum(w) == 1)
  prob = Problem(Minimize(obj),const)
  result = solve(prob)
  w = as.numeric(result$getValue(w))
  return(w)
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


Gibbs = function(X, y, sigma, grid, lambda=1e-3, eta.w = 1.0, gibbs.inner.iters = 100,
                  w.init = NULL, phi.init=NULL, burn.iters=200,
                  beta.true = NULL, bounds = c(-3, 3), max.iter=10000, print.iter=100, verbose=TRUE) {
  p = ncol(X)
  n = nrow(X)
  K = length(grid)
  if (is.null(w.init)) {
    w.init = rep(1/K,K)
  } else {
    w.init[which(w.init<1e-5)] = 1e-5
    w.init = w.init/sum(w.init)
  }
  if (is.null(phi.init)) { phi.init = rep(0,p) }
  second_moments = apply(X^2,2,sum)
  phi = phi.init
  w = w.init
  if (verbose) {
    print(sprintf("Gibbs with sigma = %f", sigma))
  }
  r = y - X%*%phi
  delta = grid[2] - grid[1]
  D = discrete.der(K,delta)
  M = sqrt(lambda/(2*delta))*D
  count = rep(0, K)
  hist = list()
  t.start = Sys.time()
  for (iter in 1:(max.iter+1)) {
    # Save history, print error in w
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
    ### update phi
    for (j in 1:p) {
      rj = r + X[,j] * phi[j]
      c1 = sum(rj * X[,j])/sigma^2
      c2 = second_moments[j]/(2*sigma^2)
      log_posteriors = log(w) + c1 * grid - c2*grid^2
      posteriors = stable.softmax(log_posteriors)
      ## generate a sample
      igrid = sample(1:K, size = 1, prob = posteriors)
      count[igrid] = 1+count[igrid]
      phi[j] = grid[igrid]
      r = rj - X[,j] * phi[j]
    }
    if (iter %% gibbs.inner.iters == 0 && iter >= burn.iters) {
      ## update w
      qbar = count/(p*gibbs.inner.iters)
      if (lambda > 0) {
        try({w.new = spline.optimize(qbar,M)
             w = w * (1-eta.w) + w.new * eta.w}, silent=TRUE)
      } else {
        w.new = qbar
        w = w * (1-eta.w) + w.new * eta.w
      }
      w[which(w<1e-5)] = 1e-5
      w = w/sum(w)
      count = rep(0, K)
    }
  }
  return(list(w = w, hist = hist, eta.w = eta.w, gibbs.inner.iters = gibbs.inner.iters))
}

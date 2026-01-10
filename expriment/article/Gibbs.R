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

Gibbs = function(X, y, sigma, Theta, lambda=1e-3, eta.w = 1.0, gibbs.inner.iters = 1,
                  w.init = NULL, phi.init=NULL, burn.iters=0,
                  w.true = NULL, max.iter=10000, save.iter=1, print.iter=100, verbose=TRUE) {
  p = ncol(X)
  n = nrow(X)
  K = length(Theta)
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
  delta = Theta[2] - Theta[1]
  D = discrete.der(K,delta)
  M = sqrt(lambda/(2*delta))*D
  count = rep(0, K)
  hist = list()
  t.start = Sys.time()
  for (iter in 1:(max.iter+1)) {
    if ((iter %% save.iter == 1) | (save.iter == 1)) {
      # Save history, print error in w
      t.curr = Sys.time()
      hist[[iter]] = list(w = w, phi = phi, time = t.curr - t.start)
    }
    if ((iter %% print.iter == 1) | (print.iter == 1)) {
      if (verbose) {
        if (is.null(w.true)) { err = NA }
        else { err = sum(abs(w-w.true))/2 }
        print(sprintf("Iteration: %d, TV from truth: %f", iter, err))
      }
    }
    ### update phi
    for (j in 1:p) {
      rj = r + X[,j] * phi[j]
      c1 = sum(rj * X[,j])/sigma^2
      c2 = second_moments[j]/(2*sigma^2)
      log_posteriors = log(w) + c1 * Theta - c2*Theta^2
      posteriors = stable.softmax(log_posteriors)
      ## generate a sample
      itheta = sample(1:K, size = 1, prob = posteriors)
      count[itheta] = 1+count[itheta]
      phi[j] = Theta[itheta]
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

VI = function(X, y, sigma, Theta, lambda=1e-3, eta.w = 1.0,
              w.init=NULL, posterior.init=NULL, w.true = NULL,
              max.iter=1000, save.iter=1, print.iter=1, verbose=TRUE, plot=FALSE){
  p = ncol(X)
  n = nrow(X)
  K = length(Theta)
  if (is.null(w.init)) {
    w.init = rep(1/K,K)
  } else {
    w.init[which(w.init<1e-5)] = 1e-5
    w.init = w.init/sum(w.init)
  }
  if (is.null(posterior.init)) { posterior.init = w.init %o% rep(1,p) }
  second_moments = apply(X^2,2,sum)
  posterior_weights = posterior.init
  w = w.init
  hist = list()
  if (verbose) { print(sprintf("VI with sigma = %f", sigma)) }
  phi = t(Theta%*%posterior_weights)
  r = y - X%*%phi
  delta = Theta[2]-Theta[1]
  D = discrete.der(K,delta)
  M = sqrt(lambda/(2*delta))*D
  t.start = Sys.time()
  for (iter in 1:(max.iter+1)) {
    if ((iter %% save.iter == 1) | (save.iter == 1)) {
      # Save history, print error in w
      t.curr = Sys.time()
      hist[[iter]] = list(w = w, phi = phi, time = t.curr - t.start)
    }
    if ((iter %% print.iter == 1) | (print.iter == 1)) {
      if (verbose) {
        if (is.null(w.true)) { err = NA }
        else { err = sum(abs(w-w.true))/2 }
        print(sprintf("Iteration: %d, TV from truth: %f", iter, err))
      }
      if (plot) {
        plot(Theta,w.true,type="l")
        lines(Theta,w,col="green")
      }
    }
    ## update phi
    for(j in 1:p) {
      rj = r + X[,j] * phi[j]
      c1 = sum(rj * X[,j])/sigma^2
      c2 = second_moments[j]/(2*sigma^2)
      log_posteriors = log(w) + c1 * Theta - c2*Theta^2
      posteriors = stable.softmax(log_posteriors)
      posterior_weights[,j] = posteriors
      phi[j] = sum(Theta * posteriors)
      r = rj - X[,j] * phi[j]
    }
    ## update w
    qbar = rowMeans(posterior_weights)
    if (lambda > 0) {
      try({w.new = spline.optimize(qbar,M)
           w = w * (1-eta.w) + w.new * eta.w}, silent=TRUE)
    } else {
      w.new = qbar
      w = w * (1-eta.w) + w.new * eta.w
    }
    w[which(w<1e-5)] = 1e-5
    w = w/sum(w)
  }
  return(list(w = w, hist = hist, eta.w = eta.w))
}


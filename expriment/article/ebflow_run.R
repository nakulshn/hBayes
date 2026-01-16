source("./EBflow.R")
library(pracma)


# ---- CLI args: --seed, --sim, --scenario ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed   <- as.integer(get_arg("--seed", 1))



cat("EBflow: seed =", seed, "\n")


lambda            <- as.numeric(get_arg("--lambda", 3e-3))
max_iter          <- as.integer(get_arg("--max_iter", 10000))


# ---------------------------
# Load data
# ---------------------------
datadir <- file.path(getwd(), "data")
datafile <- file.path(datadir, sprintf("data_seed_%03d.rds", seed))
if (!file.exists(datafile)) {
  stop("Data file not found: ", datafile, "\n",
       "Did you generate it with the same --seed and datadir?")
}

d <- readRDS(datafile)

# Expect these fields from your generator:
X <- d$X
y <- d$y

sigma <- d$sigma_true
beta_true <- d$beta_true

cat("Loaded:", datafile, "\n",
    "n =", nrow(X), "p =", ncol(X), "sigma =", sigma, "\n")

bounds <- c(-3, 3)



# ---------------------------
# EBflow adapter for "bundle"
# ---------------------------
make_ebflow_steps <- function(max.iter = 10000, burn.iter = 200,
                              eta.phi = "decay", eta.w = "decay") {

  if (is.character(eta.phi) && eta.phi == "decay") {
    eta.phi.vec <- pracma::logspace(0, -1, max.iter - burn.iter)
  } else {
    eta.phi.vec <- rep(as.numeric(eta.phi), max.iter - burn.iter)
  }
  eta.phi.vec <- c(rep(1, burn.iter), eta.phi.vec)

  if (is.character(eta.w) && eta.w == "decay") {
    eta.w.vec <- pracma::logspace(-2, -3, max.iter - burn.iter)
  } else {
    eta.w.vec <- rep(as.numeric(eta.w), max.iter - burn.iter)
  }
  eta.w.vec <- c(rep(0, burn.iter), eta.w.vec)

  list(eta.phi = eta.phi.vec, eta.w = eta.w.vec)
}

run_ebflow <- function(X, y, sigma, bounds, grid = NULL,
                       K = 65,
                       lambda = 3e-3,
                       precondition = TRUE,
                       seed = 123,
                       burn.iter = 200,
                       max.iter = 10000,
                       eta.phi = "decay",
                       eta.w = "decay",
                       predict = TRUE) {


  if (is.null(grid)) {
    a <- bounds[1]
    b <- bounds[2]
    grid <- seq(a, b, length.out = K)
  }

  steps <- make_ebflow_steps(max.iter = max.iter, burn.iter = burn.iter,
                             eta.phi = eta.phi, eta.w = eta.w)

  set.seed(seed)
  ebflow_out <- gradient.flow.EB(
    X = X, y = y, sigma = sigma, grid = grid,
    lambda = lambda,
    eta.w = steps$eta.w,
    eta.phi = steps$eta.phi,
    tausq.scale = 0.5,
    precondition = precondition,
    max.iter = max.iter,
    verbose = TRUE,
    beta.true = beta_true,
    bounds = bounds
  )

  ebflow_out$grid <- grid
  ebflow_out
}

# ---------------------------
# RUN EBflow + SAVE
# ---------------------------
set.seed(seed + 300)

eb_out <- run_ebflow(
  X, y, sigma, bounds,
  K = 65,
  lambda = lambda,
  precondition = TRUE,
  seed = seed + 300,
  burn.iter = 200,
  max.iter = max_iter,
  eta.phi = "decay",
  eta.w = "decay"
)

w = eb_out$w

# ---- Plot CDFs: final w vs empirical beta_true ----
w_final <- eb_out$w
grid_pts <- eb_out$grid

xfine <- seq(bounds[1], bounds[2], length.out = 10000)
F_true <- ecdf(beta_true)(xfine)
F_est  <- cdf_from_w_grid(w_final, grid_pts, xfine)

plot(xfine, F_true, type = "l",
     xlab = "x", ylab = "CDF",
     main = sprintf("CDFs after %d iters (seed=%d)", max_iter, seed))
lines(xfine, F_est, lty = 2)

legend("topleft",
       legend = c("Empirical CDF (true betas)", "CDF from final w (grid)"),
       lty = c(1, 2), bty = "n")

# Optional: save to file
outdir <- file.path(getwd(), "ebflow_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pngfile <- file.path(outdir, sprintf("cdf_compare_seed_%03d.png", seed))
png(pngfile, width = 900, height = 700)
plot(xfine, F_true, type = "l",
     xlab = "x", ylab = "CDF",
     main = sprintf("CDFs after %d iters (seed=%d)", max_iter, seed))
lines(xfine, F_est, lty = 2)
legend("topleft",
       legend = c("Empirical CDF (true betas)", "CDF from final w (grid)"),
       lty = c(1, 2), bty = "n")
dev.off()

cat("Saved CDF plot:", pngfile, "\n")





w_final <- eb_out$w
phi = eb_out$hist[[max_iter+1]]$phi
steps <- make_ebflow_steps(max.iter = max_iter, burn.iter = 200,
                             eta.phi = "decay", eta.w = "decay")


tmp.out <- gradient.flow.EB(
    X = X, y = y, sigma = sigma, grid = grid_pts,
    lambda = 0,
    eta.w = 0,
    w.init=w_final,
    phi.init=phi,
    eta.phi = steps$eta.phi[max_iter],
    precondition = TRUE,
    max.iter = 200000,
    verbose = TRUE,
    beta.true = beta_true,
    bounds = bounds
  )

posterior.mean.noreduce = function(w,phi,grid_pts,tau) {
  p = nrow(phi); niters = ncol(phi); p.tot = niters*p; K = length(grid_pts)
  dim(phi) = p.tot
  centered = rep(1,K) %o% phi - grid_pts %o% rep(1,p.tot)
  log.density.w = dnorm(centered, sd=tau, log=TRUE) + log(w)
  col.max = apply(log.density.w,2,max)
  density.w = t(exp(t(log.density.w)-col.max))
  theta.mean = colMeans(grid_pts*density.w)/colMeans(density.w)
  dim(theta.mean) = c(p,niters)
#   return(rowMeans(theta.mean))
    return(theta.mean)
}

posterior_mean_one_phi <- function(w, phi_vec, grid_pts, tau) {
  # phi_vec: length p
  K <- length(grid_pts)
  p <- length(phi_vec)

  # K x p
  centered <- outer(grid_pts, phi_vec, FUN = function(th, ph) ph - th)
  logdw <- dnorm(centered, sd = tau, log = TRUE) + log(w)

  # stabilize per-column
  col_max <- apply(logdw, 2, max)
  dw <- exp(sweep(logdw, 2, col_max, "-"))

  num <- colSums(dw * grid_pts)
  den <- colSums(dw)
  num / den
}


# phi = sapply(seq(1,100001,1), function(it) tmp.out$hist[[it]]$phi)
# beta.samples = posterior.mean.noreduce(w,phi,grid_pts,tmp.out$tau)
# eb_out$beta.samples = beta.samples

niters <- 200001
p <- length(tmp.out$hist[[1]]$phi)
beta.samples <- matrix(NA_real_, nrow = p, ncol = niters)

for (it in 1:niters) {
  beta.samples[, it] <- posterior_mean_one_phi(w, tmp.out$hist[[it]]$phi, grid_pts, tmp.out$tau)
  if (it %% 1000 == 0) cat("done", it, "\n")
}

eb_out$beta.samples <- beta.samples


# ---------------------------
# Save output (filename includes seed)
# ---------------------------
outdir <- file.path(getwd(), "ebflow_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("ebflow_out_seed_%03d.rds", seed))

saveRDS(
  eb_out,
  file = outfile
)

cat("Saved:", outfile, "\n")

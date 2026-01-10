source("./EBflow.R")
library(pracma)

# ---- CLI args: --seed, --runid, --sim, --scenario ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed   <- as.integer(get_arg("--seed", 1))
runid  <- as.integer(get_arg("--runid", seed))

# optional controls
sim_id <- as.integer(get_arg("--sim", 1))
scenario_id <- as.integer(get_arg("--scenario", 2))

cat("EBflow: seed =", seed, " runid =", runid,
    " sim_id =", sim_id, " scenario_id =", scenario_id, "\n")

# ---------------------------
# PATHS (match PT output convention)
# ---------------------------
run_dir <- file.path(getwd(), "oracle_section4_out", paste0("run_", runid, "_seed_", seed))


# ---------------------------
# USER SETTINGS
# ---------------------------
sigma_true <- 1.0
bounds <- c(-3, 3)



# ---------------------------
# Data generation (same as original)
# ---------------------------
n <- 8000
p <- 800

X.mat <- matrix(rnorm(n*p, mean = 0, sd = sqrt(1/n)), n, p)

betas <- matrix(0, p, 3)
betas[,1] <- c(rep(-10,100), rep(10,100), rep(0,p - 200))
betas[,2] <- rnorm(p, mean = 0, sd = 1)
#resample betas that are out of [-3, 3] so they all lie in [-3, 3]
out_of_bounds <- which(betas[,2] < -3 | betas[,2] > 3)
while(length(out_of_bounds) > 0) {
  betas[out_of_bounds,2] <- rnorm(length(out_of_bounds), mean = 0, sd = 1)
  out_of_bounds <- which(betas[,2] < -3 | betas[,2] > 3)
}
betas[,3] <- rbinom(p, 1, 0.5) * rnorm(p, mean = 7, sd = 1)

y <- as.numeric(X.mat %*% betas[,2] + sigma_true * rnorm(n))

# ---------------------------
# BUILD BUNDLE (use saved y/X/beta_true!)
# ---------------------------
bundle <- list(
  X = X.mat,
  y = y,
  beta_true = betas[,2],
  sigma_true = sigma_true,
  bounds = bounds
)


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

run_ebflow <- function(bundle,
                       Theta = NULL,
                       K = 65,
                       lambda = 3e-3,
                       precondition = TRUE,
                       seed = 123,
                       burn.iter = 200,
                       max.iter = 10000,
                       save.iter = 100,
                       eta.phi = "decay",
                       eta.w = "decay",
                       tausq.scale = 0.5,
                       predict = TRUE) {

  X <- bundle$X
  y <- bundle$y
  sigma <- bundle$sigma_true

  if (is.null(Theta)) {
    a <- bundle$bounds[1]
    b <- bundle$bounds[2]
    Theta <- seq(a, b, length.out = K)
  }

  steps <- make_ebflow_steps(max.iter = max.iter, burn.iter = burn.iter,
                             eta.phi = eta.phi, eta.w = eta.w)


gen_prior = function(K, prior) {
  w.grid = seq(from=-3, to=3, length=K)
  if (prior == 'gaussian') {
    w.true = dnorm(w.grid, mean=0, sd=1)
  } else if (prior == 'bimodal') {
    w.true = dnorm(w.grid, mean=-1.5, sd=0.5) + dnorm(w.grid, mean=1.5, sd=0.5)
  } else if (prior == 'cauchy') {
    w.true = dcauchy(w.grid,0,0.6)
  } else if (prior == 'skew') {
    w.true = dnorm(w.grid, mean=-2, sd=0.5) + dnorm(w.grid, mean=-1.5, sd=1) + dnorm(w.grid, mean=0, sd=2)
  }
  w.true = w.true/sum(w.true)
  return(list(w.grid=w.grid, w.true=w.true))
}


  set.seed(seed)
  ebflow_out <- gradient.flow.EB(
    X = X, y = y, sigma = sigma, Theta = Theta,
    lambda = lambda,
    eta.w = steps$eta.w,
    eta.phi = steps$eta.phi,
    tausq.scale = tausq.scale,
    precondition = precondition,
    save.iter = save.iter,
    max.iter = max.iter,
    verbose = TRUE,
    w.true = gen_prior(K, 'gaussian')$w.true
  )

  ebflow_out$bundle_meta <- bundle$meta
  ebflow_out$Theta <- Theta
  ebflow_out
}

# ---------------------------
# RUN EBflow + SAVE
# ---------------------------
set.seed(seed + runid)

eb_out <- run_ebflow(
  bundle,
  K = 65,
  lambda = 1e-3,
  precondition = TRUE,
  seed = 123,
  burn.iter = 200,
  max.iter = 10000,
  save.iter = 100,
  eta.phi = "decay",
  eta.w = "decay"
)

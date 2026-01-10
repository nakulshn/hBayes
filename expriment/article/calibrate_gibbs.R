# run_Gibbs_on_PT.R

library(pracma)   # only if you want logspace etc (not required here)
library(CVXR)     # required by spline.optimize()

source("./Gibbs.R")  # <-- contains Gibbs(), stable.softmax(), etc.

# ---- CLI args: --seed, --runid, --sim, --scenario ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed   <- as.integer(get_arg("--seed", 1))
runid  <- as.integer(get_arg("--runid", seed))

sim_id      <- as.integer(get_arg("--sim", 1))
scenario_id <- as.integer(get_arg("--scenario", 2))

# ---- Gibbs tuning args (optional) ----
K_grid            <- as.integer(get_arg("--K", 65))
lambda            <- as.numeric(get_arg("--lambda", 3e-3))
eta_w             <- as.numeric(get_arg("--eta_w", 1.0))
gibbs_inner_iters <- as.integer(get_arg("--gibbs_inner", 100))
burn_iters        <- as.integer(get_arg("--burn", 200))
max_iter          <- as.integer(get_arg("--max_iter", 10000))
save_iter         <- as.integer(get_arg("--save_iter", 100))
print_iter         <- as.integer(get_arg("--print_iter", 100))

cat("Gibbs runner:",
    "seed =", seed,
    "runid =", runid,
    "sim =", sim_id,
    "scenario =", scenario_id, "\n")


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

# If you ever change sigma in PT, you can store it in the scenario list.
# For now PT uses sigma_true <- 1 (hard-coded).
sigma <- 1

# Use same bounds as PT/EBflow (or let CLI override if you want)
bounds <- c(-3, 3)

# ---------------------------
# Theta grid
# ---------------------------
Theta <- seq(bounds[1], bounds[2], length.out = K_grid)

# Optional init:
# - w.init: start uniform (default inside Gibbs)
# - phi.init: could start at ridge or lse from PT; ridge tends to be stable
# phi_init <- scen$beta_ridge
# Ensure it lives on the grid (optional). If you want EXACT grid states:
# phi_init <- Theta[findInterval(phi_init, Theta, all.inside = TRUE)]

# ---------------------------
# Run Gibbs
# ---------------------------
set.seed(seed + runid)

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



gibbs_out <- Gibbs(
  X = X.mat,
  y = y,
  sigma = sigma,
  Theta = Theta,
  lambda = lambda,
  eta.w = eta_w,
  gibbs.inner.iters = gibbs_inner_iters,
  w.init = NULL,
  phi.init = NULL,
  burn.iters = burn_iters,
  w.true = gen_prior(K_grid, 'gaussian')$w.true,
  max.iter = max_iter,
  save.iter = save_iter,
  print.iter = print_iter,
  verbose = TRUE
)


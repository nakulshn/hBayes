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

lambda <- as.numeric(get_arg("--lambda", 3e-3))



cat("EBflow: seed =", seed, " runid =", runid,
    " sim_id =", sim_id, " scenario_id =", scenario_id, "\n")

# ---------------------------
# PATHS (match PT output convention)
# ---------------------------
run_dir <- file.path(getwd(), "oracle_section4_out", paste0("run_", runid, "_seed_", seed))
run_rds_path <- file.path(run_dir, "oracle_section4_thin_results.rds")

if (!file.exists(run_rds_path)) {
  stop("Could not find PT results RDS at:\n  ", run_rds_path,
       "\nDid you run the PT script with the same --seed/--runid?")
}

# ---------------------------
# USER SETTINGS
# ---------------------------
sigma_true <- 1.0
bounds <- c(-3, 3)

# ---------------------------
# LOAD SAVED RUN + PICK SCENARIO
# ---------------------------
obj <- readRDS(run_rds_path)

scen_name <- paste0("scenario", scenario_id)

if (sim_id < 1 || sim_id > length(obj$results)) {
  stop("sim_id out of range. Have ", length(obj$results), " sims in obj$results.")
}

scen <- obj$results[[sim_id]]$scenario[[scen_name]]
if (is.null(scen)) stop("Could not find requested sim/scenario in saved RDS.")

# ---------------------------
# BUILD BUNDLE (use saved y/X/beta_true!)
# ---------------------------
bundle <- list(
  X = obj$X,
  y = scen$y,
  beta_true = scen$beta_true,
  sigma_true = sigma_true,
  bounds = bounds,
  init = list(
    beta_ridge = scen$beta_ridge,
    beta_lasso = scen$beta_lasso,
    beta_lse   = scen$beta_lse
  ),
  meta = list(
    run_dir = run_dir,
    run_rds_path = run_rds_path,
    seed = seed,
    runid = runid,
    sim_id = sim_id,
    scenario_id = scenario_id,
    n = nrow(obj$X),
    p = ncol(obj$X)
  )
)

bundle_out_path <- file.path(
  run_dir,
  sprintf("bundle_run%s_seed%d_lambda%.3f.rds", runid, seed, lambda)
)

saveRDS(bundle, bundle_out_path)
cat("Saved bundle to:\n  ", bundle_out_path, "\n\n", sep = "")

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
    w.true = NULL
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
  lambda = lambda,
  precondition = TRUE,
  seed = 123,
  burn.iter = 200,
  max.iter = 10000,
  save.iter = 100,
  eta.phi = "decay",
  eta.w = "decay"
)

w = eb_out$w
phi = eb_out$hist[[10000+1]]$phi
X <- bundle$X
y <- bundle$y
sigma <- bundle$sigma_true
steps <- make_ebflow_steps(max.iter = 10000, burn.iter = 200,
                             eta.phi = "decay", eta.w = "decay")
a <- bundle$bounds[1]
b <- bundle$bounds[2]
Theta <- seq(a, b, length.out = 65)

tmp.out <- gradient.flow.EB(
    X = X, y = y, sigma = sigma, Theta = Theta,
    lambda = 0,
    eta.w = 0,
    w.init=w,
    phi.init=phi,
    eta.phi = steps$eta.phi[10000],
    tausq.scale = 0.5,
    precondition = TRUE,
    save.iter = 1,
    max.iter = 100000,
    verbose = TRUE,
    w.true = NULL
  )

posterior.mean.noreduce = function(w,phi,Theta,tau) {
  p = nrow(phi); niters = ncol(phi); p.tot = niters*p; K = length(Theta)
  dim(phi) = p.tot
  centered = rep(1,K) %o% phi - Theta %o% rep(1,p.tot)
  log.density.w = dnorm(centered, sd=tau, log=TRUE) + log(w)
  col.max = apply(log.density.w,2,max)
  density.w = t(exp(t(log.density.w)-col.max))
  theta.mean = colMeans(Theta*density.w)/colMeans(density.w)
  dim(theta.mean) = c(p,niters)
#   return(rowMeans(theta.mean))
    return(theta.mean)
}

posterior_mean_one_phi <- function(w, phi_vec, Theta, tau) {
  # phi_vec: length p
  K <- length(Theta)
  p <- length(phi_vec)

  # K x p
  centered <- outer(Theta, phi_vec, FUN = function(th, ph) ph - th)
  logdw <- dnorm(centered, sd = tau, log = TRUE) + log(w)

  # stabilize per-column
  col_max <- apply(logdw, 2, max)
  dw <- exp(sweep(logdw, 2, col_max, "-"))

  num <- colSums(dw * Theta)
  den <- colSums(dw)
  num / den
}


# phi = sapply(seq(1,100001,1), function(it) tmp.out$hist[[it]]$phi)
# beta.samples = posterior.mean.noreduce(w,phi,Theta,tmp.out$tau)
# eb_out$beta.samples = beta.samples

niters <- 100001
p <- length(tmp.out$hist[[1]]$phi)
beta.samples <- matrix(NA_real_, nrow = p, ncol = niters)

for (it in 1:niters) {
  beta.samples[, it] <- posterior_mean_one_phi(w, tmp.out$hist[[it]]$phi, Theta, tmp.out$tau)
  if (it %% 1000 == 0) cat("done", it, "\n")
}

eb_out$beta.samples <- beta.samples



eb_out_path <- file.path(
  run_dir,
  sprintf("ebflow_out_run%s_seed%d_lambda%.3f.rds", runid, seed, lambda)
)

saveRDS(eb_out, eb_out_path)
cat("Saved EBflow output to:\n  ", eb_out_path, "\n", sep = "")

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
# Locate PT output (matches your PT script)
# ---------------------------
run_dir <- file.path(getwd(), "oracle_section4_out", paste0("run_", runid, "_seed_", seed))
run_rds_path <- file.path(run_dir, "oracle_section4_thin_results.rds")

if (!file.exists(run_rds_path)) {
  stop("Could not find PT results at:\n  ", run_rds_path,
       "\nDid you run PT with the same --seed/--runid?")
}

obj <- readRDS(run_rds_path)

# ---------------------------
# Pull sim/scenario data
# ---------------------------
scen_name <- paste0("scenario", scenario_id)

if (sim_id < 1 || sim_id > length(obj$results)) {
  stop("sim_id out of range. obj$results has ", length(obj$results), " sims.")
}

scen <- obj$results[[sim_id]]$scenario[[scen_name]]
if (is.null(scen)) stop("Could not find requested sim/scenario in saved RDS.")

X <- obj$X
y <- scen$y


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

gibbs_out <- Gibbs(
  X = X,
  y = y,
  sigma = sigma,
  Theta = Theta,
  lambda = lambda,
  eta.w = eta_w,
  gibbs.inner.iters = gibbs_inner_iters,
  w.init = NULL,
  phi.init = NULL,
  burn.iters = burn_iters,
  w.true = NULL,
  max.iter = max_iter,
  save.iter = save_iter,
  print.iter = print_iter,
  verbose = TRUE
)

# Attach meta for provenance
gibbs_out$meta <- list(
  run_dir = run_dir,
  run_rds_path = run_rds_path,
  seed = seed,
  runid = runid,
  sim_id = sim_id,
  scenario_id = scenario_id,
  n = nrow(X),
  p = ncol(X),
  bounds = bounds,
  K = K_grid,
  lambda = lambda,
  eta_w = eta_w,
  gibbs_inner_iters = gibbs_inner_iters,
  burn_iters = burn_iters,
  max_iter = max_iter,
  save_iter = save_iter
  )



w = gibbs_out$w
phi = gibbs_out$hist[[max_iter+1]]$phi
tmp.out <- Gibbs(
    X = X,
    y = y,
    sigma = sigma,
    Theta = Theta,
    lambda = 0,
    eta.w = 1.0,
    gibbs.inner.iters = 1000000,
    w.init = w,
    phi.init = phi,
    burn.iters = 0,
    w.true = NULL,
    max.iter = 100000,
    save.iter = 1,
    print.iter = 1000,
    verbose = TRUE
  )

beta.samples = sapply(seq(1,100001,1), function(it) tmp.out$hist[[it]]$phi)
gibbs_out$beta.samples = beta.samples

# ---------------------------
# Save output in run folder
# ---------------------------
out_path <- file.path(
  run_dir,
  sprintf("gibbs_out_run%s_seed%d_lambda%.3f.rds", runid, seed, lambda)
)

saveRDS(gibbs_out, out_path)
cat("Saved Gibbs output to:\n  ", out_path, "\n", sep = "")

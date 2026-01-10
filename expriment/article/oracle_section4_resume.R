# oracle_section4_plots.R
library(NPBayes)
library(glmnet)

# ---------------------------
# Settings you asked for
# ---------------------------
n.gibbs <- 100000
burnin  <- 10
nsim    <- 1

# Which beta scenario to run? (original file has j=1..3)
# If you want all 3 scenarios, set scenario_ids <- 1:3
scenario_ids <- 2  # start with 1 to match a single-figure workflow

# Output folder
out_dir <- file.path(getwd(), "oracle_section4_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1)

# ---------------------------
# Data generation (same as original)
# ---------------------------
n <- 8000
p <- 800

X.mat <- matrix(rnorm(n*p, mean = 0, sd = sqrt(1/n)), n, p)


rds_path <- file.path(out_dir, "oracle_section4_all_results.rds")
obj <- readRDS(rds_path)

# pick which sim + scenario to extend
s <- 1
j <- obj$settings$scenario_ids[1]   # or explicitly: j <- 2

hOld <- obj$results[[s]]$scenario[[paste0("scenario", j)]]$hBayes_raw

beta_last    <- hOld$beta.gibbs[nrow(hOld$beta.gibbs), ]
betaind_last <- hOld$delta.gibbs[nrow(hOld$delta.gibbs), ]
pi_last      <- hOld$pi.gibbs[nrow(hOld$pi.gibbs), ]
a_vec        <- hOld$a.vec

# reconstruct X and Y from what you saved in the scenario
X <- X.mat  # assuming you still have it in your environment exactly as before
Y <- obj$results[[s]]$scenario[[paste0("scenario", j)]]$y
sigma <- 1  # your script uses sigma_true=1 and passes it as fixed.sigma
Xbeta_last <- as.numeric(X %*% beta_last)

n.more <- 50000

hNew <- gibbs.normal.fixed.sigma.resume(
  n.more,
  Y = Y, X = X, sigma = sigma,
  a.int = range(a_vec),   # should be c(-3, 3) in your script
  L = 6,
  init = list(beta = beta_last,
              beta.ind = betaind_last,
              pi = pi_last,
              Xbeta = Xbeta_last),
  cpp = TRUE
)

# Append, dropping the duplicated initial row from the new run (row 1)
hCombined <- hOld
hCombined$pi.gibbs    <- rbind(hOld$pi.gibbs,    hNew$pi.gibbs[-1, , drop=FALSE])
hCombined$beta.gibbs  <- rbind(hOld$beta.gibbs,  hNew$beta.gibbs[-1, , drop=FALSE])
hCombined$delta.gibbs <- rbind(hOld$delta.gibbs, hNew$delta.gibbs[-1, , drop=FALSE])

# Replace in your results object
obj$results[[s]]$scenario[[paste0("scenario", j)]]$hBayes_raw <- hCombined

# Update settings if you want it to reflect the new total
obj$settings$n.gibbs <- nrow(hCombined$pi.gibbs)

# Save a new RDS (don’t overwrite the old one until you’ve verified)
saveRDS(obj, file = file.path(out_dir, "oracle_section4_all_results_RESUMED.rds"))

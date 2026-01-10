library(NPBayes)
library(glmnet)

# ---- CLI args: --seed, --runid, (optional) --scenario ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed   <- as.integer(get_arg("--seed", 1))
runid  <- as.integer(get_arg("--runid", seed))
# optional: allow overriding scenario_ids from CLI
scenario_cli <- get_arg("--scenario", NA)
if (!is.na(scenario_cli)) {
  scenario_ids <- as.integer(strsplit(scenario_cli, ",")[[1]])
} else {
  scenario_ids <- 2
}

set.seed(seed)

# Unique output folder per run
out_dir <- file.path(getwd(), "bimodal_oracle_section4_out", paste0("run_", runid, "_seed_", seed))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Running with seed =", seed, " runid =", runid, "\n")
cat("Output dir:", out_dir, "\n")


# ---------------------------
# Settings you asked for
# ---------------------------
n.gibbs <- 200000
burnin  <- 10
nsim    <- 1

# ---------------------------
# Data generation (same as original)
# ---------------------------
n <- 8000
p <- 800

X.mat <- matrix(rnorm(n*p, mean = 0, sd = sqrt(1/n)), n, p)

betas <- matrix(0, p, 3)
betas[,1] <- c(rep(-10,100), rep(10,100), rep(0,p - 200))
# betas[,2] <- rnorm(p, mean = 0, sd = 1)
# #resample betas that are out of [-3, 3] so they all lie in [-3, 3]
# out_of_bounds <- which(betas[,2] < -3 | betas[,2] > 3)
# while(length(out_of_bounds) > 0) {
#   betas[out_of_bounds,2] <- rnorm(length(out_of_bounds), mean = 0, sd = 1)
#   out_of_bounds <- which(betas[,2] < -3 | betas[,2] > 3)
# }

# initial draw from mixture
mix <- rbinom(p, 1, 0.5)
betas[,2] <- rnorm(p, mean = ifelse(mix == 1, 1.5, -1.5), sd = 0.5)

# optional: reject and resample values outside [-3, 3]
out_of_bounds <- which(betas[,2] < -3 | betas[,2] > 3)
while(length(out_of_bounds) > 0) {
  mix2 <- rbinom(length(out_of_bounds), 1, 0.5)
  betas[out_of_bounds,2] <- rnorm(
    length(out_of_bounds),
    mean = ifelse(mix2 == 1, 1.5, -1.5),
    sd = 0.5
  )
  out_of_bounds <- which(betas[,2] < -3 | betas[,2] > 3)
}


betas[,3] <- rbinom(p, 1, 0.5) * rnorm(p, mean = 7, sd = 1)

inv.log.odds <- function(Xbeta) 1/(1 + exp(-Xbeta))

# ---------------------------
# Helper: make CDF bands from pi.gibbs draws
# pi.gibbs is stored in log space in NPBayes scripts (they use exp(pi.gibbs))
# ---------------------------
cdf_bands_from_pi <- function(pi_gibbs_log, a_vec, burnin) {
  post_idx <- (burnin + 1):nrow(pi_gibbs_log)
  pi_post  <- exp(pi_gibbs_log[post_idx, , drop = FALSE])   # T_post x (K-1)

  # Normalize each draw defensively (should already sum to ~1, but better safe)
  pi_post <- pi_post / rowSums(pi_post)

  # CDF for each draw on support a_vec[-1] with a leading 0 at a_vec[1]
  cdf_post <- t(apply(pi_post, 1, cumsum))                  # T_post x (K-1)
  cdf_post <- cbind(0, cdf_post)                            # T_post x K

  # pointwise quantiles across MCMC draws
  qs <- apply(cdf_post, 2, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  # qs is 3 x K
  list(cdf_draws = cdf_post, q025 = qs[1,], q50 = qs[2,], q975 = qs[3,])
}

# ---------------------------
# Storage
# ---------------------------
results <- vector("list", length = nsim)
names(results) <- paste0("sim", seq_len(nsim))

mse_mat <- array(NA_real_, dim = c(4, nsim, length(scenario_ids)),
                 dimnames = list(c("lse","hBayes","lasso","ridge"), #"oracle",
                                 paste0("sim",1:nsim),
                                 paste0("scenario",scenario_ids)))

# ---------------------------
# Main loop
# ---------------------------
for (s in 1:nsim) {
  cat("=== sim", s, "of", nsim, "===\n")

  sim_out <- list()
  sim_out$scenario <- list()

  for (jj in seq_along(scenario_ids)) {
    j <- scenario_ids[jj]
    cat("  scenario", j, "\n")

    sigma_true <- 1.0
    y <- as.numeric(X.mat %*% betas[,j] + sigma_true * rnorm(n))


    # # oracle (as in original script)
    # beta.oracle <- gibbs.logistic.permute.oracle(n.gibbs, y, X.mat, betas[,j], FALSE)
    # beta.oracle <- colMeans(beta.oracle[(burnin+1):n.gibbs, , drop=FALSE])

    # MLE
fit <- lm(y ~ X.mat + 0)
beta.lse <- coef(fit)

# lasso
fit.glmnet <- cv.glmnet(X.mat, y, alpha = 1, intercept = FALSE, family = "gaussian")
beta.lasso <- as.numeric(coef(fit.glmnet, s="lambda.min"))[-1]

# ridge
fit.glmnet <- cv.glmnet(X.mat, y, alpha = 0, intercept = FALSE, family = "gaussian")
beta.ridge <- as.numeric(coef(fit.glmnet, s="lambda.min"))[-1]

    # prior support bounds (same logic as original)
a.min <- -3 #min(beta.lse) - 1
a.max <- 3 #max(beta.lse) + 1

set.seed(seed + runid)

    # hBayes Gibbs (C++ backend)
hBeta <- gibbs.normal.fixed.sigma(
  n.gibbs, y, X.mat,
  sigma_true,
  c(a.min, a.max),
  6,
  as.vector(beta.ridge),
  cpp = TRUE
)


    beta.hBeta <- colMeans(hBeta$beta.gibbs[(burnin+1):n.gibbs, , drop=FALSE])

    # MSEs
    mse_mat["lse",   s, jj] <- mean((beta.lse    - betas[,j])^2)
    # mse_mat["oracle",s, jj] <- mean((beta.oracle - betas[,j])^2)
    mse_mat["hBayes",s, jj] <- mean((beta.hBeta  - betas[,j])^2)
    mse_mat["lasso", s, jj] <- mean((beta.lasso  - betas[,j])^2)
    mse_mat["ridge", s, jj] <- mean((beta.ridge  - betas[,j])^2)

    # CDF bands for the learned prior on beta
    # NOTE: hBeta$a.vec is the grid; exp(hBeta$pi.gibbs[t,]) is mass on intervals a.vec[-1]
    cdf_info <- cdf_bands_from_pi(hBeta$pi.gibbs, hBeta$a.vec, burnin)

    # Save per-scenario outputs (this is "all data" for later debugging/analysis)
    sim_out$scenario[[paste0("scenario",j)]] <- list(
      y = y,
      beta_true = betas[,j],
      beta_lse = beta.lse,
    #   beta_oracle = beta.oracle,
      beta_hBayes = beta.hBeta,
      beta_lasso = beta.lasso,
      beta_ridge = beta.ridge,
      a_vec = hBeta$a.vec,
      pi_gibbs_log = hBeta$pi.gibbs,
      cdf_q025 = cdf_info$q025,
      cdf_q50  = cdf_info$q50,
      cdf_q975 = cdf_info$q975,
      hBayes_raw = hBeta
    )

    # ---- Plot per-sim CDF bands ----
    fig_file <- file.path(out_dir, sprintf("cdf_sim%02d_scenario%d.png", s, j))
    png(fig_file, width = 900, height = 650)
    plot(hBeta$a.vec, cdf_info$q50, type="l", lwd=2,
         xlab=expression(beta), ylab="CDF",
         main=sprintf("hBayes prior CDF (sim %d, scenario %d)", s, j))
    lines(hBeta$a.vec, cdf_info$q025, lty=2, lwd=2)
    lines(hBeta$a.vec, cdf_info$q975, lty=2, lwd=2)
    dev.off()
  }

  results[[s]] <- sim_out
}

# ---------------------------
# Pooled / average CDF plot across sims
# (You said: “average/concatenation of all sims” — this is the clean version:
#   average the posterior-median CDF curves pointwise.)
# ---------------------------
for (jj in seq_along(scenario_ids)) {
  j <- scenario_ids[jj]

  # assume a_vec is identical across sims for a fixed scenario (it should be, unless a.min/a.max vary).
  # If you see mismatches, we can switch to an interpolation to a common grid.
  a0 <- results[[1]]$scenario[[paste0("scenario",j)]]$a_vec

  cdf_mat <- do.call(rbind, lapply(results, function(sim_out) {
    sim_out$scenario[[paste0("scenario",j)]]$cdf_q50
  }))

  avg_cdf <- colMeans(cdf_mat)

  pooled_file <- file.path(out_dir, sprintf("cdf_pooled_avg_scenario%d.png", j))
  png(pooled_file, width = 900, height = 650)
  plot(a0, avg_cdf, type="l", lwd=2,
       xlab=expression(beta), ylab="CDF",
       main=sprintf("Pooled average CDF across %d sims (scenario %d)", nsim, j))
  dev.off()
}

# ---------------------------
# Save everything
# ---------------------------
saveRDS(list(
  settings = list(n.gibbs=n.gibbs, burnin=burnin, nsim=nsim, scenario_ids=scenario_ids),
  mse_mat = mse_mat,
  results = results,
  X = X.mat
), file = file.path(out_dir, "oracle_section4_all_results.rds"))

cat("\nDone.\nOutputs in:", out_dir, "\n")
print(apply(mse_mat, c(1,3), mean, na.rm=TRUE))

# compute_mse_across_runs.R
# Only include folders run_X_seed_X
# Compute MSE of posterior-mean beta with burnin=100000, max_points=100000

root_dir <- file.path(getwd(), "pt_output")

sim_id     <- 1
scenario   <- 2
burnin     <- 100000
max_points <- 100000

all_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)

# keep only run_X_seed_X
run_dirs <- all_dirs[grepl("^run_[0-9]+_seed_[0-9]+$", basename(all_dirs))]
run_dirs <- run_dirs[
  sapply(basename(run_dirs), function(d) {
    x <- sub("^run_([0-9]+)_seed_([0-9]+)$", "\\1,\\2", d)
    a <- as.integer(strsplit(x, ",")[[1]])
    a[1] == a[2]
  })
]

mse_one_run <- function(rds_file, sim_id, scenario, burnin, max_points) {
  obj <- readRDS(rds_file)
  sc  <- obj$results[[sim_id]]$scenario[[paste0("scenario", scenario)]]

  beta_true <- sc$beta_true

  # --- posterior mean beta (as you already do) ---
  beta_draws <- sc$hBayes_raw$beta.gibbs
  Ttot <- nrow(beta_draws)
  if (burnin >= Ttot) return(c(mse_postmean = NA_real_, mse_ols = NA_real_))

  post  <- beta_draws[(burnin + 1):Ttot, , drop = FALSE]
  Tpost <- nrow(post)
  if (Tpost <= 0) return(c(mse_postmean = NA_real_, mse_ols = NA_real_))

  use_T <- min(max_points, Tpost)
  beta_hat <- colMeans(post[1:use_T, , drop = FALSE])
  mse_postmean <- mean((beta_hat - beta_true)^2)

  # --- OLS (LSE) beta ---
  if (!is.null(sc$beta_lse)) {
    beta_ols <- as.numeric(sc$beta_lse)
    mse_ols  <- mean((beta_ols - beta_true)^2)
  } else {
    mse_ols <- NA_real_
  }

  c(mse_postmean = mse_postmean, mse_ols = mse_ols)
}


res <- data.frame(run = basename(run_dirs),
                  mse_postmean = NA_real_,
                  mse_ols      = NA_real_)

for (i in seq_along(run_dirs)) {
  rds_file <- file.path(run_dirs[i], "oracle_section4_all_results.rds")
  if (!file.exists(rds_file)) next

  out <- mse_one_run(rds_file, sim_id, scenario, burnin, max_points)
  res$mse_postmean[i] <- out["mse_postmean"]
  res$mse_ols[i]      <- out["mse_ols"]

  cat(sprintf("[%d/%d] %s  MSE(postmean)=%.6g  MSE(OLS)=%.6g\n",
              i, length(run_dirs), res$run[i],
              res$mse_postmean[i], res$mse_ols[i]))
}

ok <- is.finite(res$mse)

cat("\n---- Summary (run_X_seed_X only) ----\n")

ok_post <- is.finite(res$mse_postmean)
cat("Runs used (postmean):", sum(ok_post), "\n")
cat(sprintf("Mean MSE (postmean): %.6g\n", mean(res$mse_postmean[ok_post])))
cat(sprintf("SD   MSE (postmean): %.6g\n", sd(res$mse_postmean[ok_post])))

ok_ols <- is.finite(res$mse_ols)
cat("\nRuns used (OLS):", sum(ok_ols), "\n")
cat(sprintf("Mean MSE (OLS): %.6g\n", mean(res$mse_ols[ok_ols])))
cat(sprintf("SD   MSE (OLS): %.6g\n", sd(res$mse_ols[ok_ols])))

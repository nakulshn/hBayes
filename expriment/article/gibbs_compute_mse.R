# compute_gibbs_postmean_mse_across_runs.R
# For gaussian_gibbs_output/run_X_seed_X, compute posterior-mean beta MSE vs beta_true
# using burnin = 0 and (up to) 100000 draws per run. Then summarize mean/sd across runs.

root_dir <- file.path(getwd(), "gaussian_gibbs_output")

# where beta_true lives (same structure as your plotting code)
pt_root_dir <- file.path(getwd(), "pt_output")

sim_id   <- 1
scenario <- 2

burnin     <- 0
use_T_max  <- 100000

# ---- find run dirs: only run_X_seed_X with X==seed ----
all_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- all_dirs[grepl("^run_[0-9]+_seed_[0-9]+$", basename(all_dirs))]
run_dirs <- run_dirs[
  sapply(basename(run_dirs), function(d) {
    x <- sub("^run_([0-9]+)_seed_([0-9]+)$", "\\1,\\2", d)
    a <- as.integer(strsplit(x, ",")[[1]])
    a[1] == a[2]
  })
]

# helper: get beta_true for a run from pt_output (mirrors your plot script)
get_beta_true <- function(run_basename, sim_id, scenario) {
  pt_rds <- file.path(pt_root_dir, run_basename, "oracle_section4_all_results.rds")
  if (!file.exists(pt_rds)) return(NULL)
  pt_obj <- readRDS(pt_rds)
  sc <- pt_obj$results[[sim_id]]$scenario[[paste0("scenario", scenario)]]
  sc$beta_true
}

get_gibbs_rds <- function(run_dir) {
  # dir: run_1_seed_1  --> file: gibbs_out_run1_seed1_lambda0.001.rds
  run_name <- basename(run_dir)

  run_compact <- sub("^run_([0-9]+)_seed_([0-9]+)$",
                     "run\\1_seed\\2",
                     run_name)

  fname <- paste0("gibbs_out_", run_compact, "_lambda0.001.rds")
  fpath <- file.path(run_dir, fname)

  if (!file.exists(fpath)) {
    cat("  [DEBUG] Expected gibbs file not found:\n")
    cat("          ", fpath, "\n")
    cat("  [DEBUG] Files present:\n")
    print(basename(list.files(run_dir)))
    return(NA_character_)
  }

  fpath
}

mse_one_run <- function(run_dir, sim_id, scenario, burnin, use_T_max) {
  run_name <- basename(run_dir)

  beta_true <- get_beta_true(run_name, sim_id, scenario)
  if (is.null(beta_true)) return(NA_real_)

  gibbs_rds <- get_gibbs_rds(run_dir)
  if (!is.character(gibbs_rds) || is.na(gibbs_rds) || !file.exists(gibbs_rds)) {
    return(NA_real_)
  }

  obj <- readRDS(gibbs_rds)
  if (is.null(obj$beta.samples)) return(NA_real_)

  beta_draws <- obj$beta.samples  # expected: n.gibbs x p
  if (!is.matrix(beta_draws)) beta_draws <- as.matrix(beta_draws)

  # transpose if needed (same logic as your plot script)
  if (ncol(beta_draws) != length(beta_true) && nrow(beta_draws) == length(beta_true)) {
    beta_draws <- t(beta_draws)
  }

  Ttot <- nrow(beta_draws)
  p    <- ncol(beta_draws)
  if (p != length(beta_true)) return(NA_real_)
  if (burnin >= Ttot) return(NA_real_)

  post  <- beta_draws[(burnin + 1):Ttot, , drop = FALSE]
  Tpost <- nrow(post)
  if (Tpost <= 0) return(NA_real_)

  use_T <- min(use_T_max, Tpost)
  beta_hat <- colMeans(post[1:use_T, , drop = FALSE])

  mean((beta_hat - beta_true)^2)
}

# ---- run over dirs ----
res <- data.frame(
  run = basename(run_dirs),
  mse_postmean = NA_real_
)

for (i in seq_along(run_dirs)) {
  res$mse_postmean[i] <- mse_one_run(
    run_dirs[i], sim_id, scenario, burnin, use_T_max
  )

  cat(sprintf("[%d/%d] %s  MSE(postmean, burn=%d, T<=%d)=%.6g\n",
              i, length(run_dirs), res$run[i], burnin, use_T_max, res$mse_postmean[i]))
}

ok <- is.finite(res$mse_postmean)

cat("\n---- Summary (gaussian_gibbs_output, run_X_seed_X only) ----\n")
cat("Runs used:", sum(ok), "out of", nrow(res), "\n")
cat(sprintf("Mean MSE: %.6g\n", mean(res$mse_postmean[ok])))
cat(sprintf("SD   MSE: %.6g\n", sd(res$mse_postmean[ok])))

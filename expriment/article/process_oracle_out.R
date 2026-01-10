# paperstyle_cdf_from_rds.R
# Recreate Figure-1-style CDF plots from oracle_section4_all_results.rds
# with configurable burn-in and number of iterations used in the CDF bands.

library(NPBayes)

in_rds  <- file.path(getwd(), "pt_out", "run_2_seed_2", "oracle_section4_all_results.rds") #oracle_section4_out_gaussian_40000_sigma1
out_dir <- file.path(getwd(), "pt_out", "run_2_seed_2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# User-configurable knobs
# ---------------------------
sim_id <- 1          # which simulation to plot
burnin <- 50000        # burn-in iterations to discard
n_use  <- 150000       # how many post-burnin draws to use (NULL = all)
# scenario_ids_override <- NULL  # e.g., 1:3; keep NULL to use settings

# ---------------------------
# Helper: bands from stored log(pi) draws
# ---------------------------
cdf_bands_from_pi <- function(pi_gibbs_log, burnin, n_use = NULL) {
  T <- nrow(pi_gibbs_log)
  if (burnin >= T) stop("burnin must be < number of MCMC iterations saved")

  post_idx <- (burnin + 1):T
  if (!is.null(n_use)) {
    n_use <- as.integer(n_use)
    if (n_use <= 0) stop("n_use must be positive or NULL")
    post_idx <- post_idx[seq_len(min(n_use, length(post_idx)))]
  }

  pi_post <- exp(pi_gibbs_log[post_idx, , drop = FALSE])
  pi_post <- pi_post / rowSums(pi_post)

  cdf_post <- t(apply(pi_post, 1, cumsum))
  cdf_post <- cbind(0, cdf_post)

  # posterior mean CDF (pointwise)
  mean_cdf <- colMeans(cdf_post)

  # pointwise quantiles (for bands)
  qs <- apply(cdf_post, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

  list(
    mean = mean_cdf,
    q025 = qs[1,], q50 = qs[2,], q975 = qs[3,],
    n_draws = nrow(cdf_post), burnin = burnin
  )
}

# ---------------------------
# Load results
# ---------------------------
obj <- readRDS(in_rds)
results  <- obj$results
settings <- obj$settings

scenario_ids <- settings$scenario_ids
# if (!is.null(scenario_ids_override)) scenario_ids <- scenario_ids_override

stopifnot(sim_id >= 1, sim_id <= length(results))

for (j in scenario_ids) {
  scen_name <- paste0("scenario", j)
  sc <- results[[sim_id]]$scenario[[scen_name]]

  a <- sc$a_vec

  # Empirical CDFs
  cdf_true <- ecdf(sc$beta_true)
  cdf_mle  <- ecdf(sc$beta_ridge)

  # Recompute hBayes bands from stored pi draws (log-space)
  cdf_info <- cdf_bands_from_pi(sc$pi_gibbs_log, burnin = burnin, n_use = n_use)
  q025 <- cdf_info$q025
  q50  <- cdf_info$q50
  q975 <- cdf_info$q975

  out_file <- file.path(
    out_dir,
    sprintf("CDF_paperstyle_sim%02d_scenario%d_burn%d_nuse%s.png",
            sim_id, j, burnin, ifelse(is.null(n_use), "ALL", as.character(n_use)))
  )
  png(out_file, width = 900, height = 650)

  par(mfrow = c(1,1))

  plot(a, cdf_mle(a),
       type="l", lwd=2, col="black",
       xlab=expression(beta), ylab="CDF", ylim=c(0,1))

  lines(a, cdf_true(a), lwd=2, col="orange", lty=1)

  if (j == 2) { mu_true <- 0; sd_true <- 1
  } else {     mu_true <- mean(sc$beta_true); sd_true <- sd(sc$beta_true) }
  lines(a, pnorm(a, mean=mu_true, sd=sd_true), lwd=2, col="orange", lty=3)

lines(a, cdf_info$mean, lwd=2, col="darkgreen")
#   lines(a, q025, lwd=2, col="darkgreen", lty=2)
#   lines(a, q975, lwd=2, col="darkgreen", lty=2)

  legend("topleft",
         legend=c("MLE empirical",
                  "True empirical",
                  sprintf("True Gaussian N(%.2f, %.2f^2)", mu_true, sd_true),
                  sprintf("hBayes median (post draws=%d)", cdf_info$n_draws)#,
                #   "hBayes 95% band"
                  ),
         col=c("black","orange","orange","darkgreen"),#,"darkgreen"),
         lty=c(1,1,3,1),#,2),
         lwd=c(2,2,2,2),#,2),
         bty="n")

  dev.off()
  cat("Wrote:", out_file, "\n")
}

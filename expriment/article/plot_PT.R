# pt_overlay_10_cdfs.R

# ---- knobs ----
root_dir    <- file.path(getwd(), "bimodal_pt_output")
pt_rds_name <- "oracle_section4_all_results.rds"
sim_id      <- 1
scenario_id <- 2
burnin      <- 50000
n_use       <- 100000  # set NULL to use all post-burnin
out_png     <- file.path(root_dir, sprintf("PT_overlay_sim%02d_scen%d.png", sim_id, scenario_id))

# ---- helpers ----
truncN01_cdf <- function(x, a = -1, b = 1) {
  Z <- pnorm(b) - pnorm(a)
  pmin(1, pmax(0, (pnorm(x) - pnorm(a)) / Z))
}

trunc_bimodal_cdf <- function(x,
                              a = -1, b = 1,
                              mu = c(-1.5, 1.5),
                              sigma = c(0.5, 0.5),
                              w = c(0.5, 0.5)) {

  stopifnot(length(mu) == length(sigma),
            length(mu) == length(w))

  # mixture CDF at x
  F_mix <- function(z) {
    out <- 0
    for (k in seq_along(mu)) {
      out <- out + w[k] * pnorm(z, mean = mu[k], sd = sigma[k])
    }
    out
  }

  Za <- F_mix(a)
  Zb <- F_mix(b)
  Z  <- Zb - Za

  pmin(1, pmax(0, (F_mix(x) - Za) / Z))
}


post_mean_cdf_from_logpi <- function(pi_gibbs_log, burnin, n_use = NULL) {
  T <- nrow(pi_gibbs_log)
  stopifnot(burnin < T)
  idx <- (burnin + 1):T
  if (!is.null(n_use)) idx <- idx[seq_len(min(length(idx), n_use))]

  pi_post <- exp(pi_gibbs_log[idx, , drop = FALSE])
  pi_post <- pi_post / rowSums(pi_post)

  F_mean <- c(0, cumsum(colMeans(pi_post)))  # length 65
  list(F_mean = F_mean, n_draws = length(idx))
}

get_scenario <- function(obj, sim_id, scenario_id) {
  sc_name <- paste0("scenario", scenario_id)
  obj$results[[sim_id]]$scenario[[sc_name]]
}

keep_run_eq_seed <- function(dir) {
  nm <- basename(dir)
  parts <- strsplit(nm, "_")[[1]]
  if (length(parts) < 4) return(FALSE)
  # expecting "run_X_seed_Y"
  run_num  <- suppressWarnings(as.integer(parts[2]))
  seed_num <- suppressWarnings(as.integer(parts[4]))
  is.finite(run_num) && is.finite(seed_num) && (run_num == seed_num)
}

# ---- gather CDFs ----
run_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[sapply(run_dirs, keep_run_eq_seed)]

a_common <- NULL
F_list   <- list()
labels   <- character()

for (d in run_dirs) {
  pt_path <- file.path(d, pt_rds_name)
  if (!file.exists(pt_path)) next

  obj <- readRDS(pt_path)
  sc  <- get_scenario(obj, sim_id, scenario_id)

  a <- sc$a_vec
  info <- post_mean_cdf_from_logpi(sc$pi_gibbs_log, burnin, n_use)

  if (is.null(a_common)) a_common <- a
  if (length(a) != length(a_common) || any(abs(a - a_common) > 1e-12)) {
    stop("Not all runs share the same a_vec grid. (Need a common grid to overlay.)")
  }

  F_list[[length(F_list) + 1]] <- info$F_mean
  labels <- c(labels, basename(d))
}

if (length(F_list) == 0) stop("No runs found / no matching rds files.")

# ---- plot: all PT lines + true truncated normal ----
# ---- plot: all PT lines + mean PT + true truncated normal ----
F_true <- trunc_bimodal_cdf(a_common, a = -3, b = 3)

F_mat  <- do.call(rbind, F_list)    # n_runs x 65
F_mean <- colMeans(F_mat)

png(out_png, width = 900, height = 650)

plot(a_common, F_true, type = "n", ylim = c(0, 1),
     xlab = expression(beta), ylab = "CDF",
     main = sprintf("PT posterior mean CDFs (n=%d) | sim %d, scenario %d",
                    nrow(F_mat), sim_id, scenario_id))

# individual PT runs (thin, semi-transparent)
pt_col <- grDevices::adjustcolor("gray30", alpha.f = 0.35)
for (i in seq_len(nrow(F_mat))) {
  lines(a_common, F_mat[i, ], col = pt_col, lwd = 1.3)
}

# mean PT CDF (thick, black)
lines(a_common, F_mean, col = "black", lwd = 3)

# true truncated normal (thick, blue)
lines(a_common, F_true, col = "dodgerblue3", lwd = 3, lty = 2)

legend("topleft",
       legend = c("PT posterior mean (each run)",
                  "Mean PT CDF across runs",
                  "Trunc 1/2*N(-1.5,0.5^2) + 1/2*N(1.5,0.5^2) on [-3,3]"),
       col = c(pt_col, "black", "dodgerblue3"),
       lwd = c(2, 3, 3),
       lty = c(1, 1, 2),
       bty = "n")

dev.off()
message("Wrote: ", out_png)

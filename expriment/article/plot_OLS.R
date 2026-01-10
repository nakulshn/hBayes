# lse_overlay_empirical_cdfs.R

# ---- knobs ----
root_dir    <- file.path(getwd(), "bimodal_pt_output")
rds_name    <- "oracle_section4_all_results.rds"
sim_id      <- 1
scenario_id <- 2
out_png     <- file.path(root_dir, "OLS_empCDF.png")

# ---- helpers ----
trunc_bimodal_cdf <- function(x,
                              a = -3, b = 3,
                              mu = c(-1.5, 1.5),
                              sigma = c(0.5, 0.5),
                              w = c(0.5, 0.5)) {
  stopifnot(length(mu) == length(sigma), length(mu) == length(w))

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

get_scenario <- function(obj, sim_id, scenario_id) {
  sc_name <- paste0("scenario", scenario_id)
  obj$results[[sim_id]]$scenario[[sc_name]]
}

keep_run_eq_seed <- function(dir) {
  nm <- basename(dir)
  parts <- strsplit(nm, "_")[[1]]
  if (length(parts) < 4) return(FALSE)
  run_num  <- suppressWarnings(as.integer(parts[2]))
  seed_num <- suppressWarnings(as.integer(parts[4]))
  is.finite(run_num) && is.finite(seed_num) && (run_num == seed_num)
}

# empirical CDF evaluated on grid x_grid
emp_cdf_on_grid <- function(samples, x_grid) {
  samples <- samples[is.finite(samples)]
  if (length(samples) == 0) stop("beta_lse has no finite values.")
  # F(x) = P(sample <= x)
  sapply(x_grid, function(x) mean(samples <= x))
}

# ---- gather LSE empirical CDFs ----
run_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[sapply(run_dirs, keep_run_eq_seed)]

a_common <- NULL
F_list   <- list()
labels   <- character()

for (d in run_dirs) {
  path <- file.path(d, rds_name)
  if (!file.exists(path)) next

  obj <- readRDS(path)
  sc  <- get_scenario(obj, sim_id, scenario_id)

  if (is.null(sc$a_vec)) stop("Missing sc$a_vec in: ", path)
  if (is.null(sc$beta_lse)) stop("Missing sc$beta_lse in: ", path)

  a <- sc$a_vec
  beta_lse <- sc$beta_lse

  if (is.null(a_common)) a_common <- a
  if (length(a) != length(a_common) || any(abs(a - a_common) > 1e-12)) {
    stop("Not all runs share the same a_vec grid. (Need a common grid to overlay.)")
  }

  F_emp <- emp_cdf_on_grid(beta_lse, a_common)

  F_list[[length(F_list) + 1]] <- F_emp
  labels <- c(labels, basename(d))
}

if (length(F_list) == 0) stop("No runs found / no matching rds files.")

# ---- plot ----
F_true <- trunc_bimodal_cdf(a_common, a = -3, b = 3)#, mu = 0, sigma = 1, w = 1)

F_mat  <- do.call(rbind, F_list)  # n_runs x length(a_common)
F_mean <- colMeans(F_mat)

png(out_png, width = 900, height = 650)

plot(a_common, F_true, type = "n", ylim = c(0, 1),
     xlab = expression(beta), ylab = "Empirical CDF",
     main = sprintf("LSE empirical CDFs (n=%d runs) | sim %d, scenario %d",
                    nrow(F_mat), sim_id, scenario_id))

# individual runs (thin)
lse_col <- grDevices::adjustcolor("gray30", alpha.f = 0.35)
for (i in seq_len(nrow(F_mat))) {
  lines(a_common, F_mat[i, ], col = lse_col, lwd = 1.3)
}

# mean empirical CDF across runs (optional)
lines(a_common, F_mean, col = "black", lwd = 3)

# true prior CDF
lines(a_common, F_true, col = "dodgerblue3", lwd = 3, lty = 2)

legend("topleft",
       legend = c("Empirical CDF of beta_OLS (each run)",
                  "Mean beta_OLS empirical CDF across runs",
                  "True prior CDF (truncated bimodal on [-3,3])"),
       col = c(lse_col, "black", "dodgerblue3"),
       lwd = c(2, 3, 3),
       lty = c(1, 1, 2),
       bty = "n")

dev.off()
message("Wrote: ", out_png)

# pt_overlay_10_pdfs.R

# ---- knobs ----
root_dir    <- file.path(getwd(), "pt_output")
pt_rds_name <- "oracle_section4_all_results.rds"
sim_id      <- 1
scenario_id <- 2
burnin      <- 100000
n_use       <- 100000  # set NULL to use all post-burnin
out_png     <- file.path(root_dir, sprintf("PT_overlay_PDF_sim%02d_scen%d.png", sim_id, scenario_id))

# ---- helpers ----
# Truncated bimodal *density* on [a,b]
trunc_bimodal_pdf <- function(x,
                              a = -3, b = 3,
                              mu = c(-1.5, 1.5),
                              sigma = c(0.5, 0.5),
                              w = c(0.5, 0.5)) {
  stopifnot(length(mu) == length(sigma), length(mu) == length(w))

  # mixture pdf
  f_mix <- function(z) {
    out <- 0
    for (k in seq_along(mu)) {
      out <- out + w[k] * dnorm(z, mean = mu[k], sd = sigma[k])
    }
    out
  }

  # mixture CDF for normalization constant on [a,b]
  F_mix <- function(z) {
    out <- 0
    for (k in seq_along(mu)) {
      out <- out + w[k] * pnorm(z, mean = mu[k], sd = sigma[k])
    }
    out
  }

  Z <- F_mix(b) - F_mix(a)
  f <- f_mix(x) / Z
  # zero outside [a,b]
  f * (x >= a & x <= b)
}

post_mean_pi_from_logpi <- function(pi_gibbs_log, burnin, n_use = NULL) {
  T <- nrow(pi_gibbs_log)
  stopifnot(burnin < T)
  idx <- (burnin + 1):T
  if (!is.null(n_use)) idx <- idx[seq_len(min(length(idx), n_use))]

  pi_post <- exp(pi_gibbs_log[idx, , drop = FALSE])
  pi_post <- pi_post / rowSums(pi_post)

  pi_mean <- colMeans(pi_post)  # length K = (#bins)
  list(pi_mean = pi_mean, n_draws = length(idx))
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

# Convert bin masses -> piecewise-constant density at bin midpoints
pi_to_density <- function(a_vec, pi) {
  da <- diff(a_vec)                       # length K
  mid <- (a_vec[-1] + a_vec[-length(a_vec)]) / 2
  f <- pi / da
  list(x = mid, f = f)
}

# ---- gather PDFs ----
run_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[sapply(run_dirs, keep_run_eq_seed)]

a_common <- NULL
dens_list <- list()
labels    <- character()

for (d in run_dirs) {
  pt_path <- file.path(d, pt_rds_name)
  if (!file.exists(pt_path)) next

  obj <- readRDS(pt_path)
  sc  <- get_scenario(obj, sim_id, scenario_id)

  a <- sc$a_vec
  info <- post_mean_pi_from_logpi(sc$pi_gibbs_log, burnin, n_use)

  if (is.null(a_common)) a_common <- a
  if (length(a) != length(a_common) || any(abs(a - a_common) > 1e-12)) {
    stop("Not all runs share the same a_vec grid. (Need a common grid to overlay.)")
  }

  dens <- pi_to_density(a, info$pi_mean)
  dens_list[[length(dens_list) + 1]] <- dens$f
  labels <- c(labels, basename(d))
}

if (length(dens_list) == 0) stop("No runs found / no matching rds files.")

# Common x-grid for densities
mid_common <- (a_common[-1] + a_common[-length(a_common)]) / 2

# True truncated bimodal PDF on the same x grid
f_true <- trunc_bimodal_pdf(mid_common, a = -3, b = 3, mu = 0, sigma = 1, w = 1)


# Matrix of densities: n_runs x K
f_mat  <- do.call(rbind, dens_list)
f_mean <- colMeans(f_mat)

# y-limit
ymax <- max(c(f_mat, f_true), finite = TRUE)

png(out_png, width = 900, height = 650)

plot(mid_common, f_true, type = "n", ylim = c(0, ymax),
     xlab = expression(beta), ylab = "PDF",
     main = sprintf("PT posterior mean PDFs (n=%d) | sim %d, scenario %d",
                    nrow(f_mat), sim_id, scenario_id))

# individual PT runs (thin, semi-transparent)
pt_col <- grDevices::adjustcolor("gray30", alpha.f = 0.35)
for (i in seq_len(nrow(f_mat))) {
  lines(mid_common, f_mat[i, ], col = pt_col, lwd = 1.3)
}

# mean PDF across runs (thick, black)
lines(mid_common, f_mean, col = "black", lwd = 3)

# true truncated bimodal PDF (thick, blue)
lines(mid_common, f_true, col = "dodgerblue3", lwd = 3, lty = 2)

legend("topright",
       legend = c("PT posterior mean (each run)",
                  "Mean PT PDF across runs",
                  "Trunc 1/2*N(-1.5,0.5^2)+1/2*N(1.5,0.5^2) on [-3,3]"),
       col = c(pt_col, "black", "dodgerblue3"),
       lwd = c(2, 3, 3),
       lty = c(1, 1, 2),
       bty = "n")

dev.off()
message("Wrote: ", out_png)

# plot_beta_mse_vs_iters.R
# Compute/plot MSE of posterior-mean beta vs true beta as function of post-burnin iters.

rds_file    <- file.path(getwd(), "pt_output", "run_2_seed_2", "oracle_section4_all_results.rds")
out_file    <- file.path(getwd(), "pt_output", "run_2_seed_2", "PT_plot_mse_vs_iters.png")

sim_id     <- 1
scenario   <- 2

burnins    <- c(50000, 100000)  # toggle as you like
step       <- 100                                 # evaluate every 'step' iters (keeps it fast)
max_points <- 50000                                # cap #points per curve (optional)

obj <- readRDS(rds_file)
sc  <- obj$results[[sim_id]]$scenario[[paste0("scenario", scenario)]]

beta_true  <- sc$beta_true
beta_draws <- sc$hBayes_raw$beta.gibbs            # matrix: n.gibbs x p

Ttot <- nrow(beta_draws)
p    <- ncol(beta_draws)

mse_curve <- function(burn, step, max_points) {
  if (burn >= Ttot - 1) return(NULL)
  post <- beta_draws[(burn + 1):Ttot, , drop = FALSE]
  Tpost <- nrow(post)

  # cumulative posterior means at selected iteration counts
  idx <- seq(step, Tpost, by = step)
  if (length(idx) > max_points) idx <- idx[seq_len(max_points)]

  cs  <- apply(post, 2, cumsum)                   # Tpost x p
  mu  <- cs[idx, , drop = FALSE] / idx            # posterior mean using first t draws
  mse <- rowMeans((mu - matrix(beta_true, nrow(mu), p, byrow = TRUE))^2)

  data.frame(t = idx, mse = mse, burn = burn)
}

curves <- do.call(rbind, lapply(burnins, mse_curve, step = step, max_points = max_points))
stopifnot(!is.null(curves) && nrow(curves) > 0)

# ---- Plot (multiple burn-ins as multiple curves) ----
png(out_file, width = 900, height = 650)

plot(NA, xlim = range(curves$t), ylim = range(curves$mse),
     xlab = "Iterations after burn-in (t)",
     ylab = "MSE of posterior mean beta",
     main = sprintf("MSE vs iterations (sim %d, scenario %d)", sim_id, scenario))

for (b in burnins) {
  d <- curves[curves$burn == b, ]
  if (nrow(d) == 0) next
  lines(d$t, d$mse, lwd = 2)
}

legend("topright",
       legend = paste0("burn=", burnins),
       lwd = 2, bty = "n")
dev.off()

cat("Saved plot to:", out_file, "\n")

# plot_polyatree_plus_ebflow_cdf.R

library(NPBayes)

# --- paths (match your current workflow) ---
in_rds    <- file.path(getwd(), "pt_out", "run_1_seed_1", "oracle_section4_all_results.rds")  # Polya-tree results
eb_out_rds <- file.path(getwd(), "pt_out", "run_1_seed_1", "ebflow_out_run1_seed1_sim01_scenario2.rds") # EBflow output

out_dir <- file.path(getwd(), "pt_out", "run_1_seed_1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# knobs (same as your plot script)
# ---------------------------
sim_id <- 1
burnin <- 50000
n_use  <- 150000


# ---------------------------
# Helper: 1D KS distance from CDFs on a common grid
# ---------------------------
KS_from_cdf_grid <- function(x, F, G) {
  stopifnot(length(x) == length(F), length(x) == length(G), length(x) >= 2)
  o <- order(x)
  x <- x[o]; F <- F[o]; G <- G[o]

  D <- abs(F - G)
  sum(0.5 * (D[-1] + D[-length(D)]) * diff(x))

#   stopifnot(length(x) == length(F), length(x) == length(G), length(x) >= 1)
#   o <- order(x)
#   F <- F[o]; G <- G[o]
#   max(abs(F - G), na.rm = TRUE)

}


# ---------------------------
# Helper: bands from stored log(pi) draws (your existing logic)
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

  mean_cdf <- colMeans(cdf_post)
  qs <- apply(cdf_post, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

  list(mean = mean_cdf, q025 = qs[1,], q50 = qs[2,], q975 = qs[3,], n_draws = nrow(cdf_post))
}

# ---------------------------
# Helper: robustly extract last w-hat and Theta from ebflow_out
# (since EBflow.R isn't in this chat, we defensively check common field names)
# ---------------------------
extract_final_w_theta <- function(eb) {
  # Theta grid
  Theta <- NULL
  for (nm in c("Theta", "theta", "grid", "Theta.grid", "theta.grid")) {
    if (!is.null(eb[[nm]]) && is.numeric(eb[[nm]])) { Theta <- eb[[nm]]; break }
  }
  if (is.null(Theta)) stop("Couldn't find Theta grid in ebflow_out. Try str(eb_out) and adjust field name.")

  # final w
  w <- NULL
  # common “current” fields
  for (nm in c("w", "w.hat", "w_hat", "w.npmle", "w_npmle")) {
    if (!is.null(eb[[nm]]) && is.numeric(eb[[nm]])) { w <- eb[[nm]]; break }
  }
  # common “saved path” fields (matrix/list)
  if (is.null(w)) {
    for (nm in c("w.save", "w_saved", "w.hist", "w_history", "w.path", "w_path")) {
      obj <- eb[[nm]]
      if (is.null(obj)) next
      if (is.matrix(obj) && ncol(obj) == length(Theta)) { w <- obj[nrow(obj), ]; break }
      if (is.list(obj) && length(obj) > 0 && is.numeric(obj[[length(obj)]])) { w <- obj[[length(obj)]]; break }
    }
  }
  if (is.null(w)) stop("Couldn't find final w in ebflow_out. Try str(eb_out) and adjust field name.")

  if (length(w) != length(Theta)) stop(sprintf("length(w)=%d but length(Theta)=%d", length(w), length(Theta)))

  # normalize (NPMLE weights should already be nonneg + sum 1, but be safe)
  w <- pmax(w, 0)
  if (sum(w) <= 0) stop("Final w sums to 0 after truncation; something is wrong.")
  w <- w / sum(w)

  list(Theta = Theta, w = w)
}

smooth_pt_pdf_qp <- function(p, a, lambda = 1e-2,
                             boundary_knots = c(-1.5, 0, 1.5),
                             boundary_weight = 0, boundary_sigma_bins = 1) {
  # p: length-64 bin masses (sum to 1)
  # a: length-65 grid edges (monotone), defines 64 bins. For you, a = sc$a_vec
  # lambda: overall roughness strength
  # boundary_weight: extra penalty mass near coarse boundaries (0 = off)
  # boundary_sigma_bins: spread (in bins) of extra penalty around those boundaries

  stopifnot(length(a) == length(p) + 1)
  n <- length(p)

  # bin widths and pdf heights
  dx <- diff(a)
  stopifnot(all(dx > 0))
  # treat as uniform if it's intended to be dyadic equal bins
  # but we keep dx vector so this works even if slightly non-uniform
  # heights h = q / dx

  # Second-difference operator on heights: (h_{i-1} - 2 h_i + h_{i+1})
  # Size: (n-2) x n
  D2 <- matrix(0, nrow = n - 2, ncol = n)
  for (i in 2:(n - 1)) {
    row <- i - 1
    D2[row, i - 1] <- 1
    D2[row, i]     <- -2
    D2[row, i + 1] <- 1
  }

  # Convert heights to masses: h = q / dx  =>  h = diag(1/dx) q
  S <- diag(1 / dx, n, n)

  # Roughness term: || D2 * h ||^2 = || D2 * S * q ||^2
  R <- D2 %*% S  # (n-2) x n

  # Optional: extra penalty near coarse boundaries (in *bin index* space)
  w <- rep(1, n - 2)
  if (boundary_weight > 0) {
    mids <- 0.5 * (a[-1] + a[-length(a)])   # bin midpoints length n
    # map each second-diff row (centered at bin i) to midpoint i
    # row r corresponds to i = r+1 (since i runs 2..n-1)
    row_mids <- mids[2:(n - 1)]
    for (bk in boundary_knots) {
      # distance in bins, approximated by distance / mean(dx)
      dist_bins <- abs(row_mids - bk) / mean(dx)
      w <- w + boundary_weight * exp(-0.5 * (dist_bins / boundary_sigma_bins)^2)
    }
  }
  W <- diag(w, n - 2, n - 2)

  # Objective:
  # 0.5 * (q-p)'(q-p) + 0.5*lambda * (R q)' W (R q)
  # => 0.5 q' (I + lambda R' W R) q - p' q + const
  Dmat <- diag(1, n) + lambda * t(R) %*% W %*% R
  dvec <- p

  # Constraints for quadprog:
  # quadprog solves: min 0.5 b' D b - d' b  s.t. A' b >= b0
  # We want:
  #   q_i >= 0
  #   sum(q) = 1  -> implement as two inequalities: sum(q) >= 1 and -sum(q) >= -1
  Amat <- cbind(diag(1, n),
                rep(1, n),
                -rep(1, n))
  bvec <- c(rep(0, n), 1, -1)

  # Solve
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Please install quadprog: install.packages('quadprog')")
  }
  sol <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 0)

  q <- sol$solution
  # numerical cleanup
  q[q < 0] <- 0
  q <- q / sum(q)

  # Return everything you need
  list(
    q = q,                         # smoothed bin masses (len 64)
    pdf = q / dx,                  # smoothed pdf heights (len 64)
    cdf = c(0, cumsum(q)),          # smoothed CDF on a grid (len 65)
    value = sol$value,
    weights = w
  )
}


# ---------------------------
# Load results
# ---------------------------
obj <- readRDS(in_rds)
results  <- obj$results
settings <- obj$settings
scenario_ids <- settings$scenario_ids

eb_out <- readRDS(eb_out_rds)
eb <- extract_final_w_theta(eb_out)
Theta <- eb$Theta
cdf_eb <- cumsum(eb$w)

# step-function CDF evaluator on any x (so we can draw on the Polya-tree 'a' grid)
cdf_eb_fun <- approxfun(Theta, cdf_eb, method = "constant", f = 0, yleft = 0, yright = 1)

stopifnot(sim_id >= 1, sim_id <= length(results))

for (j in scenario_ids) {
  scen_name <- paste0("scenario", j)
  sc <- results[[sim_id]]$scenario[[scen_name]]
  a <- sc$a_vec

  cdf_true <- ecdf(sc$beta_true)
  cdf_mle  <- ecdf(sc$beta_lse)

  cdf_info <- cdf_bands_from_pi(sc$pi_gibbs_log, burnin = burnin, n_use = n_use)


F_pt <- cdf_info$mean           # length 65
# p_pt <- diff(F_pt)              # length 64, bin masses

# sm <- smooth_pt_pdf_qp(
#   p = p_pt,
#   a = a,
#   lambda = 1e-2,
#   boundary_knots = c(-1.5, 0, 1.5),
#   boundary_weight = 5,          # try 0 (off), 1, 5, 10
#   boundary_sigma_bins = 1       # spread (in bins)
# )

# F_pt <- sm$cdf


F_eb <- cdf_eb_fun(a)

W1_pt <- KS_from_cdf_grid(a, F_pt, cdf_true(a))
W1_eb <- KS_from_cdf_grid(a, F_eb, cdf_true(a))
cat(sprintf("sim %d scenario %d: KS(PT, True Empirical CDF) = %.6f\n", sim_id, j, W1_pt))
cat(sprintf("sim %d scenario %d: KS(EBflow, True Empirical CDF) = %.6f\n", sim_id, j, W1_eb))

cat("out_dir:", out_dir, "\n")

  out_file <- file.path(out_dir,
                        sprintf("CDF_plus_EBflow_sim%02d_scenario%d.png", sim_id, j))
  png(out_file, width = 900, height = 650)

  plot(a, cdf_mle(a),
       type="l", lwd=2, col="black",
       xlab=expression(beta), ylab="CDF", ylim=c(0,1))

  lines(a, cdf_true(a), lwd=2, col="orange", lty=1)

  if (j == 2) { mu_true <- 0; sd_true <- 1
  } else {     mu_true <- mean(sc$beta_true); sd_true <- sd(sc$beta_true) }
#   lines(a, pnorm(a, mean=mu_true, sd=sd_true), lwd=2, col="orange", lty=3)

lines(a, F_pt, lwd=2, col="darkgreen")
lines(a, cdf_info$q025, lwd=2, col="darkgreen", lty=2)
lines(a, cdf_info$q975, lwd=2, col="darkgreen", lty=2)

  # --- EBflow NPMLE CDF overlay ---
  lines(a, cdf_eb_fun(a), lwd=2, col="blue")

  legend("topleft",
         legend=c("MLE empirical",
                  "True empirical",
                #   sprintf("True Gaussian N(%.2f, %.2f^2)", mu_true, sd_true),
                  sprintf("hBayes mean (post draws=%d)", cdf_info$n_draws),
                  "hBayes 95% band",
                  "EBflow NPMLE CDF (final w)"),
         col=c("black","orange","darkgreen","darkgreen","blue"),
         lty=c(1,1,1,2,1),
         lwd=c(2,2,2,2,2),
         bty="n")

  dev.off()
  cat("Wrote:", out_file, "\n")
}

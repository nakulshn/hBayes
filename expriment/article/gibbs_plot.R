# gibbs_overlay_10_cdfs.R

# ---- knobs ----
root_dir      <- file.path(getwd(), "bimodal_pt_output")  # <- changed
sim_id        <- 1
scenario_id   <- 2

# If your Gibbs files have consistent naming, set a pattern:
# e.g. "gibbs_out_run*_seed*_lambda*.rds"
gibbs_pattern <- "\\.rds$"

out_png <- file.path(root_dir, sprintf("Gibbs_overlay_sim%02d_scen%d.png", sim_id, scenario_id))

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


keep_run_eq_seed <- function(dir) {
  nm <- basename(dir)
  parts <- strsplit(nm, "_")[[1]]
  if (length(parts) < 4) return(FALSE)
  # expecting "run_X_seed_Y"
  run_num  <- suppressWarnings(as.integer(parts[2]))
  seed_num <- suppressWarnings(as.integer(parts[4]))
  is.finite(run_num) && is.finite(seed_num) && (run_num == seed_num)
}

parse_run_seed <- function(dir) {
  nm <- basename(dir)                 # "run_2_seed_2"
  parts <- strsplit(nm, "_")[[1]]     # c("run","2","seed","2")
  if (length(parts) < 4) return(NULL)
  run  <- suppressWarnings(as.integer(parts[2]))
  seed <- suppressWarnings(as.integer(parts[4]))
  if (!is.finite(run) || !is.finite(seed)) return(NULL)
  list(run = run, seed = seed)
}

make_gibbs_filename <- function(run, seed, lambda = 0.001) {
  # IMPORTANT: match the exact formatting of your saved files.
  # Your example: gibbs_out_run2_seed2_lambda0.001.rds
  sprintf("gibbs_out_run%d_seed%d_lambda%.3f.rds", run, seed, lambda)
}


# Grab iter w
get_final_w <- function(g) {
  if (is.null(g$hist) || length(g$hist) == 0) stop("No hist in gibbs_out.")
  #print length(g$hist)
  cat("Number of Gibbs iterations stored in hist:", length(g$hist), "\n")
  w <- g$hist[[10001]]$w
  w <- as.numeric(w)
  w / sum(w)
}

w_to_cdf <- function(w) {
    cumsum(w)
    #TODO implement the alternate version projecting discrete pdf onto "uniform on bins" pdf
}

# w: length 65 masses at boundary points Theta[1:65]
# returns bin masses p: length 64, sum(p)=1
disc_to_hist_bins <- function(w) {
  w <- as.numeric(w)
  w <- w / sum(w)
  Kp1 <- length(w)        # 65
  K   <- Kp1 - 1          # 64

  p <- numeric(K)
  p[1]   <- w[1] + 0.5 * w[2]
  if (K > 2) {
    for (i in 2:(K-1)) {
      p[i] <- 0.5 * w[i] + 0.5 * w[i+1]
    }
  }
  p[K]   <- 0.5 * w[K] + w[Kp1]

  # numerical safety
  p[p < 0] <- 0
  p / sum(p)
}

# Given bin masses p (length 64) on intervals [Theta[i], Theta[i+1]],
# produce CDF values at the 65 boundary points (piecewise-linear CDF evaluated at boundaries)
hist_bins_to_cdf_on_boundaries <- function(p) {
  p <- as.numeric(p)
  p <- p / sum(p)
  c(0, cumsum(p))  # length 65, F(Theta[1])=0, F(Theta[65])=1
}

# Your old: cumsum(w) is the discrete step CDF (right-continuous) at boundary points.
# New: project to histogram first, then evaluate histogram CDF at boundary points.
w_to_cdf_histproj <- function(w) {
  p <- disc_to_hist_bins(w)
  hist_bins_to_cdf_on_boundaries(p)
}



# ---- Gibbs CDF evaluators + W1 vs true ----

# step CDF (right-continuous) of point-mass weights w at support points Theta
# returns F(x) for vector x
disc_cdf_eval <- function(x, Theta, w) {
  Theta <- as.numeric(Theta)
  w <- as.numeric(w); w <- w / sum(w)

  o <- order(Theta); Theta <- Theta[o]; w <- w[o]
  cw <- cumsum(w)

  # right-continuous: F(x)=P(Theta<=x)
  idx <- findInterval(x, Theta, rightmost.closed = TRUE, all.inside = FALSE)
  # idx in 0..length(Theta)
  out <- ifelse(idx <= 0, 0,
                ifelse(idx >= length(cw), 1, cw[idx]))
  pmin(1, pmax(0, out))
}

# piecewise-linear CDF for histogram with bin masses p on intervals [Theta[i], Theta[i+1]]
hist_cdf_eval <- function(x, Theta, p) {
  Theta <- as.numeric(Theta)
  p <- as.numeric(p); p <- p / sum(p)

  o <- order(Theta); Theta <- Theta[o]
  K <- length(Theta) - 1
  stopifnot(length(p) == K)

  # clamp x to [min,max]
  x0 <- pmin(max(Theta), pmax(min(Theta), x))
  i  <- findInterval(x0, Theta, rightmost.closed = TRUE, all.inside = TRUE)
  i  <- pmin(i, K)

  L  <- Theta[i]; R <- Theta[i + 1]
  t  <- ifelse(R > L, (x0 - L) / (R - L), 0)

  F_left <- c(0, cumsum(p))[i]
  F_left + t * p[i]
}

# Toggle: "disc" (step) or "hist" (disc->hist projection then piecewise-linear)
gibbs_cdf_eval <- function(x, Theta, w, mode = c("disc", "hist")) {
  mode <- match.arg(mode)
  if (mode == "disc") {
    disc_cdf_eval(x, Theta, w)
  } else {
    p <- disc_to_hist_bins(w)  # you already defined this
    hist_cdf_eval(x, Theta, p)
  }
}

# W1 on a fixed grid: ∫ |F-G| dx via trapezoid (you already have W1_on_grid in PT script)
W1_on_grid <- function(x, F, G) {
  D <- abs(F - G)
  sum(0.5 * (D[-1] + D[-length(D)]) * diff(x))
}

# Compute W1(Gibbs, true) for a single w (length 65)
W1_gibbs_vs_true <- function(w, Theta, F_true_fun, a, b, n_grid = 10000,
                             mode = c("disc", "hist")) {
  mode <- match.arg(mode)
  xg <- seq(a, b, length.out = n_grid)
  F_true_g <- F_true_fun(xg)
  F_gibbs_g <- gibbs_cdf_eval(xg, Theta, w, mode = mode)
  W1_on_grid(xg, F_gibbs_g, F_true_g)
}


# build truncated CDF from pdf via grid + cumulative trapezoid; returns a function F(x)
make_trunc_cdf_from_pdf <- function(pdf, a, b, n_grid = 20001) {
  xg <- seq(a, b, length.out = n_grid)
  fg <- pdf(xg)
  if (any(!is.finite(fg)) || any(fg < 0)) stop("pdf must be finite and nonnegative on [a,b].")

  dx  <- diff(xg)
  inc <- 0.5 * (fg[-1] + fg[-length(fg)]) * dx
  cg  <- c(0, cumsum(inc))               # unnormalized CDF at xg
  Z   <- cg[length(cg)]
  if (!is.finite(Z) || Z <= 0) stop("Bad normalization constant; check pdf or bounds.")

  Fg <- pmin(1, pmax(0, cg / Z))
  approxfun(xg, Fg, method = "linear", yleft = 0, yright = 1, rule = 2)
}


# ---- gather CDFs ----
run_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[sapply(run_dirs, keep_run_eq_seed)]

F_list       <- list()
labels       <- character()

lambda_target <- 0.001

# ---- W1 settings ----
a_trunc <- -3
b_trunc <-  3
n_grid_w1 <- 10000

# Prefer: reuse PT's true-CDF builder for consistency
# true_pdf <- function(x) dnorm(x, 0, 1)

#true_pdf bimodal
true_pdf <- function(x) { 0.5 * dnorm(x, -1.5, 0.5) + 0.5 * dnorm(x, 1.5, 0.5) }


F_true_fun <- make_trunc_cdf_from_pdf(true_pdf, a_trunc, b_trunc, n_grid_w1)

w1_gibbs_disc_vals <- numeric(0)
w1_gibbs_hist_vals <- numeric(0)


for (d in run_dirs) {
#gibbs_path <- #gibbs_out_run2_seed2_lambda0.001.rds

  rs <- parse_run_seed(d)
  if (is.null(rs)) next

  file_name  <- make_gibbs_filename(rs$run, rs$seed, lambda_target)
  gibbs_path <- file.path(d, file_name)
  g <- readRDS(gibbs_path)

  # Optional: enforce sim/scenario match if present in meta
  if (!is.null(g$meta$sim_id) && g$meta$sim_id != sim_id) next
  if (!is.null(g$meta$scenario_id) && g$meta$scenario_id != scenario_id) next

  Theta <- seq(-3, 3, length.out = 65)
  w_fin <- get_final_w(g)

  if (length(w_fin) != length(Theta)) {
    stop(sprintf("Length mismatch in %s: length(w)=%d, length(Theta)=%d",
                 gibbs_path, length(w_fin), length(Theta)))
  }


w1_disc <- W1_gibbs_vs_true(w_fin, Theta, F_true_fun, a_trunc, b_trunc,
                            n_grid = n_grid_w1, mode = "disc")
w1_hist <- W1_gibbs_vs_true(w_fin, Theta, F_true_fun, a_trunc, b_trunc,
                            n_grid = n_grid_w1, mode = "hist")
w1_gibbs_disc_vals <- c(w1_gibbs_disc_vals, w1_disc)
w1_gibbs_hist_vals <- c(w1_gibbs_hist_vals, w1_hist)

message(sprintf("%s | W1(Gibbs disc,true)=%.6g | W1(Gibbs hist,true)=%.6g",
                basename(d), w1_disc, w1_hist))



  F_list[[length(F_list) + 1]] <- w_to_cdf_histproj(w_fin) #w_to_cdf(w_fin)  # length K+1
  labels <- c(labels, basename(d))
}

if (length(F_list) == 0) stop("No runs found / no matching Gibbs rds files.")


cat(sprintf("\nGibbs W1 (disc) over %d runs: mean = %.6g, sd = %.6g\n",
            length(w1_gibbs_disc_vals),
            mean(w1_gibbs_disc_vals), sd(w1_gibbs_disc_vals)))

cat(sprintf("Gibbs W1 (disc->hist) over %d runs: mean = %.6g, sd = %.6g\n",
            length(w1_gibbs_hist_vals),
            mean(w1_gibbs_hist_vals), sd(w1_gibbs_hist_vals)))


# ---- plot: all Gibbs final-draw CDFs + true truncated normal ----
# Note: CDF vector is length K+1; x-grid needs same length.
# We’ll plot CDF at c(Theta[1], Theta) so it aligns with the leading 0.
x_cdf <- Theta

# F_true <- truncN01_cdf(x_cdf, a = -3, b = 3)
F_true <- trunc_bimodal_cdf(x_cdf, a = -3, b = 3)

F_mat  <- do.call(rbind, F_list)    # n_runs x (K+1)
F_mean <- colMeans(F_mat)

png(out_png, width = 900, height = 650)

plot(x_cdf, F_true, type = "n", ylim = c(0, 1),
     xlab = expression(beta), ylab = "CDF",
     main = sprintf("Gibbs final-draw CDFs (n=%d) | sim %d, scenario %d",
                    nrow(F_mat), sim_id, scenario_id))

# individual Gibbs runs (thin, semi-transparent)
gibbs_col <- grDevices::adjustcolor("gray30", alpha.f = 0.35)
for (i in seq_len(nrow(F_mat))) {
  lines(x_cdf, F_mat[i, ], col = gibbs_col, lwd = 1.3)
}

# mean across runs (optional but handy)
lines(x_cdf, F_mean, col = "black", lwd = 3)

# true truncated normal (thick, blue)
lines(x_cdf, F_true, col = "dodgerblue3", lwd = 3, lty = 2)

legend("topleft",
       legend = c("Gibbs final draw (each run)",
                  "Mean CDF across runs",
                  "Trunc N(0,1) on [-3,3]"),
       col = c(gibbs_col, "black", "dodgerblue3"),
       lwd = c(2, 3, 3),
       lty = c(1, 1, 2),
       bty = "n")

dev.off()
message("Wrote: ", out_png)

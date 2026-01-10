# pt_w1_vs_trueprior.R

# ---- knobs ----
root_dir    <- file.path(getwd(), "pt_output")
pt_rds_name <- "oracle_section4_all_results.rds"
sim_id      <- 1
scenario_id <- 2
burnin      <- 100000
n_use       <- 100000   # set NULL to use all post-burnin
n_grid      <- 10000    # grid size for W1
a_trunc     <- -3
b_trunc     <-  3

# ---- user-supplied true prior pdf on R (can be unnormalized) ----
true_pdf <- function(x) dnorm(x, 0, 1)  # N(0,1) (will be truncated+renormalized to [-3,3])

# ---- helpers ----

# posterior mean CDF on bin boundaries
post_mean_cdf_from_logpi <- function(pi_gibbs_log, burnin, n_use = NULL) {
  T <- nrow(pi_gibbs_log)
  stopifnot(burnin < T)
  idx <- (burnin + 1):T
  if (!is.null(n_use)) idx <- idx[seq_len(min(length(idx), n_use))]

  pi_post <- exp(pi_gibbs_log[idx, , drop = FALSE])
  pi_post <- pi_post / rowSums(pi_post)

  p_mean <- colMeans(pi_post)     # bin masses, length 64
  F_mean <- c(0, cumsum(p_mean))  # CDF at boundaries, length 65
  F_mean
}

# evaluate PT CDF at arbitrary x (piecewise linear between boundaries)
pt_cdf_eval <- function(x, a, F_a) {
  stopifnot(length(a) == length(F_a), length(a) >= 2)
  o <- order(a); a <- a[o]; F_a <- F_a[o]

  x0 <- pmin(max(a), pmax(min(a), x))
  i  <- findInterval(x0, a, rightmost.closed = TRUE, all.inside = TRUE)
  i  <- pmin(i, length(a) - 1)

  L  <- a[i];     R  <- a[i + 1]
  FL <- F_a[i];   FR <- F_a[i + 1]
  t  <- ifelse(R > L, (x0 - L) / (R - L), 0)

  FL + t * (FR - FL)
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

# W1 on a grid via trapezoid: ∫ |F-G| dx
W1_on_grid <- function(x, F, G) {
  D <- abs(F - G)
  sum(0.5 * (D[-1] + D[-length(D)]) * diff(x))
}

get_scenario <- function(obj, sim_id, scenario_id) {
  sc_name <- paste0("scenario", scenario_id)
  obj$results[[sim_id]]$scenario[[sc_name]]
}

# ---- main ----

# true CDF function (build once)
F_true <- make_trunc_cdf_from_pdf(true_pdf, a_trunc, b_trunc)

# grid for W1 (fixed once)
xg <- seq(a_trunc, b_trunc, length.out = n_grid)
F_true_g <- F_true(xg)


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

w1_vals <- numeric(0)

for (d in run_dirs) {
  pt_path <- file.path(d, pt_rds_name)
  if (!file.exists(pt_path)) next

  obj <- readRDS(pt_path)
  sc  <- get_scenario(obj, sim_id, scenario_id)

  a <- as.numeric(sc$a_vec)
  F_pt_a <- post_mean_cdf_from_logpi(sc$pi_gibbs_log, burnin, n_use)

  # evaluate PT CDF on the common grid
  F_pt_g <- pt_cdf_eval(xg, a, F_pt_a)

  w1 <- W1_on_grid(xg, F_pt_g, F_true_g)
  w1_vals <- c(w1_vals, w1)

  message(sprintf("%s | W1(PT, true prior) = %.6g", basename(d), w1))
}

if (length(w1_vals) == 0) stop("No runs found / no matching rds files.")

cat(sprintf("\nW1 over %d runs: mean = %.6g, sd = %.6g\n",
            length(w1_vals), mean(w1_vals), sd(w1_vals)))
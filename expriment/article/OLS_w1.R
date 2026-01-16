#!/usr/bin/env Rscript

# Computes OLS beta_hat from seed-001 data, clips to [-3,3],
# then computes grid-based W1 via integral of |ecdf - ecdf| on 10,000-point grid.

w1_cdf_grid_emp_vs_emp <- function(a, b, bounds = c(-3, 3), nbins = 10000) {
  xfine <- seq(bounds[1], bounds[2], length.out = nbins)
  dx <- xfine[2] - xfine[1]
  Fa <- ecdf(a)(xfine)
  Fb <- ecdf(b)(xfine)
  sum(abs(Fa - Fb)) * dx
}

clip_to_bounds <- function(z, bounds = c(-3, 3)) {
  pmin(pmax(z, bounds[1]), bounds[2])
}

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

# Defaults assume your generator saved to ./data/data_seed_001.rds
infile <- get_arg("--infile", file.path(getwd(), "data", "data_seed_001.rds"))
nbins  <- as.integer(get_arg("--nbins", 10000))
lo     <- as.numeric(get_arg("--lo", -3))
hi     <- as.numeric(get_arg("--hi",  3))
bounds <- c(lo, hi)

if (!file.exists(infile)) {
  stop("Input file not found: ", infile, "\n",
       "Pass --infile /path/to/data_seed_001.rds if needed.")
}

dat <- readRDS(infile)

X <- dat$X
y <- dat$y
beta_true <- dat$beta_true

# OLS: beta_hat = (X'X)^{-1} X'y (use solve on symmetric p x p system)
beta_ols <- solve(crossprod(X), crossprod(X, y))

# Clip OLS to [-3, 3] for comparability
beta_ols_clip  <- clip_to_bounds(beta_ols, bounds)
beta_true_clip <- clip_to_bounds(beta_true, bounds)

w1 <- w1_cdf_grid_emp_vs_emp(beta_true_clip, beta_ols_clip, bounds = bounds, nbins = nbins)

cat(sprintf("Seed: %s\n", dat$seed))
cat(sprintf("Input: %s\n", infile))
cat(sprintf("n=%d, p=%d, nbins=%d, bounds=[%.3f, %.3f]\n", dat$n, dat$p, nbins, bounds[1], bounds[2]))
cat(sprintf("W1_grid_ecdf (clipped to bounds): %.10f\n", w1))

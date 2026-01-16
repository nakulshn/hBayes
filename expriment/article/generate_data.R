#!/usr/bin/env Rscript

library(truncnorm)
library(MASS)   # for mvrnorm


# Outputs an .rds file labeled by seed.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed  <- as.integer(get_arg("--seed", 1))
outdir <- get_arg("--outdir", file.path(getwd(), "data"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

cat("Generating data with seed =", seed, "\n")
cat("Output dir:", outdir, "\n")

# ---------------------------
# Settings (same as your script)
# ---------------------------
n <- 4000
p <- 800
sigma_true <- 1.0

# ---------------------------
# Generate X -- iid N(0, 1/n) entries
# ---------------------------
X <- matrix(rnorm(n * p, mean = 0, sd = sqrt(1 / n)), n, p)



# # ---------------------------
# # Generate X -- iid multivariate Gaussian rows, with correlated coordinates.
# # ---------------------------

# rho <- 0.7  # correlation level

# # Toeplitz covariance
# Sigma <- toeplitz(rho^(0:(p - 1)))

# # Generate iid rows
# X <- MASS::mvrnorm(
#   n = n,
#   mu = rep(0, p),
#   Sigma = Sigma
# ) / sqrt(n)


# # ---------------------------
# # Generate betas: Basic Gaussian N(0, 1) truncated at [-3, 3].
# # ---------------------------

# beta <- rtruncnorm(
#   p,
#   a = -3, b = 3,
#   mean = 0,
#   sd = 1.0
# )

# # ---------------------------
# # Generate betas: truncated continuous spike-and-slab Seeds 1, 2, 3
# # ---------------------------
# pi_slab <- 0.5
# tau0 <- 0.04
# tau1 <- 1.0

# mix <- rbinom(p, 1, pi_slab)
# beta <- numeric(p)

# beta[mix == 0] <- rtruncnorm(
#   sum(mix == 0),
#   a = -3, b = 3,
#   mean = 0,
#   sd = tau0
# )

# beta[mix == 1] <- rtruncnorm(
#   sum(mix == 1),
#   a = -3, b = 3,
#   mean = 0,
#   sd = tau1
# )


# ---------------------------
# Generate betas: bimodal Seed 4,5,6
# ---------------------------
pi_left_mode <- 0.5
tau0 <- 0.2
tau1 <- 0.2

mix <- rbinom(p, 1, pi_left_mode)
beta <- numeric(p)

beta[mix == 0] <- rtruncnorm(
  sum(mix == 0),
  a = -3, b = 3,
  mean = -1.5,
  sd = tau0
)

beta[mix == 1] <- rtruncnorm(
  sum(mix == 1),
  a = -3, b = 3,
  mean = 1.5,
  sd = tau1
)


# ---------------------------
# Generate y
# ---------------------------
y <- as.numeric(X %*% beta + sigma_true * rnorm(n))

# ---------------------------
# Save (filename includes seed)
# ---------------------------
outfile <- file.path(outdir, sprintf("data_seed_%03d.rds", seed))

saveRDS(
  list(
    seed = seed,
    n = n,
    p = p,
    sigma_true = sigma_true,
    X = X,
    beta_true = beta,
    y = y
  ),
  file = outfile
)

cat("Saved:", outfile, "\n")

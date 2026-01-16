# run_Gibbs_on_PT.R

library(pracma)   # only if you want logspace etc (not required here)
library(CVXR)     # required by spline.optimize()

source("./Gibbs.R")  # <-- contains Gibbs(), stable.softmax(), etc.

# ---- CLI args: --seed, --sim, --scenario ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed   <- as.integer(get_arg("--seed", 1))


# ---- Gibbs tuning args (optional) ----
K_grid            <- as.integer(get_arg("--K", 65))
lambda            <- as.numeric(get_arg("--lambda", 1e-3))
eta_w             <- as.numeric(get_arg("--eta_w", 1.0))
max_iter          <- as.integer(get_arg("--max_iter", 10000))
gibbs_inner_iters <- as.integer(get_arg("--gibbs_inner_iters", 100))
cat("Gibbs runner:",
    "seed =", seed, "\n")


# ---------------------------
# Load data
# ---------------------------
datadir <- file.path(getwd(), "data")
datafile <- file.path(datadir, sprintf("data_seed_%03d.rds", seed))
if (!file.exists(datafile)) {
  stop("Data file not found: ", datafile, "\n",
       "Did you generate it with the same --seed and datadir?")
}

d <- readRDS(datafile)

# Expect these fields from your generator:
X <- d$X
y <- d$y

sigma <- d$sigma_true
beta_true <- d$beta_true

cat("Loaded:", datafile, "\n",
    "n =", nrow(X), "p =", ncol(X), "sigma =", sigma, "\n")

bounds <- c(-3, 3)

# ---------------------------
# grid
# ---------------------------
grid <- seq(bounds[1], bounds[2], length.out = K_grid)


# ---------------------------
# Run Gibbs
# ---------------------------
set.seed(seed + 100)

gibbs_out <- Gibbs(
  X = X,
  y = y,
  sigma = sigma,
  grid = grid,
  lambda = lambda,
  eta.w = eta_w,
  beta.true = beta_true,   # <-- NEW
  bounds = bounds,         # <-- NEW (for W1 CDF range)
  max.iter = max_iter,
  verbose = TRUE,
  gibbs.inner.iters = gibbs_inner_iters
)

# ---- Plot CDFs: final w vs empirical beta_true ----
w_final <- gibbs_out$w

xfine <- seq(bounds[1], bounds[2], length.out = 10000)
F_true <- ecdf(beta_true)(xfine)
F_est  <- cdf_from_w_grid(w_final, grid, xfine)

plot(xfine, F_true, type = "l",
     xlab = "x", ylab = "CDF",
     main = sprintf("CDFs after %d iters (seed=%d)", max_iter, seed))
lines(xfine, F_est, lty = 2)

legend("topleft",
       legend = c("Empirical CDF (true betas)", "CDF from final w (grid)"),
       lty = c(1, 2), bty = "n")

# Optional: save to file
outdir <- file.path(getwd(), "gibbs_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pngfile <- file.path(outdir, sprintf("cdf_compare_seed_%03d.png", seed))
png(pngfile, width = 900, height = 700)
plot(xfine, F_true, type = "l",
     xlab = "x", ylab = "CDF",
     main = sprintf("CDFs after %d iters (seed=%d)", max_iter, seed))
lines(xfine, F_est, lty = 2)
legend("topleft",
       legend = c("Empirical CDF (true betas)", "CDF from final w (grid)"),
       lty = c(1, 2), bty = "n")
dev.off()

cat("Saved CDF plot:", pngfile, "\n")



w = gibbs_out$w
phi = gibbs_out$hist[[max_iter+1]]$phi
tmp.out <- Gibbs(
    X = X,
    y = y,
    sigma = sigma,
    grid = grid,
    lambda = 0,
    eta.w = 1.0,
    gibbs.inner.iters = 1000000,
    w.init = w,
    phi.init = phi,
    burn.iters = 0,
    beta.true = beta_true,   # <-- NEW
    bounds = bounds,         # <-- NEW (for W1 CDF range)
    max.iter = 200000,
    verbose = TRUE
  )

beta.samples = sapply(seq(1,200001,1), function(it) tmp.out$hist[[it]]$phi)
gibbs_out$beta.samples = beta.samples



# ---------------------------
# Save output (filename includes seed)
# ---------------------------
outdir <- file.path(getwd(), "gibbs_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("gibbs_out_seed_%03d.rds", seed))

saveRDS(
  gibbs_out,
  file = outfile
)

cat("Saved:", outfile, "\n")

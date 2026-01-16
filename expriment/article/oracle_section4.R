library(NPBayes)
library(glmnet)
source("C:/Users/nakul/Downloads/hBayes/R/normal.gibbs.R")

# ---- CLI args: --seed, --sim, --scenario ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed   <- as.integer(get_arg("--seed", 1))

set.seed(seed + 200)

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

# ---------------------------
# Settings you asked for
# ---------------------------
n.gibbs <- 200000
burnin  <- 10


# ridge
fit.glmnet <- cv.glmnet(X, y, alpha = 0, intercept = FALSE, family = "gaussian")
beta.ridge <- as.numeric(coef(fit.glmnet, s="lambda.min"))[-1]

# prior support bounds
a.min <- -3
a.max <- 3


# hBayes Gibbs (C++ backend)
hBeta <- gibbs.normal.fixed.sigma(
  n.gibbs, y, X,
  sigma,
  c(a.min, a.max),
  6,
  as.vector(beta.ridge),
  cpp = TRUE,
  beta.true = beta_true
)



# Build plot grid
grid65 <- seq(a.min, a.max, length.out = 65)
w_grid <- project_hist_to_grid(hBeta$a.vec, hBeta$w_mean, grid65)

xfine <- seq(a.min, a.max, length.out = 10000)
F_true <- ecdf(beta_true)(xfine)
F_est  <- cdf_from_w_grid(w_grid, grid65, xfine)

plot(xfine, F_true, type="l", xlab="x", ylab="CDF",
     main="CDF: empirical true betas vs projected PT prior (65-grid)")
lines(xfine, F_est, lty=2)
legend("topleft", c("Empirical CDF (true betas)", "Projected PT -> grid CDF"),
       lty=c(1,2), bty="n")


outdir <- file.path(getwd(), "pt_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pngfile <- file.path(outdir, sprintf("pt_cdf_compare_seed_%03d.png", seed))
png(pngfile, width = 900, height = 700)

plot(xfine, F_true, type = "l",
     xlab = "x", ylab = "CDF",
     main = sprintf("Polya-tree mean prior CDF vs true (seed=%d)", seed))
lines(xfine, F_est, lty = 2)
legend("topleft",
       legend = c("Empirical CDF (true betas)", "CDF of mean prior (PT)"),
       lty = c(1, 2), bty = "n")

dev.off()
cat("Saved CDF plot:", pngfile, "\n")



# ---------------------------
# Save everything
# ---------------------------
outdir <- file.path(getwd(), "pt_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

saveRDS(list(
  settings = list(n.gibbs=n.gibbs, burnin=burnin),
  hBeta = hBeta
), file = file.path(outdir, sprintf("pt_out_seed_%03d.rds", seed)))

cat("\nDone.\nOutputs in:", outdir, "\n")

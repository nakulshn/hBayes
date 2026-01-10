## ======================================================================
##  Libraries
## ======================================================================

library(NPBayes)
library(glmnet)
library(cmdstanr)
library(brms)
options(brms.backend = "cmdstanr")
library(parallel)

## ======================================================================
##  Helper to silence functions
## ======================================================================

quiet <- function(expr) {
    res <- NULL
    suppressMessages(
        suppressWarnings(
            capture.output(
                res <- eval.parent(substitute(expr)),
                file = NULL
            )
        )
    )
    res
}

## ======================================================================
##  Global settings
## ======================================================================

n.gibbs <- 200
n.brms  <- 200
sim     <- 30L        # total number of simulations
burnin  <- ceiling(0.2 * n.gibbs)
n       <- 400
p       <- 10 * 8
verbose = TRUE
## ---- read sim_start and sim_end from command line ---------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
    sim_start <- as.integer(args[1])
    sim_end   <- as.integer(args[2])
} else {
    sim_start <- 1L
    sim_end   <- 8
}

if (is.na(sim_start) || is.na(sim_end)) {
    stop("sim_start and sim_end must be integers.")
}
if (sim_start < 1L || sim_end > sim || sim_start > sim_end) {
    stop("sim_start and sim_end must satisfy 1 <= sim_start <= sim_end <= sim.")
}

sims_to_run <- sim_start:sim_end

cat("Running sims", sim_start, "to", sim_end, "of", sim, "\n")

## ---- priors -----------------------------------------------------------

prior.dense <- prior(
    R2D2(
        mean_R2  = 0.9,
        prec_R2  = 1
    ),
    class = "b"
)
prior.sparse <- prior(
    R2D2(
        mean_R2  = 0.1,
        prec_R2  = 1
    ),
    class = "b"
)

## ======================================================================
##  Data setup
## ======================================================================

set.seed(1)
X.mat <- array(
    rnorm(n * p, mean = 0, sd = sqrt(1 / n)),
    dim = c(n, p)
)

set.seed(1)
betas      <- array(dim = c(p, 3))
betas[, 1] <- c(rep(-10, p/8), rep(10, p/8), rep(0, p*6/8))
betas[, 2] <- rnorm(p, mean = 3, sd = 4)
betas[, 3] <- rbinom(p, 1, 0.5) * rnorm(p, mean = 7, sd = 1)

inv.log.odds <- function(Xbeta) {
    1 / (1 + exp(-Xbeta))
}

## ---- standardize X only for Stan/brms ---------------------------------

X.mu <- colMeans(X.mat)
X.sd <- apply(X.mat, 2, sd)
X.sd[ X.sd == 0 ] <- 1  # avoid division by zero just in case

X.std <- sweep(X.mat, 2, X.mu, "-")
X.std <- sweep(X.std, 2, X.sd, "/")

## ======================================================================
##  Core allocation logic (Lunarc-friendly)
## ======================================================================

cores_per_brms <- 2L

## As requested:
cores_total <- min(cores_per_brms * length(sims_to_run)*3, parallel::detectCores())
n_jobs      <- max(1L, floor(cores_total / cores_per_brms))

options(mc.cores = n_jobs)

Sys.setenv(OMP_NUM_THREADS = "1")    # keep BLAS/OpenMP single-threaded
cat("Detected cores:", parallel::detectCores(),
    "| cores_total used:", cores_total,
    "| mclapply jobs:", n_jobs,
    "| cores_per_brms:", cores_per_brms, "\n")

## ======================================================================
##  Method names (for column order in printing)
## ======================================================================

method_names <- c(
    "MLE",
    "Oracle",
    "hBeta",
    "Lasso",
    "Ridge",
    "MLE_adj",
    "R2D2_sparse",
    "R2D2_dense"
)
sim.mse.mat <- array(
    NA_real_,
    dim = c(length(method_names), sim, 3L),
    dimnames = list(
        method = method_names,
        sim    = paste0("sim", seq_len(sim)),
        beta   = paste0("beta_set", 1:3)
    )
)
## Print CSV header once (results lines come from run_one_ij)
cat(
    paste(c("i", "j", method_names), collapse = ","),
    "\n",
    sep = ""
)

## ======================================================================
##  One (i, j) run — prints its own result line
## ======================================================================

run_one_ij <- function(i, j,
                       X.mat, X.std, X.sd, betas,
                       n.gibbs, burnin,
                       prior.sparse, prior.dense,
                       n.brms,
                       verbose=FALSE) {

    ## Make results reproducible per (i,j)
    set.seed(100000 * i + 100 * j + 1L)
    cat('starting (',i,j,')\n')
    n <- nrow(X.mat)

    ## Simulate response using ORIGINAL X
    y <- rbinom(n, 1, inv.log.odds(X.mat %*% betas[, j]))

    if(verbose){
        ## Oracle Gibbs (original scale)
        cat('starting oracle gibbs for (',i,j,')\n')
        beta.oracle.draws <-
            gibbs.logistic.permute.oracle(
                n.gibbs, y, X.mat, betas[, j], FALSE
            )
    }else{
        ## Oracle Gibbs (original scale)
        beta.oracle.draws <- quiet(
            gibbs.logistic.permute.oracle(
                n.gibbs, y, X.mat, betas[, j], FALSE
            )
        )
    }

    beta.oracle <- colMeans(
        beta.oracle.draws[burnin:n.gibbs, , drop = FALSE]
    )

    ## MLE (original scale)
    fit      <- quiet(glm(y ~ X.mat + 0, family = binomial))
    beta.mle <- coef(fit)
    beta.mle.adj <- beta.mle / 1.499

    ## Lasso (original scale)
    fit.glmnet <- quiet(
        cv.glmnet(
            X.mat, y,
            alpha        = 1,
            intercept    = FALSE,
            family       = "binomial",
            type.measure = "class"
        )
    )
    beta.lasso <- as.numeric(coef(fit.glmnet, s = "lambda.min"))[-1]

    ## Ridge
    fit.glmnet <- quiet(
        cv.glmnet(
            X.mat, y,
            alpha        = 0,
            intercept    = FALSE,
            family       = "binomial",
            type.measure = "class"
        )
    )
    beta.ridge <- as.numeric(coef(fit.glmnet, s = "lambda.min"))[-1]

    ## Hierarchical Beta (NPBayes, original scale)
    a.min <- min(c(-24, beta.mle - 0.5))
    a.max <- max(c( 24, beta.mle + 0.5))


    if(verbose){

        cat('starting gibbs for (',i,j,')\n')
        hBeta <-
            gibbs.logistic(
                n.gibbs,
                y,
                X.mat,
                c(a.min, a.max),
                6,
                beta = as.vector(beta.ridge),
                cpp  = TRUE
            )

    }else{
        hBeta <-
            gibbs.logistic(
                n.gibbs,
                y,
                X.mat,
                c(a.min, a.max),
                6,
                beta = as.vector(beta.ridge),
                cpp  = TRUE
            )
    }
    beta.hBeta <- colMeans(
        hBeta$beta.gibbs[burnin:n.gibbs, , drop = FALSE]
    )

    ## brms (Stan) – use STANDARDIZED X, then transform back
    dat1 <- data.frame(y = y, X.std)
    if(verbose){

        cat('starting brm sparse for (',i,j,')\n')
        fit.R2D2.sparse <- brm(
                formula = y ~ -1 + .,
                data    = dat1,
                family  = bernoulli(link = "logit"),
                prior   = prior.sparse,
                chains  = 4,
                cores   = cores_per_brms,
                iter    = n.brms,
                seed    = 1,
                refresh = 0
            )
    }else{
        fit.R2D2.sparse <- quiet(
            brm(
                formula = y ~ -1 + .,
                data    = dat1,
                family  = bernoulli(link = "logit"),
                prior   = prior.sparse,
                chains  = 4,
                cores   = cores_per_brms,
                iter    = n.brms,
                seed    = 1,
                refresh = 0
            )
        )
    }
    beta.R2D2.sparse.std <- fixef(fit.R2D2.sparse)[, 1]
    rm(fit.R2D2.sparse)
    beta.R2D2.sparse <- as.numeric(beta.R2D2.sparse.std) / X.sd

    if(verbose){

        cat('starting brm dense for (',i,j,')\n')
        fit.R2D2.dense <- brm(
            formula = y ~ -1 + .,
            data    = dat1,
            family  = bernoulli(link = "logit"),
            prior   = prior.dense,
            chains  = 4,
            cores   = cores_per_brms,
            iter    = n.brms,
            seed    = 1,
            refresh = 0
        )
    }else{
        fit.R2D2.dense <- quiet(
            brm(
                formula = y ~ -1 + .,
                data    = dat1,
                family  = bernoulli(link = "logit"),
                prior   = prior.dense,
                chains  = 4,
                cores   = cores_per_brms,
                iter    = n.brms,
                seed    = 1,
                refresh = 0
            )
        )
    }
    beta.R2D2.dense.std <- fixef(fit.R2D2.dense)[, 1]
    beta.R2D2.dense <- as.numeric(beta.R2D2.dense.std) / X.sd
    rm(fit.R2D2.dense)

    ## MSEs (original scale)
    mse_vec <- c(
        sqrt(mean((beta.mle         - betas[, j])^2)),
        sqrt(mean((beta.oracle      - betas[, j])^2)),
        sqrt(mean((beta.hBeta       - betas[, j])^2)),
        sqrt(mean((beta.lasso       - betas[, j])^2)),
        sqrt(mean((beta.ridge       - betas[, j])^2)),
        sqrt(mean((beta.mle.adj     - betas[, j])^2)),
        sqrt(mean((beta.R2D2.sparse - betas[, j])^2)),
        sqrt(mean((beta.R2D2.dense  - betas[, j])^2))
    )

    ## ---- PRINT RESULT LINE HERE -----------------------------------------
    cat(
        paste(
            c(
                i,
                j,
                format(mse_vec, digits = 6)
            ),
            collapse = ","
        ),
        "\n",
        sep = ""
    )
    flush.console()


    list(i = i, j = j, mse = mse_vec)
}

## ======================================================================
##  Run sims: parallel over i and j
## ======================================================================

combo_ij <- expand.grid(
    i = sims_to_run,
    j = 1:3
)

## We don't need the returned list; printing is done inside run_one_ij

res_list <- mclapply(
    X = seq_len(nrow(combo_ij)),
    FUN = function(k) {
        i <- combo_ij$i[k]
        j <- combo_ij$j[k]

        run_one_ij(
            i = i, j = j,
            X.mat   = X.mat,
            X.std   = X.std,
            X.sd    = X.sd,
            betas   = betas,
            n.gibbs = n.gibbs,
            burnin  = burnin,
            prior.sparse = prior.sparse,
            prior.dense  = prior.dense,
            n.brms  = n.brms,
            verbose=verbose
        )
    }
)


for (res in res_list) {
    sim.mse.mat[, res$i, res$j] <- res$mse
}

save(
    sim.mse.mat,
    file = sprintf("sim_mse_mat_%02d_%02d.RData", sim_start, sim_end)
)
cat("Done.\n")

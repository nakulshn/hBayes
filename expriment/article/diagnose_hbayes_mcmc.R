# diagnose_hbayes_mcmc.R
# Usage:
#   Rscript diagnose_hbayes_mcmc.R path/to/oracle_section4_all_results.rds

suppressPackageStartupMessages({
  library(coda)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide path to the .rds file as the first argument.")
rds_path <- args[1]

obj <- readRDS(rds_path)

# ---- helpers ----

safe_mcmc <- function(x) {
  # accept vector or matrix; return coda::mcmc
  if (is.null(dim(x))) return(mcmc(as.numeric(x)))
  mcmc(x)
}

thin_matrix <- function(mat, thin = 1) {
  if (thin <= 1) return(mat)
  mat[seq(1, nrow(mat), by = thin), , drop = FALSE]
}

pick_beta_indices <- function(beta_true = NULL, beta_lse = NULL, k = 20) {
  # choose a representative set of coefficients to diagnose
  # preference: largest |beta_true|, else largest |beta_lse|, else random
  p <- if (!is.null(beta_true)) length(beta_true) else if (!is.null(beta_lse)) length(beta_lse) else NA_integer_
  if (is.na(p)) return(integer(0))

  if (!is.null(beta_true)) {
    ord <- order(abs(beta_true), decreasing = TRUE)
    return(unique(c(ord[1:min(k, p)], sample.int(p, min(k, p)))))
  }
  if (!is.null(beta_lse)) {
    ord <- order(abs(beta_lse), decreasing = TRUE)
    return(unique(c(ord[1:min(k, p)], sample.int(p, min(k, p)))))
  }
  unique(sample.int(p, min(2*k, p)))
}

running_mean <- function(x) cumsum(x) / seq_along(x)

geweke_grid <- function(x_mcmc, burnins) {
  # returns a data.frame with burnin and frac_bad (|z|>1.96) over selected params
  # x_mcmc: matrix iterations x params
  out <- data.frame(burnin = burnins, frac_bad = NA_real_)
  for (i in seq_along(burnins)) {
    b <- burnins[i]
    if (b >= nrow(x_mcmc) - 50) { out$frac_bad[i] <- NA; next }
    xm <- mcmc(x_mcmc[(b+1):nrow(x_mcmc), , drop = FALSE])
    z <- tryCatch(geweke.diag(xm)$z, error = function(e) rep(NA_real_, ncol(x_mcmc)))
    out$frac_bad[i] <- mean(abs(z) > 1.96, na.rm = TRUE)
  }
  out
}

mean_stability_grid <- function(x_mcmc, burnins) {
  # compares posterior mean after burnin b to posterior mean after a large burnin baseline
  # returns avg relative L2 difference across params
  T <- nrow(x_mcmc)
  b0 <- burnins[length(burnins)]  # baseline: last burnin in grid
  if (b0 >= T - 50) b0 <- floor(0.5 * T)
  mu0 <- colMeans(x_mcmc[(b0+1):T, , drop = FALSE])

  out <- data.frame(burnin = burnins, rel_l2 = NA_real_)
  denom <- sqrt(sum(mu0^2)) + 1e-12
  for (i in seq_along(burnins)) {
    b <- burnins[i]
    if (b >= T - 50) { out$rel_l2[i] <- NA; next }
    mu <- colMeans(x_mcmc[(b+1):T, , drop = FALSE])
    out$rel_l2[i] <- sqrt(sum((mu - mu0)^2)) / denom
  }
  out
}

make_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

# ---- locate draws inside your saved structure ----
# Expected structure from your earlier script:
# obj$results[[s]]$scenario[["scenario2"]]$hBayes_raw (or hBayes_raw itself)
# We'll diagnose *all sims* and scenarios found.

out_dir <- make_dir(file.path(dirname(rds_path), "mcmc_diagnostics"))
cat("Diagnostics output dir:", out_dir, "\n")

all_sims <- obj$results
if (is.null(all_sims) || length(all_sims) == 0) stop("Couldn't find obj$results in the RDS.")

# Gather available (sim,scenario) combos
pairs <- list()
for (s in seq_along(all_sims)) {
  sc_list <- all_sims[[s]]$scenario
  if (is.null(sc_list)) next
  for (nm in names(sc_list)) {
    if (!is.list(sc_list[[nm]])) next
    hb <- sc_list[[nm]]$hBayes_raw
    if (!is.null(hb) && is.list(hb)) {
      pairs[[length(pairs)+1]] <- list(sim = s, scenario_name = nm, entry = sc_list[[nm]])
    }
  }
}
if (length(pairs) == 0) stop("No hBayes_raw found under obj$results[[s]]$scenario[[...]]$hBayes_raw")

cat("Found", length(pairs), "sim×scenario runs with hBayes_raw.\n")

# ---- main diagnostics per run ----
for (pp in pairs) {
  s  <- pp$sim
  nm <- pp$scenario_name
  entry <- pp$entry
  hb <- entry$hBayes_raw

  beta_draws <- hb$beta.gibbs
  pi_draws   <- hb$pi.gibbs
  sigma_draws <- hb$sigma.gibbs  # may be NULL (fixed sigma run)

  if (is.null(beta_draws) || is.null(pi_draws)) next

  T <- nrow(beta_draws)
  p <- ncol(beta_draws)
  I <- ncol(pi_draws)

  # Choose a small set of betas to inspect closely
  beta_idx <- pick_beta_indices(beta_true = entry$beta_true, beta_lse = entry$beta_lse, k = 20)
  beta_sel <- beta_draws[, beta_idx, drop = FALSE]

  # For pi, inspect a handful of bins (e.g. quartiles)
  pi_idx <- unique(pmax(1, pmin(I, round(c(1, I*0.25, I*0.5, I*0.75, I)))))
  pi_sel <- pi_draws[, pi_idx, drop = FALSE]

  # Candidate burn-ins to evaluate
  burn_grid <- unique(pmax(0, round(seq(0, 0.5*T, length.out = 11))))
  burn_grid <- burn_grid[burn_grid < (T - 50)]

  tag <- sprintf("sim%02d_%s_T%d", s, nm, T)
  run_dir <- make_dir(file.path(out_dir, tag))

  # ----- ESS / autocorr on selected betas -----
  m_beta <- safe_mcmc(beta_sel)
  ess_beta <- effectiveSize(m_beta)

  # crude rule: if median ESS is low, you probably need more iterations
  ess_summary <- summary(as.numeric(ess_beta))
  writeLines(capture.output(print(ess_summary)), con = file.path(run_dir, "beta_ess_summary.txt"))

  # ----- Geweke vs burn-in (selected betas) -----
  g_beta <- geweke_grid(beta_sel, burn_grid)
  ms_beta <- mean_stability_grid(beta_sel, burn_grid)

  write.csv(g_beta,  file.path(run_dir, "geweke_beta_by_burnin.csv"), row.names = FALSE)
  write.csv(ms_beta, file.path(run_dir, "mean_stability_beta_by_burnin.csv"), row.names = FALSE)

  # ----- also do it for selected pi components -----
  g_pi <- geweke_grid(pi_sel, burn_grid)
  ms_pi <- mean_stability_grid(pi_sel, burn_grid)
  write.csv(g_pi,  file.path(run_dir, "geweke_pi_by_burnin.csv"), row.names = FALSE)
  write.csv(ms_pi, file.path(run_dir, "mean_stability_pi_by_burnin.csv"), row.names = FALSE)

  # ----- plots -----
  png(file.path(run_dir, "trace_beta_selected.png"), width = 1200, height = 800)
  par(mfrow = c(4, 5), mar = c(3, 3, 2, 1))
  for (j in seq_len(ncol(beta_sel))) {
    plot(beta_sel[, j], type = "l", xlab = "iter", ylab = "", main = paste0("beta[", beta_idx[j], "]"))
  }
  dev.off()

  png(file.path(run_dir, "running_mean_beta_selected.png"), width = 1200, height = 800)
  par(mfrow = c(4, 5), mar = c(3, 3, 2, 1))
  for (j in seq_len(ncol(beta_sel))) {
    rm <- running_mean(beta_sel[, j])
    plot(rm, type = "l", xlab = "iter", ylab = "", main = paste0("run mean beta[", beta_idx[j], "]"))
  }
  dev.off()

  png(file.path(run_dir, "burnin_geweke_beta.png"), width = 900, height = 600)
  plot(g_beta$burnin, g_beta$frac_bad, type = "b", xlab = "burn-in", ylab = "fraction |z|>1.96",
       main = paste0(tag, "  Geweke (betas)"))
  abline(h = 0.05, lty = 2)
  dev.off()

  png(file.path(run_dir, "burnin_mean_stability_beta.png"), width = 900, height = 600)
  plot(ms_beta$burnin, ms_beta$rel_l2, type = "b", xlab = "burn-in", ylab = "rel L2 diff vs baseline",
       main = paste0(tag, "  Mean stability (betas)"))
  dev.off()

  png(file.path(run_dir, "burnin_geweke_pi.png"), width = 900, height = 600)
  plot(g_pi$burnin, g_pi$frac_bad, type = "b", xlab = "burn-in", ylab = "fraction |z|>1.96",
       main = paste0(tag, "  Geweke (pi selected bins)"))
  abline(h = 0.05, lty = 2)
  dev.off()

  png(file.path(run_dir, "burnin_mean_stability_pi.png"), width = 900, height = 600)
  plot(ms_pi$burnin, ms_pi$rel_l2, type = "b", xlab = "burn-in", ylab = "rel L2 diff vs baseline",
       main = paste0(tag, "  Mean stability (pi selected bins)"))
  dev.off()

  # sigma diagnostics if present
  if (!is.null(sigma_draws)) {
    png(file.path(run_dir, "trace_sigma.png"), width = 900, height = 600)
    plot(sigma_draws, type = "l", xlab = "iter", ylab = "sigma", main = paste0(tag, " sigma trace"))
    dev.off()

    m_sig <- safe_mcmc(sigma_draws)
    ess_sig <- effectiveSize(m_sig)
    writeLines(paste("ESS(sigma):", as.numeric(ess_sig)), con = file.path(run_dir, "sigma_ess.txt"))
  }

  # ----- quick textual suggestion -----
  # heuristic: pick smallest burnin where geweke frac_bad <= 0.05 AND mean stability <= 0.02
  pick_burn <- function(g, ms) {
    ok <- which(!is.na(g$frac_bad) & !is.na(ms$rel_l2) & g$frac_bad <= 0.05 & ms$rel_l2 <= 0.02)
    if (length(ok) == 0) return(NA_integer_)
    g$burnin[min(ok)]
  }
  b_suggest_beta <- pick_burn(g_beta, ms_beta)
  b_suggest_pi   <- pick_burn(g_pi, ms_pi)

  msg <- c(
    paste0("Run: ", tag),
    paste0("T = ", T, ", p = ", p, ", I = ", I),
    paste0("Median ESS(beta selected) = ", round(median(as.numeric(ess_beta)), 1),
           " (min=", round(min(as.numeric(ess_beta)), 1), ", max=", round(max(as.numeric(ess_beta)), 1), ")"),
    paste0("Suggested burn-in (beta): ", b_suggest_beta),
    paste0("Suggested burn-in (pi):   ", b_suggest_pi),
    "Heuristic notes:",
    "- If ESS is low (<~200 for key params), consider running longer or thinning only for storage (not for inference).",
    "- If Geweke frac_bad stays high even after large burn-in, you likely need more iterations / better mixing."
  )
  writeLines(msg, con = file.path(run_dir, "summary.txt"))
}

cat("\nDone. See:", out_dir, "\n")

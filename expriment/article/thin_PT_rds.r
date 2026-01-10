# thin_PT_rds.R
# Usage:
# Rscript thin_PT_rds.R --seed 1 --runid 1

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

seed  <- as.integer(get_arg("--seed", 1))
runid <- as.integer(get_arg("--runid", seed))

run_dir <- file.path(getwd(), "bimodal_oracle_section4_out", paste0("run_", runid, "_seed_", seed))
big_path  <- file.path(run_dir, "oracle_section4_all_results.rds")
thin_path <- file.path(run_dir, "oracle_section4_thin_results.rds")

if (!file.exists(big_path)) stop("Missing: ", big_path)

obj <- readRDS(big_path)

# keep only what downstream needs
thin <- list(
  X = obj$X,
  results = lapply(obj$results, function(sim) {
    # sim$scenario is (apparently) a named list like scenario1, scenario2, ...
    list(
      scenario = lapply(sim$scenario, function(scen) {
        list(
          y = scen$y,
          beta_true  = scen$beta_true,
          beta_ridge = scen$beta_ridge,
          beta_lasso = scen$beta_lasso,
          beta_lse   = scen$beta_lse
        )
      })
    )
  })
)

saveRDS(thin, thin_path)

cat("Wrote thin RDS:\n  ", thin_path, "\n", sep = "")
cat("Original size:", file.info(big_path)$size, "bytes\n")
cat("Thin size:    ", file.info(thin_path)$size, "bytes\n")

rm(obj, thin)
invisible(gc())

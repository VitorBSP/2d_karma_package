# Summarize the Monte Carlo results and compare them with Tables 1-3.
#
# Usage: Rscript analysis/simulations/summarize_monte_carlo.R
#          [--input=analysis/output/simulations] [--scenario_set=paper]
#
# For each parameter: mean of the estimates, bias (mean - truth), MSE
# (mean of squared errors), sample skewness and kurtosis (moment
# estimators, kurtosis 3 under normality), and coverage of the nominal 95%
# Wald intervals estimate +/- z_0.975 * se for the three standard error
# types. Replications whose fit failed or did not converge are excluded and
# counted. MSE here is mean((estimate - truth)^2); the original scripts used
# bias^2 + sd^2 with the n - 1 denominator, which differs by O(1/R).

if (!file.exists("analysis/R/common.R")) {
  stop("Run this script from the repository root.")
}
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "simulations", "config.R"))
opts <- parse_args(list(input = file.path("analysis", "output",
                                          "simulations"),
                        scenario_set = "paper"))
dir <- file.path(opts$input, opts$scenario_set)
files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0L) stop("No results found in ", dir)
z <- stats::qnorm(1 - (1 - confidence_level) / 2)

summaries <- lapply(files, function(f) {
  res <- readRDS(f)
  ok <- is.na(res$error) & res$convergence %in% 0L &
    stats::complete.cases(res$estimate)
  est <- res$estimate[ok, , drop = FALSE]
  coverage <- function(se) {
    se <- se[ok, , drop = FALSE]
    lo <- est - z * se
    hi <- est + z * se
    truth <- matrix(res$truth, nrow(est), ncol(est), byrow = TRUE)
    colMeans(truth >= lo & truth <= hi, na.rm = TRUE)
  }
  scenario <- sub("_.*$", "", basename(f))
  data.frame(
    scenario = scenario, size = res$config$size,
    parameter = res$parameter_names, truth = res$truth,
    replications = length(ok), used = sum(ok),
    failed = sum(!is.na(res$error)),
    not_converged = sum(is.na(res$error) & !res$convergence %in% 0L),
    mean = colMeans(est),
    bias = colMeans(est) - res$truth,
    mse = colMeans(sweep(est, 2, res$truth)^2),
    skewness = apply(est, 2, sample_skewness),
    kurtosis = apply(est, 2, sample_kurtosis),
    coverage_expected = coverage(res$se_expected),
    coverage_observed = coverage(res$se_observed),
    coverage_legacy = coverage(res$se_legacy),
    row.names = NULL
  )
})
summary_table <- do.call(rbind, summaries)
summary_table <- summary_table[order(summary_table$scenario,
                                     summary_table$size), ]
utils::write.csv(summary_table, file.path(dir, "summary.csv"),
                 row.names = FALSE)

if (opts$scenario_set == "paper") {
  paper <- utils::read.csv(file.path("analysis", "simulations",
                                     "paper-tables.csv"))
  stats_cols <- setdiff(names(paper), c("scenario", "size", "parameter"))
  names(paper)[match(stats_cols, names(paper))] <- paste0(stats_cols,
                                                          "_paper")
  cmp <- merge(summary_table, paper, by = c("scenario", "size", "parameter"))
  # Monte Carlo standard errors of the computed mean and coverage, to judge
  # whether a difference from the article is beyond simulation noise.
  sd_est <- sqrt(pmax(cmp$mse - cmp$bias^2, 0))
  cmp$mc_se_mean <- sd_est / sqrt(cmp$used)
  cmp$mean_diff_in_mc_se <- (cmp$mean - cmp$mean_paper) / cmp$mc_se_mean
  cmp$mc_se_coverage <- sqrt(0.95 * 0.05 / cmp$used)
  utils::write.csv(cmp, file.path(dir, "comparison-with-paper.csv"),
                   row.names = FALSE)
  print(cmp[, c("scenario", "size", "parameter", "used", "mean",
                "mean_paper", "mean_diff_in_mc_se", "mse", "mse_paper",
                "coverage_expected", "coverage_legacy", "coverage_paper")],
        digits = 4, row.names = FALSE)
} else {
  print(summary_table, digits = 4, row.names = FALSE)
}

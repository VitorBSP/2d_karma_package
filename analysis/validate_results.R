# Validate that the reproduction pipelines still give the results obtained
# so far.
#
# Usage (from the repository root, after installing karma2d):
#   Rscript analysis/validate_results.R                 # rerun + check all
#   Rscript analysis/validate_results.R --check=regression
#   Rscript analysis/validate_results.R --check=paper
#   Rscript analysis/validate_results.R --rerun=FALSE --output=DIR
#   Rscript analysis/validate_results.R --rerun=FALSE \
#       --mc_dir=analysis/output/simulations_r1000 --check=paper
#   Rscript analysis/validate_results.R --update_snapshot=TRUE
#
# Two checks:
#
#   regression  Every number produced by the pipelines is compared with the
#               snapshot in analysis/reference/ (results already obtained
#               with this code). Tolerance: relative 1e-5 (estimates may
#               differ in the last digits across BLAS/platforms; anything
#               larger means the results changed). Monte Carlo summaries are
#               compared only with a snapshot of the same configuration
#               (number of replications and sizes), since the streams are
#               deterministic given the seed.
#
#   paper       The results are compared with the values printed in the
#               article: NDVI Tables 4-5 within one unit of the last printed
#               decimal; Monte Carlo means within 3 Monte Carlo standard
#               errors and coverages within 3 * sqrt(p (1 - p) / R) + 0.01,
#               with p the article's coverage.
#               Divergences already documented in docs/technical-audit.md
#               are listed in analysis/reference/known-divergences.csv and
#               reported as "known" instead of failing.
#
# Steps run with --rerun=TRUE (default): the NDVI pipeline and the quick
# Monte Carlo (20 replications, sizes 10 and 25) into --output. Use
# --rerun=FALSE to validate outputs already produced (e.g. a full Monte Carlo
# run written with --output=DIR). The script exits with status 1 when any
# unexpected difference is found and writes validation-report.csv to
# --output.

if (!file.exists("analysis/R/common.R")) {
  stop("Run this script from the repository root.")
}
source(file.path("analysis", "R", "common.R"))
opts <- parse_args(list(check = "all", rerun = "TRUE",
                        output = file.path("analysis", "output",
                                           "validation"),
                        ndvi_dir = "", mc_dir = "",
                        rel_tol = "1e-5", update_snapshot = "FALSE"))
reference_dir <- file.path("analysis", "reference")
out <- opts$output
ndvi_dir <- if (nzchar(opts$ndvi_dir)) opts$ndvi_dir else file.path(out, "ndvi")
mc_dir <- if (nzchar(opts$mc_dir)) opts$mc_dir else file.path(out, "mc")
rel_tol <- as.numeric(opts$rel_tol)
checks <- if (opts$check == "all") c("regression", "paper") else opts$check
stopifnot(all(checks %in% c("regression", "paper")))

run_script <- function(script, args) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c(script, args), stdout = "", stderr = "")
  if (status != 0) stop(script, " failed with status ", status)
}

if (as.logical(opts$rerun)) {
  message("Running the NDVI pipeline ...")
  run_script(file.path("analysis", "ndvi", "run_ndvi.R"),
             paste0("--output=", ndvi_dir))
  message("Running the quick Monte Carlo ...")
  run_script(file.path("analysis", "simulations", "run_monte_carlo.R"),
             c("--mode=quick", paste0("--output=", mc_dir)))
  run_script(file.path("analysis", "simulations", "summarize_monte_carlo.R"),
             paste0("--input=", mc_dir))
}

# Result tables produced by the pipelines and the columns identifying a row.
result_tables <- list(
  ndvi_table4 = list(file = file.path(ndvi_dir, "table4-estimates.csv"),
                     keys = c("model", "parameter")),
  ndvi_table5 = list(file = file.path(ndvi_dir, "table5-metrics.csv"),
                     keys = c("model", "region")),
  ndvi_detection = list(file = file.path(ndvi_dir,
                                         "anomaly-detection.csv"),
                        keys = c("element", "model")),
  ndvi_model_selection = list(file = file.path(ndvi_dir,
                                               "model-selection.csv"),
                              keys = c("p", "q")),
  monte_carlo = list(file = file.path(mc_dir, "paper", "summary.csv"),
                     keys = c("scenario", "size", "parameter"))
)

# Snapshot file for a table; Monte Carlo snapshots are per configuration.
snapshot_name <- function(name, table) {
  if (name != "monte_carlo") return(paste0(name, ".csv"))
  sprintf("monte_carlo_R%d_sizes%s.csv", max(table$replications),
          paste(sort(unique(table$size)), collapse = "-"))
}

read_result <- function(name) {
  path <- result_tables[[name]]$file
  if (file.exists(path)) utils::read.csv(path, check.names = FALSE) else NULL
}

if (as.logical(opts$update_snapshot)) {
  dir.create(reference_dir, showWarnings = FALSE)
  for (name in names(result_tables)) {
    table <- read_result(name)
    if (is.null(table)) next
    dest <- file.path(reference_dir, snapshot_name(name, table))
    utils::write.csv(table, dest, row.names = FALSE)
    message("Snapshot written: ", dest)
  }
  quit(save = "no", status = 0)
}

# Compare two tables cell by cell on their numeric columns.
compare_tables <- function(new, ref, keys, check, label) {
  key_new <- do.call(paste, c(new[keys], sep = " | "))
  key_ref <- do.call(paste, c(ref[keys], sep = " | "))
  rows <- list()
  missing <- setdiff(key_ref, key_new)
  if (length(missing)) {
    rows[[length(rows) + 1]] <- data.frame(
      check = check, table = label, row = missing, column = "(row)",
      reference = NA, computed = NA, difference = NA, tolerance = NA,
      status = "FAIL: row missing")
  }
  numeric_cols <- setdiff(names(ref)[vapply(ref, is.numeric, logical(1))],
                          keys)
  for (col in intersect(numeric_cols, names(new))) {
    a <- ref[[col]][match(intersect(key_ref, key_new), key_ref)]
    b <- new[[col]][match(intersect(key_ref, key_new), key_new)]
    tol <- rel_tol * pmax(abs(a), 1e-3)
    both_na <- is.na(a) & is.na(b)
    ok <- both_na | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
    rows[[length(rows) + 1]] <- data.frame(
      check = check, table = label, row = intersect(key_ref, key_new),
      column = col, reference = a, computed = b, difference = b - a,
      tolerance = tol, status = ifelse(ok, "ok", "FAIL"))
  }
  # Non-numeric columns (e.g. inference status) must match exactly.
  for (col in setdiff(intersect(names(ref), names(new)),
                      c(keys, numeric_cols))) {
    a <- as.character(ref[[col]][match(intersect(key_ref, key_new),
                                       key_ref)])
    b <- as.character(new[[col]][match(intersect(key_ref, key_new),
                                       key_new)])
    ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    rows[[length(rows) + 1]] <- data.frame(
      check = check, table = label, row = intersect(key_ref, key_new),
      column = col, reference = NA, computed = NA, difference = NA,
      tolerance = NA, status = ifelse(ok, "ok", "FAIL"))
  }
  do.call(rbind, rows)
}

report <- list()

if ("regression" %in% checks) {
  for (name in names(result_tables)) {
    new <- read_result(name)
    if (is.null(new)) {
      message("regression: no output for ", name, " (skipped)")
      next
    }
    ref_path <- file.path(reference_dir, snapshot_name(name, new))
    if (!file.exists(ref_path)) {
      message("regression: no snapshot ", basename(ref_path),
              " for this configuration (skipped)")
      next
    }
    ref <- utils::read.csv(ref_path, check.names = FALSE)
    report[[length(report) + 1]] <- compare_tables(
      new, ref, result_tables[[name]]$keys, "regression", name)
  }
}

if ("paper" %in% checks) {
  known <- utils::read.csv(file.path(reference_dir, "known-divergences.csv"))
  is_known <- function(table, row) {
    paste(table, row) %in% paste(known$table, known$row)
  }
  ndvi_cmp_path <- file.path(ndvi_dir, "comparison-with-paper.csv")
  if (file.exists(ndvi_cmp_path)) {
    cmp <- utils::read.csv(ndvi_cmp_path)
    tol <- 10^(-cmp$decimals)
    status <- ifelse(cmp$agrees, "ok",
                     ifelse(is_known("ndvi_paper", cmp$quantity), "known",
                            "FAIL"))
    report[[length(report) + 1]] <- data.frame(
      check = "paper", table = "ndvi_paper", row = cmp$quantity,
      column = "value", reference = cmp$published, computed = cmp$computed,
      difference = cmp$computed - cmp$published, tolerance = tol,
      status = status)
  }
  mc_cmp_path <- file.path(mc_dir, "paper", "comparison-with-paper.csv")
  if (file.exists(mc_cmp_path)) {
    cmp <- utils::read.csv(mc_cmp_path)
    row <- paste(cmp$scenario, cmp$size, cmp$parameter)
    tol_mean <- 3 * cmp$mc_se_mean
    p_cov <- pmin(pmax(cmp$coverage_paper, 0.01), 0.99)
    tol_cov <- 3 * sqrt(p_cov * (1 - p_cov) / cmp$used) + 0.01
    mk <- function(column, ref, val, tol) {
      ok <- abs(val - ref) <= tol
      key <- paste(row, column)
      data.frame(check = "paper", table = "monte_carlo_paper", row = key,
                 column = column, reference = ref, computed = val,
                 difference = val - ref, tolerance = tol,
                 status = ifelse(ok, "ok",
                                 ifelse(is_known("monte_carlo_paper", key),
                                        "known", "FAIL")))
    }
    report[[length(report) + 1]] <- rbind(
      mk("mean", cmp$mean_paper, cmp$mean, tol_mean),
      mk("coverage", cmp$coverage_paper, cmp$coverage_expected, tol_cov)
    )
    if (max(cmp$used) < 1000) {
      message("paper: the Monte Carlo has only ", max(cmp$used),
              " replications; its comparison with the article is",
              " indicative only (wide tolerances).")
    }
  }
}

if (length(report) == 0L) stop("Nothing was validated; check --output.")
report <- do.call(rbind, report)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(report, file.path(out, "validation-report.csv"),
                 row.names = FALSE)

summary_table <- as.data.frame.matrix(table(
  paste(report$check, report$table), report$status))
print(summary_table)
failures <- report[startsWith(report$status, "FAIL"), ]
if (nrow(failures) > 0) {
  cat("\nUnexpected differences:\n")
  print(utils::head(failures[, c("check", "table", "row", "column",
                                 "reference", "computed", "tolerance")],
                    50), row.names = FALSE)
  cat(sprintf("\n%d unexpected difference(s); full report in %s\n",
              nrow(failures), file.path(out, "validation-report.csv")))
  quit(save = "no", status = 1)
}
cat(sprintf("\nAll %d comparisons passed (%d known divergences).\n",
            nrow(report), sum(report$status == "known")))

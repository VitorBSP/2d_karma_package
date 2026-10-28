# Monte Carlo study of the CMLE of the 2-D KARMA model (Section 4.1).
#
# Usage (from the repository root, after installing karma2d):
#   Rscript analysis/simulations/run_monte_carlo.R --mode=quick
#   Rscript analysis/simulations/run_monte_carlo.R --mode=full --cores=8
# Options:
#   --mode=quick|full          replications and image sizes (config.R)
#   --replications=N           override the number of replications
#   --scenario_set=paper|legacy_scripts
#   --scenarios=karma10,karma01,karma11
#   --cores=N                  parallel workers (forked; use 1 on Windows)
#   --seed=N                   master seed
#   --output=DIR               results directory (outside the package)
#
# Reproducibility: every replication r of every (scenario, size) cell gets
# its own L'Ecuyer-CMRG stream, derived deterministically from the master
# seed, so results do not depend on the number of cores or on the order in
# which cells are run. Each cell is saved as soon as it finishes
# (checkpoint); an existing cell file with the same configuration is
# reused, so an interrupted run can be resumed with the same command.
#
# Each replication records the estimates, the convergence code, and
# standard errors from (i) the expected information (package default),
# (ii) the observed information, and (iii) the information used by the
# original scripts (expected blocks with the observed shape-shape entry).
# The original scripts did not compute coverage rates, so it is not known
# which matrix produced the coverage rates printed in the article.

if (!file.exists("analysis/R/common.R")) {
  stop("Run this script from the repository root.")
}
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "simulations", "config.R"))

opts <- parse_args(list(mode = "quick", replications = "",
                        scenario_set = "paper",
                        scenarios = "karma10,karma01,karma11", cores = "1",
                        seed = "20260101",
                        output = file.path("analysis", "output",
                                           "simulations")))
mode <- modes[[opts$mode]]
if (is.null(mode)) stop("--mode must be 'quick' or 'full'.")
replications <- if (nzchar(opts$replications))
  as.integer(opts$replications) else mode$replications
scenarios <- scenario_sets[[opts$scenario_set]]
if (is.null(scenarios)) stop("Unknown --scenario_set.")
scenario_names <- strsplit(opts$scenarios, ",")[[1]]
stopifnot(all(scenario_names %in% names(scenarios)))
cores <- as.integer(opts$cores)
master_seed <- as.integer(opts$seed)
out_dir <- file.path(opts$output, opts$scenario_set)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

one_replication <- function(stream, s, size) {
  assign(".Random.seed", stream, envir = globalenv())
  order <- c((sqrt(length(s$ar) + 1) - 1), (sqrt(length(s$ma) + 1) - 1))
  truth <- c(s$alpha, s$ar, s$ma, s$shape)
  n_par <- length(truth)
  empty <- rep(NA_real_, n_par)
  result <- list(estimate = empty, se_expected = empty, se_observed = empty,
                 se_legacy = empty, convergence = NA_integer_,
                 error = NA_character_)
  y <- tryCatch(karma2d_simulate(size, size, alpha = s$alpha, ar = s$ar,
                                 ma = s$ma, shape = s$shape),
                error = function(e) e)
  if (inherits(y, "error")) {
    result$error <- paste("simulation:", conditionMessage(y))
    return(result)
  }
  fit <- tryCatch(suppressWarnings(karma2d_fit(y, order)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    result$error <- paste("fit:", conditionMessage(fit))
    return(result)
  }
  est <- coef(fit)
  se_from <- function(info) {
    v <- tryCatch(diag(solve(info)), error = function(e) empty)
    ifelse(is.finite(v) & v > 0, sqrt(pmax(v, 0)), NA_real_)
  }
  obs <- tryCatch(karma2d_information(y, est, order, type = "observed"),
                  error = function(e) NULL)
  legacy <- fit$information
  if (!is.null(obs)) legacy[n_par, n_par] <- obs[n_par, n_par]
  result$estimate <- unname(est)
  result$se_expected <- unname(sqrt(diag(fit$vcov)))
  result$se_observed <- if (is.null(obs)) empty else unname(se_from(obs))
  result$se_legacy <- unname(se_from(legacy))
  result$convergence <- fit$convergence
  result
}

# Independent stream for replication r of cell (i, j): advance the master
# L'Ecuyer-CMRG stream once per (cell, replication), in a fixed order.
cell_streams <- function(cell_index, n_rep) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(master_seed + 1000L * cell_index)
  s <- .Random.seed
  streams <- vector("list", n_rep)
  for (r in seq_len(n_rep)) {
    streams[[r]] <- s
    s <- parallel::nextRNGStream(s)
  }
  streams
}

cells <- expand.grid(size = mode$sizes, scenario = scenario_names,
                     stringsAsFactors = FALSE)
for (k in seq_len(nrow(cells))) {
  name <- cells$scenario[k]
  size <- cells$size[k]
  s <- scenarios[[name]]
  cell_index <- match(name, names(scenario_sets[[1]])) * 100L + size
  file <- file.path(out_dir, sprintf("%s_%dx%d.rds", name, size, size))
  config <- list(scenario = s, size = size, replications = replications,
                 seed = master_seed, cell_index = cell_index,
                 scenario_set = opts$scenario_set,
                 karma2d_version = as.character(
                   utils::packageVersion("karma2d")))
  if (file.exists(file)) {
    old <- readRDS(file)
    if (identical(old$config[names(config) != "karma2d_version"],
                  config[names(config) != "karma2d_version"])) {
      message("Reusing checkpoint ", file)
      next
    }
  }
  message(sprintf("%s %dx%d: %d replications on %d core(s)", name, size,
                  size, replications, cores))
  streams <- cell_streams(cell_index, replications)
  started <- Sys.time()
  runs <- if (cores > 1L) {
    parallel::mclapply(streams, one_replication, s = s, size = size,
                       mc.cores = cores, mc.preschedule = TRUE)
  } else {
    lapply(streams, one_replication, s = s, size = size)
  }
  pick <- function(field) do.call(rbind, lapply(runs, `[[`, field))
  saveRDS(list(
    config = config,
    parameter_names = karma2d:::karma2d_coef_names(
      c(sqrt(length(s$ar) + 1) - 1, sqrt(length(s$ma) + 1) - 1)),
    truth = c(s$alpha, s$ar, s$ma, s$shape),
    estimate = pick("estimate"), se_expected = pick("se_expected"),
    se_observed = pick("se_observed"), se_legacy = pick("se_legacy"),
    convergence = vapply(runs, `[[`, integer(1), "convergence"),
    error = vapply(runs, `[[`, character(1), "error"),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started,
                                          units = "secs"))
  ), file)
}
write_session_info(out_dir)
message("Done. Summarize with: Rscript analysis/simulations/",
        "summarize_monte_carlo.R --input=", opts$output)

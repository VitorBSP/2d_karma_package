# Shared helpers for the reproduction pipelines in analysis/.
# These scripts are NOT part of the installed package; they depend on the
# installed karma2d package and use a few of its internal functions through
# `karma2d:::` (documented where used).

suppressPackageStartupMessages(library(karma2d))

# Parse --key=value command line arguments into a named list.
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  out <- defaults
  for (a in args) {
    kv <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("Arguments must look like --key=value: ", a)
    key <- gsub("-", "_", kv[2])
    if (!key %in% names(defaults)) stop("Unknown argument --", kv[2])
    out[[key]] <- methods::as(kv[3], class(defaults[[key]]))
  }
  out
}

# Repository root, found from the location of this file or the working dir.
repo_root <- function() {
  candidates <- c(getwd(), dirname(getwd()), dirname(dirname(getwd())))
  hit <- candidates[file.exists(file.path(candidates, "DESCRIPTION")) &
                      dir.exists(file.path(candidates, "analysis"))]
  if (length(hit) == 0L) {
    stop("Run the pipelines from the repository root (or a subfolder).")
  }
  normalizePath(hit[1])
}

check_sha256 <- function(path, expected) {
  if (nzchar(Sys.which("sha256sum"))) {
    got <- strsplit(system2("sha256sum", shQuote(path), stdout = TRUE),
                    " ")[[1]][1]
  } else if (requireNamespace("digest", quietly = TRUE)) {
    got <- digest::digest(path, algo = "sha256", file = TRUE)
  } else {
    warning("Cannot verify checksum of ", path, " (no sha256sum/digest).")
    return(invisible(NA))
  }
  if (!identical(got, expected)) {
    stop("Checksum mismatch for ", path, ":\n  expected ", expected,
         "\n  got      ", got)
  }
  invisible(TRUE)
}

# Write session information next to the results.
write_session_info <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    paste("Date (UTC):", format(Sys.time(), tz = "UTC")),
    paste("karma2d:", as.character(utils::packageVersion("karma2d"))),
    utils::capture.output(utils::sessionInfo()),
    if (requireNamespace("terra", quietly = TRUE))
      paste("GDAL:", paste(terra::gdal(lib = "all"), collapse = " / "))
  ), file.path(dir, "session-info.txt"))
}

# Sample skewness and kurtosis (moment estimators, kurtosis = 3 under
# normality), as computed by moments::skewness/kurtosis in the original
# scripts.
sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- x - mean(x)
  mean(m^3) / mean(m^2)^1.5
}
sample_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  m <- x - mean(x)
  mean(m^4) / mean(m^2)^2
}

# Information matrix used by the original scripts (and therefore by the
# published standard errors): expected gamma-gamma and gamma-shape blocks
# combined with the OBSERVED shape-shape entry. Kept only to reproduce the
# published numbers; see docs/technical-audit.md.
legacy_information <- function(y, parameters, order) {
  k <- karma2d_information(y, parameters, order, type = "expected")
  o <- karma2d_information(y, parameters, order, type = "observed")
  n <- nrow(k)
  k[n, n] <- o[n, n]
  k
}

image_png <- function(file, m, main, zlim = range(m, na.rm = TRUE),
                      col = grDevices::hcl.colors(64, "YlGn", rev = TRUE)) {
  grDevices::png(file, width = 900, height = 900, res = 150)
  on.exit(grDevices::dev.off())
  graphics::image(t(m[nrow(m):1, , drop = FALSE]), zlim = zlim, col = col,
                  main = main, axes = FALSE, asp = 1)
}

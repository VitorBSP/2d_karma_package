# GRDI application of Pereira, Palm and Bayer (2026), Section 4.2.2.
#
# Usage (from the repository root, after installing karma2d):
#   Rscript analysis/grdi/run_grdi.R --input=PATH/nasa_corte_india_expandido.csv
#       [--output=analysis/output/grdi] [--karma_order=0,1]
#       [--gaussian_order=0,1] [--select=TRUE]
#
# DATA AVAILABILITY. The original scripts read
# `grdi/data/nasa_corte_india_expandido.csv`, a CSV export of a crop of the
# Global Gridded Relative Deprivation Index (GRDI) v1 raster
# (povmap-grdi-v1.tif, NASA SEDAC/CIESIN). Neither that CSV nor the raster
# was committed to the repository (both were listed in .gitignore), and the
# crop extent used to create the CSV was not recorded. The CSV that IS in
# the repository (analysis/data/grdi/nasa_corte_india.csv, 500 x 1001) is
# not read by any original script and is too small for the window used
# below. Therefore this pipeline cannot reproduce the article until the
# expanded CSV is supplied; see docs/data-provenance.md.
#
# Steps (as in the original scripts, analysis/legacy/grdi/*_v2.R):
#   1. Read the CSV (header row = column names), take rows 1200:1300 and
#      columns 2200:2300 of the data matrix, keep the first 100 x 100 block,
#      divide by 100, and replace missing values by the image mean.
#   2. Phase I window: rows/columns 51:100. Synthetic anomaly: pixels
#      [19:23, 41:50] set to 0.1.
#   3.-7. As for NDVI, with a 3 x 5 rectangular structuring element.

local({
  if (!file.exists("analysis/R/common.R")) {
    stop("Run this script from the repository root.")
  }
  source(file.path("analysis", "R", "common.R"))
  source(file.path("analysis", "R", "gaussian_arma2d.R"))
  source(file.path("analysis", "R", "application.R"))
})

opts <- parse_args(list(input = "", output = file.path("analysis", "output",
                                                         "grdi"),
                        karma_order = "0,1", gaussian_order = "0,1",
                        select = "TRUE"))
if (!nzchar(opts$input) || !file.exists(opts$input)) {
  stop(paste(
    "The GRDI input CSV is required: --input=PATH/nasa_corte_india_expandido.csv.",
    "It is not distributed with the repository; see the header of this",
    "script and docs/data-provenance.md."
  ), call. = FALSE)
}

raw <- as.matrix(utils::read.csv(opts$input, check.names = FALSE))
if (nrow(raw) < 1300 || ncol(raw) < 2300) {
  stop(sprintf(paste(
    "The input has %d x %d cells; the original window (rows 1200:1300,",
    "columns 2200:2300) needs at least 1300 x 2300."
  ), nrow(raw), ncol(raw)), call. = FALSE)
}
grdi <- raw[1200:1300, 2200:2300][1:100, 1:100] / 100
n_missing <- sum(is.na(grdi))
grdi[is.na(grdi)] <- mean(grdi, na.rm = TRUE) # as in the original scripts
message(sprintf("GRDI window: %d missing values replaced by the mean.",
                n_missing))
if (any(grdi <= 0 | grdi >= 1)) {
  stop("GRDI values on the unit scale must lie strictly inside (0, 1).")
}

run_application(list(
  name = "GRDI",
  y = grdi,
  phase1 = list(rows = 51:100, cols = 51:100),
  anomaly = list(rows = 19:23, cols = 41:50, value = 0.1),
  karma_order = as.integer(strsplit(opts$karma_order, ",")[[1]]),
  gaussian_order = as.integer(strsplit(opts$gaussian_order, ",")[[1]]),
  structuring_element = karma2d_structuring_element(3, 5),
  select = as.logical(opts$select),
  article_choice = "2-D KARMA(0, 1)",
  published = utils::read.csv(file.path("analysis", "grdi",
                                        "published-values.csv")),
  to_display = function(v) 100 * v,
  display_range = c(0, 100),
  out_dir = opts$output
))

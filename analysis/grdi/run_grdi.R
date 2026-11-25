# GRDI application of Pereira, Palm and Bayer (2026), Section 4.2.2.
#
# Usage (from the repository root, after installing karma2d):
#   Rscript analysis/grdi/run_grdi.R --input=PATH/nasa_corte_india_expandido.csv
#       [--output=analysis/output/grdi] [--karma_order=0,1]
#       [--gaussian_order=0,1] [--select=TRUE] [--procedure=article]
#       [--window=1200,2200] [--missing_value=0.01]
#
# Functional test on the committed extraction (NOT the article's window):
#   Rscript analysis/grdi/run_grdi.R \
#     --input=analysis/data/grdi/nasa_corte_india.csv \
#     --window=246,674 --missing_value=0.01
#
# DATA AVAILABILITY. The original scripts read
# `grdi/data/nasa_corte_india_expandido.csv`, a CSV export of a crop of the
# Global Gridded Relative Deprivation Index (GRDI) v1 raster
# (povmap-grdi-v1.tif, NASA SEDAC/CIESIN). Neither that CSV nor the raster
# was committed to the repository (both were listed in .gitignore), and the
# crop extent used to create the CSV was not recorded. The CSV that IS in
# the repository (analysis/data/grdi/nasa_corte_india.csv, 501 x 1001) is
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
#
# PROCEDURE. The authors' archive that reproduces the NDVI results (Tables
# 4-5) did not include the GRDI scripts. The anomaly value (0.1) and the
# window above are therefore those of the scripts first committed to the
# repository (analysis/legacy/grdi/), not verified against the article's
# own computation; the NDVI procedure "article" (effective-region
# calibration metrics, unpadded recursion cropped by 2w, opening then
# closing, Gaussian residuals divided by sigma_hat) is applied for
# consistency. --procedure=padded gives the earlier pipeline.

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
                        select = "TRUE", procedure = "article",
                        window = "1200,2200", missing_value = ""))
if (!nzchar(opts$input) || !file.exists(opts$input)) {
  stop(paste(
    "The GRDI input CSV is required: --input=PATH/nasa_corte_india_expandido.csv.",
    "It is not distributed with the repository; see the header of this",
    "script and docs/data-provenance.md."
  ), call. = FALSE)
}

# --window=ROW,COL is the top-left cell of the 100 x 100 window in the data
# matrix. The default is the window of the original scripts
# (raw[1200:1300, 2200:2300][1:100, 1:100]). --missing_value=X treats cells
# equal to X (before dividing by 100) as missing; the committed extraction
# analysis/data/grdi/nasa_corte_india.csv codes them as 0.01.
corner <- as.integer(strsplit(opts$window, ",")[[1]])
article_window <- identical(corner, c(1200L, 2200L))
raw <- as.matrix(utils::read.csv(opts$input, check.names = FALSE))
rows <- corner[1] + 0:99
cols <- corner[2] + 0:99
if (max(rows) > nrow(raw) || max(cols) > ncol(raw)) {
  stop(sprintf(paste(
    "The input has %d x %d cells; the window starting at (%d, %d) needs at",
    "least %d x %d. The committed extraction does not contain the article's",
    "window; use --window to choose another one (see docs/data-provenance.md)."
  ), nrow(raw), ncol(raw), corner[1], corner[2], max(rows), max(cols)),
  call. = FALSE)
}
grdi <- raw[rows, cols]
if (nzchar(opts$missing_value)) {
  grdi[grdi == as.numeric(opts$missing_value)] <- NA
}
grdi <- grdi / 100
n_missing <- sum(is.na(grdi))
grdi[is.na(grdi)] <- mean(grdi, na.rm = TRUE) # as in the original scripts
message(sprintf("GRDI window: %d missing values replaced by the mean.",
                n_missing))
if (any(grdi <= 0 | grdi >= 1)) {
  stop("GRDI values on the unit scale must lie strictly inside (0, 1).")
}
if (!article_window) {
  message("This is NOT the article's window: no comparison with Tables 6-7.")
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
  procedure = opts$procedure,
  article_choice = "2-D KARMA(0, 1)",
  published = if (article_window) {
    utils::read.csv(file.path("analysis", "grdi", "published-values.csv"))
  },
  to_display = function(v) 100 * v,
  display_range = c(0, 100),
  out_dir = opts$output
))

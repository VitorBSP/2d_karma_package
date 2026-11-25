# NDVI application of Pereira, Palm and Bayer (2026), Section 4.2.1.
#
# Usage (from the repository root, after installing karma2d):
#   Rscript analysis/ndvi/run_ndvi.R [--output=analysis/output/ndvi]
#       [--karma_order=1,0] [--gaussian_order=0,1] [--select=TRUE]
#       [--procedure=article]
#
# Steps (--procedure=article, the procedure of the authors' scripts that
# produced Tables 4-5 and Figures 4-7; see analysis/R/application.R)
#   1. Read the MODIS NDVI scene (checksum-verified) and rescale it to (0, 1)
#      by (x + 1) / 2.
#   2. Phase I (calibration) window: rows/columns 51:100. Synthetic anomaly
#      (as in the authors' scripts): pixels [5:20, 15:30] set to 0.01.
#   3. Model selection over p, q in {0, ..., 3} on the Phase I window
#      (optional, --select): table of significance and accuracy metrics.
#   4. Fit 2-D KARMA(--karma_order) and the Gaussian 2-D ARMA
#      (--gaussian_order) benchmark on the Phase I window (Table 4).
#   5. Calibration metrics on the effective region of the Phase I fit; full
#      scene: recursion (Eq. 8) on the unpadded image, first 2w rows and
#      columns dropped; MAE, MSE, SMAPE, RMSLE (Table 5), plus the Gaussian
#      rows "as published" (with the defects of the printed Table 5).
#   6. Quantile residuals (KARMA) and standardized residuals (Gaussian,
#      divided by sigma_hat), control limits +/-3, opening then closing with
#      a 3 x 3 square (Figures 4-7).
#   7. Write tables, figures, the comparison with the published values and
#      session information to --output.
# --procedure=padded runs the earlier pipeline of this repository (duplicate
# padding, opening only, Gaussian residuals divided by their sd).

local({
  root <- if (file.exists("analysis/R/common.R")) "." else
    stop("Run this script from the repository root.")
  source(file.path(root, "analysis", "R", "common.R"))
  source(file.path(root, "analysis", "R", "gaussian_arma2d.R"))
  source(file.path(root, "analysis", "R", "application.R"))
})

opts <- parse_args(list(output = file.path("analysis", "output", "ndvi"),
                        karma_order = "1,0", gaussian_order = "0,1",
                        select = "TRUE", procedure = "article"))
karma_order <- as.integer(strsplit(opts$karma_order, ",")[[1]])
gaussian_order <- as.integer(strsplit(opts$gaussian_order, ",")[[1]])
out_dir <- opts$output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Data -------------------------------------------------------------------
tif <- file.path("analysis", "data", "ndvi", "NDVI_2019_10.tif")
check_sha256(tif, paste0("aac7580d93b6a36c43e9aeaa394ea194",
                         "c7600bc94aea8eb72a5be3c7c8862173"))
ndvi <- terra::as.matrix(terra::rast(tif), wide = TRUE)
dimnames(ndvi) <- NULL
y <- (ndvi + 1) / 2

# 2.-7. Phase I window (rows/columns 51:100) and synthetic anomaly
# (pixels [5:20, 15:30] set to 0.01, as in the authors' scripts) -----------
run_application(list(
  name = "NDVI",
  y = y,
  phase1 = list(rows = 51:100, cols = 51:100),
  anomaly = list(rows = 5:20, cols = 15:30, value = 0.01),
  karma_order = karma_order,
  gaussian_order = gaussian_order,
  structuring_element = karma2d_structuring_element(3, 3),
  select = as.logical(opts$select),
  procedure = opts$procedure,
  article_choice = "2-D KARMA(1, 0)",
  published = utils::read.csv(file.path("analysis", "ndvi",
                                        "published-values.csv")),
  to_display = function(v) 2 * v - 1,
  display_range = c(-1, 1),
  out_dir = out_dir
))

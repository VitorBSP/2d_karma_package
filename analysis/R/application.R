# Shared pipeline of the NDVI and GRDI applications (Section 4.2).
#
# config: list with
#   name, y (image on the unit scale, no NA), phase1 (list rows, cols),
#   anomaly (list rows, cols, value), karma_order, gaussian_order,
#   structuring_element, select (logical), published (data frame or NULL),
#   out_dir, and optionally
#   procedure        "article" (default) or "padded", see below;
#   post_processing  "opening_closing" (default for "article") or "opening"
#                    (default for "padded"), applied to the |r| > 3 map with
#                    structuring_element.
#
# procedure = "article" follows the authors' scripts that produced Tables
# 4-5 and Figures 4-7 of the article:
#   * calibration metrics compare the Phase I window with the fitted values
#     of the same fit on its effective region (rows/columns > w, w = max(p,
#     q)), without padding;
#   * full scene: the recursion (Eq. 8) runs on the unpadded image with the
#     anomaly and the first 2w rows and columns are dropped from the image,
#     the prediction and the residuals (98 x 98 for w = 1);
#   * |r| > 3, then opening followed by closing with structuring_element;
#   * Gaussian standardized residuals (y - prediction) / sigma_hat;
#   * the Gaussian rows printed in the article are also reproduced, labelled
#     "as published", with their two defects (see gaussian_arma2d.R).
# procedure = "padded" is the earlier pipeline of this repository: duplicate
# padding, metrics on the full N x M image and on the Phase I window of the
# padded prediction, opening only, Gaussian residuals standardized by their
# standard deviation over the scene.
run_application <- function(config) {
  y <- config$y
  phase1 <- config$phase1
  anomaly <- config$anomaly
  karma_order <- config$karma_order
  gaussian_order <- config$gaussian_order
  out_dir <- config$out_dir
  procedure <- if (is.null(config$procedure)) "article" else config$procedure
  procedure <- match.arg(procedure, c("article", "padded"))
  post_processing <- config$post_processing
  if (is.null(post_processing)) {
    post_processing <- if (procedure == "article") "opening_closing" else
      "opening"
  }
  post_processing <- match.arg(post_processing,
                               c("opening_closing", "opening"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  y_train <- y[phase1$rows, phase1$cols]
  y_anomaly <- y
  y_anomaly[anomaly$rows, anomaly$cols] <- anomaly$value
  to_display <- config$to_display
  display_range <- config$display_range

  # Effective region of the Phase I window and crop of the full scene.
  effective <- function(m, w) {
    m[seq(w + 1, nrow(m)), seq(w + 1, ncol(m)), drop = FALSE]
  }
  crop <- function(m, w) {
    if (w == 0) return(m)
    m[-seq_len(2 * w), -seq_len(2 * w), drop = FALSE]
  }
  karma_calibration_metrics <- function(fit) {
    w <- max(fit$order)
    karma2d_metrics(effective(y_train, w), effective(fitted(fit), w))
  }

  # 3. Model selection --------------------------------------------------------
  select_models <- function(y_train, y_full) {
    grid <- expand.grid(p = 0:3, q = 0:3)
    grid <- grid[!(grid$p == 0 & grid$q == 0), ]
    rows <- lapply(seq_len(nrow(grid)), function(i) {
      ord <- c(grid$p[i], grid$q[i])
      fit <- tryCatch(suppressWarnings(karma2d_fit(y_train, ord)),
                      error = function(e) NULL)
      if (is.null(fit)) {
        return(data.frame(p = ord[1], q = ord[2], n_par = NA,
                          converged = FALSE, all_significant = NA,
                          inference = "fit failed", MAE = NA, MSE = NA,
                          SMAPE = NA, RMSLE = NA))
      }
      tab <- summary(fit)$coefficients
      m <- if (procedure == "article") {
        karma_calibration_metrics(fit)
      } else {
        pred <- predict(fit, newdata = y_full, padding = "duplicate")
        karma2d_metrics(y_full[phase1$rows, phase1$cols],
                        pred[phase1$rows, phase1$cols])
      }
      data.frame(p = ord[1], q = ord[2], n_par = nrow(tab),
                 converged = fit$convergence == 0L,
                 all_significant = all(tab[, "Pr(>|z|)"] < 0.05),
                 inference = fit$inference$status,
                 MAE = m[["MAE"]], MSE = m[["MSE"]], SMAPE = m[["SMAPE"]],
                 RMSLE = m[["RMSLE"]])
    })
    do.call(rbind, rows)
  }

  if (isTRUE(config$select)) {
    selection <- select_models(y_train, y_anomaly)
    utils::write.csv(selection, file.path(out_dir, "model-selection.csv"),
                     row.names = FALSE)
    # Operationalization of the article's criteria (identification by
    # significance, then parsimony, ties broken by MSE). The article also
    # used graphical residual checks, which are not automated here.
    cand <- selection[selection$converged & selection$all_significant %in%
                        TRUE & selection$inference == "ok", ]
    cand <- cand[order(cand$n_par, cand$MSE), ]
    message("Automated selection: 2-D KARMA(", cand$p[1], ", ", cand$q[1],
            "); article: ", config$article_choice)
  }

  # 4. Fits on the Phase I window ---------------------------------------------
  fit <- karma2d_fit(y_train, karma_order)
  legacy_se <- sqrt(diag(solve(legacy_information(y_train, coef(fit),
                                                   karma_order))))
  karma_table <- data.frame(
    model = sprintf("2-D KARMA(%d,%d)", karma_order[1], karma_order[2]),
    parameter = names(coef(fit)),
    summary(fit)$coefficients,
    se_legacy_information = legacy_se, check.names = FALSE, row.names = NULL
  )
  gfit <- gaussian_arma2d_fit(y_train, gaussian_order)
  gauss_table <- data.frame(
    model = sprintf("2-D ARMA(%d,%d)", gaussian_order[1], gaussian_order[2]),
    parameter = names(gfit$coefficients),
    Estimate = gfit$coefficients, `Std. Error` = gfit$se,
    `z value` = gfit$coefficients / gfit$se,
    `Pr(>|z|)` = 2 * stats::pnorm(-abs(gfit$coefficients / gfit$se)),
    se_legacy_information = NA, check.names = FALSE, row.names = NULL
  )
  table4 <- rbind(karma_table, gauss_table)
  utils::write.csv(table4, file.path(out_dir, "table4-estimates.csv"),
                   row.names = FALSE)
  karma_model <- karma_table$model[1]
  gauss_model <- gauss_table$model[1]
  gauss_published_model <- paste(gauss_model, "as published")

  # 5. Prediction and accuracy ------------------------------------------------
  # Each model is described by the observed image, its prediction and its
  # residuals on the region used for the full scene, and by the observed and
  # predicted values of the calibration area.
  truth <- matrix(FALSE, nrow(y), ncol(y))
  truth[anomaly$rows, anomaly$cols] <- TRUE
  wk <- max(karma_order)
  wg <- max(gaussian_order)
  sigma_hat <- gfit$coefficients[["sigma"]]
  if (procedure == "article") {
    gauss_calib <- effective(gaussian_arma2d_predict(gfit, y_train,
                                                     padding = "none"), wg)
    obs_gauss_calib <- effective(y_train, wg)
    y_gauss <- crop(y_anomaly, wg)
    pred_gauss <- crop(gaussian_arma2d_predict(gfit, y_anomaly,
                                               padding = "none"), wg)
    pred_gauss_pub <- crop(gaussian_arma2d_published(gfit, y_anomaly), wg)
    stopifnot(nrow(gauss_calib) == ncol(gauss_calib))
    models <- list(
      karma = list(
        model = karma_model, variant = "KARMA",
        calib_obs = effective(y_train, wk),
        calib_pred = effective(fitted(fit), wk),
        y = crop(y_anomaly, wk),
        pred = crop(predict(fit, newdata = y_anomaly, padding = "none"), wk),
        resid = crop(karma2d_residuals(fit, newdata = y_anomaly,
                                       padding = "none"), wk),
        truth = crop(truth, wk), offset = 2 * wk
      ),
      gauss = list(
        model = gauss_model, variant = "Gaussian corrected",
        calib_obs = obs_gauss_calib, calib_pred = gauss_calib,
        y = y_gauss, pred = pred_gauss,
        resid = (y_gauss - pred_gauss) / sigma_hat,
        truth = crop(truth, wg), offset = 2 * wg
      ),
      # Defects of the published Gaussian rows, reproduced for traceability:
      # (a) calibration prediction compared with its transpose (a vector in
      # row-major order against a column-major matrix); (b) full scene with
      # the identity-scale coefficients pushed through the logit recursion.
      gauss_published = list(
        model = gauss_published_model, variant = "Gaussian as published",
        calib_obs = obs_gauss_calib, calib_pred = t(gauss_calib),
        y = y_gauss, pred = pred_gauss_pub,
        resid = (y_gauss - pred_gauss_pub) / sigma_hat,
        truth = crop(truth, wg), offset = 2 * wg
      )
    )
  } else {
    pred_gauss <- gaussian_arma2d_predict(gfit, y_anomaly,
                                          padding = "duplicate")
    pred_karma <- predict(fit, newdata = y_anomaly, padding = "duplicate")
    r_gauss_raw <- y_anomaly - pred_gauss
    models <- list(
      karma = list(
        model = karma_model, variant = "KARMA",
        calib_obs = y_anomaly[phase1$rows, phase1$cols],
        calib_pred = pred_karma[phase1$rows, phase1$cols],
        y = y_anomaly, pred = pred_karma,
        resid = karma2d_residuals(fit, newdata = y_anomaly,
                                  padding = "duplicate"),
        truth = truth, offset = 0
      ),
      gauss = list(
        model = gauss_model, variant = "Gaussian",
        calib_obs = y_anomaly[phase1$rows, phase1$cols],
        calib_pred = pred_gauss[phase1$rows, phase1$cols],
        y = y_anomaly, pred = pred_gauss,
        # as in the scripts first committed to the repository:
        # standardized by the standard deviation over the scene
        resid = r_gauss_raw / stats::sd(as.vector(r_gauss_raw)),
        truth = truth, offset = 0
      )
    )
  }
  metric_rows <- function(m) {
    rbind(
      data.frame(model = m$model, region = "Calibration Area",
                 t(karma2d_metrics(m$calib_obs, m$calib_pred))),
      data.frame(model = m$model, region = "Full Scene",
                 t(karma2d_metrics(m$y, m$pred)))
    )
  }
  table5 <- do.call(rbind, lapply(unname(models), metric_rows))
  utils::write.csv(table5, file.path(out_dir, "table5-metrics.csv"),
                   row.names = FALSE)

  # 6. Residuals and anomaly detection ----------------------------------------
  post_process <- function(mask, se) {
    switch(post_processing,
           opening = morph_open(mask, se),
           opening_closing = morph_close(morph_open(mask, se), se))
  }
  detect <- function(r, se) {
    det <- karma2d_detect_anomalies(r, 3, se, morphology = "none")
    list(threshold_mask = det$threshold_mask,
         detected_mask = post_process(det$threshold_mask, se))
  }
  detection_summary <- function(m, det, variant = m$variant) {
    d <- det$detected_mask
    truth <- m$truth
    data.frame(variant = variant, model = m$model,
               flagged_threshold = sum(det$threshold_mask),
               flagged_after_post_processing = sum(d),
               true_positive = sum(d & truth),
               false_positive = sum(d & !truth),
               recall = sum(d & truth) / sum(truth),
               precision = if (sum(d) > 0) sum(d & truth) / sum(d) else NA)
  }
  se3 <- config$structuring_element
  dets <- lapply(models, function(m) detect(m$resid, se3))
  detection <- do.call(rbind, lapply(names(models), function(k) {
    detection_summary(models[[k]], dets[[k]])
  }))
  if (procedure == "padded") {
    # The scripts first committed to the repository called EBImage::erode()/
    # dilate() without `kern`, i.e. with EBImage's default
    # makeBrush(5, shape = "diamond"); reported for traceability.
    diamond <- abs(row(matrix(0, 5, 5)) - 3) +
      abs(col(matrix(0, 5, 5)) - 3) <= 2
    detection <- rbind(detection, do.call(rbind, lapply(models, function(m) {
      detection_summary(m, detect(m$resid, diamond),
                        paste(m$variant, "(EBImage 5x5 diamond)"))
    })))
  }
  rownames(detection) <- NULL
  utils::write.csv(detection, file.path(out_dir, "anomaly-detection.csv"),
                   row.names = FALSE)

  # Figures
  mk <- models$karma
  mg <- models$gauss
  image_png(file.path(out_dir, "fig-a-image-with-anomaly.png"),
            to_display(mk$y), paste(config$name, "with synthetic anomaly"),
            zlim = display_range)
  image_png(file.path(out_dir, "fig-b-karma-prediction.png"),
            to_display(mk$pred), "2-D KARMA prediction",
            zlim = display_range)
  image_png(file.path(out_dir, "fig-c-gaussian-prediction.png"),
            to_display(mg$pred), "Gaussian 2-D ARMA prediction",
            zlim = display_range)
  div <- grDevices::hcl.colors(64, "Blue-Red 3")
  image_png(file.path(out_dir, "fig-residuals-karma.png"), mk$resid,
            "2-D KARMA quantile residuals", col = div)
  image_png(file.path(out_dir, "fig-residuals-gaussian.png"), mg$resid,
            "Gaussian standardized residuals", col = div)
  bw <- c("white", "black")
  image_png(file.path(out_dir, "fig-detection-karma-threshold.png"),
            dets$karma$threshold_mask * 1, "KARMA |r| > 3", zlim = 0:1,
            col = bw)
  image_png(file.path(out_dir, "fig-detection-karma-opened.png"),
            dets$karma$detected_mask * 1, "KARMA detected", zlim = 0:1,
            col = bw)
  image_png(file.path(out_dir, "fig-detection-gaussian-threshold.png"),
            dets$gauss$threshold_mask * 1, "Gaussian |r| > 3", zlim = 0:1,
            col = bw)
  image_png(file.path(out_dir, "fig-detection-gaussian-opened.png"),
            dets$gauss$detected_mask * 1, "Gaussian detected", zlim = 0:1,
            col = bw)
  if (!is.null(dets$gauss_published)) {
    grDevices::png(file.path(out_dir,
                             "fig-detection-gaussian-as-published.png"),
                   width = 1400, height = 700, res = 150)
    graphics::par(mfrow = c(1, 2))
    for (k in c("threshold_mask", "detected_mask")) {
      m <- dets$gauss_published[[k]] * 1
      graphics::image(t(m[rev(seq_len(nrow(m))), , drop = FALSE]),
                      zlim = 0:1, col = bw, axes = FALSE, asp = 1,
                      main = paste("Gaussian as published,",
                                   if (k == "threshold_mask") "|r| > 3" else
                                     "detected"))
    }
    grDevices::dev.off()
  }
  grDevices::png(file.path(out_dir, "fig-residual-histograms.png"),
                 width = 1400, height = 700, res = 150)
  graphics::par(mfrow = c(1, 2))
  phase1_part <- function(m) {
    m$resid[phase1$rows - m$offset, phase1$cols - m$offset]
  }
  graphics::hist(phase1_part(mk), breaks = 40,
                 main = "2-D KARMA quantile residuals (Phase I)", xlab = "")
  graphics::hist(phase1_part(mg), breaks = 40,
                 main = "Gaussian standardized residuals (Phase I)", xlab = "")
  grDevices::dev.off()

  # 7. Comparison with the published values ------------------------------------
  published <- config$published
  if (!is.null(published)) {
    computed <- c(
      stats::setNames(table4$Estimate, paste(table4$model, table4$parameter,
                                             "estimate")),
      stats::setNames(table4$`Std. Error`, paste(table4$model, table4$parameter,
                                                 "se")),
      stats::setNames(karma_table$se_legacy_information,
                      paste(karma_table$model, karma_table$parameter,
                            "se_legacy")),
      unlist(lapply(seq_len(nrow(table5)), function(i) {
        stats::setNames(unlist(table5[i, c("MAE", "MSE", "SMAPE", "RMSLE")]),
                        paste(table5$model[i], table5$region[i],
                              c("MAE", "MSE", "SMAPE", "RMSLE")))
      }))
    )
    published$computed <- unname(computed[published$quantity])
    published$abs_diff <- abs(published$computed - published$published)
    # Agreement: within one unit of the last decimal reported in the article.
    published$agrees <- published$abs_diff <= 10^(-published$decimals) + 1e-12
    utils::write.csv(published, file.path(out_dir, "comparison-with-paper.csv"),
                     row.names = FALSE)
  }
  write_session_info(out_dir)

  cat("Procedure:", procedure, "; post-processing:", post_processing, "\n")
  print(table4, digits = 4)
  print(table5, digits = 4)
  print(detection, digits = 3)
  if (!is.null(published)) {
    cat("\nComparison with the article:\n")
    print(published[, c("quantity", "published", "computed", "agrees")],
          digits = 5)
  }
  invisible(list(fit = fit, gaussian_fit = gfit, table4 = table4,
                 table5 = table5, detection = detection,
                 comparison = published, procedure = procedure))
}

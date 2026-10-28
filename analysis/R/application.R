# Shared pipeline of the NDVI and GRDI applications (Section 4.2).
#
# config: list with
#   name, y (image on the unit scale, no NA), phase1 (list rows, cols),
#   anomaly (list rows, cols, value), karma_order, gaussian_order,
#   structuring_element, select (logical), published (data frame or NULL),
#   out_dir.
run_application <- function(config) {
  y <- config$y
  phase1 <- config$phase1
  anomaly <- config$anomaly
  karma_order <- config$karma_order
  gaussian_order <- config$gaussian_order
  out_dir <- config$out_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  y_train <- y[phase1$rows, phase1$cols]
  y_anomaly <- y
  y_anomaly[anomaly$rows, anomaly$cols] <- anomaly$value
  to_display <- config$to_display
  display_range <- config$display_range

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
      pred <- predict(fit, newdata = y_full, padding = "duplicate")
      m <- karma2d_metrics(y_full[phase1$rows, phase1$cols],
                           pred[phase1$rows, phase1$cols])
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

  # 5. Prediction and accuracy --------------------------------------------------
  pred_karma <- predict(fit, newdata = y_anomaly, padding = "duplicate")
  pred_gauss <- gaussian_arma2d_predict(gfit, y_anomaly, padding = "duplicate")
  metric_rows <- function(model, pred) {
    rbind(
      data.frame(model = model, region = "Calibration Area",
                 t(karma2d_metrics(y_anomaly[phase1$rows, phase1$cols],
                                   pred[phase1$rows, phase1$cols]))),
      data.frame(model = model, region = "Full Scene",
                 t(karma2d_metrics(y_anomaly, pred)))
    )
  }
  table5 <- rbind(metric_rows(gauss_table$model[1], pred_gauss),
                  metric_rows(karma_table$model[1], pred_karma))
  utils::write.csv(table5, file.path(out_dir, "table5-metrics.csv"),
                   row.names = FALSE)

  # 6. Residuals and anomaly detection ----------------------------------------
  r_karma <- karma2d_residuals(fit, newdata = y_anomaly, padding = "duplicate")
  r_gauss_raw <- y_anomaly - pred_gauss
  r_gauss <- r_gauss_raw / stats::sd(as.vector(r_gauss_raw)) # as in the
  # original scripts: standardized by the standard deviation over the scene
  se3 <- config$structuring_element
  det_karma <- karma2d_detect_anomalies(r_karma, 3, se3)
  det_gauss <- karma2d_detect_anomalies(r_gauss, 3, se3)
  truth <- matrix(FALSE, nrow(y), ncol(y))
  truth[anomaly$rows, anomaly$cols] <- TRUE
  detection_summary <- function(model, det) {
    d <- det$detected_mask
    data.frame(model = model,
               flagged_threshold = sum(det$threshold_mask),
               flagged_after_opening = sum(d),
               true_positive = sum(d & truth),
               false_positive = sum(d & !truth),
               recall = sum(d & truth) / sum(truth),
               precision = if (sum(d) > 0) sum(d & truth) / sum(d) else NA)
  }
  # The original scripts called EBImage::erode()/dilate() without `kern`,
  # i.e. with EBImage's default makeBrush(5, shape = "diamond"), not the
  # element stated in the article; both are reported.
  diamond <- abs(row(matrix(0, 5, 5)) - 3) + abs(col(matrix(0, 5, 5)) - 3) <= 2
  det_karma_d <- karma2d_detect_anomalies(r_karma, 3, diamond)
  det_gauss_d <- karma2d_detect_anomalies(r_gauss, 3, diamond)
  detection <- rbind(
    cbind(element = "article", rbind(
      detection_summary(karma_table$model[1], det_karma),
      detection_summary(gauss_table$model[1], det_gauss))),
    cbind(element = "EBImage default 5x5 diamond", rbind(
      detection_summary(karma_table$model[1], det_karma_d),
      detection_summary(gauss_table$model[1], det_gauss_d)))
  )
  utils::write.csv(detection, file.path(out_dir, "anomaly-detection.csv"),
                   row.names = FALSE)

  # Figures
  image_png(file.path(out_dir, "fig-a-image-with-anomaly.png"),
            to_display(y_anomaly), paste(config$name, "with synthetic anomaly"),
            zlim = display_range)
  image_png(file.path(out_dir, "fig-b-karma-prediction.png"),
            to_display(pred_karma), "2-D KARMA prediction",
            zlim = display_range)
  image_png(file.path(out_dir, "fig-c-gaussian-prediction.png"),
            to_display(pred_gauss), "Gaussian 2-D ARMA prediction",
            zlim = display_range)
  div <- grDevices::hcl.colors(64, "Blue-Red 3")
  image_png(file.path(out_dir, "fig-residuals-karma.png"), r_karma,
            "2-D KARMA quantile residuals", col = div)
  image_png(file.path(out_dir, "fig-residuals-gaussian.png"), r_gauss,
            "Gaussian standardized residuals", col = div)
  bw <- c("white", "black")
  image_png(file.path(out_dir, "fig-detection-karma-threshold.png"),
            det_karma$threshold_mask * 1, "KARMA |r| > 3", zlim = 0:1, col = bw)
  image_png(file.path(out_dir, "fig-detection-karma-opened.png"),
            det_karma$detected_mask * 1, "KARMA detected", zlim = 0:1, col = bw)
  image_png(file.path(out_dir, "fig-detection-gaussian-threshold.png"),
            det_gauss$threshold_mask * 1, "Gaussian |r| > 3", zlim = 0:1,
            col = bw)
  image_png(file.path(out_dir, "fig-detection-gaussian-opened.png"),
            det_gauss$detected_mask * 1, "Gaussian detected", zlim = 0:1,
            col = bw)
  grDevices::png(file.path(out_dir, "fig-residual-histograms.png"),
                 width = 1400, height = 700, res = 150)
  graphics::par(mfrow = c(1, 2))
  graphics::hist(r_karma[phase1$rows, phase1$cols], breaks = 40,
                 main = "2-D KARMA quantile residuals (Phase I)", xlab = "")
  graphics::hist(r_gauss[phase1$rows, phase1$cols], breaks = 40,
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
                 comparison = published))
}

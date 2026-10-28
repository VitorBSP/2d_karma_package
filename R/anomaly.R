#' Binary mathematical morphology
#'
#' Erosion, dilation, opening and closing of a binary image by a structuring
#' element. The origin of the structuring element is its central pixel
#' (for an even dimension, the pixel `floor(d / 2) + 1`).
#'
#' Border convention: pixels outside the image never remove pixels during
#' erosion and never add pixels during dilation (equivalently, the image is
#' padded with `TRUE` for erosion and `FALSE` for dilation), so that the
#' opening of an object touching the border is not truncated artificially.
#'
#' @param mask Logical (or 0/1) matrix.
#' @param structuring_element Logical (or 0/1) matrix; see
#'   [karma2d_structuring_element()].
#' @return A logical matrix with the dimension of `mask`.
#' @examples
#' m <- matrix(FALSE, 9, 9)
#' m[3:6, 3:7] <- TRUE
#' m[9, 1] <- TRUE # isolated pixel
#' se <- karma2d_structuring_element(3, 3)
#' morph_open(m, se) # removes the isolated pixel, keeps the block
#' @name morphology
NULL

check_mask <- function(mask, name) {
  if (!is.matrix(mask) || anyNA(mask) ||
      !(is.logical(mask) || all(mask %in% c(0, 1)))) {
    stop(sprintf("`%s` must be a logical or 0/1 matrix without NA.", name),
         call. = FALSE)
  }
  m <- matrix(as.logical(mask), nrow(mask), ncol(mask))
  m
}

# Row/column offsets of the "on" pixels of a structuring element.
se_offsets <- function(se) {
  se <- check_mask(se, "structuring_element")
  if (!any(se)) stop("The structuring element is empty.", call. = FALSE)
  on <- which(se, arr.ind = TRUE)
  origin <- floor(dim(se) / 2) + 1
  cbind(on[, 1] - origin[1], on[, 2] - origin[2])
}

# mask shifted so that out[i, j] = mask[i + dr, j + dc] (fill outside).
shift_mask <- function(mask, dr, dc, fill) {
  n <- nrow(mask)
  k <- ncol(mask)
  out <- matrix(fill, n, k)
  rows <- seq_len(n)
  cols <- seq_len(k)
  src_r <- rows + dr
  src_c <- cols + dc
  ok_r <- src_r >= 1 & src_r <= n
  ok_c <- src_c >= 1 & src_c <= k
  out[rows[ok_r], cols[ok_c]] <- mask[src_r[ok_r], src_c[ok_c]]
  out
}

#' @rdname morphology
#' @export
morph_erode <- function(mask, structuring_element) {
  mask <- check_mask(mask, "mask")
  off <- se_offsets(structuring_element)
  out <- matrix(TRUE, nrow(mask), ncol(mask))
  for (s in seq_len(nrow(off))) {
    out <- out & shift_mask(mask, off[s, 1], off[s, 2], fill = TRUE)
  }
  out
}

#' @rdname morphology
#' @export
morph_dilate <- function(mask, structuring_element) {
  mask <- check_mask(mask, "mask")
  off <- se_offsets(structuring_element)
  out <- matrix(FALSE, nrow(mask), ncol(mask))
  for (s in seq_len(nrow(off))) {
    out <- out | shift_mask(mask, -off[s, 1], -off[s, 2], fill = FALSE)
  }
  out
}

#' @rdname morphology
#' @export
morph_open <- function(mask, structuring_element) {
  morph_dilate(morph_erode(mask, structuring_element), structuring_element)
}

#' @rdname morphology
#' @export
morph_close <- function(mask, structuring_element) {
  morph_erode(morph_dilate(mask, structuring_element), structuring_element)
}

#' Rectangular structuring element
#'
#' @param nrow,ncol Dimensions. The article uses a 3 x 3 square for the NDVI
#'   application and a 3 x 5 rectangle for the GRDI application.
#' @return A logical matrix of `TRUE` values.
#' @examples
#' karma2d_structuring_element(3, 5)
#' @export
karma2d_structuring_element <- function(nrow = 3L, ncol = nrow) {
  check_scalar_number(nrow, "nrow", positive = TRUE)
  check_scalar_number(ncol, "ncol", positive = TRUE)
  matrix(TRUE, as.integer(nrow), as.integer(ncol))
}

#' Residual-based anomaly detection
#'
#' Implements the detector of Section 3.2 of the article: quantile residuals
#' are compared with Shewhart-type control limits \eqn{\pm L} (default
#' \eqn{L = 3}), and the resulting binary map is optionally post-processed by
#' mathematical morphology (erosion followed by dilation, i.e. an opening).
#'
#' @param residuals A matrix of residuals (typically quantile residuals from
#'   [karma2d_residuals()]), or a `karma2d_fit` object, whose quantile
#'   residuals on its own data are used.
#' @param limit Control limit \eqn{L > 0}; pixels with \eqn{|r| > L} are
#'   flagged. Alternatively, a vector `c(lower, upper)` of asymmetric limits.
#' @param structuring_element Structuring element for the post-processing,
#'   see [karma2d_structuring_element()].
#' @param morphology `"opening"` (erosion then dilation, as in the article),
#'   `"closing"`, `"erosion"`, `"dilation"` or `"none"`.
#' @return An object of class `karma2d_anomaly`: a list with the residual
#'   matrix, the limits, the thresholded map `threshold_mask`, the
#'   post-processed map `detected_mask` (logical matrices; `NA` residuals are
#'   never flagged), the structuring element and the morphology used.
#' @examples
#' y <- karma2d_simulate(40, 40, alpha = 0.66, ar = c(0.28, 0.19, 0.17),
#'                       shape = 19.6, seed = 11)
#' fit <- karma2d_fit(y[21:40, 21:40], order = c(1, 0))
#' y_anom <- y
#' y_anom[5:9, 5:10] <- 0.3 # synthetic anomaly
#' r <- karma2d_residuals(fit, newdata = y_anom, padding = "duplicate")
#' det <- karma2d_detect_anomalies(r, limit = 3)
#' det
#' which(det$detected_mask, arr.ind = TRUE)[1:5, ]
#' @export
karma2d_detect_anomalies <- function(residuals, limit = 3,
                                     structuring_element =
                                       karma2d_structuring_element(3, 3),
                                     morphology = c("opening", "closing",
                                                    "erosion", "dilation",
                                                    "none")) {
  morphology <- match.arg(morphology)
  if (inherits(residuals, "karma2d_fit")) residuals <- residuals$residuals
  check_image(residuals, "residuals")
  if (!is.numeric(limit) || !length(limit) %in% 1:2 || anyNA(limit)) {
    stop("`limit` must be a positive number or a vector c(lower, upper).",
         call. = FALSE)
  }
  limits <- if (length(limit) == 1L) {
    if (limit <= 0) stop("`limit` must be positive.", call. = FALSE)
    c(-limit, limit)
  } else {
    sort(limit)
  }
  flagged <- !is.na(residuals) &
    (residuals < limits[1] | residuals > limits[2])
  threshold_mask <- matrix(flagged, nrow(residuals), ncol(residuals))
  detected <- switch(morphology,
    opening = morph_open(threshold_mask, structuring_element),
    closing = morph_close(threshold_mask, structuring_element),
    erosion = morph_erode(threshold_mask, structuring_element),
    dilation = morph_dilate(threshold_mask, structuring_element),
    none = threshold_mask
  )
  structure(
    list(residuals = residuals, limits = limits,
         threshold_mask = threshold_mask, detected_mask = detected,
         structuring_element = check_mask(structuring_element,
                                          "structuring_element"),
         morphology = morphology),
    class = "karma2d_anomaly"
  )
}

#' @export
print.karma2d_anomaly <- function(x, ...) {
  cat("2-D KARMA residual-based anomaly detection\n")
  cat(sprintf("Control limits: [%g, %g]\n", x$limits[1], x$limits[2]))
  cat(sprintf("Pixels outside the limits: %d of %d\n",
              sum(x$threshold_mask), length(x$threshold_mask)))
  cat(sprintf("Detected after %s (%d x %d element): %d\n", x$morphology,
              nrow(x$structuring_element), ncol(x$structuring_element),
              sum(x$detected_mask)))
  invisible(x)
}

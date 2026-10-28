# Internal helpers. None of these are exported.

# Recycle arguments to a common length following base R conventions:
# zero-length input gives zero-length output.
recycle_args <- function(...) {
  args <- list(...)
  lens <- lengths(args)
  if (any(lens == 0L)) {
    return(lapply(args, function(a) a[0]))
  }
  n <- max(lens)
  lapply(args, rep_len, length.out = n)
}

# log(1 - exp(a)) for a <= 0, accurate over the whole range
# (Maechler, 2012, "Accurately computing log(1 - exp(-|a|))").
log1mexp <- function(a) {
  out <- a
  small <- !is.na(a) & a > -log(2)
  out[small] <- log(-expm1(a[small]))
  out[!small] <- log1p(-exp(a[!small]))
  out
}

# u / log(1 - u) for u in [0, 1), with the limit -1 at u = 0.
# Used for mu^shape when it is tiny or underflows to zero.
ratio_u_log1mu <- function(u) {
  out <- u / log1p(-u)
  tiny <- !is.na(u) & u < 1e-10
  out[tiny] <- -1 + u[tiny] / 2
  out
}

check_scalar_number <- function(x, name, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("`%s` must be a single finite number.", name), call. = FALSE)
  }
  if (positive && x <= 0) {
    stop(sprintf("`%s` must be positive.", name), call. = FALSE)
  }
  invisible(x)
}

check_order <- function(order) {
  if (!is.numeric(order) || length(order) != 2L || anyNA(order) ||
        any(order < 0) || any(order != round(order))) {
    stop("`order` must be a vector c(p, q) of two non-negative integers.",
         call. = FALSE)
  }
  as.integer(order)
}

check_support <- function(lower, upper) {
  check_scalar_number(lower, "lower")
  check_scalar_number(upper, "upper")
  if (lower >= upper) {
    stop("`lower` must be smaller than `upper`.", call. = FALSE)
  }
  invisible(TRUE)
}

check_image <- function(y, name = "y") {
  if (!is.matrix(y) || !is.numeric(y)) {
    stop(sprintf("`%s` must be a numeric matrix.", name), call. = FALSE)
  }
  invisible(y)
}

# Rescale an image to the unit interval and check that every pixel lies
# strictly inside the support.
to_unit_scale <- function(y, lower, upper, name = "y") {
  check_image(y, name)
  if (anyNA(y) || any(!is.finite(y))) {
    stop(sprintf(
      "`%s` contains missing or non-finite values; impute or crop them first.",
      name
    ), call. = FALSE)
  }
  z <- (y - lower) / (upper - lower)
  if (any(z <= 0 | z >= 1)) {
    stop(sprintf(
      "All values of `%s` must lie strictly inside (lower, upper) = (%s, %s).",
      name, format(lower), format(upper)
    ), call. = FALSE)
  }
  z
}

supported_links <- c("logit", "probit", "cloglog", "cauchit")

make_karma2d_link <- function(link) {
  if (!is.character(link) || length(link) != 1L ||
        !link %in% supported_links) {
    stop(sprintf(
      "`link` must be one of %s.",
      paste0('"', supported_links, '"', collapse = ", ")
    ), call. = FALSE)
  }
  stats::make.link(link)
}

# Evaluate `code` with the given seed and restore the caller's RNG state on
# exit, so that functions taking a `seed` argument have no side effect on the
# global random number stream.
with_local_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code)
  }
  check_scalar_number(seed, "seed")
  global <- globalenv()
  had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = global, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = global) # nolint
    } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  })
  set.seed(seed)
  code
}

# Search the committed GRDI extraction for the window used in the article.
#
# Usage (from the repository root, after installing karma2d):
#   Rscript analysis/grdi/search_article_window.R [--cores=3]
#
# The authors state that analysis/data/grdi/nasa_corte_india.csv is an
# extraction of nasa_corte_india_expandido.csv, the file read by the original
# GRDI scripts (window raw[1200:1300, 2200:2300][1:100, 1:100]). The offset of
# the extraction is unknown, so every 100 x 100 window is tried:
#   1. Preprocess as the original scripts: divide by 100 and replace missing
#      values (coded 0.01 in this file) by the mean of the window.
#   2. Screen all windows by the statistics of the Phase I block
#      [51:100, 51:100] implied by Table 6: mean 0.742 (Gaussian intercept)
#      and standard deviation 0.0733 * sqrt(1 + sum(theta^2)) = 0.0865.
#   3. Rank the survivors by the spatial autocorrelations implied by the
#      published Gaussian MA(1): 0.380, 0.347, 0.187 at lags (0,1), (1,0),
#      (1,1).
#   4. Fit the Gaussian ARMA(0,1) and 2-D KARMA(0,1) to the 120 best
#      candidates and measure the distance to Table 6 in published standard
#      errors. The article's window would give a distance near zero.

if (!file.exists("analysis/R/common.R")) {
  stop("Run this script from the repository root.")
}
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "R", "gaussian_arma2d.R"))
opts <- parse_args(list(cores = "3"))

raw <- as.matrix(utils::read.csv(file.path("analysis", "data", "grdi",
                                           "nasa_corte_india.csv"),
                                 check.names = FALSE))
dimnames(raw) <- NULL
missing <- raw == 0.01
values <- ifelse(missing, 0, raw / 100)

# Sums over all h x h blocks from 2-D cumulative sums: block rows r+1..r+h.
cumsum2 <- function(a) {
  a <- rbind(0, cbind(0, a))
  t(apply(apply(a, 2, cumsum), 1, cumsum))
}
block <- function(s, r, c, h) {
  s[cbind(r + h + 1, c + h + 1)] - s[cbind(r + 1, c + h + 1)] -
    s[cbind(r + h + 1, c + 1)] + s[cbind(r + 1, c + 1)]
}
s1 <- cumsum2(values)
s2 <- cumsum2(values^2)
sm <- cumsum2(missing * 1)
r <- rep(seq_len(nrow(raw) - 99) - 1, ncol(raw) - 99)
c <- rep(seq_len(ncol(raw) - 99) - 1, each = nrow(raw) - 99)
fill <- block(s1, r, c, 100) / (1e4 - block(sm, r, c, 100))
n_fill <- block(sm, r + 50, c + 50, 50)
sum1 <- block(s1, r + 50, c + 50, 50) + n_fill * fill
sum2 <- block(s2, r + 50, c + 50, 50) + n_fill * fill^2
mean1 <- sum1 / 2500
sd1 <- sqrt(pmax(sum2 - 2500 * mean1^2, 0) / 2499)
windows <- data.frame(row = r + 1, col = c + 1, mean = mean1, sd = sd1)
cat("Windows tried:", nrow(windows), "\n")
cand <- windows[abs(windows$mean - 0.742) < 0.005 &
                  abs(windows$sd - 0.0865) < 0.01, ]
cat("Passing the mean/sd screen:", nrow(cand), "\n")

phase1 <- function(row, col) {
  w <- raw[row:(row + 99), col:(col + 99)]
  miss <- w == 0.01
  w <- w / 100
  w[miss] <- mean(w[!miss])
  w[51:100, 51:100]
}
acf3 <- function(p) {
  z <- p - mean(p)
  v <- mean(z^2)
  c(mean(z[, -1] * z[, -50]), mean(z[-1, ] * z[-50, ]),
    mean(z[-1, -1] * z[-50, -50])) / v
}
theta <- c(0.4329, 0.3709, 0.2606)
target <- c(theta[1] + theta[2] * theta[3], theta[2] + theta[1] * theta[3],
            theta[3]) / (1 + sum(theta^2))
acfs <- t(mapply(function(a, b) acf3(phase1(a, b)), cand$row, cand$col))
cand$acf_distance <- sqrt(rowSums(sweep(acfs, 2, target)^2))
cand$screen_score <- abs(cand$mean - 0.742) / 0.005 +
  abs(cand$sd - 0.0865) / 0.01
pick <- unique(rbind(head(cand[order(cand$acf_distance), ], 60),
                     head(cand[order(cand$screen_score), ], 60)))

pub_gauss <- c(0.7420, 0.4329, 0.3709, 0.2606, 0.0733)
se_gauss <- c(0.0030, 0.0166, 0.0178, 0.0190, 0.0011)
pub_karma <- c(1.1114, 0.3467, 0.2819, 0.1848, 12.3857)
se_karma <- c(0.0139, 0.0161, 0.0163, 0.0165, 0.2221)
fits <- parallel::mclapply(seq_len(nrow(pick)), function(i) {
  p <- phase1(pick$row[i], pick$col[i])
  g <- gaussian_arma2d_fit(p, c(0L, 1L))$coefficients
  k <- coef(suppressWarnings(karma2d_fit(p, c(0, 1))))
  c(gauss_se_units = max(abs(g - pub_gauss) / se_gauss),
    karma_se_units = max(abs(k - pub_karma) / se_karma))
}, mc.cores = as.integer(opts$cores))
pick <- cbind(pick, do.call(rbind, fits))
pick$distance <- pmax(pick$gauss_se_units, pick$karma_se_units)
pick <- pick[order(pick$distance), ]
cat("\nBest candidates (distance = largest deviation from Table 6, in",
    "published standard errors):\n")
print(head(pick[, c("row", "col", "mean", "sd", "gauss_se_units",
                    "karma_se_units", "distance")], 5),
      digits = 3, row.names = FALSE)
cat(sprintf(paste0("\nSmallest distance: %.1f standard errors. The article's",
                   " window would give about 0.\n"), pick$distance[1]))

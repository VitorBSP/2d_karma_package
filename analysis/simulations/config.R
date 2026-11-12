# Monte Carlo configuration (Section 4.1 of the article).
#
# Two scenario sets are provided:
#   "paper":  true parameters printed in Tables 1-3 of the article. The
#             article states that the KARMA(1,0) and KARMA(0,1) values were
#             taken from the NDVI and GRDI fits (Tables 4 and 6), which is
#             consistent with the printed values.
#   "legacy_scripts": the values hard-coded in the original scripts
#             analysis/legacy/simulations/script_2dkarma_{10,01,11}.R. They
#             differ from Tables 1-2 for KARMA(1,0) and KARMA(0,1); KARMA(1,1)
#             is identical. Kept so the discrepancy can be studied; the
#             article's tables are the reference.

scenario_sets <- list(
  paper = list(
    karma10 = list(alpha = 0.6609, ar = c(0.2819, 0.1897, 0.1653),
                   ma = numeric(0), shape = 19.6065, table = "Table 1"),
    karma01 = list(alpha = 1.1114, ar = numeric(0),
                   ma = c(0.3467, 0.2819, 0.1848), shape = 12.3857,
                   table = "Table 2"),
    karma11 = list(alpha = 0.3, ar = c(0.3, -0.15, 0.23),
                   ma = c(0.15, -0.1, 0.25), shape = 17, table = "Table 3")
  ),
  legacy_scripts = list(
    karma10 = list(alpha = 0.67, ar = c(0.31, 0.21, -0.15), ma = numeric(0),
                   shape = 23, table = NA),
    karma01 = list(alpha = 0.12, ar = numeric(0), ma = c(0.1, 0.02, -0.12),
                   shape = 11, table = NA),
    karma11 = list(alpha = 0.3, ar = c(0.3, -0.15, 0.23),
                   ma = c(0.15, -0.1, 0.25), shape = 17, table = NA)
  )
)

# quick: development and CI (minutes); full: configuration of the article.
modes <- list(
  quick = list(replications = 20L, sizes = c(10L, 25L)),
  full = list(replications = 10000L, sizes = c(10L, 25L, 50L))
)

confidence_level <- 0.95

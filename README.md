# karma2d: Two-Dimensional Kumaraswamy ARMA Models for Bounded Image Data

`karma2d` implements the two-dimensional Kumaraswamy autoregressive moving
average (2-D KARMA) model of Pereira, Palm and Bayer (2026) for images whose
pixel values are bounded in an interval, such as vegetation indices,
reflectances or rates. Conditionally on a strongly causal neighbourhood, each
pixel follows a Kumaraswamy distribution with median
$\mu[n,m]$ and shape $\varphi$, and the median follows a spatial ARMA
structure through a link function $g$:

$$
g(\mu[n,m]) = \alpha + \sum_{(i,j)} \phi_{(i,j)}\, g(y[n-i,m-j])
  + \sum_{(k,l)} \theta_{(k,l)}\, e[n-k,m-l].
$$

The package provides conditional maximum likelihood estimation with analytic
score and information, simulation, prediction, quantile residuals, and
residual-based anomaly detection with morphological post-processing. It is
written in pure R and depends only on `stats`.

## Development status

This is a pre-release version (0.1.0.9000). It is **not on CRAN**, the API may
change without notice, and CRAN compatibility is not claimed. The accompanying
article is an unpublished manuscript.

## Installation

```r
# install.packages("remotes")
remotes::install_github("VitorBSP/2d_karma_package")
```

## Example

```r
library(karma2d)

# Simulate a KARMA(1,0) image and fit the model
y <- karma2d_simulate(60, 60, alpha = 0.66, ar = c(0.28, 0.19, 0.17),
                      shape = 19.6, seed = 42)
fit <- karma2d_fit(y, order = c(1, 0))
summary(fit)

# Fitted medians (first max(p, q) rows/columns are NA with padding = "none")
mu_hat <- predict(fit)

# Quantile residuals
res <- karma2d_residuals(fit)

# Anomaly detection on the NDVI image (natural scale -1 to 1):
# fit on the Phase I window, inject a synthetic anomaly, monitor the scene
data(ndvi_paraiba)
fit_ndvi <- karma2d_fit(ndvi_paraiba[51:100, 51:100], order = c(1, 0),
                        lower = -1, upper = 1)
scene <- ndvi_paraiba
scene[5:20, 15:30] <- -0.4
res_full <- karma2d_residuals(fit_ndvi, newdata = scene,
                              padding = "duplicate")
det <- karma2d_detect_anomalies(res_full, limit = 3,
                                structuring_element =
                                  karma2d_structuring_element(3, 3))
det
```

The summary reports estimates close to the true values (for this seed,
`alpha` = 0.659, shape = 19.29) with standard errors from the conditional
expected information. Run the code to see the full output.

## What the package provides

| Function(s) | Purpose |
|---|---|
| `dkuma()`, `pkuma()`, `qkuma()`, `rkuma()` | Kumaraswamy distribution in the median parametrization |
| `karma2d_lags()` | Ordering and labels of the lag coefficients for a given order |
| `karma2d_simulate()` | Simulate an image from a 2-D KARMA model |
| `karma2d_fit()`, `karma2d_control()`, `karma2d_start()` | Conditional maximum likelihood fit, optimizer settings, starting values |
| `print()`, `summary()`, `coef()`, `vcov()`, `logLik()`, `nobs()`, `confint()` | Standard S3 methods for `karma2d_fit` (`confint()` through the default method) |
| `predict()`, `fitted()` | Filtered medians or linear predictor; border handling `padding = "none"` or `"duplicate"` |
| `residuals()`, `karma2d_residuals()` | Quantile residuals, optionally on new data |
| `karma2d_fisher()` | Conditional expected or observed information of a fit |
| `karma2d_loglik()`, `karma2d_score()`, `karma2d_information()` | Log-likelihood, analytic score and information at given parameters |
| `karma2d_diagnostics()` | Residual diagnostics |
| `karma2d_detect_anomalies()` | Control-chart anomaly detection with morphological post-processing |
| `karma2d_structuring_element()`, `morph_erode()`, `morph_dilate()`, `morph_open()`, `morph_close()` | Base-R binary morphology |
| `karma2d_metrics()` | Image quality metrics between observed and predicted images |
| `ndvi_paraiba` | NDVI image used in the article |

## Corrections relative to the article

1. **Fisher information (Appendix B).** Appendix B of the article has sign
   errors in the conditional Fisher information: $K_{\mu\varphi}$ has the
   opposite sign, and two terms of $K_{\varphi\varphi}$ have the wrong sign,
   which can make $K_{\varphi\varphi}$ negative. The package implements a
   re-derived expected information, verified against numerical integration.
2. **Published standard errors.** The information matrix in the original
   scripts used the observed $\varphi\varphi$ entry (mixed with expected
   entries elsewhere). This is the matrix that produced the published standard
   errors.
3. **Numerical precision of the Kumaraswamy functions.** The original
   functions lost precision when $\mu^\varphi$ was tiny (for example,
   `qkuma(0.5)` returned 0.200012 instead of the median 0.2). The package uses
   `log1p`/`expm1` to fix this.

## Reproducing the article

See `vignette("reproducing-the-paper")` and the `analysis/` folder of the
repository (excluded from the package build). In summary:

- **NDVI, Table 4:** reproduced, except the shape standard error. The published
  value (0.2855) comes from the original mixed information matrix; the
  corrected value is 0.356.
- **NDVI, Table 5:** only partly reproduced; the cause is unresolved.
- **GRDI:** not reproducible, because the input data are not in the repository.
- **Monte Carlo study:** a pipeline is provided with quick and full modes.

## Limitations

- Only isotropic orders $(p, q)$ are supported; there is no subset or sparse
  lag selection.
- No stationarity or invertibility conditions are proven for the model.
- The first $\max(p, q)$ rows and columns (borders) are conditioned on, not
  modelled.
- Large images are slow in pure R; a 50 x 50 ARMA(1,1) fit takes about 1.5 s.

## Citation

```r
citation("karma2d")
```

## License

MIT, see `LICENSE`.
<!-- The license choice must be confirmed by the authors. -->

## References

- Pereira, V. B. S., Palm, B. G. and Bayer, F. M. (2026). A Kumaraswamy-Based
  Spatial Model for Bounded Image Data. Manuscript.
- Bayer, F. M., Bayer, D. M. and Pumi, G. (2017). Kumaraswamy autoregressive
  moving average models for double bounded environmental data. *Journal of
  Hydrology*, 555, 385-396.
- Palm, B. G., Bayer, F. M. and Cintra, R. J. (2022). 2-D Rayleigh
  autoregressive moving average model for SAR image modeling. *Computational
  Statistics & Data Analysis*, 171, 107453.
- Mitnik, P. A. and Baek, S. (2013). The Kumaraswamy distribution:
  median-dispersion re-parameterizations for regression modeling and
  simulation-based estimation. *Statistical Papers*, 54, 177-192.
- Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals. *Journal
  of Computational and Graphical Statistics*, 5, 236-244.

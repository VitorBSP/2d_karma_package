# karma2d 0.1.0.9000

First pre-release.

* Kumaraswamy distribution in the median parametrization: `dkuma()`,
  `pkuma()`, `qkuma()`, `rkuma()`.
* 2-D KARMA simulation (`karma2d_simulate()`) and conditional maximum
  likelihood fitting (`karma2d_fit()`, `karma2d_control()`, `karma2d_start()`)
  with analytic score, for isotropic orders `c(p, q)` and several links.
* S3 methods for fits: `print()`, `summary()`, `coef()`, `vcov()`, `logLik()`,
  `nobs()`, `predict()` (with `padding = "none"` or `"duplicate"`), `fitted()`
  and `residuals()`; `confint()` works through the default method.
* Quantile residuals (`karma2d_residuals()`), conditional expected and observed
  information (`karma2d_fisher()`, `karma2d_information()`), log-likelihood and
  score (`karma2d_loglik()`, `karma2d_score()`), and diagnostics
  (`karma2d_diagnostics()`).
* Residual-based anomaly detection with base-R morphology
  (`karma2d_detect_anomalies()`, `karma2d_structuring_element()`,
  `morph_erode()`, `morph_dilate()`, `morph_open()`, `morph_close()`) and image
  quality metrics (`karma2d_metrics()`).
* Dataset `ndvi_paraiba`.
* Corrections relative to the article:
    * The conditional Fisher information of Appendix B has sign errors
      ($K_{\mu\varphi}$ has the opposite sign; two terms of $K_{\varphi\varphi}$
      have the wrong sign). A re-derived expected information, verified against
      numerical integration, is implemented.
    * The original scripts used the observed $\varphi\varphi$ entry in the
      information matrix; this produced the published standard errors.
    * The original Kumaraswamy functions lost precision when $\mu^\varphi$ was
      tiny (e.g. `qkuma(0.5)` returned 0.200012 instead of 0.2); fixed with
      `log1p`/`expm1`.

## Changes relative to the original research scripts

* The legacy scripts are kept unchanged in `analysis/legacy`.
* `rgdal`/`raster` were replaced by `terra` (the NDVI matrix is identical).
* `EBImage` was replaced by base-R morphology. EBImage's default structuring
  element is a 5x5 diamond, whereas the article states a 3x3 square / 3x5
  rectangle; the pipelines use the article's elements.
* The subset-lag feature of the original simulator was removed because its
  index arithmetic selected the wrong lags.
* The shape parameter is optimized on the log scale by default.
* The default `reltol` is 1e-10.

# Technical audit of the 2-D KARMA code base

This document records the audit of the original research code, the
verification of its mathematics against the article, the decisions taken when
building the `karma2d` package, and the remaining open questions. Every
numerical statement below was produced by code in this repository; the
commands are given so that each can be re-checked.

* Article: V. B. S. Pereira, B. G. Palm and F. M. Bayer, *A Kumaraswamy-Based
  Spatial Model for Bounded Image Data* (manuscript, 34 pages, file
  `2_D_KARMA.pdf`). The PDF was available during the audit; the text was
  extracted with `pdftotext`, and the pages of Appendix B were also rendered
  as images and read directly, because the formula layout does not survive
  text extraction. The manuscript still contains co-author review notes
  (e.g. "{Bayer: FICOU REPETIDO AQUIIII.}"), so it is a draft.
* Original code: commit `7aa688c` of this repository. It is preserved,
  unmodified, under `analysis/legacy/` (relocated by a pure `git mv`).
* Environment of the audit: Ubuntu 24.04, R 4.3.3, numDeriv 2016.8-1.1,
  terra 1.7-65, raster 3.6, testthat 3.2.1.

## 1. Inventory of the original repository

### 1.1 Layout

| Original path | Content | New location |
|---|---|---|
| `2DKARMA/functions/` | 12 R files with the model functions | `analysis/legacy/functions/` |
| `2DKARMA/src/functions.R` | Kumaraswamy d/p/q/r functions | `analysis/legacy/src/` |
| `2DKARMA/src/kumaraswamy-distribution.cpp` | C++ Kumaraswamy, never compiled | `analysis/legacy/src/` |
| `2DKARMA/ndvi/ndvi/*.tif` (+ `.aux.xml`) | 33 MODIS NDVI rasters, 100 x 100 | `analysis/data/ndvi/` |
| `experimentos_numericos/` | Monte Carlo scripts, `testes.R`, 4 PNGs | `analysis/legacy/simulations/` |
| `ndvi/` | NDVI application, KARMA and Gaussian | `analysis/legacy/ndvi/` |
| `grdi/*.R` | GRDI application, 5 variants | `analysis/legacy/grdi/` |
| `grdi/data/nasa_corte_india.csv` | 500 x 1001 GRDI crop (unused by any script) | `analysis/data/grdi/` |

There was no `DESCRIPTION`, `NAMESPACE`, test, documentation or licence. All
scripts used `source()` with paths relative to the repository root and
`library()` calls at the top; several relied on objects defined in the global
environment (for example `m`, `y`, `link` inside `escore()`).

### 1.2 Functions

| Function (file) | Responsibility | Status after audit |
|---|---|---|
| `dkuma`, `pkuma`, `qkuma`, `rkuma` (`src/functions.R`) | Kumaraswamy, median parametrization | Correct on (0, 1) except precision loss (Sec. 3.1); `a`, `b` handled inconsistently |
| C++ `cpp_dkumar` etc. (`src/kumaraswamy-distribution.cpp`) | Kumaraswamy in the **standard** (a, b) parametrization | Not usable: `#include "shared.h"` (from the `extraDistr` package) is absent, and `GETV`, `VALID_PROB`, `rng_unif` are defined there, so the file cannot compile. It is also a different parametrization. Not carried over. |
| `loglik` (`log_vero_2dkarma.R`) | Conditional log-likelihood | Correct (Sec. 3.3); link hard-coded to logit |
| `build_2dkarma_quantities` (`2dkarma_base.R`) | Recursion for eta, errors and d eta / d gamma | Correct; reference for the new code |
| `escore_reformed` (`new_escore_2dkarma.R`) | Score from `build_2dkarma_quantities` | Correct (matches numerical derivatives) |
| `escore` (`escore_reformed_2dkarma.R`) | Score, explicit loops | Correct, identical to the above |
| `escore` (`escore_reformed_2dkarma (cópia).R`) | Older copy | Uses undefined `m`, `ar`/`ma`; column-major flattening of `deta.dalpha` mismatched with row-major elsewhere. Obsolete. |
| `escore` (`escore.R`, `escore_2dkarma.R`) | Early scores | Broken: undefined objects (`m`, `y`, `link`, `Ualpha`), derivative recursion averages (`apply(..., mean)`) over flattened rows, Rayleigh leftovers (`dmu` of the RARMA model). Not usable. |
| `zerar_fora_diagonal` (`escore.R`) | Zero the off-diagonal of a matrix | Unused |
| `fisher_information_2dkarma` (`fisher_2dkarma.R`) | "Fisher" information | gamma-gamma and gamma-shape blocks are the expected information with the correct sign; the shape-shape entry is the **observed** information (Sec. 3.5) |
| `fisher_information_2dkarma_vectorized` | Same, block-wise | Identical to the above to 3e-12; its `log_prec_scale` option multiplies by the Jacobian correctly |
| `simulacao` (`simu_2dkarma.R` and copies in the MC scripts) | Simulation | Correct for full orders; subset-order feature wrong (Sec. 3.2); prints every pixel; calls `source()` |
| `loglik_gaussiana`, `escore_gaussiana` (`log_vero_arma.R`, `escore_arma.R`) | Gaussian 2-D ARMA benchmark | The gradient uses the lag window `(i-p):(i-1), (j-p):(j-1)`, a different neighbourhood from the log-likelihood (`(i-p):i, (j-p):j` minus the centre), so it is not the gradient of the objective |
| `fit_with_fisher_polishing`, `fit_with_fisher_puro_core` (`testes.R`) | Fisher scoring experiments | Add a ridge to the information silently; not carried over |
| `arma_mle_inference` (application scripts) | Standard errors from `-optim()$hessian` | On failure adds a ridge (`-min(eigen) + 1e-8`) without telling the user; not carried over |
| `preencher_com_media`, `mse`, `mae`, `mape`, `r_2`, `rmsle`, plotting helpers (application scripts) | Utilities | `mape` is MAPE, not the SMAPE reported in the article; in `aplicacao_grdi_arma_gaussiana_v2.R` it is `abs(verdadeiro - predito / verdadeiro)` (operator precedence bug) |

### 1.3 Duplicated or conflicting implementations

* Five definitions named `escore` (three files plus copies inside each Monte
  Carlo script); two correct, three broken.
* Two Fisher implementations (equivalent).
* `simulacao`, `dkuma`, `loglik`, `escore` are pasted into every Monte Carlo
  script, so the scripts do not use the files in `2DKARMA/functions/`.
* Each application has a "v1" and a "v2" script differing only by window
  sizes; the GRDI Gaussian model has three variants.

### 1.4 External dependencies of the original code

`numDeriv`, `tidyverse`, `moments`, `fastrep` (Monte Carlo); `abind`,
`raster`, `rgdal`, `EBImage`, `tidyverse`, `reshape2`, `readr`, `MASS`
(applications). `rgdal` was archived on CRAN in October 2023. It was loaded but
not called directly; `raster::stack()` read the GeoTIFFs. The package now
uses `terra`, and `terra::as.matrix(terra::rast(f), wide = TRUE)` returns a
matrix **identical** to `raster::as.matrix(raster::stack(files)[[32]])`
(maximum absolute difference 0, checked with raster 3.6 and terra 1.7-65).
`EBImage` (Bioconductor) was replaced by a base-R binary morphology
(Sec. 6). The installed package depends only on `stats`.

## 2. Package design

* Single package root at the repository root (`DESCRIPTION`, `R/`, `man/`,
  `tests/`, `vignettes/`, `data/`, `data-raw/`, `inst/`). `analysis/`, `docs/`
  and `data-raw/` are excluded from the build.
* Pure R. The recursion was profiled before deciding on C++: a 50 x 50
  ARMA(1,1) fit takes about 1.5 s. That is fast enough for the package. The
  only use case where C++ would matter is the full Monte Carlo study, which
  runs offline and in parallel. The C++ file in the original repository was
  not usable (Sec. 1.2).
* Public API: `karma2d_fit()` (class `karma2d_fit` with `print`, `summary`,
  `coef`, `vcov`, `logLik`, `nobs`, `predict`, `fitted`, `residuals`;
  `confint`, `AIC`, `BIC` through the defaults), `karma2d_simulate()`,
  `karma2d_residuals()` (residuals on new data), `karma2d_detect_anomalies()`,
  `karma2d_fisher()`, `karma2d_diagnostics()`, `karma2d_metrics()`, the
  Kumaraswamy functions, `karma2d_lags()`, low-level `karma2d_loglik()`,
  `karma2d_score()`, `karma2d_information()`, morphology helpers and the
  dataset `ndvi_paraiba`. A separate `karma2d_predict()` was not added because
  `predict()` covers it.

## 3. Mathematical audit

### 3.1 Kumaraswamy distribution

Article Eq. (1): density
`phi * delta * y^(phi-1) * (1-y^phi)^(delta-1)` with
`delta = log(0.5) / log(1 - mu^phi)`, CDF `1 - (1-y^phi)^delta`.
The legacy formulas agree. Problems found:

* `log(1 - mu^phi)` loses precision when `mu^phi` is tiny. For `mu = 0.2,
  phi = 19.6` (`mu^phi = 2e-14`) the legacy `qkuma(0.5)` returns 0.200012
  instead of the median 0.2, and the density is off by 0.11 %. The new code
  uses `log1p()`, and `log(1 - y^phi)` is computed as `log1mexp(phi log y)`,
  which is also accurate for `y` within 1e-16 of 1.
* Legacy `dkuma(x, mu, phi, a, b)` divides by `b - a` but does not rescale
  `x`; `pkuma` checks `0 < x < 1` whatever `a`, `b`. In the new interface
  `x`/`q` **and** `median` are on the original scale `(lower, upper)`.
* Legacy `dkuma`/`pkuma` `stop()` if any value is outside (0, 1). The new
  functions follow R conventions (0 outside the support, `NaN` with a warning
  for invalid parameters, `NA` propagation, recycling, `log`, `lower.tail`,
  `log.p`, dimensions of matrix arguments kept).

Validated by tests: normalization (quadrature), CDF = integral of the
density, `F(median) = 0.5`, `F(F^-1(p)) = p` down to `p = 1e-10`,
monotonicity, tail identities, location-scale relation, KS test of `rkuma()`,
and agreement with the legacy functions wherever `mu^phi > 1e-8`.

### 3.2 Spatial structure

Article Eq. (4): `eta[n,m] = alpha + sum phi(i,j) g(y[n-i,m-j]) + sum
theta(k,l) e[n-k,m-l]`, `(i,j) in {0..p}^2 \ (0,0)`, coefficients ordered
`phi(0,1), phi(0,2), ..., phi(p,p)`.

The legacy code builds the lag vector as
`rev(as.vector(t(y[(i-p):i, (j-p):j])))[-1]`, a row-major flattening of the
window reversed. For `p = 1` this gives `(y[i,j-1], y[i-1,j], y[i-1,j-1])`,
i.e. `phi(0,1), phi(1,0), phi(1,1)` with **i the row lag and j the column
lag**; for `p = 2` it gives `(0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1),
(2,2)`. This matches the article's ordering and the derivative loops
(`for ii in 0:p, for jj in 0:p`). `karma2d_lags()` documents it, and unit
tests check single-coefficient predictors and a hand-written AR(1) predictor.
The counts are `(p + 1)^2 - 1` and `(q + 1)^2 - 1`.

Border convention (unchanged): pixels in the first `w = max(p, q)` rows or
columns are not modelled; their errors and derivatives are 0, and AR lags of
modelled pixels use the observed border values.

Subset orders (`AR = c(1, 3)` etc.) in the legacy simulator and predictor are
**incorrect**: they assume the coefficients of "ring" `k` occupy positions
`k^2 ... (k+1)^2 - 1` (simulator) or `2*sum(1:(k-1)) + k ...` (predictor), but
in the lexicographic ordering the lags of a given ring are not contiguous
(for `p = 2`, ring 1 is at positions 1, 3, 4). The parameter-count check
`sum(AR * 2 + 1)` is also wrong. The feature was removed; only full isotropic
orders are supported.

Simulation: the legacy scheme (grid extended by `2w`, first `w` rows and
columns initialized at `g(y) = alpha`, `e = 0`, sequential row-major draws by
inversion, first `2w` rows/columns discarded) is kept as the default of
`karma2d_simulate()` (`burn_in = 2w`). With the same seed the new simulator
reproduces the legacy images to 6e-13 (test `test-legacy-baseline.R`). The
legacy `print()` of every pixel and the `source()` call were removed. `seed`
is applied locally and the global RNG state is restored.

### 3.3 Conditional likelihood and optimization

Eq. (5) sums `log phi + log delta + (phi-1) log y + (delta-1) log(1-y^phi)`
over `n, m > w`. The legacy `loglik` and the new `karma2d_loglik()` agree to
3e-14. Sign conventions: the legacy code maximized with
`control = list(fnscale = -1)`; the package minimizes the negative
log-likelihood with the negative score. Invalid trial points get an infinite
objective (BFGS backtracks).

The shape parameter is optimized on the log scale by default
(`karma2d_control(shape_transform = "log")`): the gradient component is
multiplied by `phi` (chain rule) and estimates are returned on the original
scale. The information matrices and covariance are always computed for the
original parametrization, so no Jacobian is needed there. The identity
transform of the article remains available and reaches the same maximum
(test: log-likelihoods equal to 1e-6).

Starting values follow Sec. 2.3 of the article: `alpha0 = g(ybar)`, zero AR
and MA coefficients, `phi0 = 0.3 * (ybar (1 - ybar) / S^2 - 1)`. The Monte
Carlo scripts used exactly this; the application scripts used
`mean(g(y))` instead of `g(ybar)` for `alpha0`. This only changes the path,
not the maximum.

`optim()`'s default `reltol` (about 1.5e-8) stopped about 1e-4 away from the
maximum on the NDVI window (log-likelihood 4132). The package default is
`reltol = 1e-10`; with it the NDVI estimates reproduce the article to the
last printed digit.

### 3.4 Score vector

The two correct legacy implementations were compared with `numDeriv::grad()`
at the true parameters on simulated images: relative errors of 1e-8 (AR, MA)
and 3e-6 (ARMA, numDeriv error). The new `karma2d_score()` agrees with them
to 5e-11. Tests repeat the numDeriv comparison for AR(1), MA(1), ARMA(1,1),
ARMA(2,1), parameters away from the truth, and the probit, cloglog and
cauchit links. The recursions match Appendix A, including the MA
derivative `e[n-k,m-l] - sum theta(s,t) d eta[n-s,m-t] / d theta(k,l)` and
the shape score `1/phi + log y + c mu log mu - (delta-1) y^phi log y /
(1-y^phi)`.

### 3.5 Fisher information

Notation per pixel: `u = mu^phi`, `L = log(1-u)`, `delta = log(0.5)/L`,
`A = phi mu^(phi-1) / ((1-u) L)`, `B = u log(mu) / ((1-u) L)`,
`h(delta) = (psi(2) - psi(delta+1)) / (delta-1)`,
`t(delta) = [(psi(delta) - psi(2))^2 + psi'(2) - psi'(delta)] / (delta-2)`.

Derivation. Under the model, `E = -delta log(1 - Y^phi) ~ Exp(1)` and
`X = Y^phi ~ Beta(1, delta)`. The pixel scores are
`d l / d mu = A (1 - E)` and
`d l / d phi = 1/phi + B (1 - E) + [log X - (delta-1) X log X / (1-X)] / phi`.
Using `E[g(X)(1 - E)] = delta d/d delta E_delta[g(X)]` and the Beta moments
`E[X log X / (1-X)] = h(delta)` and
`(delta-1) E[X log^2 X / (1-X)^2] = delta t(delta)`, the pixel expected
information is

```
K_mu,mu   = A^2
K_mu,phi  = A (B + delta h / phi)
K_phi,phi = 1/phi^2 + B^2 + 2 delta B h / phi + delta t / phi^2
```

and `K = sum D (A^2 mu_eta^2) D'` (gamma block), `sum D (K_mu,phi mu_eta)`
(cross block), `sum K_phi,phi`, over the effective region, with
`D = d eta / d gamma` and `mu_eta = d mu / d eta`.

Numerical check (script in `tests/testthat/helper-karma2d.R` and
`test-likelihood.R`): the pixel information was integrated numerically, with
the substitution `t = E` on `(0, Inf)`, for
`(mu, phi) = (0.5, 3), (0.8, 17), (0.66, 19.6), (0.3, 1.5), (0.75, 12.4),
(0.9, 2), (0.95, 1.2)`, i.e. `delta` from 0.25 to 2400. Agreement is to 8
significant digits for all three entries. At the removable singularities
`delta = 1` and `delta = 2` the closed form is evaluated as the average at
`delta +/- 1e-4`; agreement with integration remains at 1e-8.

Comparison with Appendix B of the article (writing
`lambda_k = mu^(k phi - 2) / ((1-u)^k L^k)`, so `A = phi mu lambda_1` and
`A^2 = phi^2 lambda_2`):

| Entry | Article | Verified | Effect |
|---|---|---|---|
| gamma-gamma | `phi^2 lambda_2 / g'(mu)^2 D D'` | same | none |
| gamma-phi | `d / g'(mu) D` with `d = -phi mu [delta lambda_1 (1-psi(delta+1)-kappa)/((delta-1) phi) + log(mu) lambda_2]` | `-d / g'(mu) D` | Wrong sign of the cross covariances. Standard errors are unaffected (flipping the sign of one row/column of `K` is the reparametrization `phi -> -phi`). Joint Wald tests involving `phi` are affected. |
| phi-phi | `1/phi^2 - mu^2 lambda_2 log^2 mu + 2 delta mu^2 log mu lambda_1 (...) - delta [psi(delta)(psi(delta)+2(kappa-1)) - psi'(delta) + k_0] / ((delta-2) phi^2)` | `+ mu^2 lambda_2 log^2 mu` and `+ delta [...] / ((delta-2) phi^2)` | Two sign errors. The article's expression is negative at every point checked (e.g. -1.552 instead of 0.161 at `mu = 0.5, phi = 3`; -0.165 instead of 0.00686 at `mu = 0.8, phi = 17`), so it is not an information. |

The article also sums `K` over `n = 1..N, m = 1..M`, while the likelihood
sums over `n, m > w`. The package sums over the effective region, consistent
with the likelihood.

The legacy `fisher_information_2dkarma` implements the gamma blocks with the
**correct** sign, but its shape-shape entry uses the observed `y`
(terms `l2`...`l5`). It equals `-d^2 l / d phi^2` exactly (test
`legacy shape-shape entry was the observed one`). The resulting matrix is
neither the expected nor the observed information. It is, however, the matrix
that produced the standard errors of Table 4 (Sec. 5.2).

What the package provides:

* `information = "expected"` (default): the verified expected information.
* `information = "observed"`: minus the Hessian, by central differences of
  the analytic score; agrees with `numDeriv::hessian()` to 1e-9 relative.
* `vcov` is `chol2inv(chol(K))` only if `K` is finite and numerically
  positive definite. Otherwise `vcov` is `NA` with a warning. A reciprocal
  condition number below 1e-12 is flagged as "ill-conditioned" in the fit,
  in `summary()` and in `karma2d_diagnostics()`. No ridge or other
  regularization is ever applied.

## 4. Residuals and prediction

Quantile residuals `Phi^-1(F(y | mu_hat, phi_hat))` are computed from the
smaller tail on the log scale, so extreme residuals stay finite. The legacy
code replaced infinite residuals by the maximum and `NaN` by
`min(x[is.nan(x)])`, which is itself `NaN`.

Prediction (Eq. 8) with the border strategy of Sec. 3.1: the article calls it
"mirror padding" but describes, and the code implements, **duplication**
without reversal of the first `2w` rows (on top), then of the first `2w`
columns (on the left). `predict(padding = "duplicate")` reproduces it;
`padding = "none"` leaves the first `w` rows and columns `NA`.

## 5. Reproduction of the article

Commands are in `vignette("reproducing-the-paper")` and `analysis/README.md`.

### 5.1 Monte Carlo (Tables 1-3)

The true values in the original scripts differ from the article for two of
the three models:

| Model | Article (Tables 1-2) | Original script |
|---|---|---|
| KARMA(1,0) | 0.6609, (0.2819, 0.1897, 0.1653), 19.6065 | 0.67, (0.31, 0.21, -0.15), 23 |
| KARMA(0,1) | 1.1114, (0.3467, 0.2819, 0.1848), 12.3857 | 0.12, (0.1, 0.02, -0.12), 11 |
| KARMA(1,1) | 0.3, (0.3, -0.15, 0.23), (0.15, -0.1, 0.25), 17 | identical |

The article says the first two were taken from the NDVI and GRDI fits, and the
Table 1 and Table 2 values are indeed those of Tables 4 and 6. Therefore the
scripts in the repository are not the ones that produced Tables 1 and 2. The
pipeline uses the article's values (`--scenario_set=paper`) and keeps the
script values (`--scenario_set=legacy_scripts`). The original scripts did not
compute coverage rates, so which information matrix produced the published
coverage is unknown. The pipeline reports coverage with the expected, the
observed and the legacy matrix.

Results: see Sec. 5.4.

### 5.2 NDVI (Tables 4-5, Figures 4-7)

Data: `NDVI_2019_10.tif` (`stack_all[[32]]` of the original script), rescaled
by `(x + 1)/2`, Phase I `[51:100, 51:100]`, synthetic anomaly
`[5:20, 15:30] = 0.3`. Note that the committed script fits KARMA(1,1)
(`p = 1; q = 1`), while the article reports KARMA(1,0); the pipeline uses the
article's order. Output of `Rscript analysis/ndvi/run_ndvi.R`:

* **Model selection.** The automated rule (all parameters significant at 5 %,
  then fewest parameters, then smallest MSE) selects **2-D KARMA(1,0)**, as in
  the article. The article's graphical diagnostics are not automated.
* **Table 4, KARMA(1,0).** Estimates 0.6609, 0.2819, 0.1897, 0.1653, 19.6064
  (article: 19.6065). All agree within one unit of the last decimal.
  Standard errors for alpha and phi agree (0.0278/0.0277, 0.0191, 0.0210,
  0.0222). The shape standard error is **0.3564** with the verified expected
  information and 0.3145 with the observed information. The article's
  **0.2855** is reproduced exactly (0.2855) by the legacy mixed matrix.
* **Table 4, Gaussian ARMA(0,1).** Estimates (0.8623, 0.4519, 0.4348, 0.3919,
  sigma 0.0439) and standard errors reproduce the article.
* **Table 5.** For KARMA(1,0) in the calibration area, MAE 0.0290 (0.0291),
  MSE 0.00182 (0.0018) and RMSLE 0.0236 (0.0236) reproduce. SMAPE 3.458
  (3.4711) does not, nor do the full-scene metrics, e.g. SMAPE 4.58 (7.69).
  For the Gaussian model none of the metrics reproduce: MAE 0.0299 (0.0612)
  calibration and 0.0327 (0.0620) full scene. Hypotheses tested and
  rejected: (i) the metric is MAPE, as computed by the `mape()` of the
  scripts: 3.55/5.34 for KARMA and 3.60/5.38 for the Gaussian model, which
  does not match; (ii) a Gaussian ARMA(1,1), as committed in the script:
  calibration MAE 0.0240. The committed NDVI KARMA script also compares
  `y_anomaly` on the unit scale with `estimated_ndvi` on the NDVI scale
  (`* 2 - 1`), which would give MAE of order 0.3, so it cannot be the source
  of the table either. **Cause unresolved.** With the computed metrics, the
  prediction-accuracy advantage of KARMA over the Gaussian model is much
  smaller than reported (calibration MAE 0.0290 vs 0.0299; on the full scene
  the Gaussian model has slightly lower MSE and RMSLE).
* **Anomaly detection** (|r| > 3, opening with the 3 x 3 square): KARMA
  detects 286 pixels, 256 of the 256 anomalous pixels (recall 1.00,
  precision 0.90); Gaussian detects 206, 202 of them (recall 0.79, precision
  0.98). The original scripts called `EBImage::erode()`/`dilate()` without a
  `kern` argument. They therefore used EBImage's default
  `makeBrush(5, shape = "diamond")` (checked in the EBImage source), not the
  3 x 3 square stated in the article, so the published figures were probably
  produced with the diamond. The pipeline reports both.

### 5.3 GRDI (Tables 6-7, Figures 8-10)

**Not reproducible from the repository.** Every GRDI script reads
`grdi/data/nasa_corte_india_expandido.csv`, which was listed in `.gitignore`
and never committed, as was the source raster `povmap-grdi-v1.tif`. The crop
extent used to create the CSV was not recorded. The committed
`nasa_corte_india.csv` (500 x 1001) is read by no script and is smaller than
the window used (rows 1200:1300, columns 2200:2300). `analysis/grdi/run_grdi.R`
implements the pipeline (same steps, anomaly `[19:23, 41:50] = 0.1`, 3 x 5
element) and stops with an explanatory message without the input. The
committed GRDI script also fits KARMA(1,0) (`p = 1; q = 0`), while the article
reports KARMA(0,1). About 24 % of the values of the committed CSV are exactly
0.01. This is probably a fill value, but its meaning is undocumented.

### 5.4 Monte Carlo results

Full table: `docs/monte-carlo-r1000.md`. It is based on 1,000 replications
per cell; the article's 10,000-replication run was not executed in this
environment and remains pending (about 4 h on 3 cores). All 9,000 fits
converged. 51 of 54 means agree with Tables 1-3 within two Monte Carlo
standard errors. MSEs agree except for the 10 x 10 images (e.g. KARMA(1,1)
alpha 0.133 vs 0.017), which suggests extreme small-image estimates not
reflected in the article. Coverage rates reproduce the article's pattern,
including the under-coverage of phi(1,1) and theta(1,1).

## 6. Anomaly-detection post-processing

Binary erosion/dilation are implemented in base R with an arbitrary
structuring element (origin at its centre). Border convention: outside pixels
never remove pixels in an erosion and never add pixels in a dilation. EBImage's
border handling was not compared, because EBImage is not installed in the
audit environment.

## 7. Validation summary

`devtools::test()`: 10 test files, about 190 expectations, all passing. They
cover the distribution, lag structure, likelihood, score vs numDeriv,
information vs numerical integration, observed vs numDeriv Hessian,
equivalence with the legacy baseline (`data-raw/legacy-baseline.R` →
`tests/testthat/fixtures/legacy-baseline.rds`), simulation reproducibility and
global-RNG safety, fitting, S3 methods, prediction, residuals, morphology,
detection and metrics. `R CMD check --as-cran` results are reported in the
pull request.

## 8. Open questions for the authors

1. Appendix B: confirm the corrected signs of `K_mu,phi` and `K_phi,phi`.
2. Which information matrix should the article report? The published shape
   standard errors come from the mixed matrix; with the expected information
   SE(phi) increases by about 25 % on the NDVI window.
3. How were Table 5 and Table 7 computed (code version, metric definitions,
   compared images)?
4. Which code and seeds produced Tables 1-3 (the committed scripts use
   different parameters and do not compute coverage)?
5. Provide `nasa_corte_india_expandido.csv` (or the crop extent and GRDI
   version) so the GRDI application can be reproduced.
6. Which structuring element produced Figures 7 and 10?
7. Licence of the package (MIT was set provisionally) and the roles of the
   co-authors in `DESCRIPTION`.

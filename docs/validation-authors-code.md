# Validation of the authors' simulation and NDVI code

This document validates the code in the three archives supplied on
2026-11-12 (`simu.zip`, `simu_new2.7z`, `apli_ndvi.7z`), which contain the
scripts that produced Tables 1-5 of the article, and compares them with the
`karma2d` package. Each statement was obtained by running the authors' code
unmodified next to the package on the same inputs, or by reproducing a
published number to the printed precision. The archives are not committed to
this repository.

## 1. What the archives contain

| Path | Content |
|---|---|
| `simu/simu_paper/` | `fit_2dkarma.R`, `simu_2dkarma.R`, `functions.R`, Monte Carlo scripts and their outputs for the article |
| `simu/simu_new/`, `simu_new2/simu_new2/` | Same code with one variant of `fit_2dkarma.R`; Monte Carlo runs for other parameter sets |
| `apli_ndvi/apli ndvi/` (duplicated inside `simu_new2.7z`) | NDVI application: `apli_ndvi.R`, `apli_ndvi_gaussian.R`, `fitting.R`, `metrics.R`, Gaussian fit `2D-GARMA-function/fit_2dgarma.r`, the 33 NDVI rasters |

`fit_2dkarma.R` is byte-identical in `simu_paper`, `simu_new2` and
`apli ndvi`; `functions.R` and `simu_2dkarma.R` are identical everywhere.

## 2. Provenance of the published tables

| Article | Produced by | Check |
|---|---|---|
| Table 1, KARMA(1,0) | `simu_paper/simu_MC_10.R` → `resultado_2dkarma10_paper_2.txt` | all printed values identical |
| Table 2, KARMA(0,1) | `simu_paper/simu_MC_01.R` → `resultado_2dkarma01_paper.txt` | identical |
| Table 3, KARMA(1,1) | `simu_new/simu_MC_11.R` → `teste_2dkarma11.txt` (10, 25) and `simu_MC_2dkarma11_50.txt` (50) | identical |
| Table 4 | `apli_ndvi.R`, `apli_ndvi_gaussian.R` | estimates and standard errors reproduced (Sec. 4) |
| Table 5 | same | reproduced to 4 decimals, including the defects of Sec. 5 |

All Monte Carlo runs use 10,000 replications and `set.seed(1)` once at the
start. This settles question 4 of the technical audit: the scripts first
committed to this repository were not the ones behind Tables 1-3.

## 3. Core code: what is correct

Tested on 25 x 25 images simulated from the three scenarios of the article
(script: authors' `simu.2dkarma` and `fit.2dkarma` vs the package).

| Component | Result |
|---|---|
| `dkum`, `pkum`, `qkum` | Agree with `dkuma`, `pkuma`, `qkuma` to 3e-11 on a grid away from tiny `mu^phi` (same precision caveat as the legacy functions for `mu^phi < 1e-8`) |
| `simu.2dkarma` | Identical to `karma2d_simulate(burn_in = max(p, q))` with the same seed (max difference 3e-12). It discards only the `m` initialization rows and columns; the scripts first committed here discarded `2m`. The Monte Carlo pipeline now uses the authors' value. |
| Log-likelihood | Equal to `karma2d_loglik()` at the same parameters |
| Analytic score | Matches `numDeriv::grad` (3e-8 for AR and MA, 7e-6 for ARMA, which is numDeriv's own error) |
| Fitted medians | Equal to the package's at the same coefficients (5e-6 or better) |
| Lag ordering | `phi(0,1), phi(1,0), phi(1,1)`: row lag first, as in the article and the package |

## 4. Core code: findings

### 4.1 The information matrix is the mixed one (affects every standard error)

`z$K` is assembled as minus a block matrix whose gamma blocks use the
expected weights and whose shape-shape entry is `sum(L)` with `L` built from
the observed `y1` (`L2`...`L5`). The result equals, to 1e-13 relative, the
"mixed" matrix described in the technical audit (Sec. 3.5): expected
information for the gamma-gamma and gamma-shape blocks, **observed**
information for the shape-shape entry. The sign of the cross block in the
code is the correct one, opposite to Appendix B of the article.

Consequences:

* On data simulated from the model the mixed, expected and observed standard
  errors are close (e.g. shape: 0.768 / 0.766 / 0.767 on a 25 x 25 AR image),
  so the coverage rates of Tables 1-3 are barely affected.
* On the NDVI window, where the model is only an approximation, they differ:
  the published SE of the shape, 0.2855, becomes 0.356 with the expected
  information and 0.315 with the observed information.
* The matrix is neither estimator of the asymptotic covariance and is not
  guaranteed to be positive definite.

Recommendation: use the expected information (corrected formulas) or the
observed information, and state which in the article.

### 4.2 The optimizer stops before the maximum

`optim()` is called with its default `reltol`. At the returned estimates the
score is not negligible (largest component 0.125 on a 25 x 25 AR image), and
the estimates differ from the maximum by up to 0.004, about 0.13 standard
errors. This adds noise but no visible bias to the Monte Carlo summaries.
The package uses `reltol = 1e-10`.

### 4.3 Monte Carlo loop

* Replications whose fit throws an error or does not converge are discarded
  and redrawn, without being counted. Measured with the authors' loop on
  10 x 10 images (300 kept replications): KARMA(1,1) 3 discarded (1.0 %),
  KARMA(1,0) none. The effect on the tables is small, but the number of
  discarded replications should be reported.
* The ARMA branch runs `optim()` a second time without the gradient and with
  `hessian = TRUE` (`opt2`), whose result is never used, and every fit calls
  `numDeriv::grad`. Both only cost time.
* Coverage uses 1.96 and the mixed-matrix standard errors.

### 4.4 The `simu_new` variant of `fit_2dkarma.R` breaks the pure AR model

In `simu_new/fit_2dkarma.R` the block that builds the lag matrix `P` is
commented out, but the AR score still uses `t(P)`. Every AR fit fails with
`object 'P' not found`; because the Monte Carlo loop discards failed fits
and tries again, an AR simulation with this file never ends. The ARMA and MA
branches are unaffected (Table 3 was produced with this file; Table 1 with
the `simu_paper` version, which works).

### 4.5 Multiple local maxima on small images (package and authors' code)

On 400 images of size 10 x 10 simulated from the Table 3 model (authors'
simulator, seed 1), no single-start strategy always finds the best maximum.
Counting fits whose log-likelihood is more than 0.01 below the best found by
any of the strategies:

| Strategy | Below best | MSE of alpha |
|---|---|---|
| Authors' settings (shape on its own scale, default `reltol`) | 29 of 400 | 0.0176 |
| Package before the fix (log shape, unscaled objective) | 23 of 400 | 0.657 (one fit ended at alpha = 16.3) |
| Package now (log shape, objective scaled by the pixel count) | 24 of 400 | 0.0208 |
| Best of two shape transforms | 11 of 400 | 0.0182 |

The catastrophic fit of the earlier package version came from the first
BFGS step, which is proportional to the unscaled gradient. The objective is
now divided by the number of modelled pixels (`fnscale`), which removes it
and halves the number of function evaluations. This also explains the larger
10 x 10 MSE reported in the first version of `docs/monte-carlo-r1000.md`:
it was an artefact of that fit, not a property of the article's numbers.

On 50 x 50 images the problem does not appear: starting at the true
parameters or at the default starting values gave the same maximum in 80 of
80 replications.

### 4.6 The bias of phi(1,1) and theta(1,1) in Table 3 does not vanish relative to the standard error

Table 3 shows bias 0.095, 0.042, 0.022 for `phi(1,1)` at sizes 10, 25, 50
and coverage 0.746, 0.855, 0.858. The check above shows this is the actual
maximum of the conditional likelihood, not an optimizer artefact. The bias
decreases roughly like 1/N, the same rate as the standard error of an
N x N image, so the ratio bias/SE and the coverage stay about constant. This
is the behaviour expected from edge effects of conditional likelihoods on
two-dimensional lattices (the errors of the first rows and columns are set
to zero); see Guyon (1982, *Biometrika* 69, 95-105) for the Gaussian case.
The statement of the article that the coverage rates "gradually converge to
the nominal value" is therefore not supported for these two parameters, and
the usual Wald intervals for them are too optimistic at every size studied.

## 5. NDVI application

### 5.1 Procedure actually used (differs from the scripts first committed)

* Synthetic anomaly: pixels `[5:20, 15:30]` set to **0.01** (the committed
  script used 0.3).
* Calibration metrics: `y_train` against the fitted medians of the Phase I
  fit, effective region only.
* Full scene: `fitting()` is called with the **unpadded** image; the padded
  matrix `y_new` is built but never used. The first `2m` rows and columns of
  the prediction are dropped. The border strategy described in Sec. 3.1 of
  the article is therefore not what produced the results.
* Post-processing: `mmand::opening` followed by `mmand::closing` with a
  3 x 3 box. The article says "erosion followed by dilation".
* Gaussian residuals: `(y - prediction) / sigma_hat`.

With this procedure the package reproduces the KARMA rows of Table 5 to the
printed precision: calibration 0.0291, 0.0018, 3.4711, 0.0236; full scene
0.0305, 0.0037, 7.6934, 0.0407.

### 5.2 Two defects in the Gaussian benchmark (they change the conclusions)

**(a) Calibration area: prediction compared with its transpose.**
`fit.2dgarma()` returns `fitted = mu` as a vector in row-major order
(`as.vector(t(eta[...]))`), and `apli_ndvi_gaussian.R` computes
`mae(real_cal, fit$fitted)` with `real_cal` a matrix, which R recycles in
column-major order. Each pixel `[i, j]` is compared with the prediction for
`[j, i]`. (`fit.2dkarma()` returns a matrix, so the KARMA rows are
unaffected.)

**(b) Full scene: Gaussian coefficients pushed through the logit link.**
`fitting()` hard-codes `linkfun = logit` and `linkinv = plogis`, and is
called with the Gaussian coefficients, which were estimated on the identity
scale. The "Gaussian prediction" is `plogis(alpha + theta' e)` with
`e = logit(y) - eta`: with `alpha = 0.8623` it is centred near 0.70 for
pixels near 0.85.

Both are confirmed by reproducing the published numbers exactly:

| Gaussian ARMA(0,1) | MAE | MSE | SMAPE | RMSLE |
|---|---|---|---|---|
| Calibration, article | 0.0612 | 0.0061 | 7.2126 | 0.0425 |
| Calibration, prediction transposed (defect a) | 0.0612 | 0.0061 | 7.2125 | 0.0425 |
| **Calibration, corrected** | **0.0303** | **0.0019** | **3.5984** | **0.0242** |
| Full scene, article | 0.0620 | 0.0057 | 11.4619 | 0.0453 |
| Full scene, logit link (defect b) | 0.0620 | 0.0057 | 11.4619 | 0.0453 |
| **Full scene, corrected** | **0.0371** | **0.0060** | **8.2569** | **0.0587** |
| 2-D KARMA(1,0), calibration | 0.0291 | 0.0018 | 3.4711 | 0.0236 |
| 2-D KARMA(1,0), full scene | 0.0305 | 0.0037 | 7.6934 | 0.0407 |

Effect on the article's claims:

* Calibration area: the article reports reductions of about 54 % in MAE and
  68 % in MSE. Corrected, they are about **4 %** and **6 %**.
* Full scene: reported 54 % (MAE) and 36 % (MSE); corrected, about **18 %**
  and **38 %**. KARMA remains better on the full scene.
* Anomaly detection (threshold 3, opening then closing, 98 x 98 region,
  256 anomalous pixels):

| Model | Flagged | After post-processing | Recall | Precision |
|---|---|---|---|---|
| 2-D KARMA(1,0) | 372 | 289 | 1.00 | 0.89 |
| Gaussian, as published (defect b) | 337 | 0 | 0.00 | – |
| **Gaussian, corrected** | 418 | 297 | 0.98 | 0.85 |

  The statement that the Gaussian model "fails to clearly delineate the
  anomaly region" (Figure 7(d)) is a consequence of defect (b). With the
  correct Gaussian prediction the anomaly is detected almost as well as with
  2-D KARMA.

The GRDI scripts were not in the archives. If they use the same
`fitting()` and `fit.2dgarma()`, Tables 6-7 and Figure 10 have the same two
defects; this should be checked by the authors.

## 6. Summary for the authors

1. Correct and consistent with the package: distribution functions,
   simulator, likelihood, score, lag ordering, the KARMA rows of Tables 4-5
   and the point estimates of Tables 1-3.
2. To fix before publication: the Gaussian benchmark (Sec. 5.2), which
   changes the size of the reported gains and reverses the detection
   comparison; the standard errors (Sec. 4.1); Appendix B (technical audit,
   Sec. 3.5); the description of the border handling, anomaly value and
   morphological operations (Sec. 5.1); the coverage statement for the ARMA
   model (Sec. 4.6).
3. To report: the number of discarded Monte Carlo replications, and the
   optimizer tolerance.

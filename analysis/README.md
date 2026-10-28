# Reproduction of the article

This folder reproduces the numerical results of Pereira, Palm and Bayer
(2026). It is **not part of the installed package**: it is excluded from the
package build, may use extra packages (`terra`), and writes results to
`analysis/output/` (git-ignored).

Run everything from the repository root after installing the package
(`R CMD INSTALL .` or `remotes::install_github("VitorBSP/2d_karma_package")`).

| Step | Command | Output |
|---|---|---|
| Monte Carlo, quick check (20 replications, sizes 10 and 25) | `Rscript analysis/simulations/run_monte_carlo.R --mode=quick` | `analysis/output/simulations/paper/*.rds` |
| Monte Carlo, article configuration (10,000 replications, sizes 10, 25, 50) | `Rscript analysis/simulations/run_monte_carlo.R --mode=full --cores=8` | same; resumable per cell |
| Monte Carlo summary and comparison with Tables 1-3 | `Rscript analysis/simulations/summarize_monte_carlo.R` | `summary.csv`, `comparison-with-paper.csv` |
| NDVI application (Tables 4-5, Figures 4-7) | `Rscript analysis/ndvi/run_ndvi.R` | `analysis/output/ndvi/` |
| GRDI application (Tables 6-7, Figures 8-10) | `Rscript analysis/grdi/run_grdi.R --input=PATH/nasa_corte_india_expandido.csv` | `analysis/output/grdi/`; the input file is not in the repository |
| Validate that the pipelines still give the results already obtained, and compare with the article | `Rscript analysis/validate_results.R` | `analysis/output/validation/validation-report.csv`; exit status 1 on any unexpected difference |

### Validating the results

`analysis/validate_results.R` reruns the NDVI pipeline and the quick Monte
Carlo (about 1.5 min) and performs two checks:

* **regression**: every number is compared with the snapshot of the results
  already obtained, stored in `analysis/reference/` (relative tolerance
  1e-5). Monte Carlo summaries are compared with the snapshot of the same
  configuration: `monte_carlo_R20_sizes10-25.csv` (quick mode) or
  `monte_carlo_R1000_sizes10-25-50.csv` (the 1,000-replication run of
  `docs/monte-carlo-r1000.md`).
* **paper**: comparison with the article (NDVI within one unit of the last
  printed decimal; Monte Carlo means within 3 Monte Carlo standard errors,
  coverages within `3 * sqrt(p (1 - p) / R) + 0.01`). Divergences already
  documented are listed with their reason in
  `analysis/reference/known-divergences.csv` and reported as `known`.

Validate outputs already produced without rerunning, e.g. the
1,000-replication run or a full run:

```
Rscript analysis/validate_results.R --rerun=FALSE \
  --mc_dir=analysis/output/simulations_r1000 --ndvi_dir=analysis/output/ndvi
```

After an intended change of results (for example a correction), regenerate
the snapshot explicitly with `--update_snapshot=TRUE` and commit it, so that
the change is visible in review. The `validate-results` GitHub workflow runs
the validation on every pull request.

Other options: `--replications=N`, `--scenario_set=legacy_scripts` (the
parameter values hard-coded in the original scripts), `--scenarios=karma11`,
`--seed=N`, `--output=DIR`; for the applications `--karma_order=p,q`,
`--gaussian_order=p,q`, `--select=FALSE` (skip model selection).

Each output directory contains `session-info.txt` (R, package and GDAL
versions). Monte Carlo replications use independent L'Ecuyer-CMRG streams
derived from `--seed`, so results do not depend on the number of cores.
On Windows use `--cores=1` (forking is unavailable).

For the full experiments, recording the package versions with `renv` is
recommended (`renv::init()` then `renv::snapshot()`); it is not required for
ordinary users of the package.

Layout:

* `R/`: shared helpers, the Gaussian 2-D ARMA benchmark, the application
  pipeline.
* `simulations/`, `ndvi/`, `grdi/`: entry points, configuration, published
  values used for comparison.
* `data/`: raw data and `SHA256SUMS` (see `docs/data-provenance.md`).
* `legacy/`: the original scripts, unchanged, kept for traceability. They do
  not run as they are (missing files, archived packages, absolute layout).

Status of the reproduction and all known divergences:
`docs/technical-audit.md`, Section 5.

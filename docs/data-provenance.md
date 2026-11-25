# Data provenance

SHA-256 checksums of every raw file are listed in `analysis/data/SHA256SUMS`
(check with `sha256sum -c analysis/data/SHA256SUMS` from the repository root).
Fields marked *unknown* were not recorded in the original project. They
should be completed by the authors before publication.

## NDVI rasters (`analysis/data/ndvi/`)

| Field | Value |
|---|---|
| Content | 33 GeoTIFF files `NDVI_YYYY_KK.tif` (2018: composites 01-22; 2019: 01-11) with GDAL `.aux.xml` statistics side-cars |
| Product | MODIS/Terra Vegetation Indices 16-Day L3 Global 250 m (MOD13Q1), NASA LP DAAC, as cited in the article; product collection/version *unknown* |
| URL | https://modis.gsfc.nasa.gov/data/dataprod/mod13.php (product page); the original download source *unknown* |
| Acquisition date | *unknown*; `KK` is presumably the 16-day composite index within the year |
| Processing | Reprojected and cropped by the authors (procedure *unknown*) |
| Dimensions / resolution | 100 x 100 pixels; 0.0024129 x 0.0024130 degrees (about 250 m) |
| Extent / CRS | lon -37.17953 to -36.93824, lat -7.14512 to -6.90382; WGS 84 (EPSG:4326); semi-arid region of Paraiba, Brazil |
| Units / scale | NDVI, dimensionless, natural range (-1, 1); observed range -0.20 to 1.00 across files |
| Missing values | None in any of the 33 files |
| Used in the article | `NDVI_2019_10.tif` (the 32nd file in alphabetical order, `stack_all[[32]]` in the original script), the "May 2019" scene; SHA-256 `aac7580d93b6a36c43e9aeaa394ea194c7600bc94aea8eb72a5be3c7c8862173` |
| Licence | NASA EOSDIS data are distributed without restrictions on reuse (NASA Earth Science Data and Information Policy); citation of the product is requested |
| Transformations in the analysis | `y = (NDVI + 1) / 2`; Phase I window `[51:100, 51:100]`; synthetic anomaly `[5:20, 15:30] <- 0.3` on the unit scale (NDVI -0.4) |

The package dataset `ndvi_paraiba` (`data/ndvi_paraiba.rda`) is
`NDVI_2019_10.tif` read by `data-raw/ndvi_paraiba.R` with `terra`. The script
verifies the MD5 checksum of the source; the dataset keeps the NDVI scale and
stores the extent and CRS as attributes. `terra` and the `raster` reader of
the original scripts give identical matrices.

## GRDI (`analysis/data/grdi/`)

| Field | Value |
|---|---|
| Product | Global Gridded Relative Deprivation Index (GRDI), v1, NASA SEDAC / CIESIN (`povmap-grdi-v1`) |
| URL | https://sedac.ciesin.columbia.edu/data/set/povmap-grdi-v1 (now distributed through NASA Earthdata) |
| Licence | Creative Commons Attribution 4.0 (CC BY 4.0) according to SEDAC; to be confirmed for the exact version used |
| File in the repository | `nasa_corte_india.csv`: 500 rows x 1001 columns (header `V1..V1001`), values 0.01-89.29, no missing values, 24 % of the cells exactly 0.01 (meaning *unknown*, possibly a fill value); SHA-256 `abfcf86e383b2b037714a4a306f09d1823fee47241bd3208abb6c590b8063fe6`. It is **not read by any original script**. |
| File used by the article | `nasa_corte_india_expandido.csv`, **not in the repository** (122 MB). Published by the authors on 2026-11-25 in the private Hugging Face dataset `VitorBSP/nasa_corte_india_expandido` (https://huggingface.co/datasets/VitorBSP/nasa_corte_india_expandido); SHA-256 `96e53f58a27cbbff874fed6caa2c70f4bde13c882f99fb5a37691234659b6458`. 4,501 rows x 5,001 columns (header `V1..V5001`), 82.6 % missing (`NA`), values 3.7-99.8. CRS and resolution of the source raster still not recorded. |
| Relation between the two files | The committed extraction is exactly `expandido[1501:2001, 2501:3501]` with `NA` written as 0.01 (all cells agree). |
| Transformations in the original scripts | rows 1200:1300, columns 2200:2300 of the data matrix, first 100 x 100 block, divided by 100, missing values replaced by the mean; Phase I `[51:100, 51:100]`; synthetic anomaly `[19:23, 41:50] <- 0.1` |

**Reproduction with the expanded file (2026-11-25).** Download (the dataset
is private; a read token is needed):

```
curl -L -H "Authorization: Bearer $HF_TOKEN" -o nasa_corte_india_expandido.csv \
  https://huggingface.co/datasets/VitorBSP/nasa_corte_india_expandido/resolve/main/nasa_corte_india_expandido.csv
sha256sum nasa_corte_india_expandido.csv   # 96e53f58...b6458
Rscript analysis/grdi/run_grdi.R --input=nasa_corte_india_expandido.csv
```

The window of the original scripts (rows 1200:1299, columns 2200:2299, 42
missing cells, 7 of them in the Phase I block) gives estimates close to, but
not equal to, Table 6:

| Parameter | KARMA(0,1) computed | Table 6 | Gaussian computed | Table 6 |
|---|---|---|---|---|
| alpha | 1.1341 | 1.1114 | 0.7466 | 0.7420 |
| theta(0,1) | 0.3452 | 0.3467 | 0.4447 | 0.4329 |
| theta(1,0) | 0.2947 | 0.2819 | 0.3985 | 0.3709 |
| theta(1,1) | 0.2136 | 0.1848 | 0.2858 | 0.2606 |
| shape / sigma | 11.833 | 12.3857 | 0.0754 | 0.0733 |

The largest deviation is 2.5 published standard errors (KARMA theta(1,1)).
Checked and ruled out as causes:

* **Window.** `analysis/grdi/search_article_window.R --input=... --missing_value=`
  tries every 100 x 100 window with at least half of its cells valid (4.1
  million). The original window is the closest to Table 6 in the whole file;
  the next best is 4.7 standard errors away.
* **Missing-value treatment.** Window mean (the original scripts), Phase I
  mean or whole-image mean change the estimates by less than 0.5 standard
  errors; the 0.01 code of the extraction moves them by 20.
* **Off-by-one cropping.** A 101 x 101 crop and a shifted Phase I give the
  same distance.

The most likely explanation is that the published file is not exactly the
version analysed for the article (for example, a later export of the
raster). Only the authors can confirm this.

Table 7 with the article's procedure (calibration area):

| Model | MAE | MSE | SMAPE | RMSLE |
|---|---|---|---|---|
| 2-D KARMA(0,1) | 0.0594 (0.0557) | 0.0060 (0.0056) | 8.45 (7.98) | 0.0466 (0.0456) |
| Gaussian, as published | 0.0806 (0.0769) | 0.0116 (0.0106) | 11.37 (10.92) | 0.0645 (0.0620) |
| Gaussian, corrected | 0.0598 | 0.0057 | 8.47 | 0.0451 |

(Values in parentheses: Table 7 of the article.) The Gaussian "as published"
variant, with the two defects found in the NDVI scripts, lands close to the
article's Gaussian rows, and the corrected benchmark lands next to KARMA.
This supports, but does not prove, that Tables 6-7 have the same defects.
With the corrected benchmark the two models predict equally well on GRDI;
in detection (3 x 5 element), KARMA has recall 1.00 and precision 0.93, the
corrected Gaussian recall 0.72 and precision 1.00, and the Gaussian "as
published" detects nothing.

**Earlier search in the extraction.** Before the expanded file was
available, the same search on the committed extraction found no window
closer than 4 standard errors (`Rscript analysis/grdi/search_article_window.R`),
consistent with the extraction not containing rows 1200:1299. The pipeline
also runs on the extraction as a functional test:

```
Rscript analysis/grdi/run_grdi.R --input=analysis/data/grdi/nasa_corte_india.csv \
  --window=246,674 --missing_value=0.01
```

## Distribution policy

* `data/`: only `ndvi_paraiba` (19 KB), small and unrestricted.
* `analysis/data/`: raw files for the reproduction pipelines. They stay in
  the repository but are excluded from the package build.
* Nothing is downloaded during installation, tests, examples or vignettes.

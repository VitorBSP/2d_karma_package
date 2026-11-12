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
| File used by the article | `nasa_corte_india_expandido.csv`, **not in the repository** (git-ignored), nor is the source raster `povmap-grdi-v1.tif`. Its crop extent, CRS and resolution were not recorded. |
| Transformations in the original scripts | rows 1200:1300, columns 2200:2300 of the data matrix, first 100 x 100 block, divided by 100, missing values replaced by the mean; Phase I `[51:100, 51:100]`; synthetic anomaly `[19:23, 41:50] <- 0.1` |

To make the GRDI application reproducible, the authors should publish the
CSV used (licence permitting, with attribution) or a script that derives it
from the official raster, stating the GRDI version, download date and crop
extent.

## Distribution policy

* `data/`: only `ndvi_paraiba` (19 KB), small and unrestricted.
* `analysis/data/`: raw files for the reproduction pipelines. They stay in
  the repository but are excluded from the package build.
* Nothing is downloaded during installation, tests, examples or vignettes.

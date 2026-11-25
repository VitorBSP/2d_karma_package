# Build the `ndvi_paraiba` dataset shipped with the package.
#
# Source file: analysis/data/ndvi/NDVI_2019_10.tif, the raster selected as
# `stack_all[[32]]` by the original NDVI application script
# (analysis/legacy/ndvi/aplicacao_ndvi_karma.R). See
# docs/data-provenance.md for provenance and checksums.
#
# terra replaces the archived rgdal/raster packages of the original script;
# both readers give identical matrices (checked with raster 3.x and
# terra 1.7, maximum absolute difference 0).
#
# Run from the repository root:  Rscript data-raw/ndvi_paraiba.R

source_file <- file.path("analysis", "data", "ndvi", "NDVI_2019_10.tif")
stopifnot(unname(tools::md5sum(source_file)) ==
            "82e6ae2b4ebf8c2a9732ffb1ffd5c317")

r <- terra::rast(source_file)
ndvi_paraiba <- terra::as.matrix(r, wide = TRUE) # rows run north to south
dimnames(ndvi_paraiba) <- NULL
attr(ndvi_paraiba, "extent") <- as.vector(terra::ext(r))
attr(ndvi_paraiba, "crs") <- "EPSG:4326"

stopifnot(identical(dim(ndvi_paraiba), c(100L, 100L)),
          !anyNA(ndvi_paraiba),
          all(ndvi_paraiba > -1 & ndvi_paraiba < 1))

save(ndvi_paraiba, file = file.path("data", "ndvi_paraiba.rda"),
     compress = "xz", version = 2)

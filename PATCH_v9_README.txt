Patch v9 - Leaflet raster fix + faster installs

What changed:
- reports/haz_leaflet.Rmd
  * Added robust terra->raster conversion (to_raster_layer()) so leaflet::addRasterImage works on GitHub Actions.
  * Explicitly uses leaflet::addRasterImage(..., project = FALSE) to avoid slow/fragile reprojection at knit time.
  * Loads raster + htmltools.

- .github/workflows/wmo_tracks_daily.yml
  * Installs raster, viridisLite, htmltools (in addition to existing deps).
  * Adds output_options(lib_dir='site_libs', self_contained=FALSE) + knit_root_dir=getwd() so dependencies land in docs/site_libs.
  * Creates docs/.nojekyll so GitHub Pages serves site_libs correctly.

How to apply:
- Unzip over the repo root (preserving paths) and commit.

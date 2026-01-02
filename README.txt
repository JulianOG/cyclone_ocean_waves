This zip updates your repo with:

1) scripts/wmo_tracks_daily.R
   - Writes per-storm NetCDF hazard files: outputs/latest/haz_###.nc (default ON)
   - Optional multi-layer GeoTIFFs (Hs0 + Pr) if you set SAVE_TIF=1
   - Always writes outputs/latest/index.csv and outputs/latest/sessionInfo.txt (even if no storms)

2) reports/haz_leaflet.Rmd
   - Builds interactive Leaflet maps from NetCDF hazards (preferred), or GeoTIFFs (fallback)
   - Exits cleanly if index.csv is missing (so the render step won't fail)

3) .github/workflows/wmo_tracks_daily.yml
   - Uses use-public-rspm:true + dependency caching to speed installs
   - Ensures /docs is created and always produces docs/index.html
   - Commits outputs/latest and docs for GitHub Pages (/docs)

What’s in this zip

1) scripts/wmo_tracks_daily.R
   - Updated to save the FULL hazard rasters (multi-layer GeoTIFFs):
       outputs/latest/haz_Hs0_###.tif, haz_Pr_###.tif, haz_Uw_###.tif, etc.
   - Adds hs0_peak_time and pr_min_time to outputs/latest/index.csv
   - You can disable saving full rasters by setting SAVE_HAZ=0.

2) reports/haz_leaflet.Rmd
   - Reads outputs/latest/index.csv and the hazard GeoTIFFs.
   - Creates interactive Leaflet maps (via leafem) showing:
       Hs0 at peak time + Hs0 max over time
       Pressure at min time + Pressure min over time
     plus observed/forecast tracks.

3) .github/workflows/wmo_tracks_daily.yml
   - Updated to install leaflet/leafem/rmarkdown and render docs/index.html daily.

How to use

- Unzip into the repo root, commit, push.
- In GitHub: Settings → Pages → Build from /docs (optional) to publish docs/index.html.


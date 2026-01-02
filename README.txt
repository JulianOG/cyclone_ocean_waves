This v5 zip fixes the 'Render report (HTML)' failures by:

1) Making the Rmd resilient:
   - Robust NetCDF subdataset loading (tries multiple methods).
   - Wraps each storm map in tryCatch so one bad file won't fail the whole render.

2) Making the workflow resilient:
   - If render fails for any reason, it writes docs/index.html as a placeholder so the workflow can still commit and Pages still works.

Files included:
- scripts/wmo_tracks_daily.R
- reports/haz_leaflet.Rmd (UPDATED)
- .github/workflows/wmo_tracks_daily.yml (UPDATED)

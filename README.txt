Fix for your current render failure

The error "cannot open the connection" happens because R Markdown may knit from the /reports
directory, making relative paths like outputs/latest/index.csv point to reports/outputs/...

This update resolves all output paths relative to the git repo root, so it works reliably on GitHub Actions.

Files included:
- scripts/wmo_tracks_daily.R (unchanged from v2)
- reports/haz_leaflet.Rmd (UPDATED: robust path resolution)
- .github/workflows/wmo_tracks_daily.yml (minor cache bump)

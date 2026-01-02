Patch v7: Fix blank Leaflet maps when knitted (especially inside {.tabset}) and on GitHub Pages.

What it changes:
- Rmd: adds JS to trigger window resize when tabs change (Leaflet redraws)
- Rmd: sets leaflet(height=650)
- Rmd: self_contained:false + lib_dir:site_libs (ensures widget assets are saved in docs/site_libs)
- Workflow: touches docs/.nojekyll (safer for Pages) and removes leafem dependency

Apply by overwriting:
- reports/haz_leaflet.Rmd
- .github/workflows/wmo_tracks_daily.yml

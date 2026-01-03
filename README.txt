Patch v8: Fixes the rmarkdown 'relativeTo(...) descendant of docs' error.

Cause:
- Setting self_contained:false + lib_dir in the Rmd header can make dependencies write to reports/site_libs
  while output_dir is docs/, and rmarkdown then errors.

Fix:
- Remove self_contained/lib_dir from the header (default is self-contained HTML).
- Keep the tabset Leaflet resize JS + explicit leaflet(height=650).

Apply by overwriting:
- reports/haz_leaflet.Rmd

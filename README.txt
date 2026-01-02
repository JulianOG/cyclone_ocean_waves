Fix for: values must be type 'character' but FUN result is type 'logical'

Cause:
- read.csv() can infer columns that are entirely NA as logical.
- vapply(..., FUN, character(1)) then fails if FUN returns logical NA.

Fix:
- In reports/haz_leaflet.Rmd, path columns are coerced with as.character() before vapply().

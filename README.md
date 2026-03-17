# cyclone_ocean_waves

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/JulianOG/cyclone_ocean_waves)

Automated tropical cyclone coastal hazard forecasting and reporting using WMO track data, the `TCHazaRds` R package, and GitHub Pages.

## Overview

`cyclone_ocean_waves` is a fully automated R-based workflow that:

- fetches active tropical cyclone track data from the WMO Severe Weather API
- computes coastal hazard fields from tropical cyclone forcing
- generates wave height and inundation proxy products
- writes quick-look visualisations and GIS-ready outputs
- publishes a daily interactive HTML report via GitHub Pages

The repository is designed to support rapid, repeatable generation of tropical cyclone ocean hazard products with minimal manual intervention.

## What the workflow produces

For each active tropical cyclone, the workflow can generate:

- **Track vectors** as GeoPackages
- **Wave hazard outputs (`Hs0`)**
  - maximum significant wave height
  - time to peak wave height
- **Inundation / surge proxy outputs (`IB`)**
  - maximum inundation proxy
  - time to peak inundation proxy
- **Quick-look JPEGs** for visual inspection
- **Multi-layer GeoTIFFs** for time-varying hazard fields
- **Optional NetCDF archives** for full gridded outputs
- **An `index.csv` catalog** listing all products and metadata
- **An interactive HTML report** with Leaflet maps

## Hazard variables

### `Hs0`
Significant wave height in metres, derived from the tropical cyclone wind field.

### `IB`
A simple inundation / barometric surge proxy, computed as:

`IB = (1010 - Pr) / 100`

where `Pr` is atmospheric pressure in hPa.

> Note: `IB` is a simplified pressure-deficit proxy and does not explicitly resolve full storm surge dynamics, wind setup, bathymetry, or local coastal geometry.

## Automated workflow

The repository uses **GitHub Actions** to run the workflow daily.

High-level steps:

1. Trigger a scheduled workflow
2. Install required system libraries and R packages
3. Fetch active WMO tropical cyclone tracks
4. Run hazard calculations using `TCHazaRds`
5. Generate raster, vector, and image outputs
6. Update the master output catalog
7. Render the interactive report from R Markdown
8. Commit outputs and publish the updated site

## Repository structure

```text
.github/workflows/       GitHub Actions automation
scripts/                 Main processing scripts
reports/                 R Markdown report templates
outputs/latest/          Latest generated hazard products
docs/                    Published GitHub Pages site

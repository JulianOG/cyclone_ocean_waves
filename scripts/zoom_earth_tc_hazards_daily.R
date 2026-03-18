suppressPackageStartupMessages({
  library(rjson)
  library(terra)
  library(TCHazaRds)
})

Sys.setenv(TZ = "UTC")

# Output folder inside the repo (works on GitHub runners)
outdir <- Sys.getenv("OUTDIR", unset = file.path("outputs", "latest"))
unlink(outdir, recursive = TRUE, force = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Controls (set via GitHub Actions env)
save_tif <- Sys.getenv("SAVE_TIF", unset = "1") == "1"      # default on (not too big)
save_nc  <- Sys.getenv("SAVE_NC", unset = "0") == "1"       # default off (NetCDF can exceed GitHub limit)
save_sum <- Sys.getenv("SAVE_SUMMARY", unset = "1") == "1"  # retained for compatibility

# Zoom Earth controls
zoom_to <- suppressWarnings(as.integer(Sys.getenv("ZOOM_TO", unset = 12)))
if (!is.finite(zoom_to)) zoom_to <- 12L
zoom_wind_units <- Sys.getenv("ZOOM_WIND_UNITS", unset = "kmh") # one of: kmh, kts, ms
default_rmw_km <- suppressWarnings(as.numeric(Sys.getenv("DEFAULT_RMW_KM", unset = NA)))

# --- helpers ----

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

CpFromVmax <- function(vMax, Ep = 1010) {
  dP <- (vMax / 0.6252)^2 / 100
  Ep - dP
}

safe_jpeg <- function(filename, plot_fun, width = 1600, height = 1200, res = 150) {
  ok <- TRUE
  tryCatch({
    grDevices::jpeg(filename, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_fun()
  }, error = function(e) {
    ok <<- FALSE
    message("JPEG failed for ", filename, ": ", conditionMessage(e))
  })
  ok
}

write_haz_raster <- function(x, filename) {
  terra::writeRaster(
    x, filename, overwrite = TRUE,
    wopt = list(
      gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6", "TILED=YES"),
      datatype = "FLT4S"
    )
  )
}

# Write compact single-layer GeoTIFFs (compressed + tiled) suitable for committing to GitHub
write_summary_tif <- function(r, filename, datatype = "FLT4S") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(
    r, filename, overwrite = TRUE, datatype = datatype,
    gdal = c("COMPRESS=DEFLATE", "ZLEVEL=9", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
}

pick_one <- function(x, nms, default = NA) {
  for (nm in nms) {
    val <- x[[nm]]
    if (!is.null(val) && length(val) > 0) {
      if (is.list(val)) val <- unlist(val, use.names = FALSE)
      val <- val[1]
      if (is.character(val)) {
        if (!is.na(val) && nzchar(val) && val != "NULL") return(val)
      } else {
        if (!is.na(val)) return(val)
      }
    }
  }
  default
}

zoom_date_string <- function(x = Sys.time(), step_hours = 6L) {
  x <- as.POSIXct(x, tz = "UTC")
  hr <- as.integer(format(x, "%H", tz = "UTC"))
  hr0 <- (hr %/% step_hours) * step_hours
  sprintf("%sT%02d:00Z", format(x, "%Y-%m-%d", tz = "UTC"), hr0)
}

zoom_as_iso <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "NULL") return(NA_character_)
  x <- as.character(x)[1]
  x <- gsub("T", " ", x, fixed = TRUE)
  x <- sub("Z$", "", x)
  if (nchar(x) == 16) x <- paste0(x, ":00")
  if (nchar(x) == 10) x <- paste(x, "00:00:00")
  x
}

zoom_list_ids <- function(date_utc = zoom_date_string(), to = 12L, include_disturbances = FALSE) {
  url <- paste0(
    "https://zoom.earth/data/storms/?date=",
    utils::URLencode(date_utc, reserved = TRUE),
    "&to=", as.integer(to)
  )

  a <- rjson::fromJSON(file = url)

  ids <- unique(as.character(unlist(a$storms %||% character())))
  if (isTRUE(include_disturbances)) {
    ids <- unique(c(ids, as.character(unlist(a$disturbances %||% character()))))
  }
  ids <- ids[!is.na(ids) & ids != "NULL" & nzchar(ids)]
  ids
}

zoom_guess_name <- function(a, storm_id) {
  nm <- pick_one(a, c("name", "title"), default = NA_character_)
  if (is.na(nm)) {
    nm <- sub("-[0-9]{4}$", "", storm_id)
    nm <- gsub("-", " ", nm)
    nm <- tools::toTitleCase(nm)
  }
  nm
}

zoom_json2spatVect <- function(a, storm_id,
                               wind_units = c("kmh", "kts", "ms"),
                               default_rmw_km = NA_real_) {

  wind_units <- match.arg(wind_units)

  tr <- a$track
  if (is.null(tr) || length(tr) < 2) stop("Not enough points in track for ", storm_id)

  pts <- lapply(tr, function(pt) {
    coords <- pt$coordinates %||% pt$coord %||% pt$coords
    if (is.null(coords) || length(coords) < 2) return(NULL)

    lon <- suppressWarnings(as.numeric(coords[[1]]))
    lat <- suppressWarnings(as.numeric(coords[[2]]))
    if (!is.finite(lon) || !is.finite(lat)) return(NULL)

    iso_time <- zoom_as_iso(pick_one(pt, c("date", "datetime", "time", "timestamp"), NA_character_))

    wind_raw <- suppressWarnings(as.numeric(
      pick_one(pt, c("wind", "windSpeed", "max_wind_speed", "vmax"), NA_real_)
    ))

    pressure <- suppressWarnings(as.numeric(
      pick_one(pt, c("pressure", "mslp", "centralPressure"), NA_real_)
    ))

    forecast_flag <- FALSE
    fc_val <- pt$forecast
    if (!is.null(fc_val)) {
      if (is.logical(fc_val)) forecast_flag <- isTRUE(fc_val[1])
      if (is.character(fc_val)) forecast_flag <- tolower(fc_val[1]) %in% c("true", "1", "yes")
      if (is.numeric(fc_val)) forecast_flag <- isTRUE(fc_val[1] != 0)
    }

    rmw_km <- suppressWarnings(as.numeric(
      pick_one(pt, c("rmw", "radiusMaxWind", "rmax"), default_rmw_km)
    ))

    data.frame(
      lon = lon,
      lat = lat,
      ISO_TIME = iso_time,
      wind_raw = wind_raw,
      pressure = pressure,
      forecast = forecast_flag,
      rmw_km = rmw_km,
      stringsAsFactors = FALSE
    )
  })

  pts <- pts[!vapply(pts, is.null, logical(1))]
  if (length(pts) < 2) stop("Not enough valid points in track for ", storm_id)

  track <- do.call(rbind, pts)
  track <- track[is.finite(track$lon) & is.finite(track$lat), , drop = FALSE]
  if (NROW(track) < 2) stop("Not enough valid lon/lat points in track for ", storm_id)

  # wind -> m/s for TCHazaRds
  vms <- suppressWarnings(as.numeric(track$wind_raw))
  if (wind_units == "kmh") vms <- vms / 3.6
  if (wind_units == "kts") vms <- vms * 0.5144
  if (wind_units == "ms")  vms <- vms

  pres <- suppressWarnings(as.numeric(track$pressure))
  bad_pres <- !is.finite(pres)
  if (any(bad_pres)) pres[bad_pres] <- CpFromVmax(vms[bad_pres], Ep = 1010)

  dataset <- ifelse(isTRUE(track$forecast), "forecast", "track")
  nm <- zoom_guess_name(a, storm_id)

  # Build line segments between consecutive points
  nt <- NROW(track)
  segs <- lapply(1:(nt - 1), function(i) {
    terra::vect(
      cbind(
        c(track$lon[i], track$lon[i + 1]),
        c(track$lat[i], track$lat[i + 1])
      ),
      type = "lines",
      crs = "EPSG:4326"
    )
  })

  track_v <- do.call(rbind, segs)

  # Segment timing from point-to-point intervals
  pt_time <- as.POSIXct(track$ISO_TIME, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  dt <- diff(as.numeric(pt_time))
  dt[!is.finite(dt) | dt <= 0] <- NA_real_

  # Zoom Earth commonly appears to be 3-hourly, so use that if timestamps are missing
  fallback_dt <- suppressWarnings(stats::median(dt[is.finite(dt) & dt > 0], na.rm = TRUE))
  if (!is.finite(fallback_dt)) fallback_dt <- 3 * 3600
  dt[!is.finite(dt) | dt <= 0] <- fallback_dt

  # Attributes expected by downstream code / TCHazaRds
  track_v$DATASET <- dataset[1:(nt - 1)]
  track_v$WID <- storm_id
  track_v$NAME <- nm
  track_v$PRES <- pres[1:(nt - 1)]

  # Zoom Earth does not expose the same WMO-style wind-radii string.
  # If TCHazaRds in your setup needs a non-NA radius proxy, set DEFAULT_RMW_KM.
  rmw_seg <- track$rmw_km[1:(nt - 1)]
  track_v$RMW_STR <- ifelse(is.finite(rmw_seg), as.character(rmw_seg), NA_character_)

  track_v$ISO_TIME <- track$ISO_TIME[1:(nt - 1)]
  track_v$forecast_time <- ifelse(
    track_v$DATASET == "forecast",
    track_v$ISO_TIME,
    NA_character_
  )
  track_v$dt <- dt

  gg <- geom(track_v)
  ggs <- seq(1, NROW(gg), 2)
  track_v$LON <- gg[, 3][ggs]
  track_v$LAT <- gg[, 4][ggs]

  track_v$max_wind_speed <- vms[1:(nt - 1)]
  track_v$STORM_SPD <- perim(track_v) / track_v$dt
  track_v$thetaFm <- 90 - TCHazaRds::returnBearing(track_v)

  track_v
}

ZOOM_TC_vl <- function(TCid,
                       wind_units = zoom_wind_units,
                       default_rmw_km = default_rmw_km) {
  url <- paste0("https://zoom.earth/data/storms/?id=", TCid)
  a <- rjson::fromJSON(file = url)
  zoom_json2spatVect(a, storm_id = TCid, wind_units = wind_units, default_rmw_km = default_rmw_km)
}

centre_for_tcid <- function(tcid) "Zoom Earth"

# --- download Zoom Earth in-force list ----
zoom_date <- Sys.getenv("ZOOM_DATE", unset = zoom_date_string())
include_disturbances <- Sys.getenv("ZOOM_INCLUDE_DISTURBANCES", unset = "0") == "1"
TCids <- zoom_list_ids(date_utc = zoom_date, to = zoom_to, include_disturbances = include_disturbances)

message("Using Zoom Earth date: ", zoom_date)
if (length(TCids) == 0) message("No active storms returned by Zoom Earth.")

# Read params once
paramsTable <- read.csv(system.file("extdata/tuningParams/defult_params.csv", package = "TCHazaRds"))

# Optionally limit how many storms to process (useful for testing)
max_tc <- suppressWarnings(as.integer(Sys.getenv("MAX_TC", unset = NA)))
if (!is.na(max_tc) && max_tc > 0) TCids <- head(TCids, max_tc)

index_rows <- list()

for (i in seq_along(TCids)) {
  TCid <- as.character(TCids[i])

  tc_gpkg <- file.path(outdir, sprintf("TC_%03d.gpkg", i))
  hs0_jpg <- file.path(outdir, sprintf("Hs0_%03d.jpg", i))
  ib_jpg  <- file.path(outdir, sprintf("IB_%03d.jpg",  i))

  # Compact summary rasters for the report (single-layer GeoTIFFs)
  summary_dir <- file.path(outdir)
  hs0_max_tif <- file.path(summary_dir, sprintf("hs0_max_%03d.tif", i))
  hs0_tpk_tif <- file.path(summary_dir, sprintf("hs0_tpkHr_%03d.tif", i))
  ibs_max_tif <- file.path(summary_dir, sprintf("ibs_max_%03d.tif", i))
  ibs_tpk_tif <- file.path(summary_dir, sprintf("ibs_tpkHr_%03d.tif", i))
  haz_nc <- file.path(outdir, sprintf("haz_%03d.nc", i))

  # Optional GeoTIFFs (multi-layer)
  haz_hs0_tif <- file.path(outdir, sprintf("haz_Hs0_%03d.tif", i))
  haz_ib_tif  <- file.path(outdir, sprintf("haz_IB_%03d.tif",  i))

  row <- data.frame(
    idx = i,
    TCid = TCid,
    name = NA_character_,
    centre = centre_for_tcid(TCid),
    gpkg = tc_gpkg,
    hs0_jpg = hs0_jpg,
    ib_jpg = ib_jpg,
    haz_nc = if (save_nc) haz_nc else NA_character_,
    haz_hs0_tif = if (save_tif) haz_hs0_tif else NA_character_,
    haz_ib_tif  = if (save_tif) haz_ib_tif else NA_character_,
    hs0_max_tif = hs0_max_tif,
    hs0_tpk_tif = hs0_tpk_tif,
    ibs_max_tif = ibs_max_tif,
    ibs_tpk_tif = ibs_tpk_tif,
    hs0_peak_time = NA_character_,
    ibs_peak_time = NA_character_,
    status = "ok",
    message = "",
    stringsAsFactors = FALSE
  )

  message("Processing ", TCid, " (", i, "/", length(TCids), ")")

  tryCatch({
    TC <- ZOOM_TC_vl(TCid)

    nm <- "NONAME"
    if (any(TC$NAME != "" & !is.na(TC$NAME))) nm <- TC$NAME[TC$NAME != "" & !is.na(TC$NAME)][1]
    row$name <- nm

    TC$TC_CENTRE <- row$centre
    file.remove(tc_gpkg)
    Sys.sleep(2)
    terra::writeVector(TC, tc_gpkg, overwrite = TRUE)

    # Raster grid around track
    r <- rast(buffer(TC, 200000), res = 0.2)
    values(r) <- 0

    land_v <- vect("osm_land_polygons_simplifyGeom_0point005_areaGT1e6_aggregated/")
    land_r <- rasterize(land_v, r, touches = TRUE, background = 0)
    inland_proximity <- terra::costDist(land_r, target = 0, scale = 1)
    GEO_land <- land_geometry(land_r, inland_proximity)

    # If there is only one line segment, duplicate it so downstream code has 2 rows to work with
    TC_for_haz <- TC
    if (nrow(TC_for_haz) == 1) {
      TC_for_haz <- rbind(TC_for_haz, TC_for_haz)
      t0 <- as.POSIXct(TC_for_haz$ISO_TIME[1], tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
      if (!is.na(t0)) {
        TC_for_haz$ISO_TIME[2] <- format(t0 + 10 * 60, "%Y-%m-%d %H:%M:%S", tz = "UTC")
      }
      if (!is.na(TC_for_haz$dt[1])) TC_for_haz$dt[2] <- TC_for_haz$dt[1]
    }

    # Time range
    t_all <- as.POSIXct(TC_for_haz$ISO_TIME, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
    tmin <- min(t_all, na.rm = TRUE)
    tmax <- max(t_all, na.rm = TRUE)
    outdate <- sort(unique(c(tmin, seq(tmin, tmax, by = 1 * 3600), tmax)))
    outdate_chr <- format(outdate, "%Y-%m-%dT%H:%M:%SZ")

    # Write NetCDF directly (TCHazaRds supports outfile)
    haz <- TCHazaRdsWindFields(
      outdate = outdate,
      TC = TC_for_haz,
      GEO_land = GEO_land,
      paramsTable = paramsTable,
      return_vars = c("Pr", "Hs0"),
      outfile = if (save_nc) haz_nc else NULL,
      overwrite = TRUE
    )

    IB <- (1010 - haz$Pr) / 100
    IB <- IB + (haz$Hs0 / haz$Hs0 - 1)

    # Name layers with time strings (helps later)
    for (vn in names(haz)) {
      try(names(haz[[vn]]) <- outdate_chr, silent = TRUE)
    }

    forecast_start <- strptime(rev(TC$ISO_TIME[TC$DATASET == "track"])[1], "%Y-%m-%d %H:%M:%S", tz = "UTC")
    fsi <- which(outdate >= forecast_start)[1]
    if (is.na(fsi)) fsi <- 1

    # Peak times (spatial extrema per time step)
    hs0_by_t <- as.numeric(terra::global(haz$Hs0, "max", na.rm = TRUE)[, 1])
    ibs_by_t <- as.numeric(terra::global(IB, "max", na.rm = TRUE)[, 1])

    hs0_peak_i <- if (all(is.na(hs0_by_t))) NA_integer_ else which.max(hs0_by_t)
    ibs_peak_i <- if (all(is.na(ibs_by_t))) NA_integer_ else which.max(ibs_by_t)

    if (!is.na(hs0_peak_i)) row$hs0_peak_time <- outdate_chr[hs0_peak_i]
    if (!is.na(ibs_peak_i)) row$ibs_peak_time <- outdate_chr[ibs_peak_i]

    # Quick-look maps
    Hs0_max <- terra::app(haz$Hs0, fun = max, na.rm = TRUE)
    IB_max  <- terra::app(IB,      fun = max, na.rm = TRUE)

    Hs0_which_max <- terra::app(haz$Hs0, fun = which.max, na.rm = TRUE) - fsi
    IB_which_max  <- terra::app(IB,      fun = which.max, na.rm = TRUE) - fsi

    # 200 km buffer outline for plots (buffer in metres in EPSG:3857)
    buf200 <- try(terra::buffer(terra::project(TC, "EPSG:3857"), 200000), silent = TRUE)
    if (!inherits(buf200, "try-error")) {
      buf200 <- try(terra::project(buf200, crs(Hs0_max)), silent = TRUE)
      if (inherits(buf200, "try-error")) buf200 <- NULL
    } else {
      buf200 <- NULL
    }

    file.remove(hs0_jpg)
    Sys.sleep(2)
    safe_jpeg(hs0_jpg, function() {
      plot(Hs0_max, main = paste0("Hs0 max (", nm, " / ", TCid, ")"))
      if (!is.null(buf200)) plot(terra::aggregate(buf200), add = TRUE, border = "grey40", lwd = 2)
      plot(TC, add = TRUE)
    })

    file.remove(ib_jpg)
    Sys.sleep(2)
    safe_jpeg(ib_jpg, function() {
      plot(IB_max, main = paste0("IB surge max (", nm, " / ", TCid, ")"))
      if (!is.null(buf200)) plot(terra::aggregate(buf200), add = TRUE, border = "grey40", lwd = 2)
      plot(TC, add = TRUE)
    })

    # Optional GeoTIFFs (Hs0 + IB only; multi-layer)
    if (save_tif) {
      file.remove(haz_hs0_tif)
      file.remove(haz_ib_tif)
      file.remove(hs0_tpk_tif)
      file.remove(ibs_tpk_tif)

      Sys.sleep(2)

      write_haz_raster(Hs0_max, haz_hs0_tif)
      write_haz_raster(IB_max, haz_ib_tif)

      write_haz_raster(Hs0_which_max, hs0_tpk_tif)
      write_haz_raster(IB_which_max, ibs_tpk_tif)
    }

  }, error = function(e) {
    row$status <- "error"
    row$message <- conditionMessage(e)
    message("ERROR for ", TCid, ": ", row$message)
  })

  index_rows[[i]] <- row
}

# Always write index/sessionInfo even if TCids is empty
index_df <- if (length(index_rows) > 0) do.call(rbind, index_rows) else data.frame()
write.csv(index_df, file.path(outdir, "index.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

message("Done. Wrote outputs to: ", normalizePath(outdir))

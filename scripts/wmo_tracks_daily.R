suppressPackageStartupMessages({
  library(rjson)
  library(terra)
  library(TCHazaRds)
})

Sys.setenv(TZ = "UTC")

# Output folder inside the repo (works on GitHub runners)
outdir <- Sys.getenv("OUTDIR", unset = file.path("outputs", "latest"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Save full hazard rasters? (set SAVE_HAZ=0 to skip)
save_haz <- Sys.getenv("SAVE_HAZ", unset = "1") == "1"

# --- helpers ----

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
  # Write a multi-layer GeoTIFF with compression + tiling
  terra::writeRaster(
    x, filename, overwrite = TRUE,
    wopt = list(
      gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6", "TILED=YES"),
      datatype = "FLT4S"
    )
  )
}

json2spatVect <- function(a, varb, previous_point = NULL) {

  track <- t(sapply(a[[varb]], function(x) x))
  nms <- colnames(track)

  track <- array(paste(track), dim = dim(track))
  colnames(track) <- nms
  track <- data.frame(track, stringsAsFactors = FALSE)

  # Prepend last observed point to the start of the forecast so line segments connect
  if (!is.null(previous_point) && NROW(previous_point) >= 1) {
    track <- rbind(NA, track)
    track$lng[1] <- previous_point$lng
    track$lat[1] <- previous_point$lat
    track$DATASET[1] <- previous_point$DATASET
    track$pressure[1] <- previous_point$pressure
    track$max_wind_speed[1] <- previous_point$max_wind_speed
    track$gust[1] <- previous_point$gust
    track$wind_radii[1] <- previous_point$wind_radii

    t0 <- strptime(previous_point$analysis_time, "%Y-%m-%d %H:%M:%S", tz = "UTC")
    track$time_interval[1] <- 0

    hr <- floor(as.numeric(track$time_interval) / 100)
    mn <- as.numeric(track$time_interval) - hr * 100
    ts <- t0 + hr * 3600 + mn * 60
    track$analysis_time <- paste(ts)
  }

  nt <- NROW(track)
  if (nt < 2) stop("Not enough points in ", varb)

  # Build line segments between consecutive points
  track_vl <- lapply(1:(nt - 1), function(i) {
    terra::vect(
      cbind(
        c(as.numeric(track$lng[i]), as.numeric(track$lng[i + 1])),
        c(as.numeric(track$lat[i]), as.numeric(track$lat[i + 1]))
      ),
      type = "lines",
      crs = "epsg:4283"
    )
  })

  track_lastpoint <- track[nt, , drop = FALSE]
  track_v <- do.call(rbind, track_vl)

  # Segment attributes (use the "start" point for each segment)
  track_v$DATASET <- varb
  track_v$WID <- track$tc_id[1:(nt - 1)]
  track_v$NAME <- track$tc_name[1:(nt - 1)]

  pres <- suppressWarnings(as.numeric(track$pressure[1:(nt - 1)]))
  vms  <- suppressWarnings(as.numeric(track$max_wind_speed[1:(nt - 1)])) * 0.5144 # knots -> m/s

  if (all(is.na(pres))) pres <- CpFromVmax(vms, Ep = 1010)
  track_v$PRES <- pres

  track_v$RMW_STR <- track$wind_radii[1:(nt - 1)]

  iso_time <- track$analysis_time[1:(nt - 1)]
  iso_time[nchar(iso_time) == 10] <- paste(iso_time[nchar(iso_time) == 10], "00:00:00")
  track_v$ISO_TIME <- iso_time

  rtm <- strptime(track_v$ISO_TIME, "%Y-%m-%d %H:%M:%S", tz = "UTC")
  dt <- diff(as.numeric(rtm))
  if (all(is.na(dt))) {
    dt <- suppressWarnings(as.numeric(track$time_interval[2:nt])) / 100 * 3600
  }
  track_v$dt <- dt

  gg <- geom(track_v)
  ggs <- seq(1, NROW(gg), 2)
  track_v$LON <- gg[, 3][ggs]
  track_v$LAT <- gg[, 4][ggs]

  track_v$max_wind_speed <- vms
  track_v$STORM_SPD <- perim(track_v) / track_v$dt
  track_v$thetaFm <- 90 - TCHazaRds::returnBearing(track_v)

  list(track_v, track_lastpoint)
}

# TC centre names from WMO page source
TC_centre_info <- c(
  "Hong Kong","Beijing","TC RSMC Honolulu","TC RSMC Miami","TC RSMC Tokyo",
  "TC RSMC New Delhi","TC RSMC La Réunion","TCWC Melbourne","TCWC Melbourne",
  "TCWC Melbourne","TC RSMC Nadi","TCWC Wellington","JTWC","Manila","Macao"
)

# --- download WMO in-force list ----

inforce <- rjson::fromJSON(file = "https://severeweather.wmo.int/json/tc_inforce.json")
dat <- t(sapply(inforce$inforce, function(x) x))
nTC <- if (!is.null(dim(dat))) dim(dat)[1] else 1

inforceTC <- data.frame(array(paste(dat), dim = c(nTC, length(dat) / nTC)), stringsAsFactors = FALSE)
names(inforceTC) <- inforce$fields

TCids <- c(inforceTC$sysid, inforceTC$same)
TCids <- TCids[TCids != "NULL" & !is.na(TCids)]
TCids <- unique(TCids)

centre_for_tcid <- function(tcid) {
  idx <- which(inforceTC$sysid == tcid | inforceTC$same == tcid)[1]
  if (length(idx) == 0 || is.na(idx)) return(NA_character_)
  cid <- suppressWarnings(as.integer(inforceTC$centerid[idx]))
  if (is.na(cid) || cid < 1 || cid > length(TC_centre_info)) return(NA_character_)
  TC_centre_info[cid]
}

WMO_TC_vl <- function(TCid) {
  a <- rjson::fromJSON(file = paste0("https://severeweather.wmo.int/v2/json/tc_", TCid, ".json"))

  vnl  <- json2spatVect(a, "track")
  track <- vnl[[1]]

  vnl2 <- json2spatVect(a, "forecast", vnl[[2]])
  forecast <- vnl2[[1]]

  # Combine as one SpatVector
  rbind(track, forecast)
}

# Read params once
paramsTable <- read.csv(system.file("extdata/tuningParams/defult_params.csv", package = "TCHazaRds"))

# Optionally limit how many storms to process (useful for testing)
max_tc <- suppressWarnings(as.integer(Sys.getenv("MAX_TC", unset = NA)))
if (!is.na(max_tc) && max_tc > 0) TCids <- head(TCids, max_tc)

# --- main loop ----

index_rows <- list()

for (i in seq_along(TCids)) {
  TCid <- TCids[i]

  tc_gpkg <- file.path(outdir, sprintf("TC_%03d.gpkg", i))
  hs0_jpg <- file.path(outdir, sprintf("Hs0_%03d.jpg", i))
  pr_jpg  <- file.path(outdir, sprintf("Pr_%03d.jpg",  i))

  # Full hazard rasters (multi-layer)
  haz_hs0_tif <- file.path(outdir, sprintf("haz_Hs0_%03d.tif", i))
  haz_pr_tif  <- file.path(outdir, sprintf("haz_Pr_%03d.tif",  i))
  haz_uw_tif  <- file.path(outdir, sprintf("haz_Uw_%03d.tif",  i))
  haz_vw_tif  <- file.path(outdir, sprintf("haz_Vw_%03d.tif",  i))
  haz_sw_tif  <- file.path(outdir, sprintf("haz_Sw_%03d.tif",  i))
  haz_dw_tif  <- file.path(outdir, sprintf("haz_Dw_%03d.tif",  i))

  row <- data.frame(
    idx = i,
    TCid = TCid,
    name = NA_character_,
    centre = centre_for_tcid(TCid),
    gpkg = tc_gpkg,
    hs0_jpg = hs0_jpg,
    pr_jpg = pr_jpg,
    haz_hs0_tif = if (save_haz) haz_hs0_tif else NA_character_,
    haz_pr_tif  = if (save_haz) haz_pr_tif  else NA_character_,
    haz_uw_tif  = if (save_haz) haz_uw_tif  else NA_character_,
    haz_vw_tif  = if (save_haz) haz_vw_tif  else NA_character_,
    haz_sw_tif  = if (save_haz) haz_sw_tif  else NA_character_,
    haz_dw_tif  = if (save_haz) haz_dw_tif  else NA_character_,
    hs0_peak_time = NA_character_,
    pr_min_time = NA_character_,
    status = "ok",
    message = "",
    stringsAsFactors = FALSE
  )

  message("Processing ", TCid, " (", i, "/", length(TCids), ")")

  tryCatch({
    TC <- WMO_TC_vl(TCid)

    # Name
    nm <- "NONAME"
    if (any(TC$NAME != "" & !is.na(TC$NAME))) nm <- TC$NAME[TC$NAME != "" & !is.na(TC$NAME)][1]
    row$name <- nm

    # Add centre field to every feature
    TC$TC_CENTRE <- row$centre

    # Write vector
    terra::writeVector(TC, tc_gpkg, overwrite = TRUE)

    # Build a modest raster grid around the track
    n_hoz_pix <- 300
    r <- rast(TC, ncol = n_hoz_pix)
    r <- rast(TC, res = rep(res(r)[1], 2))
    values(r) <- 0

    GEO_land <- land_geometry(r, r)

    # Time range for modelling
    t_all <- as.POSIXct(TC$ISO_TIME, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
    tmin <- min(t_all, na.rm = TRUE)
    tmax <- max(t_all, na.rm = TRUE)
    outdate <- seq(tmin, tmax, by = 3 * 3600)
    outdate_chr <- format(outdate, "%Y-%m-%dT%H:%M:%SZ")

    haz <- TCHazaRdsWindFields(
      outdate = outdate,
      TC = TC,
      GEO_land = GEO_land,
      paramsTable = paramsTable,
      return_vars = c("Pr", "Uw", "Vw", "Sw", "Dw", "Hs0")
    )

    # Name the layers so the time is carried through to the GeoTIFF bands
    for (vn in names(haz)) {
      try(names(haz[[vn]]) <- outdate_chr, silent = TRUE)
    }

    # Identify "peak" times (based on spatial extrema per time step)
    hs0_by_t <- as.numeric(terra::global(haz$Hs0, "max", na.rm = TRUE)[, 1])
    pr_by_t  <- as.numeric(terra::global(haz$Pr,  "min", na.rm = TRUE)[, 1])

    hs0_peak_i <- if (all(is.na(hs0_by_t))) NA_integer_ else which.max(hs0_by_t)
    pr_min_i   <- if (all(is.na(pr_by_t)))  NA_integer_ else which.min(pr_by_t)

    if (!is.na(hs0_peak_i)) row$hs0_peak_time <- outdate_chr[hs0_peak_i]
    if (!is.na(pr_min_i))   row$pr_min_time   <- outdate_chr[pr_min_i]

    # Summary rasters
    Hs0_max <- terra::app(haz$Hs0, fun = max, na.rm = TRUE)
    Pr_min  <- terra::app(haz$Pr,  fun = min, na.rm = TRUE)

    # Plots (quick-look)
    safe_jpeg(hs0_jpg, function() {
      plot(Hs0_max, main = paste0("Hs0 max (", nm, " / ", TCid, ")"))
      plot(TC, add = TRUE)
    })

    safe_jpeg(pr_jpg, function() {
      plot(Pr_min, main = paste0("Pressure min (", nm, " / ", TCid, ")"))
      plot(TC, add = TRUE)
    })

    # Save full hazard rasters (multi-layer GeoTIFFs)
    if (save_haz) {
      write_haz_raster(haz$Hs0, haz_hs0_tif)
      write_haz_raster(haz$Pr,  haz_pr_tif)
      write_haz_raster(haz$Uw,  haz_uw_tif)
      write_haz_raster(haz$Vw,  haz_vw_tif)
      write_haz_raster(haz$Sw,  haz_sw_tif)
      write_haz_raster(haz$Dw,  haz_dw_tif)
    }

  }, error = function(e) {
    row$status <- "error"
    row$message <- conditionMessage(e)
    message("ERROR for ", TCid, ": ", row$message)
  })

  index_rows[[i]] <- row
}

index_df <- do.call(rbind, index_rows)
write.csv(index_df, file.path(outdir, "index.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

message("Done. Wrote outputs to: ", normalizePath(outdir))

suppressPackageStartupMessages({
  library(rjson)
  library(terra)
  library(TCHazaRds)
})

Sys.setenv(TZ = "UTC")

# Output folder inside the repo (works on GitHub runners)
outdir <- Sys.getenv("OUTDIR", unset = file.path("outputs", "latest"))
unlink(outdir,recursive=TRUE,force = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Controls (set via GitHub Actions env)
save_tif <- Sys.getenv("SAVE_TIF", unset = "1") == "1"   # default on (not too big)
save_nc  <- Sys.getenv("SAVE_NC",  unset = "0") == "1"   # default off (NetCDF can exceed GitHub limit)
save_sum <- Sys.getenv("SAVE_SUMMARY", unset = "1") == "1" # default on (compact summary GeoTIFFs)

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

json2spatVect <- function(a, varb, previous_point = NULL) {
  track <- t(sapply(a[[varb]], function(x) x))
  nms <- colnames(track)

  track <- array(paste(track), dim = dim(track))
  colnames(track) <- nms
  track <- data.frame(track, stringsAsFactors = FALSE)
  track <- track[!is.na(track$lng),]
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
  
  if(is.null(iso_time[1])) iso_time = format(as.POSIXct(Sys.time(), tz = "UTC")+
                                               as.numeric(track$time_interval)/100*3600,
                                             "%Y-%m-%d %H:00:00")
  iso_time[nchar(iso_time) == 10] <- paste(iso_time[nchar(iso_time) == 10], "00:00:00")
  
  track_v$ISO_TIME <- iso_time
  if(!is.null(track$forecast_time[1]))
    if(!is.na(track$forecast_time[1]))
    if(track$forecast_time[1] == "NULL") track_v$forecast_time <- iso_time
  if(!is.null(track$forecast_time[1]))
    if(is.na(track$forecast_time[1])) track_v$forecast_time <- iso_time
      

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

TCids <- as.numeric(c(inforceTC$sysid, inforceTC$same))
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
  a <- rjson::fromJSON(file = paste0("https://severeweather.wmo.int/json/tc_", TCid, ".json"))
  vnl <- NULL
  try(vnl  <- json2spatVect(a, "track"))
  track <- vnl[[1]]
  vnl2 <- NULL
  try(vnl2 <- json2spatVect(a=a, varb="forecast", previous_point = vnl[[2]]))
  if(!is.null(vnl) & !is.null(vnl2)){
    forecast <- vnl2[[1]]
    track = rbind(track, forecast)
  }
  if(is.null(vnl)) track = vnl2[[1]]
  return(track)
}

# Read params once
paramsTable <- read.csv(system.file("extdata/tuningParams/defult_params.csv", package = "TCHazaRds"))

# Optionally limit how many storms to process (useful for testing)
max_tc <- suppressWarnings(as.integer(Sys.getenv("MAX_TC", unset = NA)))
if (!is.na(max_tc) && max_tc > 0) TCids <- head(TCids, max_tc)

index_rows <- list()

for (i in seq_along(TCids)) {
  TCid <- TCids[i]

  tc_gpkg <- file.path(outdir, sprintf("TC_%03d.gpkg", i))
  hs0_jpg <- file.path(outdir, sprintf("Hs0_%03d.jpg", i))
  ib_jpg  <- file.path(outdir, sprintf("IB_%03d.jpg",  i))


    # Compact summary rasters for the report (single-layer GeoTIFFs)
    summary_dir <- file.path(outdir)
    hs0_max_tif <- file.path(summary_dir, sprintf("hs0_max_%03d.tif", i))
    hs0_tpk_tif <- file.path(summary_dir, sprintf("hs0_tpkHr_%03d.tif", i))
    ibs_max_tif <- file.path(summary_dir, sprintf("ibs_max_%03d.tif", i))
    ibs_tpk_tif <- file.path(summary_dir, sprintf("ibs_tpkHr_%03d.tif", i))
  haz_nc  <- file.path(outdir, sprintf("haz_%03d.nc", i))

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
    haz_ib_tif  = if (save_tif) haz_ib_tif  else NA_character_,
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
    TC <- WMO_TC_vl(TCid)

    nm <- "NONAME"
    if (any(TC$NAME != "" & !is.na(TC$NAME))) nm <- TC$NAME[TC$NAME != "" & !is.na(TC$NAME)][1]
    row$name <- nm

    TC$TC_CENTRE <- row$centre
    file.remove(tc_gpkg)
    Sys.sleep(2)
    terra::writeVector(TC, tc_gpkg, overwrite = TRUE)

    # Raster grid around track
    r <- rast(buffer(TC,200000), res = .2)
    #r <- rast(TC, res = rep(res(r)[1], 2))
    values(r) <- 0
    
    land_v <- vect("osm_land_polygons_simplifyGeom_0point005_areaGT1e6_aggregated/")
    land_r = rasterize(land_v,r,touches=TRUE,background=0)
    inland_proximity = terra::costDist(land_r,target = 0,scale=1)
    GEO_land = land_geometry(land_r,inland_proximity)
    
    # Time range
    t_all <- as.POSIXct(TC$ISO_TIME, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
    tmin <- min(t_all, na.rm = TRUE)
    tmax <- max(t_all, na.rm = TRUE)
    outdate <- seq(tmin, tmax, by = 1 * 3600)
    outdate_chr <- format(outdate, "%Y-%m-%dT%H:%M:%SZ")

    # Write NetCDF directly (TCHazaRds supports outfile)
    
    if(length(TC) == 1) {
      TC2 = rbind(TC,TC)
      TC2$ISO_TIME[2] = paste0(substr(TC2$ISO_TIME[2],1,17),"10")
    }
    
    haz <- TCHazaRdsWindFields(
      outdate = outdate,
      TC = TC,
      GEO_land = GEO_land,
      paramsTable = paramsTable,
      return_vars = c("Pr", "Hs0"),
      outfile = if (save_nc) haz_nc else NULL,
      overwrite = TRUE
    )
    IB = (1010-haz$Pr)/100
    IB = IB + (haz$Hs0/haz$Hs0-1)
    # Name layers with time strings (helps later)
    for (vn in names(haz)) {
      try(names(haz[[vn]]) <- outdate_chr, silent = TRUE)
    }
    forcast_start = strptime(rev(TC$ISO_TIME[TC$DATASET == "track"])[1],"%Y-%m-%d %H:%M:%S",tz = "UTC")
    fsi = which(outdate >= forcast_start)[1]
    if(is.na(fsi)) fsi = 1
    
    # Peak times (spatial extrema per time step)
    hs0_by_t <- as.numeric(terra::global(haz$Hs0, "max", na.rm = TRUE)[, 1])
    ibs_by_t <- as.numeric(terra::global(IB,      "max", na.rm = TRUE)[, 1])

    hs0_peak_i <- if (all(is.na(hs0_by_t))) NA_integer_ else which.max(hs0_by_t)
    ibs_peak_i <- if (all(is.na(ibs_by_t))) NA_integer_ else which.max(ibs_by_t)
   
    if (!is.na(hs0_peak_i)) row$hs0_peak_time <- outdate_chr[hs0_peak_i]
    if (!is.na(ibs_peak_i)) row$ibs_peak_time  <- outdate_chr[ibs_peak_i]

    # Quick-look maps
    Hs0_max <- terra::app(haz$Hs0, fun = max, na.rm = TRUE)
    IB_max <-  terra::app(IB    ,  fun = max, na.rm = TRUE)
    
    Hs0_which_max <- terra::app(haz$Hs0, fun = which.max, na.rm = TRUE)-fsi
    IB_which_max <-  terra::app(IB    ,  fun = which.max, na.rm = TRUE)-fsi
     
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
      write_haz_raster(IB_max,  haz_ib_tif)

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

require(terra)
require(TCHazaRds)

library(RColorBrewer)

cols <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


TC = vect("outputs/latest/TC_001.gpkg")


r <- rast(xmin = 142.5,xmax = 147.5,ymin = -16,ymax = -9, res = .04)
r <- rast(xmin = 111, xmax = 123, ymin = -28, ymax = -15, res = .04) 
TC <- crop(TC,buffer(as.polygons(ext(r)),2))
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
outdate <- seq(tmin, tmax, by = 0.5 * 3600)
outdate_chr <- format(outdate, "%Y-%m-%dT%H:%M:%SZ")

# Write NetCDF directly (TCHazaRds supports outfile)

paramsTable <- read.csv(system.file("extdata/tuningParams/defult_params.csv", package = "TCHazaRds"))
paramsTable[c(5,7:8),2] = 2

TC$RMW = TCHazaRds::rMax_modelsR(rMaxModel = 1,TClats = TC$LAT,cPs = TC$PRES,eP = paramsTable$value[1])
TC$VMAX = TC$max_wind_speed
TC$RMAX2 = TCHazaRds::rMax2_modelsR(rMax2Model = 2,vMax = TC$VMAX,rMax = TC$RMW,TClats = TC$LAT)

haz <- TCHazaRdsWindFields(
  outdate = outdate,
  TC = TC,
  GEO_land = GEO_land,
  paramsTable = paramsTable,
  return_vars = c("Pr", "Hs0"),
)


maxhs = max(haz$Hs0)


#my bathy data.
dem <- rast("C:/Users/ogr013/OneDrive - CSIRO/Documents/OZ250_Merge_20210725.tif")
demp = terra::resample(dem,r)

maxhs2 = min(c(-0.2*demp,maxhs))
pt = c(maxhs,maxhs2)
names(pt) = c("Hs0","Hs0 20% depth limited")
plot(pt,col = cols(14*2),range = c(0,14))
plot(maxhs2,col = cols(14*2),range = c(0,14))

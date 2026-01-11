require(ncdf4)
require(terra)
require(TideHarmonics)

track = vect("outputs/latest/TC_002.gpkg")
track = buffer(vect(cbind(115.7386,-31.8256),"points",crs = rast()),10000)


tc_ext = as.numeric(ext(track)[1:4])+0

GTSM_con = function(constituent="M2"){
  nc = nc_open(paste0("http://opendap.avi.deltares.nl/thredds/dodsC/opendap/deltares/GTSM/GTSMv4.1_tide/GTSMv4.1_tide_2014_",constituent,"_rasterized.nc"))
  lat = nc$dim$lat$vals
  lon = nc$dim$lon$vals
  
  lono = lon[which(lon >= tc_ext[1] & lon <= tc_ext[2])]
  lato = lat[which(lat >= tc_ext[3] & lat <= tc_ext[4])]

  amp =  ncvar_get(nc,"wl_amp",start = c(which(lon > tc_ext[1])[1]-1,which(lat > tc_ext[3])[1]-1),
                      count =   c(length(which(lon >= tc_ext[1] & lon <= tc_ext[2])),
                                  length(which(lat >= tc_ext[3] & lat <= tc_ext[4]))))
  
  phs = ncvar_get(nc,"wl_phs",start = c(which(lon > tc_ext[1])[1]-1,which(lat > tc_ext[3])[1]-1),
                  count =   c(length(which(lon >= tc_ext[1] & lon <= tc_ext[2])),
                              length(which(lat >= tc_ext[3] & lat <= tc_ext[4]))))
  
  nc_close(nc)
  
  r_amp = rast(t(amp)[length(lato):1,],extent = c(min(lono),max(lono),min(lato),max(lato)))
  r_phs = rast(t(phs)[length(lato):1,],extent = c(min(lono),max(lono),min(lato),max(lato)))
  out = c(r_amp,r_phs)
  names(out) = c("amp","phs")
  return(out)
}

M2 = GTSM_con("M2")
S2 = GTSM_con("S2")
K1 = GTSM_con("K1")
O1 = GTSM_con("O1")

r = c(M2,S2,K1,O1)

hfit1 <- ftide(Hillarys$SeaLevel, Hillarys$DateTime, hc4)
t1 <- as.POSIXct("2012-12-31 23:00", tz = "UTC")
t2 <- as.POSIXct("2013-01-03 14:00", tz = "UTC")

rftide = function(x){
  co = as.matrix(t(array(x,dim = c(2,4))))
  row.names(co) <- c("M2","S2","K1","O1")
  colnames(co) <- c("amplitude","phase")
  hfit1$co = co
  z = predict.tide2(hfit1, t1, t2,by = 1)
  return(z)
}
r1 = aggregate(r,4)
rt = terra::app(r1,rftide)
terra::time(rt) <- seq(t1,t2,"hour")
plot(as.numeric(rt[10,10,]),type = "l")
plot(rt)


#coef(hfit1, hc = TRUE)[, 1:2] = 1
  
i=10;j=100
co = as.matrix(rbind(M2[i,j,],S2[i,j,],K1[i,j,],O1[i,j,]))
row.names(co) = c("M2","S2","K1","O1")
colnames(co) = c("amplitude","phase")
hfit1$co = co

TideHarmonics::hc4

t1 <- as.POSIXct("2012-12-31 23:00", tz = "UTC")
t2 <- as.POSIXct("2013-01-01 14:00", tz = "UTC")
z = predict.tide2(hfit1, t1, t2)



  

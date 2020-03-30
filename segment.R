## Try watershed algorithm in R
## https://cran.r-project.org/web/packages/ForestTools/vignettes/treetopAnalysis.html
library(ForestTools)
library(raster)
library(rgdal)

## load RGB and DEM
#dsm <- (raster(x = "data/Orchard_IdahoTM_DEM_clip.tif") - 974.4) #remove ground
dsm <- raster(x = "data/Orchard_IdahoTM_DEM_clip.tif") #2019 DEM data
#dsm <- (raster(x = "data/orchard2015_dem_IDTM_clip.tif") - 955.5) #remove ground
#dsm <- raster(x = "data/orchard2015_dem_IDTM_clip.tif") #2015 DEM data

plot(dsm)

#the slope/trend is much worse for 2015, need to remove that before moving forward
dsm.df <- as.data.frame(dsm, na.rm=FALSE, xy=TRUE) #convert raster to dataframe
#trend <- lm(orchard2015_dem_IDTM_clip ~ x + y, data=dsm.df) #model XY trend
trend <- lm(Orchard_IdahoTM_DEM_clip ~ x + y, data=dsm.df) #model XY trend
summary(trend)

dsm.df[!(is.na(dsm.df[,3])),3] <- trend$residuals
dsm2 <- rasterFromXYZ(dsm.df, crs=crs(dsm))
plot(dsm2)

lin <- function(x){x * 0.54 - 0.0044} #function for relating max crown height to radius

#search for crowns or "treetops" throughout the image, this step takes the longest 
ttops <- vwf(CHM = dsm2, winFun = lin, minHeight = 0.3, maxWinDiameter = NULL) #minHeight is the parameter to pick out crowns, higher will make fewer crowns

plot(dsm2 > 0.15, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
plot(ttops, col = "blue", pch = 20, cex = 0.5, add = TRUE)

mean(ttops$height)

crowns <- mcws(treetops = ttops, CHM = dsm2, minHeight = 0.2, verbose = FALSE) #minHeight here will change which pixels are picked up as shrub, a higher number means smaller polygons around each crown/ttop
plot(crowns, col = sample(rainbow(50), length(unique(crowns[])), replace = TRUE), legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n')

# crownsPoly05 <- mcws(treetops = ttops05, CHM = dsm, format = "polygons", minHeight = 0.45, verbose = FALSE)
# there's an error with the packages "format="polygons"" algorithm, it doesn't correctly convert
# therefore, trying the rasterToPolygon function
library(rgeos)
crownsPoly <- rasterToPolygons(crowns, dissolve=TRUE) #note that some crowns are multi-part

plot(dsm2, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
plot(crownsPoly, border = "blue", lwd = 0.5, add = TRUE)

crownsPoly[["crownArea"]] <- gArea(crownsPoly, byid = TRUE) #sums area of multi-part crowns

crownsPoly[["crownDiameter"]] <- sqrt(crownsPoly[["crownArea"]] / pi) * 2
mean(crownsPoly$crownDiameter)

#save polygon shapefile of segmented crowns
writeOGR(obj=crownsPoly, dsn="tempdir", layer="crowns_2019_trendremove", drive="ESRI Shapefile", overwrite=TRUE)
#writeRaster(crowns, filename="tempdir/crowns_2019.tif", datatype="INT2U", overwrite=TRUE)




## Failed attempt to use itcSegment
# library(itcSegment)
# 
# GDALinfo("data/Orchard_IdahoTM_DEM_clip.tif")
# dsm <- raster(x = "data/Orchard_IdahoTM_DEM_clip.tif")
# 
# plot(dsm)
# 
# dsm <- projectRaster(dsm, crs=crs('+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
# crs(dsm)
# 
# itcTest <- itcIMG(dsm, epsg=26911, DIST=2, searchWinSize=9) #takes way too many resources and never finishes
# summary(itcTest)
# plot(itcTest, axes=T)
# 
# #cristina on slack, check her starting values
# 
# 
# ## load RGB and DEM
# dsm <- (raster(x = "data/Orchard_IdahoTM_DEM_clip.tif") - 974.4)*150 #remove ground, scale to 250
# plot(dsm)
# rgb <- stack("C:/Users/peter.olsoy/Downloads/Orchard_Outputs/Orchard_Outputs/Flight3/Orchard_IdahoTM_RGB.tif")
# plot(rgb)
# 
# ## Clip and mask RGB to DEM
# rgb.sub <- crop(rgb, extent(dsm))
# plotRGB(rgb.sub)
# rgb.sub2 <- aggregate(rgb.sub, fact=2, fun=mean, na.rm=TRUE)
# plotRGB(rgb.sub2)
# rgb.sub3 <- resample(rgb.sub2, dsm, 'bilinear')
# plotRGB(rgb.sub3)
# rgb.sub3 <- crop(rgb.sub3, extent(dsm))
# rgb.sub3 <- mask(rgb.sub3, dsm)
# summary(rgb.sub3)
# 
# ## Make slope raster
# slope <- terrain(dsm, opt="slope", unit="degrees", neighbors=8)
# plot(slope)
# 
# ## Add DEM and slope to raster stack with RGB
# rgbds <- stack(c(dsm, slope, rgb.sub3[[2]]))
# plotRGB(rgbds)
# 
# writeRaster(rgbds, filename="C:/Users/peter.olsoy/Desktop/UAS_Data/Orchard/GEM3/Orchard_2019_DSG.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
# 
# rgbds.df <- as.data.frame(rgbds)
# 
# library(dplyr)
# library(ggplot2)
# rgbds.df %>% 
#   na.omit() %>%
#   ggplot(aes(x=Orchard_IdahoTM_RGB.1, y=Orchard_IdahoTM_RGB.2))+geom_point()
# 
# summary(rgbds.df)
# 
# rgbds.df$vvi <- ((1-abs((rgbds.df$Orchard_IdahoTM_RGB.1-30.0)/(rgbds.df$Orchard_IdahoTM_RGB.1+30.0)))*
#                    (1-abs((rgbds.df$Orchard_IdahoTM_RGB.2-50.0)/(rgbds.df$Orchard_IdahoTM_RGB.2+50.0)))*
#                    (1-abs((rgbds.df$Orchard_IdahoTM_RGB.3-0.0)/(rgbds.df$Orchard_IdahoTM_RGB.3+0.0))))
# summary(rgbds.df$vvi)
# 
# hist(rgbds.df$vvi, na.rm=TRUE)
#FUNCTIONS
#FUNCTIONS

splt=function(x,split,num) {
  
  x<-as.character(x)
  p1=strsplit(x,split=split,fixed=T)
  p2=unlist(lapply(p1,"[",num))
  return(p2)
}

library("rgdal")
library("rgeos")
library("sp")
library("raster")
library("ForestTools")

## load DEM
#dsm <- raster(x = "data/Orchard_IdahoTM_DEM_clip.tif") #2019 DEM data
dsm <- raster(x = "Data/orchard2015_dem_IDTM_clip.tif") #2015 DEM data

plot(dsm)

#kind of a decision point, but don't really want to spend too much
#computing time on this, since trend is more likely due to 
#shitty empire unmanned stuff and not likely to appear again
groundPoly <- readOGR(dsn="Data", layer="orchard_ground_poly")
plot(groundPoly, add=TRUE)

crop <- crop(dsm, groundPoly)
dsm.ground <- mask(crop, groundPoly)
plot(dsm.ground)

#the slope/trend is much worse for 2015, need to remove that before moving forward
dsm.ground.df <- as.data.frame(dsm.ground, na.rm=TRUE, xy=TRUE) #convert raster to dataframe
names(dsm.ground.df)[3] <- "z"
#trend <- lm(orchard2015_dem_IDTM_clip ~ x + y, data=dsm.df) #model XY trend
trend <- lm(z ~ x + y, data=dsm.ground.df) #model XY trend
summary(trend)

dsm.df <- as.data.frame(dsm, na.rm=TRUE, xy=TRUE)
names(dsm.df)[3] <- "z"

dsm.df$z2 <- dsm.df$z - predict.lm(trend, newdata=dsm.df)

#dsm.df[!(is.na(dsm.df[,3])),3] <- trend$residuals
dsm2 <- rasterFromXYZ(dsm.df[,c(1,2,4)], crs=crs(dsm))
plot(dsm2)

#allometric function - improve with data from Andrii
lin <- function(x){x * 0.54 - 0.0044} #function for relating max crown height to radius

###revised function: what is wrong here?
allo_fun<-function(H) {
  return(-2.178142e+01+1.718507e+00*H
         + -5.693651e-03*H^2 + 1.224483e-05*H^3)
}


#search for crowns or "treetops" throughout the image, this step takes the longest 

minHeight_ttops=0.6
minHeight_mcws=0.2
trim_tiny_polygons=0.25



minHeight_ttops=c(0.6)
minHeight_mcws=c(0.2)
trim_tiny_polygons=c(0.5,0.25)

tries=expand.grid(minHeight_ttops,minHeight_mcws,trim_tiny_polygons)

#load points in of plants
plantpts0<-readOGR(".","orchard_survival_IDTM_adjusted")

#should this be all plants for 2015?
plantpts<-plantpts0[which(plantpts0$Status_QC=="Living" | plantpts0$Status_QC=="Recent"),] #subset to living plants
#where does this living come from (field data vs. other data)?

canopy_optim=function(minHeight_ttops,minHeight_mcws,trim_tiny_polygons){

###DECISION POINT
###DECISION POINT
ttops <- vwf(CHM = dsm2, winFun = lin, minHeight = minHeight_ttops, maxWinDiameter = NULL) #minHeight is the parameter to pick out crowns, higher will make fewer crowns
###DECISION POINT
###DECISION POINT



#plot(dsm2 > 0.2, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
#plot(ttops, col = "blue", pch = 20, cex = 0.5, add = TRUE)

mean(ttops$height)

###DECISION POINT
###DECISION POINT
crowns <- mcws(treetops = ttops, CHM = dsm2, minHeight = minHeight_mcws, verbose = FALSE) #minHeight here will change which pixels are picked up as shrub, a higher number means smaller polygons around each crown/ttop
#plot(crowns, col = sample(rainbow(50), length(unique(crowns[])), replace = TRUE), legend = FALSE, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
###DECISION POINT
###DECISION POINT



#crownsPoly <- mcws(treetops = ttops, CHM = dsm2, format = "polygons", minHeight = 0.2, verbose = FALSE)
# there's an error with the packages "format="polygons"" algorithm, it doesn't correctly convert
# therefore, trying the rasterToPolygon function
crownsPoly <- rasterToPolygons(crowns, dissolve=TRUE) #note that some crowns are multi-part

#plot(dsm2, xlab = "", ylab = "", xaxt='n', yaxt = 'n')
#plot(crownsPoly, border = "blue", lwd = 0.5, add = TRUE)

d <- disaggregate(crownsPoly)

###DECISION POINT
###DECISION POINT
d2 <- d[area(d)>trim_tiny_polygons,] #used 0.25 before
###DECISION POINT
###DECISION POINT

# d2[["crownArea"]] <- gArea(d2, byid = TRUE) #calculate area of crowns, even if multi-part
# d2[["crownDiameter"]] <- sqrt(d2[["crownArea"]] / pi) * 2 #estimate diameter from area
# mean(d2$crownDiameter)

### End Peter's code ###

#crowns2019=readOGR(dsn=".","crowns2topoly") #segmented layer
crowns2015 <- d2
plot(crowns2015)

crowns2015$crown_ID<-c(1:nrow(crowns2015)) 

coverxy=gCentroid(crowns2015,byid=T) #centroid of crowns
crowns2015@data<-cbind(crowns2015@data,coverxy@coords)

#spatial overlay of points and crowns
crown_contain=over(crowns2015,plantpts,returnList=T)

names(crown_contain)=crowns2015$crown_ID #annoyingly, the first element starts at 0

crown_dat<-do.call("rbind",crown_contain) #converting it to dataframe

ori_crowns=rownames(crown_dat) #rownames are crown ids

row_crown_ID=splt(ori_crowns,".",1) #adding a column representing crown ID

crown_dat$row_crown_ID=as.numeric(row_crown_ID)

#####################################################
##ACCURACY METRICS FOR SEGMENTATION:
##ACCURACY METRICS FOR SEGMENTATION:
#####################################################
#####################################################
#####################################################

#living plants that should have been detected
#living plants that should have been detected
living_plant_detection=1-length(which(plantpts$Tag %in% crown_dat$Tag ==F))/length(plantpts$Tag)

#how many crowns there should have been (tricky, since there are recruits)
total_crowns=nrow(crowns2015)-nrow(plantpts)

#how many crowns are shared

crown_splt=splt(ori_crowns,".",2) #which ones share a crown?
crown_sharers=which(is.na(crown_splt)==F) #points in the same crown
good_crowns=crown_dat[-crown_sharers,]
shared_crowns=crown_dat[crown_sharers,]

shared_crowns=nrow(shared_crowns)/(nrow(good_crowns)+nrow(good_crowns))

return(c(living_plant_detection,total_crowns,shared_crowns))
}


results_mat<-matrix(NA,nrow=nrow(tries),ncol=3)

for(i in 1:nrow(tries)){
  
  results_mat[i,]<-canopy_optim(tries[i,1],tries[i,2],tries[i,3])
  
}


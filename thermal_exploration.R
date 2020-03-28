library("rgdal")
library("rgeos")
library("sp")
library("raster")
library("dismo")

crowns2019=readOGR(dsn=".","orchard_crowns_2019")
plot(crowns2019)

crowns2019$ID<-as.character(crowns2019$ID)
crowns2019$Crown_area<-gArea(crowns2019,byid=T)
coverxy=gCentroid(crowns2019,byid=T)

crowns2019@data<-cbind(crowns2019@data,coverxy@coords)

thermal=raster("Orchard_IdahoTM_Thermal.tif")


crown_clip<-mask(thermal,crowns2019)

crown_clip<-trim(crown_clip)

thermal_clip=projectRaster(from=crown_clip,to=crs(crowns2019),res=res(crown_clip))

plot(crown_clip)

plantpts<-readOGR(".","orchard_points_utm")

#spatial merge
contain_points=over(plantpts,crowns2019,returnList=T)

idd<-function(x) {
  return(x$ID)}

contain_ID=lapply(contain_points,idd)

contain_ID<-contain_ID[lapply(contain_ID,length)>0]

pts_merge<-data.frame(crown_ID=unlist(contain_ID),point_ID=names(contain_ID))
pts_merge$crown_ID<-as.character(pts_merge$crown_ID)
pts_merge$point_ID<-as.character(pts_merge$point_ID)


dim(pts_merge) #298 merged
length(unique(pts_merge$crown_ID)) #229 unique crown points
#this means 69 points are in the same crown
plantpts$ID<-c(1:nrow(plantpts@data))
plantpts_w_crowns=merge(plantpts,pts_merge,by.x="ID",by.y="point_ID")

plantpts_w_crowns<-plantpts_w_crowns[-which(is.na(plantpts_w_crowns$crown_ID)),]

cidU=unique(plantpts_w_crowns$crown_ID)
cidL=rep(NA,times=length(cidU))

for(i in 1:length(cidU)){
  cidL[i]<-length(which(plantpts_w_crowns$crown_ID==cidU[i]))
}

crown_duplicates=cidU[which(cidL>1)]

#step 1: set voronoi polygons around points with buffer to ensure end gets in
#step 2: clip voronoi polygons to crown segmentation

good_to_use=plantpts_w_crowns[plantpts_w_crowns$crown_ID %in% crown_duplicates==F,]


needs_voronoi=plantpts_w_crowns[plantpts_w_crowns[plantpts_w_crowns$crown_ID %in% crown_duplicates,],]

vpolygons=voronoi(needs_voronoi)
#spatial merge
contain_points_v=over(vpolygons,needs_voronoi,returnList=T)

contain_dat=do.call("rbind",contain_points_v)

vpolygons@data<-contain_dat



single_parts=readOGR("single_parts_canopy_voronoi",dsn=".")


v_points=over(single_parts,needs_voronoi,returnList=T)

idd<-function(x) {
  return(x$ID)}

v_ID=lapply(v_points,idd)

v_polygons_single<-single_parts[lapply(v_ID,length)>0,]

keepers=over(crowns2019,good_to_use,returnList=T)
k_ID=lapply(keepers,idd)
crowns_nonvoronoi=crowns2019[lapply(k_ID,length)>0,]

v_poly_ID<-v_polygons_single
v_poly_ID@data<-data.frame(crown_ID=v_polygons_single@data$crown_ID)

c_poly_ID<-crowns_nonvoronoi
c_poly_ID@data<-data.frame(crown_ID=crowns_nonvoronoi@data$ID)

final_poly<-rbind(v_poly_ID,c_poly_ID)
#writeOGR(final_poly,layer="final_poly",dsn=".",driver="ESRI Shapefile")
# clip_shp = function(small_shp, large_shp){
#   # make sure both have the same proj
#   large_shp = spTransform(large_shp, CRSobj = CRS(proj4string(small_shp)))
#   cat("About to get the intersections, will take a while...", "\n")
#   clipped_shp = rgeos::gIntersection(small_shp, large_shp, byid = T, drop_lower_td = T)
#   cat("Intersection done", "\n")
#   x = as.character(row.names(clipped_shp))
#   # these are the data to keep, can be duplicated
#   keep = gsub(pattern = "^[0-9]{1,2} (.*)$", replacement = "\\1", x)
#   large_shp_data = as.data.frame(large_shp@data[keep,])
#   row.names(clipped_shp) = row.names(large_shp_data)
#   clipped_shp = spChFIDs(clipped_shp, row.names(large_shp_data))
#   # combine and make SpatialPolygonsDataFrame back
#   clipped_shp = SpatialPolygonsDataFrame(clipped_shp, large_shp_data)
#   clipped_shp
# }
# clip_v=clip_shp(crowns2019,vpolygons)
# clip_v_final=clip_v[-which(clip_v$crown_ID%in% crown_duplicates==F),]
# writeOGR(clip_v_final,"clip_v_test",dsn=".",driver="ESRI Shapefile",overwrite=T)


fpId=over(final_poly, plantpts,returnList=T)

fdat=do.call("rbind",fpId)
repeaters=which(unlist(lapply(fpId,nrow))>1) #how do some of these not add up!
#this is because single to multiparts didn't acknoqledge ones next to each other

projectRaster(from=thermal, crs=crs(final_poly))

rID=rasterize(final_poly,y=thermal,field="crown_ID")

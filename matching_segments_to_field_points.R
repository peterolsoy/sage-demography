#FUNCTIONS
#FUNCTIONS

caughlin_polygons<-function (x, bounding.polygon) 
{
  poly_extent=extent(bounding.polygon)
  
  
  
  if (!inherits(x, "SpatialPoints")) {
    stop("Must pass a SpatialPoints* object to voronoi.polygons.")
  }
  crds = coordinates(x)
  z = deldir::deldir(crds[, 1], crds[, 2],
                     rw=c(poly_extent@xmin,poly_extent@xmax,poly_extent@ymin,poly_extent@ymax))
  w = deldir::tile.list(z)
  polys = vector(mode = "list", length = length(w))
  for (i in seq(along = polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1, ])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID = as.character(i))
  }
  SP = SpatialPolygons(polys, proj4string = CRS(proj4string(x)))
  voronoi = SpatialPolygonsDataFrame(SP, data = data.frame(x = crds[, 
                                                                    1], y = crds[, 2], area = sapply(slot(SP, "polygons"), 
                                                                                                     slot, "area"), row.names = sapply(slot(SP, "polygons"), 
                                                                                                                                       slot, "ID")))
  if (!missing(bounding.polygon)) {
    bounding.polygon <- gUnion(bounding.polygon, bounding.polygon)
    voronoi.clipped <- gIntersection(voronoi, bounding.polygon, 
                                     byid = TRUE, id = row.names(voronoi))
    df <- data.frame(voronoi)
    df$area <- sapply(slot(voronoi.clipped, "polygons"), 
                      slot, "area")
    voronoi <- SpatialPolygonsDataFrame(voronoi.clipped, 
                                        df)
  }
  voronoi
}


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
library("SDraw")

crowns2019=readOGR(dsn=".","crowns2topoly") #segmented layer
plot(crowns2019)

crowns2019$crown_ID<-c(1:nrow(crowns2019)) 

coverxy=gCentroid(crowns2019,byid=T) #centroid of crowns
crowns2019@data<-cbind(crowns2019@data,coverxy@coords)

#load points in of plants
plantpts0<-readOGR(".","orchard_survival_IDTM_adjusted")

plantpts<-plantpts0[which(plantpts0$Status_QC=="Living"),] #subset to living plants
#where does this living come from (field data vs. other data)?

#spatial overlay of points and crowns
crown_contain=over(crowns2019,plantpts,returnList=T)

names(crown_contain)=crowns2019$crown_ID #annoyingly, the first element starts at 0

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
length(which(plantpts$Tag %in% crown_dat$Tag ==F))/length(plantpts$Tag)

#how many crowns there should have been (tricky, since there are recruits)
nrow(crowns2019) 
nrow(plantpts)



crown_splt=splt(ori_crowns,".",2) #which ones share a crown?

crown_sharers=which(is.na(crown_splt)==F) #points in the same crown

good_crowns=crown_dat[-crown_sharers,]

shared_crowns=crown_dat[crown_sharers,]

shared_points0=plantpts[plantpts$Tag %in% shared_crowns$Tag==T,]
plot(shared_points0)

sc_merge=data.frame(Tag=shared_crowns$Tag,row_crown_ID=shared_crowns$row_crown_ID)

shared_points<-merge(shared_points0,sc_merge,by.x="Tag",by.y="Tag")



crown1=crowns2019[which(crowns2019$crown_ID==shared_crowns$row_crown_ID[2]),]

pt<-function(xx){
crown1=crowns2019[which(crowns2019$crown_ID==shared_crowns$row_crown_ID[xx]),]
plot(crown1)}

uni_crowns=unique(shared_crowns$row_crown_ID)

crown_list=vector("list",length=length(uni_crowns))

for(j in 1:length(crown_list)){
#loop through with voronoi.polygons
crown_list[[j]]=caughlin_polygons(x=shared_points[which(shared_points$row_crown_ID==uni_crowns[j]),],
                 bounding.polygon=crowns2019[which(crowns2019$crown_ID==uni_crowns[j]),])

}


sp_shared<-do.call("rbind",crown_list)

spmerge=over(sp_shared,shared_points,return_list=T)

sp_shared@data<-cbind(sp_shared@data,spmerge)

#ah fuck this almost works
#writeOGR(sp_all,layer="sp_all_test",dsn=".",driver="ESRI Shapefile",overwrite=T)

crowns2019.1<-crowns2019[which(crowns2019$crown_ID %in% good_crowns$row_crown_ID==T),]

sp_single<-merge(crowns2019.1,good_crowns,by.x="crown_ID",by.y="row_crown_ID")

colnames(sp_shared@data)[which(colnames(sp_shared@data)=="area")]<-"AREA_M2"
colnames(sp_shared@data)[which(colnames(sp_shared@data)=="row_crown_ID")]<-"crown_ID"
sp_single@data<-sp_single@data[,c(1,4:27)]

sp_single@data<-sp_single@data[,order(colnames(sp_single@data))]
sp_shared@data<-sp_shared@data[,order(colnames(sp_shared@data))]

final_full=rbind(sp_single,sp_shared)

writeOGR(final_full,dsn=".",layer="final_merged_segmented",driver="ESRI Shapefile",
         overwrite=T)

library("rgdal")
library("rgeos")
library("sp")
library("raster")
library("SDraw")


final_full=readOGR(dsn=".",layer="final_merged_segmented")
length(unique(final_full$crown_ID))
length(final_full$crown_ID) #why are there some repeat crown_IDs?
length(final_full$Tag) #201 here. there's something weird going on

final_full$Tag<-as.numeric(as.character(final_full$Tag))

thermal2019<-raster("Orchard_IdahoTM_Thermal.tif")

crown_IDs=rasterize(x=final_full,y=thermal2019,field="crown_ID")

crown_IDs<-trim(crown_IDs)

plot(crown_IDs)

thermal_crop<-mask(x=thermal2019,mask=final_full)

thermal_crop<-trim(thermal_crop)

thermal_IDs<-stack(thermal_crop,crown_IDs)

thermal_dat<-data.frame(rasterToPoints(thermal_IDs))
#518246 points

thermal_merge<-merge(thermal_dat,final_full@data,by.x="layer",
                     by.y="crown_ID")

#leaf_level analysis
plot(thermal_merge$Orchard_IdahoTM_Thermal~thermal_merge$subspp)

library("ggplot2")

p <- ggplot(thermal_merge, aes(x=subspp, y=Orchard_IdahoTM_Thermal)) + 
  geom_violin(fill=c("darkred")) 

p + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")


library(lme4)
m1<-lmer(Orchard_IdahoTM_Thermal~subspp+(1|Tag),data=thermal_merge)

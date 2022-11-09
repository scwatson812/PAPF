##########################
##### Load Libraries #####
##########################
library(ape)
library(maptools)
library(parallel)
library(spatstat)
library(rgeos)
library(rgdal)
library(raster)
library(sf)
library(lwgeom)
library(spdep)

##############################
##### Construct the Grid #####
##############################
grid.raster = raster(ncol = 20,nrow = 20, xmn = 0, xmx = 0.1, ymn = 0, ymx = 0.1)
fsa.sh = rasterToPolygons(grid.raster)
fsa.sh <- spTransform(fsa.sh,CRS("+init=epsg:32610"))
fsa.sh$ID = 1
plot(fsa.sh)
nb.list = poly2nb(fsa.sh)
W.mat = nb2mat(nb.list,style = 'B')

######################
##### Parameters #####
######################
min.rad = min(dist(coordinates(fsa.sh)))+1
max.rad = 5005
n.rad = 10
radius.vec = seq(min.rad,max.rad,length.out = n.rad)

############################################
##### Data Generation: CSR probability #####
############################################
n.selected.points = 40
ind <- sample(x = seq(1,dim(fsa.sh)[1]),n.selected.points,replace = FALSE)
fsa.sh$selected<- 0
fsa.sh$selected[ind] = 1

#########################
##### Plot the Data #####
#########################
plot(fsa.sh)
plot(fsa.sh[ind,],col = 'grey',add = TRUE)
plot(fsa.sh[fsa.sh$selected ==1,],col = 'blue',add = TRUE)


### Make the shapefiles
fsa.selected = fsa.sh[fsa.sh$selected == 1,]
fsa.selected.sf = st_as_sf(fsa.selected)
fsa.selected.sf = st_buffer(fsa.selected.sf,0)
ca.sh = unionSpatialPolygons(fsa.sh,fsa.sh$ID)
ca.sf = st_as_sf(ca.sh)
ca.sf = st_buffer(ca.sf,0)
plot(ca.sh)
plot(fsa.sh,add = TRUE)
plot(fsa.selected.sf$geometry,col = 'red',add = TRUE)

###############################################
##### Compute the Observed Test Statistic #####
###############################################

###Obtain the Centriods of the Observed Units
centers.obs.points =  coordinates(fsa.selected)
centers.obs.points = SpatialPoints(cbind(centers.obs.points[,1],centers.obs.points[,2]),proj4string = CRS(proj4string(ca.sh)))
centers.obs.points = st_as_sf(centers.obs.points)
n.obs.points = dim(centers.obs.points)[1]

###Create a list of length n.rad whose ith element is a vector of circles of radius radius.vec[i] centered as each of the observed areal unit centers
obs.circles.list = list()
ct = 1
for(rad in 1:n.rad){
  radius = radius.vec[rad]
  circle.obs = st_buffer(centers.obs.points,dist = radius)
  for(i in 1:dim(circle.obs)[1]){
    obs.circles.list[[ct]] = circle.obs[i,]
    ct = ct + 1
  }
}

obs.intersection.calc<-function(circle){
  obs.intersection = st_intersection(circle,fsa.selected.sf)
  final.intersection = st_intersection(ca.sf,obs.intersection)
  area = st_area(final.intersection)
  if(length(area)==0){
    area = 0
  }
  return(sum(area))
}

circle.area.calc<-function(circle){
  area = st_area(st_intersection(ca.sf,circle))
  if(length(area)==0){
    area = 0
  }
  return(sum(area))
}

observed.area.vec = unlist(lapply(obs.circles.list,obs.intersection.calc))
circle.area.vec = unlist(lapply(obs.circles.list,circle.area.calc))
observed.area.mat = matrix(observed.area.vec/circle.area.vec,n.rad,n.obs.points,byrow = TRUE)

write.table(observed.area.mat,"ObservedTestStat.txt",row.names = FALSE, col.names = FALSE)
write.table(sum(area(fsa.selected))/sum(area(fsa.sh)),"TotalAreaProp.txt",row.names = FALSE, col.names = FALSE)



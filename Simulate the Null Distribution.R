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
nb.list = poly2nb(fsa.sh)
W.mat = nb2mat(nb.list,style = 'B')

######################
##### Parameters #####
######################
###Ensure that the choice of radius.vec and n.obs.points matches the observed data
n.mc.reps = 5
min.rad = min(dist(coordinates(fsa.sh)))+1
max.rad = 5005
n.rad = 2
n.obs.points = 4
radius.vec = seq(min.rad,max.rad,length.out = n.rad)

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

##########################################################################################################################################################
##### Estimate the Distribution of the test statistic via Monte Carlo Simulation #########################################################################
##########################################################################################################################################################

intersection.calc<-function(circle){
  obs.intersection = st_intersection(circle,fsa.selected.mc.sf)
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


mc.circle.area.mat= matrix(0,n.mc.reps,n.rad*n.obs.points)
mc.prop.area.mat= matrix(0,n.mc.reps,1)
st.mc = Sys.time()
for(g in 1:n.mc.reps){
  fsa.sh.mc = fsa.sh
  select.ind = sample(size = n.obs.points,x = seq(1,dim(fsa.sh)[1]),replace = FALSE)
  fsa.sh.mc$selected <- 0
  fsa.sh.mc$selected[select.ind] = 1
  fsa.sh.mc.selected = fsa.sh.mc[fsa.sh.mc$selected == 1, ]
  fsa.selected.mc.sf = st_as_sf(fsa.sh.mc.selected)
  fsa.selected.mc.sf = st_buffer(fsa.selected.mc.sf,0)
  ###Compute Test Statistic for the given MC Simulated Data
  
  ###Create the Circles at the Observed Simulated Data points
  
  ###Obtain the Centriods of the Observed Units
  centers.obs.points =  coordinates(fsa.sh.mc.selected)
  centers.obs.points = SpatialPoints(cbind(centers.obs.points[,1],centers.obs.points[,2]),proj4string = CRS(proj4string(ca.sh)))
  centers.obs.points = st_as_sf(centers.obs.points)
  
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
  
  sim.obs.area.vec = unlist(lapply(obs.circles.list, intersection.calc))
  sim.circle.area.vec = unlist(lapply(obs.circles.list,circle.area.calc))
  mc.circle.area.mat[g,] = sim.obs.area.vec/sim.circle.area.vec
  mc.prop.area.mat[g,] = sum(area(fsa.sh.mc.selected))/sum(area(fsa.sh.mc))
  print(g)
}
sp.mc = Sys.time()

mc.circle.area.array = array(0,c(n.mc.reps,n.rad,n.obs.points))
for(g in 1:n.mc.reps){
  for(r in 1:n.rad){
    mc.circle.area.array[g,r,] = mc.circle.area.mat[g,((r-1)*n.obs.points+1):(r*n.obs.points)]
  }
}

#########################################
##### Write Out Simulations Results #####
#########################################
write.table(mc.circle.area.mat,"SimulatedProportion.txt", row.names = FALSE, col.names = FALSE)
write.table(mc.prop.area.mat,"SimulatedTotalProportion.txt", row.names = FALSE, col.names = FALSE)


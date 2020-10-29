
# This script contains all of the data exploration and cleaning procedures that were done.
## Each section is labelled. All cleaned data files are saved into the data folder within the project so they can be read into subsequent scripts where necessary

#these are the generic package needed
library(data.table)
library(sp)
library(raster)

############################
### VBA Data Exploration ###
############################
VBA_brushtail.raw <- read.csv(file = "data/VBA/brushtail.csv", stringsAsFactors = FALSE)
brushtail.dt <- as.data.table(VBA_brushtail.raw[-nrow(VBA_brushtail.raw), ])
names(VBA_brushtail.raw)
is.numeric(VBA_brushtail.raw$Survey.Start.Date)
class(VBA_brushtail.raw$Survey.Start.Date)
sum(is.na(brushtail.dt$Latitude.GDA94))  #0
sum(is.na(brushtail.dt$Longitude.GDA94)) #0 
summary(brushtail.dt$Reliability)       #14186 Acceptable
summary(brushtail.dt$Accuracy) #15 NAs, median = 50
sum(brushtail.dt[Accuracy < 200, Reliability == "Acceptable"]) #12228
brushtail.dt$Year <- year(as.Date(brushtail.dt$Survey.Start.Date, format = "%d/%m/%Y"))
hist(brushtail.dt$Year)


VBA_egk.raw <- read.csv(file = "data/VBA/egk.csv", stringsAsFactors = FALSE)
egk.dt <- as.data.table(VBA_egk.raw[-nrow(VBA_egk.raw), ])
plot(egk.dt$Survey.Start.Date)
summary(egk.dt$Reliability) #9398 Acceptable
summary(egk.dt$Accuracy) #3 NAs
sum(egk.dt[Accuracy < 200, Reliability == "Acceptable"]) #6395
egk.dt$Year <-year(as.Date(egk.dt$Survey.Start.Date, format = "%d/%m/%Y"))
summary(egk.dt$Year > 2000)
hist(egk.dt$Year)

VBA_koala.raw <- read.csv(file = "data/VBA/koala.csv", stringsAsFactors = FALSE)
koala.dt <- as.data.table(VBA_koala.raw[-nrow(VBA_koala.raw), ])
sum(koala.dt[Accuracy < 200, Reliability == "Acceptable"]) #5579
koala.dt$Year <-year(as.Date(koala.dt$Survey.Start.Date, format = "%d/%m/%Y"))


VBA_ringtail.raw <- read.csv(file = "data/VBA/ringtail.csv", stringsAsFactors = FALSE)
ringtail.dt <- as.data.table(VBA_ringtail.raw[-nrow(VBA_ringtail.raw)])
sum(ringtail.dt[Accuracy <200, Reliability == "Acceptable"]) #6258
ringtail.dt$Year <-year(as.Date(ringtail.dt$Survey.Start.Date, format = "%d/%m/%Y"))


VBA_wallaby.raw <- read.csv(file = "data/VBA/wallaby.csv", stringsAsFactors = FALSE)
all.wallaby.dt <- as.data.table(VBA_wallaby.raw[-nrow(VBA_wallaby.raw)])
sum(all.wallaby.dt[Accuracy < 200, Reliability == "Acceptable"]) #15321 ??
unique(all.wallaby.dt$Common.Name) #TOO MANY WALLABIES
unique(all.wallaby.dt$Scientific.Name)
wallaby.dt <- as.data.table(all.wallaby.dt[Scientific.Name == "Wallabia bicolor"])
sum(wallaby.dt[Accuracy <200, Reliability == "Acceptable"]) #13856... 
wallaby.dt$Year <- year(as.Date(wallaby.dt$Survey.Start.Date, format = "%d/%m/%Y"))

VBA_wombat.raw <- read.csv(file = "data/VBA/wombat.csv", stringsAsFactors = FALSE)
wombat.dt <- as.data.table(VBA_wombat.raw[-nrow(VBA_wombat.raw)])
sum(wombat.dt[Accuracy <200, Reliability == "Acceptable"]) # 5435
wombat.dt$Year <-year(as.Date(wombat.dt$Survey.Start.Date, format = "%d/%m/%Y"))

######################################################
### Making a master csv to look at landuse in QGIS ###
######################################################

VBA_all <- rbind(brushtail.dt, egk.dt, koala.dt, ringtail.dt, wallaby.dt, wombat.dt)
VBA_all <- VBA_all[VBA_all$Scientific.Name!= "",]
unique(VBA_all$Survey.method)
sum(VBA_all$Survey.method == "Translocation" | VBA_all$Survey.method == "Koala Translocation")
sum(VBA_all$Year >= 1989)
VBA_species <- subset(VBA_all,Year >= 1989)
sum(VBA_all$Reliability == "Acceptable")
write.table(VBA_species, sep= ',', "data/VBA/VBA_species.csv", col.names = T, row.names = F)

###################################
### Final species data cleaning ###
###################################

#more data exploration
sum(is.na(VBA_species$Total.Count)) #24719

#Subset koala translocations CHECK WITH NAT (translocation reliability)

# investigating the koala translocations
koala.trans <- subset(VBA_species, Survey.method == "Koala Translocation")
unique(koala.trans$Extra.Info)
table(koala.trans$Observer)
nrow(koala.trans)
sum(koala.trans$Extra.Info == "Released/Introduced")/length(VBA_species$Scientific.Name == "Phascolarctos cinereus") #41% of all koala data

hist(VBA_species$Total.Count[VBA_species$Total.Count <= 100])

#look at species counts, NAs, 0, 1, >0, huge
count.greater.100 <- subset(VBA_species, Total.Count >100)
summary(count.greater.100$Total.Count)
sum(count.greater.100)
sum(VBA_species$Total.Count == "", na.rm = F)
## cut down to 10


#look again at accuracy 

hist(VBA_species$Accuracy)
summary(VBA_species$Accuracy)
Acc.5000 <- subset(VBA_species,Accuracy <= 5000)
hist(Acc.5000$Accuracy)
hist(Acc.5000$Accuracy, xlim= c(0,1000))
plot(density(Acc.5000$Accuracy), xlim= c(0,1000))
summary(Acc.5000$Accuracy)
Acc.200 <- subset(VBA_species,Accuracy <= 200)

#Taking out species counts that are 0
count <- subset(Acc.200, Total.Count >= 1 & Total.Count <= 15 | is.na(Total.Count))

#Remove the translocation data (for now)
VBA_final <- subset(count, !Survey.method == "Koala Translocation" | is.na(Survey.method))


#only taking the necessary columns
species.dat <- VBA_final[, c("Scientific.Name","Longitude.GDA94", "Latitude.GDA94")]

#Transforming the projection from lat/long to utm
vba.crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
crs.data <- "+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

species.spdf <- SpatialPointsDataFrame(data= as.data.frame(species.dat[,1]), coords= species.dat[,2:3], proj4string = CRS(vba.crs)) 
sp.tf <- spTransform(species.spdf, CRSobj = crs.data)

sp.df <- as.data.frame(sp.tf)
names(sp.df)<- c("Scientific.Name","x", "y")
head(sp.df)


# removing points with a distance less than 100m

min.dist <- function(data, sp.name, dist){
  sp.data <- subset(data, Scientific.Name == sp.name)
  dmat <- spDists(as.matrix(sp.data[,-1]))
  dmat[upper.tri(dmat)== FALSE & dmat <= dist & dmat != 0] <- NA
  dmat.nona <- unlist(apply(dmat, 2, function(x) which(is.na(x))))
  sp.sub <- subset(sp.data[-dmat.nona,])
  return(sp.sub)
}

sci.names <- c("Phascolarctos cinereus", "Trichosurus vulpecula", "Macropus giganteus", 
                      "Pseudocheirus peregrinus", "Wallabia bicolor", "Vombatus ursinus")

species_clean <- NULL
for (s in 1:6) {
  clean.sp <- min.dist(sp.df, sci.names[s], 100)
  
  species_clean <- rbind(species_clean, clean.sp)
}

# removing duplicate points
species_clean <- unique(species_clean, by= c("Scientific.Name","x", "y")) 

### This is the final clean species data ###

write.csv(species_clean, file= "data/species_clean.csv")


# Make appendix table for all o these things
data.cleaning <- data.frame("Raw" =table(VBA_all$Scientific.Name), "Date" = table(VBA_species$Scientific.Name), "Accuracy" = table(Acc.200$Scientific.Name), "Count" = table(count$Scientific.Name), "Minimum dist"= table(species_clean$Scientific.Name))
data.cleaning$name <- c("Eastern Grey Kangaroo", "Koala", "Common Ringtail Possum", "Common Brushtail Possum", "Common Wombat", "Swamp Wallaby")
data.cleaning <- data.cleaning[,c(1,11,2,4,6,8,10)]
View(data.cleaning)
names(data.cleaning) <- c("Scientific names", "Common Names", "Raw", "Date > 1989", "Accuracy < 200", "Count <= 15", "100m between points")
## note that the data cleaning table is not finished until the number of points at each gridded layer is included
### see the "preparing the covariate list" section

##################################################
### Preparing the SDM environmental covariates ###
##################################################

# First thing to do is too understand the correlation of variables
## reading in the data
VIC_1000m_crop <-raster("data/VIC_shape/VIC_GDA9455_GRID_STATE_1000.tif")
var.path_1000 <-list.files(path= "data/grids_1000", full.names= T, pattern= "*.tif")
var.stack_1000 <- stack(var.path_1000) 
var.crop_1000 <- mask(var.stack_1000, VIC_1000m_crop)

# looking at the correlation for the 1000m
cor.crop_1000 <- layerStats(var.crop_1000, 'pearson', na.rm=T)
corr.1000.crop <- cor.crop_1000$'pearson correlation coefficient'
write.csv(corr.1000.crop, "data/grids_1000/correlation_matrix_1000_cropped.csv")

## loading and looking at the 500m correlation
VIC_500m_crop <-raster("data/VIC_shape/VIC_GDA9455_GRID_STATE_500.tif")
var.path_500 <-list.files(path= "data/grids_500", full.names= T, pattern= "*.tif")
var.stack_500 <- stack(var.path_500) 
var.crop_500 <- mask(var.stack_500, VIC_500m_crop)

cor.crop_500 <- layerStats(var.crop_500, 'pearson', na.rm=T)
corr.500.crop <- cor.crop_500$'pearson correlation coefficient'
write.csv(corr.500.crop, "data/grids_500/correlation_matrix_500_cropped.csv")

### again for the 250m
VIC_250m_crop <-raster("data/VIC_shape/VIC_GDA9455_GRID_STATE_250.tif")
var.path_250 <-list.files(path= "data/grids_250", full.names= T, pattern= "*.tif")
var.stack_250 <- stack(var.path_250) 
var.crop_250 <- mask(var.stack_250, VIC_250m_crop)

cor.crop_250 <- layerStats(var.crop_250, 'pearson', na.rm=T)
corr.250.crop <- cor.crop_250$'pearson correlation coefficient'
write.csv(corr.250.crop, "data/grids_250/correlation_matrix_250_cropped.csv")

#retaining the non-correlated variables and creating the variable stacks
var.final <- subset(var.crop_1000, c( "GRASSPROP", "GREENNESS", "HYDRODIST", "LIGHT", "PRECAN", "ROADDIST", "SLOPE", "SOILCLAY", "SOILDEPTH","SOILSAND", "SOILSILT", "TEMPMINCP", "TEMPMNDIRNG", "TEMPSEAS", "TOWNDIST", "TREEDENS"))
var.final_1000$TREEDENS[var.final_1000[["TREEDENS"]] > 100] <- NA
writeRaster(var.final_1000, "data/var.final_1000.grd", format= "raster")

var.final <- subset(var.crop_500, c( "GRASSPROP", "GREENNESS", "HYDRODIST", "LIGHT", "PRECAN", "ROADDIST", "SLOPE", "SOILCLAY", "SOILDEPTH","SOILSAND", "SOILSILT", "TEMPMINCP", "TEMPMNDIRNG", "TEMPSEAS", "TOWNDIST", "TREEDENS"))
var.final_500$TREEDENS[var.final_500[["TREEDENS"]] > 100] <- NA
writeRaster(var.final_500, "data/var.final_500.grd", format= "raster")

var.final <- subset(var.crop_250, c( "GRASSPROP", "GREENNESS", "HYDRODIST", "LIGHT", "PRECAN", "ROADDIST", "SLOPE", "SOILCLAY", "SOILDEPTH","SOILSAND", "SOILSILT", "TEMPMINCP", "TEMPMNDIRNG", "TEMPSEAS", "TOWNDIST", "TREEDENS"))
var.final_250$TREEDENS[var.final_250[["TREEDENS"]] > 100] <- NA
writeRaster(var.final_250, "data/var.final_250.grd", format= "raster")

# This will read the individual covariate stacks back in at one resolution if needed 
var.final_250 <- stack("data/var.final_250.grd") 
var.final_500 <- stack("data/var.final_500.grd") 
var.final_1000 <- stack("data/var.final_1000.grd") 

# To run the loops, I am creating a list with the three stacks
## for consistency, the order of resolutions will always go 250-500-1000, 
var.res <- list(var.final_250, var.final_500, var.final_1000)
res.names <- c("250", "500", "1000")
names(var.res) <- res.names
saveRDS(var.res, file= "data/var.res.RData")

# I need a seperate liste of covariates with the bias layers (distance to town and roads) set to their mean values for making predictions
eval.var.res <- list()

for(r in 1:3){
  idx <- which(!is.na(var.res[[r]][[1]][])) #This gves the IDs for the raster cells with data (from the first stack)
  var.temp <- var.res[[r]]
  var.temp[["ROADDIST"]][idx] <- cellStats(var.temp[["ROADDIST"]], mean)
  var.temp[["TOWNDIST"]][idx] <- cellStats(var.temp[["TOWNDIST"]], mean)
  eval.var.res[[r]] <- var.temp
}
names(eval.var.res) <- res.names
saveRDS(eval.var.res, file= "data/eval.var.res.RData")



##### Finally I am adding the calculating the number of points removed for each species 
###### when gridded at the different environmental resolutions
var.res <- readRDS("data/var.res.RData")

thin.dat <- matrix(nrow= 6, ncol= 3)
res.names <- c("250", "500", "1000")
colnames(thin.dat) <- c("gridded at 250m", "gridded at 500m", "gridded at 1000m")
sp.order <- data.cleaning$`Scientific names`
for (r in 1:3) {
  for (s in 1:6) {
    sp.subset <- subset(species.clean, Scientific.Name == sp.order[s])
    
    thin.dat[s,r] <- nrow(gridSample(sp.subset[,c("x", "y")], var.res[[r]], n=1))
    
  }
  
}
data.cleaning <- cbind(data.cleaning, thin.dat)

write.csv(data.cleaning, file = "data/data.cleaning.csv")

#############################################################
### Creating background points with 100m minimum distance ###
#############################################################

VIC_mask_100 <- disaggregate(VIC_1000m_crop, fact= 10)
centroids_100m <- xyFromCell(VIC_mask_100, which(VIC_mask_100[]==1), spatial= TRUE)
bg_randomised <- sample(centroids_100m, size = 50000)
bg_df <- as.data.frame(bg_randomised)

write.csv(bg_df, file= "data/bg_df.csv")

##################################################
### Investigating the amount of collision data ###
##################################################
library(RPostgreSQL)
library(dplyr)
library(sf)
library(rpostgis)

#naming vectors 
res.names <- c("250", "500", "1000") 
sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom")
species <- c("Koala","Brushtail_Possum", "Kangaroo", "Ringtail_Possum", "Wallaby", "Wombat") #this is how the species are listed on the spatial server
species.plot <- c("Koala","Brushtail", "Kangaroo", "Ringtail", "Wallaby", "Wombat")


## This specifies the connection to the spatial server where the road and collision information is stored.
drv <- DBI::dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- DBI::dbConnect(drv,
                      dbname = "qaeco_spatial",
                      user = "qaeco",
                      password = "Qpostgres15",
                      host = "boab.qaeco.com",
                      port = "5432")  #Connection to database server on Boab

#collecting the statewide road IDs
roads <- data.table::as.data.table(dbGetQuery(con,paste0("
  SELECT a.uid, a.length, a.x, a.y, b.tvol, c.tspd
  FROM
  (SELECT uid, ST_Length(geom)/1000 AS length, ST_X(ST_LineInterpolatePoint(geom, 0.5)) AS x, ST_Y(ST_LineInterpolatePoint(geom, 0.5)) AS y 
    FROM gis_victoria.vic_gda9455_roads_state_orig_500) AS a, 
    gis_victoria.vic_nogeom_roads_volpreds_500 AS b, 
    gis_victoria.vic_nogeom_roads_speedpreds_500 AS c
  WHERE
  a.uid = b.uid
  AND
  a.uid = c.uid
    ")))
write.csv(roads, file= "output/collision/roadID.csv")

#collecting the regional subset of road IDs
roads_CS <- data.table::as.data.table(dbGetQuery(con,paste0("
  SELECT a.uid, a.length, a.x, a.y, b.tvol, c.tspd
  FROM
  (SELECT geom FROM gis_victoria.vic_gda9455_admin_sa2 
            WHERE 
            (sa3_name16 = 'Upper Goulburn Valley' 
            OR sa3_name16 = 'Heathcote - Castlemaine - Kyneton' 
            OR sa3_name16 = 'Macedon Ranges' 
            OR sa2_name16 = 'Daylesford') 
            AND 
            (sa2_name16 !=  'Bendigo Region - South' 
            AND sa2_name16 != 'Castlemaine Region' 
            AND sa2_name16 != 'Mansfield (Vic.)' 
            AND sa2_name16 != 'Upper Yarra Valley' 
            AND sa2_name16 != 'Castlemaine')) AS z,
  (SELECT uid, ST_Length(geom)/1000 AS length, ST_X(ST_LineInterpolatePoint(geom, 0.5)) AS x, ST_Y(ST_LineInterpolatePoint(geom, 0.5)) AS y, geom 
    FROM gis_victoria.vic_gda9455_roads_state_orig_500) AS a, 
    gis_victoria.vic_nogeom_roads_volpreds_500 AS b, 
    gis_victoria.vic_nogeom_roads_speedpreds_500 AS c
  WHERE
  a.uid = b.uid
  AND
  a.uid = c.uid
  AND
  ST_Within(a.geom, z.geom)
    ")))
write.csv(roads_CS, file= "output/collision/case_study/roadID_CS.csv")


## The goal is to extract the number of recorded collisions per year for each species. This will help us understand the best split of testing/training years 
### This loop extracts the amount of collision, it is stratified per year for the statewide records and just the sum per species within the Greater Daylesford region ###
coll.data <- matrix(NA, nrow = length(species), ncol = 11)
rownames(coll.data) <- species.plot
colnames(coll.data) <- c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "statewide_total", "regional_total")

coll.coords <- list() #A list of the collision coordinates per species
coll.coords_CS <- NULL #A datafrom for the collision coordinates for all regional records

for (s in 1:6) {
  
  species_name <- species[s] # index one species at a time
  
  #This is for the regional collisions
  all_coll <- data.table::as.data.table(dbGetQuery(con, paste0("
                                     SELECT DISTINCT ON (p.id)
                                     r.uid AS uid, p.year AS year, p.x AS x, p.y AS y, CAST(1 AS INTEGER) AS coll
                                     FROM
                                     gis_victoria.vic_gda9455_roads_state_orig_500 as r,
                                     (SELECT DISTINCT ON (geom)
                                     j.id, j.year, j.x, j.y, j.geom
                                     FROM
                                     (SELECT casenumber AS id, year, ST_X(geom) AS x, ST_Y(geom) AS y, geom
                                     FROM gis_victoria.vic_gda9455_fauna_wv_2010_2013_coll
                                     WHERE
                                     species = '", species_name ,"'
                                     UNION
                                     SELECT id, year, ST_X(geom) AS x, ST_Y(geom) AS y, geom
                                     FROM gis_victoria.vic_gda9455_fauna_wv_2014_2018_coll2
                                     WHERE
                                     species = '", species_name ,"') AS j
                                     ) AS p
                                     WHERE ST_DWithin(p.geom,r.geom, 100)
                                     ORDER BY p.id, ST_Distance(p.geom,r.geom)
                                     ")))
  
  setkey(all_coll, uid)
  
  coll.data[s, 1:9] <- table(all_coll$year)
  coll.data[s, 10] <- sum(coll.data[s, 1:9])
  
  sp.coords <- all_coll[,c("year", "x", "y")]
  
  coll.coords[[s]] <- sp.coords
  
  rm(all_coll)
  
  #and now the same for the regional area
  all_coll_CS <- data.table::as.data.table(dbGetQuery(con, paste0("
                                     SELECT DISTINCT ON (p.id)
                                     r.uid AS uid, p.year AS year, p.x AS x, p.y AS y, CAST(1 AS INTEGER) AS coll
                                     FROM 
                                     (SELECT geom FROM gis_victoria.vic_gda9455_admin_sa2 
                                      WHERE 
                                      (sa3_name16 = 'Upper Goulburn Valley' 
                                      OR sa3_name16 = 'Heathcote - Castlemaine - Kyneton' 
                                      OR sa3_name16 = 'Macedon Ranges' 
                                      OR sa2_name16 = 'Daylesford') 
                                      AND 
                                      (sa2_name16 !=  'Bendigo Region - South' 
                                      AND sa2_name16 != 'Castlemaine Region' 
                                      AND sa2_name16 != 'Mansfield (Vic.)' 
                                      AND sa2_name16 != 'Upper Yarra Valley' 
                                      AND sa2_name16 != 'Castlemaine')) AS z,
                                     gis_victoria.vic_gda9455_roads_state_orig_500 as r,
                                     (SELECT DISTINCT ON (geom)
                                     j.id, j.year, j.x, j.y, j.geom
                                     FROM
                                     (SELECT casenumber AS id, year, ST_X(geom) AS x, ST_Y(geom) AS y, geom
                                     FROM gis_victoria.vic_gda9455_fauna_wv_2010_2013_coll
                                     WHERE
                                     species = '", species_name ,"'
                                     UNION
                                     SELECT id, year, ST_X(geom) AS x, ST_Y(geom) AS y, geom
                                     FROM gis_victoria.vic_gda9455_fauna_wv_2014_2018_coll2
                                     WHERE
                                     species = '", species_name ,"') AS j
                                     ) AS p
                                     WHERE ST_DWithin(p.geom,r.geom, 100)
                                     AND ST_Within(p.geom, z.geom)
                                     ORDER BY p.id, ST_Distance(p.geom,r.geom)
                                     ")))
  
  setkey(all_coll_CS, uid)
  
  coll.data[s, 11] <- nrow(all_coll_CS)
  
  sp.coords_CS <- all_coll_CS[,c("year", "x", "y")]
  sp.coords_CS$species <- rep(species.plot[s], nrow(sp.coords_CS))
  coll.coords_CS <- rbind(coll.coords_CS, sp.coords_CS)
  
  rm(all_coll_CS)
}

names(coll.coords) <- sp.names
saveRDS(coll.coords, file= "output/collision/collision.coords.RData")

write.csv(coll.coords_CS, file= "output/collision/case_study/collision.coords_CS.csv", row.names = FALSE)

write.csv(coll.data, file= "output/collision/coll.data.csv", row.names = TRUE) #this is to show the number of collisions per species per yesr

#This loop simulates the variation in the amount of collision data for testing/training based on the number of training years selected.
## We tested the variation with 6 and 7 years of training data (3 and 2 years of testing data respectively)

n_training_years <- 7
n_testing_years <- 9 - n_training_years
n_samples <- factorial(9) / factorial(n_training_years) * factorial(9 - n_training_years)

coll.train <- matrix(NA, ncol = 6, nrow = n_samples)
colnames(coll.train) <- sp.names
for (s in 1:6) {
  for (n in 1:n_samples) {
    coll.train[n, s] <- sum(sample(coll.data[s, 1:9], size = n_training_years))
  }
}

assign(paste0("coll.train_", n_training_years), as.data.frame(coll.train))

#looking at the testing data
total.coll <- coll.data[, 10]

coll.test_x <- matrix(ncol= 6, nrow = n_samples) #change the x to reflect the number of training years chosen
for (s in 1:6) {
  coll.test_x[, s] <- total.coll[s] - coll.train[, s]
}
colnames(coll.test_x) <- sp.names
assign(paste0("coll.test_", n_training_years), coll.test_x)

# colnames(coll.df_6) <- sp.names
# colnames(coll.df_7) <- sp.names

# Visualising the variation in data
par(mfrow= c(2,2))
boxplot(coll.test_6[,1], coll.test_7[,1], names= c("6:3", "7:2"), main= "koa collision test")
boxplot(coll.test_6[,3], coll.test_7[,3], names= c("6:3", "7:2"), main= "egk collision test")

boxplot(coll.train_6[,1], coll.train_7[,1], names= c("6:3", "7:2"), main= "koa collision train")
boxplot(coll.train_6[,3], coll.train_7[,3], names= c("6:3", "7:2"), main= "egk collision train")

min(coll.test_7[,1])

test.coords <- st_as_sf(coll.coords[["egk"]], coords= c("x", "y"), crs= 28355)

st_distance(test.coords, test.coords) %>%
.[lower.tri(.)] %>% 
  as.vector() %>% 
  assign("egk.dist", .) %>% 
  .[. < 100] %>% 
  hist()

koa.dist.matrix <-  st_distance(test.coords, test.coords) 
koa.dist <-  as.vector(koa.dist.matrix[lower.tri(koa.dist.matrix)])

test.dist <- st_distance(test.coords, test.coords)
test.df <- as.vector(test.dist[lower.tri(test.dist)])

class(test.df)
hist(test.df[test.df < 100])


##################################
### Visualising the study site ###
##################################
greater_daylesford <- st_read("data/VIC_shape/Daylesford_solid.shp")

VIC_latlong <- st_read("data/VIC shape/VIC_GDA94LL_ADMIN_STATE.shp")
VIC_UTM <-  st_transform(VIC_latlong, CRS("+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs"))

aus <- getData("GADM", country = "AUS", level = 1) # download map of Australia

aus <- (aus[aus$NAME_1 != "Ashmore and Cartier Islands" &
              aus$NAME_1 != "Coral Sea Islands Territory" &
              aus$NAME_1 != "Jervis Bay Territory", ])

aus <- crop(aus, extent(c(112.5, 154.0, -44.0, -10.5)))

vic <- (aus[aus$NAME_1 == "Victoria",]) # extract State of Victoria


inset_map <- tm_shape(aus) + tm_borders() + tm_layout() +
  tm_shape(vic) + tm_fill(col= "#cccccc") + tm_borders() 

png(file="output/figs/study.site.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot() +
  geom_sf(data = VIC_UTM, colour= "#333333", fill= "#cccccc") +
  geom_sf(data= greater_daylesford, fill= "#ff9933", colour= "#ff9933") +
  theme_classic() +
  annotation_scale()

print(inset_map, vp = viewport(0.8, 0.8, width = 0.25, height = 0.25))
dev.off()
#### FIN ####
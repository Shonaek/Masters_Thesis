## THis script uses loops to create the lists of data that are needed to run the BRTs for all species at each grain.

library(gbm)
library(raster)
library(sf)
library(blockCV)

#######################
### Reading in Data ###
#######################

#environmental variables and vectors
var.res <- readRDS("data/var.res.RData")
#Vic shapefile
VIC.shape <-raster("data/VIC_GDA9455_GRID_STATE_1000.tif")

#Species data and vectors
species.clean <- read.csv("data/species_clean.csv") 
scientific.names <- unique(species.clean$Scientific.Name)
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")

#Background points and spatial dataframe
bg.point.df <- read.csv("data/bg_df.csv") 
crs.info <- proj4string(VIC.shape)
bg.spdf <- SpatialPointsDataFrame(data= bg.point.df, coords= bg.point.df[,2:3], proj4string = CRS(crs.info))

##############################
### Background points list ###
##############################

# The first list I need is a list of all the background points, with the covariate information at the three grains
## The first colum identifies the 'species' ("BG" for background), then the coordinates, "presence" is the binary column (0 signifies background), 
### "site.weights" is to downweight the background points, it is currently empty because the exact weight depends on the number of species presences, and then the covariate information.

BG.dfs_thick <- list()
for(r in 1:3){ 
  bg.extracted <- raster::extract(var.res[[r]], bg.spdf)
  
  bg.df <- data.frame("Scientific.Name" = rep("BG", nrow(bg.extracted)), bg.point.df[,2:3], "presence" = rep(0, nrow(bg.extracted)), "site.weight" = NA, bg.extracted)
  
  BG.dfs_thick[[r]] <- bg.df
}
names(BG.dfs_thick) <- res.names
saveRDS(BG.dfs_thick, file= "data/BG.dfs_thick.RData")

################################
### Ungridded spatial blocks ###
################################
# This loop saves the spatial block information for each species so that the species data can be subset by fold in the subsequent loops

eval.blocks.thick <- list()

for (s in 1:6) {
  sp.sub <- subset(species.clean, Scientific.Name == scientific.names[s]) 
  sp.pres <- sp.sub[,3:4]
  sp.presence <- cbind("presence" = 1, sp.pres)
  bg.points_250 <- cbind("presence" = 0, BG.dfs_thick[["250"]][,2:3]) 
  sp.thick.PBG_250 <- rbind(sp.presence, bg.points_250)
  
  sp.sf_250 <- st_as_sf(sp.thick.PBG_250, coords = c("x", "y"), crs = crs(crs.info))
  
  
  sb <- blockCV:: spatialBlock(speciesData = sp.sf_250,
                               species = "presence",
                               k = 6,
                               theRange = 50000,
                               rasterLayer = VIC.shape,
                               selection = "random",
                               iteration = 100,
                               seed = 32)
  
  
  
  eval.blocks.thick[[s]] <- sb
}
names(eval.blocks.thick) <- spe.nam
saveRDS(eval.blocks.thick, file= "data/eval.blocks.thick.RData")

################################
### Evaluation "test" blocks ###
################################
# This loop sets up the evaluation folds for each species. The evaluation fold is the confiduation of blocks that will be withheld from the testing data

eval.master.thick <- list()
eval.sp.df.thick <- list()
eval.res.temp <- list()

for(r in 1:3){ #for each resolution
  for(s in 1:6){ #this runs through each species
    
    sp.subset <- subset(species.clean, Scientific.Name == scientific.names[s]) 
    sp.bg.extract <- raster::extract(var.res[[r]], sp.subset[,3:4])
    
    sp.bind <- data.frame(cbind("Scientific.Name" = paste0(scientific.names[[s]]) , sp.subset[,3:4], "presence" = 1, "site.weight" = NA, sp.bg.extract))
    sp.final <- rbind(sp.bind, BG.dfs_thick[[r]])
    for(f in 1:6){ #this runs through each fold
      
      rownums <- sort(unlist(eval.blocks.thick[[s]]$folds[[f]][2])) #the 2nd vector in each is the test sets. They need to be sorted so they subset properly
      fold.subset <- sp.final[rownums,]  
      fold.subset$site.weight[fold.subset$presence == 1] <- 1
      bg.weight <- (nrow(fold.subset[fold.subset$presence==1,])/nrow(fold.subset[fold.subset$presence==0,]))
      fold.subset$site.weight[fold.subset$presence != 1]  <- bg.weight
      
      eval.sp.df.thick[[f]] <- fold.subset
    }
    eval.res.temp[[s]] <- eval.sp.df.thick
  }
  names(eval.res.temp) <- spe.nam
  eval.master.thick[[r]] <- eval.res.temp
}
names(eval.master.thick) <- res.names

saveRDS(eval.master.thick, file=  "data/data_lists/eval.master.thick.RData")

###########################
### Training data folds ###
###########################
# This subsets the species/background data to the complementary folds of the test data. 
## With one set of test/evaluation folds all of the species/background data is used

eval.blocks.thick <- readRDS("data/data_lists/eval.blocks.thick.RData")

#These are the lists
train.master.thick <- list()
train.sp.df.temp <- list()
train.blocks.thick <- list()
train.res.temp <- list()

for(r in 1:3){ #for each resolution
  for(s in 1:6){ #this runs through each species
    sp.subset <- subset(species.clean, Scientific.Name == scientific.names[s]) 
    sp.coords <- sp.subset[,3:4]
    sp.presence <- cbind("presence" = 1, sp.coords)
    bg.points_res <- cbind("presence" = 0, BG.dfs_thick[[r]][,2:3])
    sp.PBG_res <- rbind(sp.presence, bg.points_res)
    
    
    sp.sf <- st_as_sf(sp.PBG_res, coords = c("x", "y"), crs = crs(crs.info))
    
    train.fold_temp <- spatialBlock(speciesData = sp.sf,      #This functions takes the spatial block information and assigns the rows to folds 
                                    species = "presence",
                                    selection = "predefined",
                                    k= 6,
                                    blocks = eval.blocks.thick[[s]]$blocks,
                                    foldsCol = "folds",
                                    seed = 32)
    train.blocks.thick[[s]] <- train.fold_temp
    
    sp.extract <- extract(var.res[[r]], sp.coords)
    sp.bind_res <- data.frame(cbind("Scientific.Name" = paste0(scientific.names[[s]]) , sp.coords, "presence" = 1, "site.weight" = NA, sp.extract))
    sp.final_res <- rbind(sp.bind_res, BG.dfs_thick[[r]])
    
    for(f in 1:6){ #this runs through each fold
      rownums <- sort(unlist(train.blocks.thick[[s]]$folds[[f]][1])) #the first vector in the folds list is a 
      
      fold.subset <- sp.final_res[rownums,]  
      fold.subset$site.weight[fold.subset$presence == 1] <- 1
      bg.weight <- (nrow(fold.subset[fold.subset$presence==1,])/nrow(fold.subset[fold.subset$presence==0,]))
      fold.subset$site.weight[fold.subset$presence != 1]  <- bg.weight
      
      train.sp.df.temp[[f]] <-fold.subset
    }
    train.res.temp[[s]] <- train.sp.df.temp
    #print(paste0("finished_", spe.nam[s], "_", res.names[r]))
  }
  names(train.res.temp) <- spe.nam
  train.master.thick[[r]] <- train.res.temp
}
names(train.master.thick) <- res.names
saveRDS(train.master.thick, file= "data/data_lists/train.master.thick.RData")


########################
### Whole Data Lists ###
########################
# To make the spatial likelihood of species occurence layers I need a datalist with the full set of species/background locations the extracted covariate data

res.temp <- list()
whole.master.thick <- list()
for(r in 1:3){ 
  for (s in 1:6){
    sp.subset <- subset(species.clean, Scientific.Name == scientific.names[s]) 
    sp.coords <- sp.subset[,3:4]
    sp.bg<- raster::extract(var.res[[r]], sp.coords)
    
    sp.bind <- data.frame(cbind("Scientific.Name" = paste0(scientific.names[[s]]) , sp.coords, "presence" = 1, "site.weight" = NA, sp.bg))
    sp.final <- rbind(sp.bind, BG.dfs_thick[[r]])
    
    sp.final$site.weight[sp.final$presence == 1] <- 1
    bg.weight <- (nrow(sp.final[sp.final$presence==1,])/nrow(sp.final[sp.final$presence==0,]))
    sp.final$site.weight[sp.final$presence != 1]  <- bg.weight
    
    res.temp[[s]] <- sp.final
  }
  names(res.temp) <- spe.nam
  whole.master.thick[[r]] <- res.temp
}
names(whole.master.thick) <- res.names
saveRDS(whole.master.thick, file= "data/data_lists/whole.master.thick.RData")

#################################
### Gridded "thin" data lists ###
#################################
#The gridded data lists are created with the exact same workflow 
## The training data thins the species presence and background to one (each) per cell to match the covariate grain
### So that there is fair testing across all grains, the evaluation (test) folds are ONLY thinned to the 250m grain. 
#### Covariate information is ectracted at the three grains but the number of presence points isn't thinned further (to give even testing)

#### RE READING IN THE DATA (incase you need it again) ####
#Vic shapefile
VIC.shape <-raster("data/VIC_GDA9455_GRID_STATE_1000.tif")

species.clean <- read.csv("data/species_clean.csv") #This is the clean, unthinned master dataframe with all 6 species

bg.point.df <- read.csv("data/bg_df.csv") #this is the dataframe of the unthinned background points (with a minimum distance of 100m) 

var.res <- readRDS("data/var.res.RData") #this is the stack of covariate data at the three resolutions

scientific.names <- unique(species.clean$Scientific.Name)
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")

#### GRIDDED BACKGROUND POINTS ####
# This is set of background points with the covariate information at each grain. 
## Background data is also thinned so that there is one background point per cell at the respective grains
crs.info <- proj4string(VIC.shape)
bg.spdf <- SpatialPointsDataFrame(data= bg.point.df, coords= bg.point.df[,2:3], proj4string = CRS(crs.info))

BG.dfs_thin <- list()
for(r in 1:3){ 
  bg.thin <- gridSample(bg.spdf, var.res[[r]], n=1)
  bg.extracted <- extract(var.res[[r]], bg.thin)
  
  bg.df <- data.frame("Scientific.Name" = rep("BG", nrow(bg.thin)) , bg.thin, "presence" = rep(0, nrow(bg.thin)), "site.weight" = NA, bg.extracted)
  
  BG.dfs_thin[[r]] <- bg.df
}
names(BG.dfs_thin) <- res.names #like the predictor covariates, I now hove a list with the background ooints gridded to each grain
saveRDS(BG.dfs_thin, file= "data/BG.dfs_thin.RData")

####SPATIAL BLOCKS ####
# The blocking structure is based on the number of presence points at the 250m grain for all species

eval.blocks.thin <- list()

for (s in 1:6) {
  sp.subset <- subset(species.clean, Scientific.Name == scientific.names[s]) 
  
  sp.thin <- gridSample(sp.subset[,3:4], var.res[["250"]], n=1)
  sp.presence <- cbind("presence" = 1, sp.thin)
  bg.points_250 <- cbind("presence" = 0, BG.dfs_thin[["250"]][,2:3])
  sp.thin.PBG_250 <- rbind(sp.presence, bg.points_250)
  
  sp.sf_250 <- st_as_sf(sp.thin.PBG_250, coords = c("x", "y"), crs = crs(crs.info))
  
  
  sb1 <- blockCV:: spatialBlock(speciesData = sp.sf_250,
                                species = "presence",
                                k = 6,
                                theRange = 50000,
                                rasterLayer = VIC.shape,
                                selection = "random",
                                iteration = 100,
                                seed = 32)
  
  
  #assign(paste0(spe.nam[s],".eval.thin"), sb1)
  eval.blocks.thin[[s]] <- sb1
}
names(eval.blocks.thin) <- spe.nam
saveRDS(eval.blocks.thin, file= "data/eval.blocks.thin.RData")

#### BACKGROUND EVALUATION POINTS ####
# The thinning structure remains fixed at the 250m grain for the evaluation background points, 
## but the covariate information is extracted for all grains
BG.eval_thin <- list()
BG.eval.points <- BG.dfs_thin[["250"]][,2:3]

for(r in 1:3){
  bg.extracted <- raster::extract(var.res[[r]], BG.eval.points)
  
  bg.df <- data.frame(cbind("Scientific.Name" = "BG" , BG.eval.points, "presence" = 0, "site.weight" = NA, bg.extracted))
  
  BG.eval_thin[[r]] <- bg.df
}
names(BG.eval_thin) <- res.names


#### EVALUTATION LISTS ####
# Similarly here you can see that the thinning is fixed at the 250 grain, rather than changing with each grain. 
## the covariate information is still being extracted at the three grains

eval.master.thin <- list()
eval.sp.df.thin <- list()
eval.res.temp <- list()


for(r in 1:3){ #for each resolution
  for(s in 1:6){ #this runs through each species
    
    sp.subset <- subset(species.clean, Scientific.Name == scientific.names[s]) 
    sp.thin_250 <- gridSample(sp.subset[,3:4], var.res[["250"]], n=1) #this thins the species data
    sp.bg.extract_250 <- raster::extract(var.res[[r]], sp.thin_250) #this extracts the species information
    
    sp.bind_250 <- data.frame(cbind("Scientific.Name" = paste0(scientific.names[[s]]) , sp.thin_250, "presence" = 1, "site.weight" = NA, sp.bg.extract_250))
    sp.final_250 <- rbind(sp.bind_250, BG.eval_thin[[r]])
    for(f in 1:6){ #this runs through each fold
      
      rownums <- sort(unlist(eval.blocks.thin[[s]]$folds[[f]][2])) #the 2nd vector in each is the test sets. They need to be sorted so they subset properly
      fold.subset <- sp.final_250[rownums,]  
      fold.subset$site.weight[fold.subset$presence == 1] <- 1
      bg.weight <- (nrow(fold.subset[fold.subset$presence==1,])/nrow(fold.subset[fold.subset$presence==0,]))
      fold.subset$site.weight[fold.subset$presence != 1]  <- bg.weight
      
      eval.sp.df.thin[[f]] <- fold.subset
    }
    eval.res.temp[[s]] <- eval.sp.df.thin
  }
  names(eval.res.temp) <- spe.nam
  eval.master.thin[[r]] <- eval.res.temp
}
names(eval.master.thin) <- res.names

saveRDS(eval.master.thin, file=  "data/data_lists/eval.master.thin.RData")

#### TRAINING FOLDS ####
# Unline the testing data, the training blocks are thinned (gridded) at each grain.

train.master.thin <- list()
train.sp.df.temp <- list()
train.blocks.thin <- list()
train.res.temp <- list()

for(r in 1:3){ 
  for (s in 1:6){
    sp.subset <- subset(species.clean, Scientific.Name == scientific.names[s]) 
    
    sp.thin <- gridSample(sp.subset[,3:4], var.res[[r]], n=1)
    sp.presence <- cbind("presence" = 1, sp.thin)
    bg.points_res <- cbind("presence" = 0, BG.dfs_thin[[r]][,2:3])
    sp.thin.PBG_res <- rbind(sp.presence, bg.points_res)
    
    
    sp.sf <- st_as_sf(sp.thin.PBG_res, coords = c("x", "y"), crs = crs(crs.info))
    
    train.fold_res <- spatialBlock(speciesData = sp.sf,
                                   species = "presence",
                                   selection = "predefined",
                                   k= 6,
                                   blocks = eval.blocks.thin[[s]]$blocks,
                                   foldsCol = "folds",
                                   seed = 32)
    train.blocks.thin[[s]] <- train.fold_res
    
    sp.extract <- extract(var.res[[r]], sp.thin)
    sp.bind_res <- data.frame(cbind("Scientific.Name" = paste0(scientific.names[[s]]) , sp.thin, "presence" = 1, "site.weight" = NA, sp.extract))
    sp.final_res <- rbind(sp.bind_res, BG.dfs_thin[[r]])
    for(f in 1:6){
      rownums <- sort(unlist(train.blocks.thin[[s]]$folds[[f]][1])) 
      
      fold.subset <- sp.final_res[rownums,]  
      fold.subset$site.weight[fold.subset$presence == 1] <- 1
      bg.weight <- (nrow(fold.subset[fold.subset$presence==1,])/nrow(fold.subset[fold.subset$presence==0,]))
      fold.subset$site.weight[fold.subset$presence != 1]  <- bg.weight
      
      train.sp.df.temp[[f]] <-fold.subset
    }
    train.res.temp[[s]] <- train.sp.df.temp
  }
  names(train.res.temp) <- spe.nam
  train.master.thin[[r]] <- train.res.temp
}
names(train.master.thin) <- res.names
saveRDS(train.master.thin, file= "data/data_lists/train.master.thin.RData")

#### WHOLE DATA MODELS ####
# The whole data models don't use spatial blocking or testing data
## The species occurences are thinned to one record per grain, as are the background points.

res.temp <- list()
whole.master.thin <- list()
for(r in 1:3){ 
  for (s in 1:6){
    sp.subset <- subset(species.clean, Scientific.Name == scientific.names[s]) 
    sp.thin <- gridSample(sp.subset[,3:4], var.res[[r]], n=1)
    sp.bg<- raster::extract(var.res[[r]], sp.thin)
    
    sp.bind <- data.frame(cbind("Scientific.Name" = paste0(scientific.names[[s]]) , sp.thin, "presence" = 1, "site.weight" = NA, sp.bg))
    sp.final <- rbind(sp.bind, BG.dfs_thin[[r]])
    
    sp.final$site.weight[sp.final$presence == 1] <- 1
    bg.weight <- (nrow(sp.final[sp.final$presence==1,])/nrow(sp.final[sp.final$presence==0,]))
    sp.final$site.weight[sp.final$presence != 1]  <- bg.weight
    
    res.temp[[s]] <- sp.final
  }
  names(res.temp) <- spe.nam
  whole.master.thin[[r]] <- res.temp
}
names(whole.master.thin) <- res.names
saveRDS(whole.master.thin, file= "data/data_lists/whole.master.thin.RData")

###########
### FIN ###
###########
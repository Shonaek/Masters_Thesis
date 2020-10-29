#This takes the raw glm model outputs, interrogates them and prepares them for visualisation

species.plot <- c("Koala","Brushtail", "Kangaroo", "Ringtail", "Wallaby", "Wombat")
res.names <- c("250", "500", "1000")
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")
sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom")

############################
### Statewide AUC Values ###
############################

coll.glm.eval <- NULL

for (r in 1:3) {
  for (s in 1:6) {
    
    
    coll.glm.raw <- read.csv(paste0("output/collision/glm_metrics_", species[s] , "_", res.names[r], ".csv")) #this individually reads in the csv files that were saved
    
    coll.glm.sp <- cbind.data.frame("AUC_train" = coll.glm.raw$AUC_train, "AUC_test" = coll.glm.raw$AUC_test, "species" = species.plot[s], "res"= res.names[r])
    coll.glm.sp$sd_train <-  sd(coll.glm.sp$AUC_train)
    coll.glm.sp$se_train <- (coll.glm.sp$sd_train/sqrt(nrow(coll.glm.sp)))
    coll.glm.sp$sd_test <-  sd(coll.glm.sp$AUC_test)
    coll.glm.sp$se_test <- (coll.glm.sp$sd_test/sqrt(nrow(coll.glm.sp)))
    coll.glm.eval <- rbind(coll.glm.eval, coll.glm.sp, stringsAsFactors = FALSE)
  }
}

write.csv(coll.glm.eval, file= "output/collision/coll.glm.eval.csv", row.names = FALSE)
################################
### Regional (CS) AUC values ###
################################

coll.auc_CS <- NULL

for (r in 1:3) {
  for (s in 1:6) {
    
    coll.glm.CS <- readRDS(paste0("output/collision/case_study/coll_glm_", sp.names[s], "_", res.names[r], "_CS.RData"))
    glm.data.CS <- coll.data.list_CS[[r]][[s]] 
    
    AUC <- au.roc(glm.data.CS$coll, predict(coll.glm.CS, glm.data.CS, type="response"))
    
    AUC.sp <- data.frame("AUC_CS" = AUC, "species" = species.plot[s], "res" = res.names[r])
    
    coll.auc_CS <- rbind(coll.auc_CS, AUC.sp)
  }
}
write.csv(coll.auc_CS, file = "output/collision/case_study/coll.auc_CS.csv")



########################################################
### Collision/Background overlap with SDM Difference ###
########################################################
# So the idea here is to extract the collision and background points from the SDM difference layers

#First we need to read in the relevant data
diff.stack.thick <- readRDS("output/1000-250.stack.thick.RData")
coll.dat <- readRDS("output/collision/collision.coords.RData") #from the data exploration script
coll.data.list_CS <- readRDS("output/collision/case_study/coll.data.list_CS.RData")
roadID <- read.csv("output/collision/roadID.csv")

diff.stack.thick.abs <- abs(diff.stack.thick) #absolute difference of the 1000-250m layers

utm.crs <- crs(diff.stack.thick)

road.spdf <- SpatialPointsDataFrame(coords = roadID[,4:5], data= roadID, proj4string = utm.crs) 

diff.coll <- NULL
for (s in 1:6) {
  diff <- raster::extract(diff.stack.thick.abs[[s]], coll.dat[[sp.names[s]]][,2:3]) #mcapply
  diff <- na.omit(diff)
  temp <- data.frame("diff" = diff, "species"= species.plot[s], "type"= "collision", "area" = "statewide")
  diff.coll <- rbind(diff.coll, temp)
  
  road.diff <- raster::extract(diff.stack.thick.abs[[s]], road.spdf)
  road.diff <- na.omit(road.diff)
  road.temp <- data.frame("diff" = road.diff, "species"= species.plot[s], "type"= "non-collision", "area" = "statewide")
  diff.coll <- rbind(diff.coll, road.temp)
  
  coll.cs <- coll.data.list_CS[["1000"]][[s]][,c("uid", "coll")] #the resolution shouln't matter I just need the collisions and the uid
  coords.cs <- merge(coll.cs, road.spdf, by = "uid")
  diff.bg.cs <- raster::extract(diff.stack.thick.abs[[s]], coords.cs[, c("x", "y")]) 
  diff.bg.cs <- na.omit(diff.bg.cs)
  cs.bg <- data.frame("diff" = diff.bg.cs, "species"= species.plot[s], "type"= "non-collision", "area" = "regional")
  diff.coll <- rbind(diff.coll, cs.bg)
  
  coords.coll <- filter(coords.cs, coll == 1)
  diff.sp.cs <- raster::extract(diff.stack.thick.abs[[s]], coords.coll[, c("x", "y")]) 
  diff.sp.cs <- na.omit(diff.sp.cs)
  cs.sp <- data.frame("diff" = diff.sp.cs, "species"= species.plot[s], "type"= "collision", "area" = "regional")
  diff.coll <- rbind(diff.coll, cs.sp)
  
}

write.csv(diff.coll, file= "output/collision/collision.density.csv", row.names = FALSE)

#############################
### collision predictions ###
#############################

### STATEWIDE ###

roadID <- read.csv("output/collision/roadID.csv")
coll.data.list <- readRDS("output/collision/coll.data.list.RData")

#Creating a dataframe with the road coordinates and making the collision predictions
coll.preds.res <- list()
coll.preds.list <- list()

for (r in 1:3) {
  for (s in 1:6) {
    
    glm.data <- coll.data.list[[r]][[s]] 
    coll.glm <- readRDS(paste0("output/collision/coll_glm_", sp.names[s], "_", res.names[r], ".RData"))
    
    coll.pred <- predict(coll.glm, glm.data, type= "response") 
    preds <- as.data.frame(cbind(uid = glm.data$uid, species = glm.data$species, coll= glm.data$coll, coll.pred))
    
    pred.df <- merge(preds, roadID[, 2:7], by= "uid") 
    
    coll.preds.res[[s]] <- pred.df
  }
  names(coll.preds.res) <- sp.names
  coll.preds.list[[r]] <- coll.preds.res
}
names(coll.preds.list) <- res.names

saveRDS(coll.preds.list, file = "output/collision/coll.preds.list.RData")


### REGIONAL ###
roadID_CS <- read.csv(file= "output/collision/case_study/roadID_CS.csv") #this dataframe contains the coordinates of the road segments in the study area
coll.data.list_CS <- readRDS("output/collision/case_study/coll.data.list_CS.RData")
coll.preds.res <- list()
coll.preds.list_CS <- list()

#This loop takes the collision data, makes the predictions and then matches the uid of the roads to the coordinates

for (r in 1:3) {
  for (s in 1:6) {
    
    glm.data <- coll.data.list_CS[[r]][[s]] 
    coll.glm <- readRDS(paste0("output/collision/case_study/coll_glm_", sp.names[s], "_", res.names[r], "_CS.RData"))
    
    coll.pred <- predict(coll.glm, glm.data, type= "response") 
    preds <- as.data.frame(cbind(uid = glm.data$uid, species = glm.data$species, coll= glm.data$coll, coll.pred))
    
    pred.df <- merge(preds, roadID_CS[, 2:7], by= "uid") 
    
    coll.preds.res[[s]] <- pred.df
  }
  names(coll.preds.res) <- sp.names
  coll.preds.list_CS[[r]] <- coll.preds.res
}
names(coll.preds.list_CS) <- res.names
saveRDS(coll.preds.list_CS, file= "output/collision/case_study/coll.preds.list_CS.RData")


###################################################
### Standardising and exporting the predictions ###
###################################################
#This takes the raw collision predictions, and creates csv files per grain containing the standardised collision predictions
## At each segment plus the mean across all species. The uid is necessary to link the dataframes with the road segment shapefle in QGIS

coll.preds.list <- readRDS("output/collision/coll.preds.list.RData") 
coll.preds.list_CS <- readRDS("output/collision/case_study/coll.preds.list_CS.RData")

for (r in 1:3) {
  
  #starewide predictions
  coll.preds.unlist <- do.call(cbind, lapply(coll.preds.list[[res.names[r]]], FUN =  function(x) x$coll.pred))
  coll.preds.std <- apply(coll.preds.unlist, 2, FUN= function(x) (x-min(x))/(max(x)-min(x)))
  coll.preds.std <- as.data.frame(coll.preds.std)
  coll.preds.std$mean <- apply(coll.preds.std, 1, FUN = mean)
  coll.preds.std$uid <- coll.preds.list[[res.names[r]]][[1]]$uid
  assign(paste0("collision_predictions_std_", res.names[r]), coll.preds.std)
  
  #regional predictions
  cs.coll.preds.unlist <- do.call(cbind, lapply(coll.preds.list_CS[[res.names[r]]], FUN =  function(x) x$coll.pred))
  cs.coll.preds.std <- apply(cs.coll.preds.unlist, 2, FUN= function(x) (x-min(x))/(max(x)-min(x)))
  cs.coll.preds.std <- as.data.frame(cs.coll.preds.std)
  cs.coll.preds.std$mean <- apply(cs.coll.preds.std, 1, FUN = mean)
  cs.coll.preds.std$uid <- coll.preds.list_CS[[res.names[r]]][[1]]$uid
  
  #now I want to add a columns to the regional area that show the cropped statewide predictions for each grain
  statewide.crop <- merge(coll.preds.std[, c("uid", "mean")], cs.coll.preds.std[, c("uid", "mean")], by= "uid", drop.y = TRUE)
  cs.coll.preds.std$statewide <- statewide.crop$mean.x
  #difference between the mean regional prediction values and the statewide values at each segment
  cs.coll.preds.std$diff <- cs.coll.preds.std$mean - cs.coll.preds.std$statewide
  
  assign(paste0("collision_predictions_std_CS_", res.names[r]), cs.coll.preds.std)
}



write.csv(collision_predictions_std_250, file = "output/collision/collision_predictions_std_250.csv", row.names = FALSE)
write.csv(collision_predictions_std_500, file = "output/collision/collision_predictions_std_500.csv", row.names = FALSE)
write.csv(collision_predictions_std_1000, file = "output/collision/collision_predictions_std_1000.csv", row.names = FALSE)

write.csv(collision_predictions_std_CS_250, file = "output/collision/case_study/collision_predictions_CS_std_250.csv", row.names = FALSE)
write.csv(collision_predictions_std_CS_500, file = "output/collision/case_study/collision_predictions_CS_std_500.csv", row.names = FALSE)
write.csv(collision_predictions_std_CS_1000, file = "output/collision/case_study/collision_predictions_CS_std_1000.csv", row.names = FALSE)


###################################
### Looking at model parameters ###
###################################

### STATEWIDE ###

coll.coeffs_VIC <- NULL
for (s in 1:6) {
  for (r in 1:3) {
    
    coll.glm <- readRDS(paste0("output/collision/coll_glm_", sp.names[s], "_", res.names[r], ".RData"))
    coll.coeff <- t(coll.glm$coefficients)
    coll.coeff.data <- as.data.frame(coll.coeff)
    coll.coeff.data$species <- species.plot[s]
    coll.coeff.data$res <- res.names[r]
    coll.coeff.data$area <- "all_VIC"
    
    coll.coeffs_VIC <- rbind(coll.coeff.data, coll.coeffs_VIC )
    
  }}

write.csv(coll.coeffs_VIC, file= "output/collision/coll.coeffs_VIC.csv", row.names = FALSE)

### REGIONAL ###
coll.coeffs_CS <- NULL
for (s in 1:6) {
  for (r in 1:3) {
    
    coll.glm <- readRDS(paste0("output/collision/case_study/coll_glm_", sp.names[s], "_", res.names[r], "_CS.RData"))
    coll.coeff <- t(coll.glm$coefficients)
    coll.coeff.data <- as.data.frame(coll.coeff)
    coll.coeff.data$species <- species.plot[s]
    coll.coeff.data$res <- res.names[r]
    coll.coeff.data$area <- "case_study"
    
    coll.coeffs_CS <- rbind(coll.coeff.data, coll.coeffs_CS )
    
  }}
write.csv(coll.coeffs_CS, file= "output/collision/case_study/coll.coeffs_CS.csv", row.names = FALSE)

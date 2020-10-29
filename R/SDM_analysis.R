## This script follows on from the BRT_Loops. It does the analysis on model performance and spatial predictions

######################
### AUC evaluation ###
######################
res.names <- c("250", "500", "1000")
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")
sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom")
external.preds_thick <- readRDS("output/external.preds_thick.RData")

## This function takes the observation data (presence or background) and predicted likelihood of occurence within that cell and returns the area under the reciever operator curve (AUC)

au.roc <- function (obsdat, preddat) {
  if (length(obsdat) != length(preddat)) 
    stop("obs and preds must be equal lengths")
  n.x <- length(obsdat[obsdat == 0])
  n.y <- length(obsdat[obsdat == 1])
  xy <- c(preddat[obsdat == 0], preddat[obsdat == 1])
  rnk <- rank(xy)
  wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * 
                                                                     n.y)
  return(round(wilc, 4))
}

#The external.preds_thick list contains the presence/background points of the testing fold and predicted likelihoods of occurences at those sites
eval.auc.thick <- list()

#this loop calculates the AUC values for each fold, the standard deviations and means are calculated across folds for each species and resolution
## note the SD were not shown in the final figures due to low statistical power, instead the AUC value for each fold is plotted
for (r in 1:3) {
  eval.auc_res <- data.frame(matrix(NA, nrow=6, ncol=10))
  names(eval.auc_res) <- c(paste0("fold_", 1:6), "sd", "mean", "upper", "lower")
  rownames(eval.auc_res) <- spe.nam
  for (s in 1:6) {
    sp.res <- external.preds_thick[[r]][[s]]
    for (f in 1:6) {
      sp.fold <- sp.res[[f]]
      auc.fold <- au.roc(sp.fold[,1], sp.fold[,2])
      eval.auc_res[s,f] <- auc.fold
    }
    eval.auc_res[s,7] <- sd(eval.auc_res[s, 1:6])
    
  }
  eval.auc_res[,8] <- round(rowMeans(eval.auc_res[1:6]), 3)
  eval.auc_res[,9] <- round(eval.auc_res[,8]+ eval.auc_res[,7], 3)
  eval.auc_res[,10] <- round(eval.auc_res[,8]- eval.auc_res[,7], 3)
  eval.auc.thick[[r]]  <- eval.auc_res  
}
names(eval.auc.thick) <- res.names
saveRDS(eval.auc.thick, file= "output/eval.auc.thick.RData")

##############################
### Kulczynski Coefficient ###
##############################
library(raster)
library(dismo)

sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom")

distribution.layers.thick <- readRDS("output/distribution.layers.thick.RData") #This is the predicted rasters

KUL.thick <- list() #This list will contain the statistics for all species

for (s in 1:6) {
  #this resamples the three distribution grains to all be 250m
  sp_250 <- distribution.layers.thick[["250"]][[spe.nam[s]]]
  sp_500 <- resample(distribution.layers.thick[["500"]][[spe.nam[s]]], distribution.layers.thick[["250"]][[spe.nam[s]]])
  sp_1000 <- resample(distribution.layers.thick[["1000"]][[spe.nam[s]]], distribution.layers.thick[["250"]][[spe.nam[s]]])
  
  #if your data is already at a common grain, you can start here, by stacking your data and extracting sample values
  mydat <- stack(sp_1000, sp_500, sp_250) 
  
  samp <- dismo::randomPoints(sp_250, 1000000) # get a sample of points; used 1M (14%) here since there are millions of values
  vals <- raster::extract(mydat, samp)
  vals <- na.omit(vals)
  colnames(vals) <- c(paste0(sp.names[s], "_1000"), paste0(sp.names[s], "_500"), paste0(sp.names[s], "_250"))
  
  
  X <- data.frame(vals)
  OUT <- matrix(NA, nrow=ncol(X), ncol=ncol(X))
  rownames(OUT) <- colnames(OUT) <- names(X)  #add names for simplicity
  
  OUT[lower.tri(OUT)] <- combn(X, 2, function(x)
  {
    two <- sum(x[,2])
    mins <- apply(x,1,min)
    mins.sum <- sum(mins)
    round((two - mins.sum)/two,3)
  })
  
  OUT <- t(OUT)
  
  
  OUT[lower.tri(OUT)] <- combn(X, 2, function(x)
  {
    one <- sum(x[,1])
    mins <- apply(x,1,min)
    mins.sum <- sum(mins)
    round((one - mins.sum)/one,3)
  })
  
  KUL.thick[[s]] <- OUT 
}
names(KUL.thick) <- sp.names
saveRDS(KUL.thick, file = "output/KUL.thick.RData")


#################################################################
### Difference between the 1000m and 250m spatial predictions ###
#################################################################
distribution.layers.thick <- readRDS("output/distribution.layers.thick.RData")

difference.thick <- list()
for (s in 1:6) {
  sp.resample <- resample(distribution.layers.thick[["1000"]][[s]], distribution.layers.thick[["250"]][[s]])
  difference.thick[[s]] <- (sp.resample- distribution.layers.thick[["250"]][[s]])
}
names(difference.thick) <- sp.names
diff.thick <- stack(difference.thick)
diff.stack.plot <- subset(diff.thick, order(species))
saveRDS(diff.stack.thick, file= "output/1000-250.stack.thick.RData")

################################
### Thinned model evaluation ###
################################

### AUC ###
eval.auc.thin <- list()

#this loop calculates the AUC values for each fold, the standard deviations and means are calculated across folds for each species and resolution
## note the SD were not shown in the final figures due to low statistical power, instead the AUC value for each fold is plotted
for (r in 1:3) {
  eval.auc_res <- data.frame(matrix(NA, nrow=6, ncol=10))
  names(eval.auc_res) <- c(paste0("fold_", 1:6), "sd", "mean", "upper", "lower")
  rownames(eval.auc_res) <- spe.nam
  for (s in 1:6) {
    sp.res <- external.preds_thin[[r]][[s]]
    for (f in 1:6) {
      sp.fold <- sp.res[[f]]
      auc.fold <- au.roc(sp.fold[,1], sp.fold[,2])
      eval.auc_res[s,f] <- auc.fold
    }
    eval.auc_res[s,7] <- sd(eval.auc_res[s, 1:6])
    
  }
  eval.auc_res[,8] <- round(rowMeans(eval.auc_res[1:6]), 3)
  eval.auc_res[,9] <- round(eval.auc_res[,8]+ eval.auc_res[,7], 3)
  eval.auc_res[,10] <- round(eval.auc_res[,8]- eval.auc_res[,7], 3)
  eval.auc.thin[[r]]  <- eval.auc_res  
}
names(eval.auc.thin) <- res.names
saveRDS(eval.auc.thin, file= "output/eval.auc.thin.RData")

### KULCZYNSKI ###

distribution.layers.thin <- readRDS("output/distribution.layers.thin.RData") #This is the list of predicted rasters

KUL.thin <- list() #This list will contain the statistics for all species

for (s in 1:6) {
  
  sp_250 <- distribution.layers.thin[["250"]][[spe.nam[s]]]
  sp_500 <- resample(distribution.layers.thin[["500"]][[spe.nam[s]]], distribution.layers.thin[["250"]][[spe.nam[s]]])
  sp_1000 <- resample(distribution.layers.thin[["1000"]][[spe.nam[s]]], distribution.layers.thin[["250"]][[spe.nam[s]]])
  mydat <- stack(sp_1000, sp_500, sp_250)
  samp <- dismo::randomPoints(sp_250, 1000000) # get a sample of points; used 1M (14%) here since there are millions of values
  vals <- raster::extract(mydat, samp)
  vals <- na.omit(vals)
  colnames(vals) <- c(paste0(sp.names[s], "_1000"), paste0(sp.names[s], "_500"), paste0(sp.names[s], "_250"))
  
  
  X <- data.frame(vals)
  OUT <- matrix(NA, nrow=ncol(X), ncol=ncol(X))
  rownames(OUT) <- colnames(OUT) <- names(X)  #add names for simplicity
  
  OUT[lower.tri(OUT)] <- combn(X, 2, function(x)
  {
    two <- sum(x[,2])
    mins <- apply(x,1,min)
    mins.sum <- sum(mins)
    round((two - mins.sum)/two,3)
  })
  
  OUT <- t(OUT)
  
  
  OUT[lower.tri(OUT)] <- combn(X, 2, function(x)
  {
    one <- sum(x[,1])
    mins <- apply(x,1,min)
    mins.sum <- sum(mins)
    round((one - mins.sum)/one,3)
  })
  
  KUL.thin[[s]] <- OUT 
}
names(KUL.thin) <- sp.names
saveRDS(KUL.thin, file = "output/KUL.thin.RData")

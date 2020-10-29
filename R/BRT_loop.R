
#Ths script is for running all of the BRT models and collecting the model performance information

#The script is seperated into two distinct sections: first the code is run on the ungridded data (which allows for multiple species presences per cell)
## The ungridded data is presented in the manuscript but the gridded (one presence per cell)  results appear in the appendices, so the gridded loops are presented in the second half of this script

#necessary packages
library(gbm)
library(raster)
library(dismo)

############################################################
### Reading in the data for the ungridded (thick) models ###
############################################################
# see the "SDM_data_lists" script for how they were created
eval.master.thick <- readRDS("data/data_lists/eval.master.thick.RData")
train.master.thick <- readRDS("data/data_lists/train.master.thick.RData")
whole.master.thick <- readRDS("data/data_lists/whole.master.thick.RData")

#These are the name vectors
res.names <- c("250", "500", "1000")
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")

var.res <- readRDS("data/var.res.RData") #the environmental covariates
eval.var.res <- readRDS("data/eval.var.res.RData") #the covariates for making model predictions (bias layers set to mean values)

#########################################################
### Running the spatially-blocked (evaluation) models ###
#########################################################

#Setting up the lists for the evaluation stats. All lists are nested by resolution 
## there are temporary lists which collect all of the species data at one resolution. so they can then be assigned to the final list which contains all three resolutions

internal.stats_thick <- list() #the number of trees, learning rate and discrimination mean (internal AUC) were extracted from each model
int.stats.sp <- list() #This is a temporary list
external.preds_thick <- list() # This contains the evaluation data (species presence/BG and species predictions) at all resolutions, species and folds
ex.preds.folds <- list() #temporary list for the evaluation data of one species at all of the folds
ex.preds.sp <- list() #temp list for all of the species with their fold data at one resolution

#### RUNNING THE BRT PER FOLD
for (r in 1:3) {
  for (s in 1:6) {
    traindat.sp <- train.master.thick[[r]][[s]]
    evaldat.sp <- eval.master.thick[[r]][[s]]
    sp.internal.df <- data.frame(matrix(NA, nrow=6, ncol=3))
    names(sp.internal.df) <- c("lr","best.trees", "AUC")
    for (f in 1:6) {
      folddat.train <- traindat.sp[[f]]
      brt.train  <- gbm.step(data=folddat.train, 
                             gbm.x = c(6:21),
                             gbm.y = 4,
                             family = "bernoulli",
                             site.weights = folddat.train$site.weight,
                             tree.complexity = 1,
                             learning.rate = 0.05,
                             bag.fraction = 0.75,
                             plot.main=FALSE, silent=TRUE, verbose=FALSE)
      
      
      # saving the internal evaluation stats from the training fold
      sp.internal.df[f,1] <- brt.train$gbm.call$learning.rate
      sp.internal.df[f,2] <- brt.train$gbm.call$best.trees
      sp.internal.df[f,3] <- brt.train$cv.statistics$discrimination.mean
      
      #external evaluation to the test fold
      folddat.eval <- evaldat.sp[[f]]
      
      preds <- predict.gbm(brt.train, folddat.eval, n.trees=brt.train$gbm.call$best.trees, type="response")
      pred.dat <- cbind(folddat.eval$presence, preds)
      ex.preds.folds[[f]] <- pred.dat
    }
    #The data from each fold, gets saved into the relavent species list for that resolution
    ex.preds.sp[[s]] <- ex.preds.folds
    int.stats.sp[[s]] <- sp.internal.df
    print(paste0("finished.", spe.nam[s], "_", res.names[r]))
  }
  names(ex.preds.sp) <- spe.nam
  names(int.stats.sp) <- spe.nam
  external.preds_thick[[r]] <- ex.preds.sp
  internal.stats_thick[[r]] <- int.stats.sp
}
names(external.preds_thick) <- res.names
saveRDS(external.preds_thick, file= "output/external.preds_thick.RData")

names(internal.stats_thick) <- res.names
saveRDS(internal.stats_thick, file= "output/internal.stats_thick.RData")


#############################################################
### Whole data models for species predictions (ungridded) ###
#############################################################

#now that we have the evaluation metrics, we need to run a final brt odel without spatial blocking to make the likelihood of occurence layers 

#These are all of the temporary and final lists needed for the loop
internal.stats.whole.thick <- list() #saving the number of trees, learning rate, internal AUC and variable importance
distribution.layers.thick <- list() #The spatial likelihood of occurence layers
dist.res <- list()
brt.stack.res <- list()
brt.list.whole.thick <- list() #The whole BRT models

### WHOLE DATA BRT LOOP ###
for (r in 1:3){
  
  internal.stats_res <- data.frame(matrix(NA, nrow=6, ncol=19))
  names(internal.stats_res) <- c("lr", "best.trees", "AUC",  names(whole.master.thick[[1]][[1]][6:21]))
  rownames(internal.stats_res) <- spe.nam
  
  for (s in 1:6){
    
    sp.dat <- whole.master.thick[[r]][[s]]
    
    brt.whole  <- gbm.step(data=sp.dat, 
                           gbm.x = c(6:21),
                           gbm.y = 4,
                           family = "bernoulli",
                           site.weights = sp.dat$site.weight,
                           tree.complexity = 1,
                           learning.rate = 0.05,
                           bag.fraction = 0.75,
                           plot.main=FALSE, silent=TRUE, verbose=FALSE)
    
    
    #internal evaluation
    internal.stats_res[s,1] <- brt.whole$gbm.call$learning.rate
    internal.stats_res[s,2] <- brt.whole$gbm.call$best.trees
    internal.stats_res[s,3] <- brt.whole$cv.statistics$discrimination.mean
    
    match.var <- match(brt.whole$contributions$var, names(internal.stats_res)) 
    
    internal.stats_res[s,match.var] <- brt.whole$contributions$rel.inf
    
    pred.sp <- predict(eval.var.res[[r]], brt.whole, n.trees = internal.stats_res[s,2], type = "response") #making the predictions
    
    dist.res[[s]] <- pred.sp
    
    brt.stack.res[[s]] <- brt.whole
    
    print(paste0("finished.", spe.nam[s], "_", res.names[r]))
  }
  internal.stats.whole.thick[[r]] <-internal.stats_res
  
  names(dist.res) <- spe.nam
  distribution.layers.thick[[r]] <- dist.res
  
  names(brt.stack.res) <- spe.nam
  brt.list.whole.thick[[r]] <- brt.stack.res
}

#saving all of the outputs
names(internal.stats.whole.thick) <- res.names
saveRDS(internal.stats.whole.thick, file= "output/internal.stats.whole.thick.RData")

names(distribution.layers.thick) <- res.names
saveRDS(distribution.layers.thick, file= "output/distribution.layers.thick.RData")

names(brt.list.whole.thick) <- res.names
saveRDS(brt.list.whole.thick, file= "output/brt.list.whole.thick.RData")

write.csv(internal.stats.whole.thick[["250"]], file = "output/var.importance.thick_250.csv")
write.csv(internal.stats.whole.thick[["500"]], file = "output/var.importance.thick_500.csv")
write.csv(internal.stats.whole.thick[["1000"]], file = "output/var.importance.thick_1000.csv")

###############################################################################
### This section contains the same methods above, but with the gridded data ###
###############################################################################
### Reading in the gridded data ###


### THese are the lists of with the gridded species data in spatial-blocks
eval.master.thin <- readRDS("data/data_lists/eval.master.thin.RData")
train.master.thin <- readRDS("data/data_lists/train.master.thin.RData")
whole.master.thin <- readRDS("data/data_lists/whole.master.thin.RData")

# This will read in the exact same vector names and environmental covariates as the ungridded data, they're just repeated for convenience if you're only running one data type
res.names <- c("250", "500", "1000")
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")

var.res <- readRDS("data/var.res.RData") #the environmental covariates as model inputs
eval.var.res <- readRDS("data/eval.var.res.RData") #the environmental covariates with the bias-correction for model prediction



### Gridded "Thin" evaluation models ###

#Setting up the lists for the evaluation stats
internal.stats_thin <- list()
int.stats.sp <- list()
external.preds_thin <- list()
ex.preds.folds <- list()
ex.preds.sp <- list()

#running the brt models per fold
for (r in 1:3) {
  for (s in 1:6) {
    traindat.sp <- train.master.thin[[r]][[s]]
    evaldat.sp <- eval.master.thin[[r]][[s]]
    sp.internal.df <- data.frame(matrix(NA, nrow=6, ncol=3))
    names(sp.internal.df) <- c("lr","best.trees", "AUC")
    for (f in 1:6) {
      folddat.train <- traindat.sp[[f]]
      brt.train  <- gbm.step(data=folddat.train, 
                             gbm.x = c(6:21),
                             gbm.y = 4,
                             family = "bernoulli",
                             site.weights = folddat.train$site.weight,
                             tree.complexity = 1,
                             learning.rate = 0.05,
                             bag.fraction = 0.75,
                             plot.main=FALSE, silent=TRUE, verbose=FALSE)
      
      
      # saving the internal evaluation stats from the training fold
      sp.internal.df[f,1] <- brt.train$gbm.call$learning.rate
      sp.internal.df[f,2] <- brt.train$gbm.call$best.trees
      sp.internal.df[f,3] <- brt.train$cv.statistics$discrimination.mean
      
      #external evaluation to the test fold
      folddat.eval <- evaldat.sp[[f]]
      
      preds <- predict.gbm(brt.train, folddat.eval, n.trees=brt.train$gbm.call$best.trees, type="response")
      pred.dat <- cbind(folddat.eval$presence, preds)
      ex.preds.folds[[f]] <- pred.dat
    }
    #The data from each fold, gets saved into the relavent species list for that resolution
    ex.preds.sp[[s]] <- ex.preds.folds
    int.stats.sp[[s]] <- sp.internal.df
    print(paste0("finished.", spe.nam[s], "_", res.names[r]))
  }
  names(ex.preds.sp) <- spe.nam
  names(int.stats.sp) <- spe.nam
  external.preds_thin[[r]] <- ex.preds.sp
  internal.stats_thin[[r]] <- int.stats.sp
}
names(external.preds_thin) <- res.names
saveRDS(external.preds_thin, file= "output/external.preds_thin.RData")

names(internal.stats_thin) <- res.names
saveRDS(internal.stats_thin, file= "output/internal.stats_thin.RData")



### Thin whole-data models ###

internal.stats.whole.thin <- list()
distribution.layers.thin <- list()
dist.res <- list()
brt.stack.res <- list()
brt.list.whole.thin <- list()
for (r in 1:3){
  
  internal.stats_res <- data.frame(matrix(NA, nrow=6, ncol=19))
  names(internal.stats_res) <- c("lr", "best.trees", "AUC",  names(whole.master.thin[[1]][[1]][6:21]))
  rownames(internal.stats_res) <- spe.nam
  
  for (s in 1:6){
    
    sp.dat <- whole.master.thin[[r]][[s]]
    
    brt.whole  <- gbm.step(data=sp.dat, 
                           gbm.x = c(6:21),
                           gbm.y = 4,
                           family = "bernoulli",
                           site.weights = sp.dat$site.weight,
                           tree.complexity = 1,
                           learning.rate = 0.05,
                           bag.fraction = 0.75,
                           plot.main=FALSE, silent=TRUE, verbose=FALSE)
    
    
    #internal evaluation
    internal.stats_res[s,1] <- brt.whole$gbm.call$learning.rate
    internal.stats_res[s,2] <- brt.whole$gbm.call$best.trees
    internal.stats_res[s,3] <- brt.whole$cv.statistics$discrimination.mean
    
    match.var <- match(brt.whole$contributions$var, names(internal.stats_res)) 
    
    internal.stats_res[s,match.var] <- brt.whole$contributions$rel.inf
    
    pred.sp <- predict(eval.var.res[[r]], brt.whole, n.trees = internal.stats_res[s,2], type = "response")
    
    dist.res[[s]] <- pred.sp
    
    brt.stack.res[[s]] <- brt.whole
    
    print(paste0("finished.", spe.nam[s], "_", res.names[r]))
  }
  internal.stats.whole.thin[[r]] <-internal.stats_res
  
  names(dist.res) <- spe.nam
  distribution.layers.thin[[r]] <- dist.res
  
  names(brt.stack.res) <- spe.nam
  brt.list.whole.thin[[r]] <- brt.stack.res
}
names(internal.stats.whole.thin) <- res.names
saveRDS(internal.stats.whole.thin, file= "output/internal.stats.whole.thin.RData")

names(distribution.layers.thin) <- res.names
saveRDS(distribution.layers.thin, file= "output/distribution.layers.thin.RData")

names(brt.list.whole.thin) <- res.names
saveRDS(brt.list.whole.thin, file= "output/brt.list.whole.thin.RData")


internal.stats.whole.thin <- readRDS("output/internal.stats.whole.thin.RData")

write.csv(internal.stats.whole.thin[["250"]], file = "output/var.importance.thin_250.csv")
write.csv(internal.stats.whole.thin[["500"]], file = "output/var.importance.thin_500.csv")
write.csv(internal.stats.whole.thin[["1000"]], file = "output/var.importance.thin_1000.csv")
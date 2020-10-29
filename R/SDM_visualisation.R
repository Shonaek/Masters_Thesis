library(ggplot2)
library(gridExtra)
library(corrplot)
library(raster)
library(rasterVis)
library(tidyr)
library(dplyr)
library(purrr)

####################
### name vectors ###
####################

res.names <- c("250", "500", "1000")
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")
sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom")
species.plot <- c("Koala","Brushtail", "Kangaroo", "Ringtail", "Wallaby", "Wombat") #this is for tbe plots that don't fit the possum addition
cbPalette <- c("#FFBD00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00")
cb18 <- rep(c("#FFBD00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00"), each= 3)


#########################
### Plotting the AUCs ###
#########################
eval.auc.thick <- readRDS("output/eval.auc.thick.RData")

auc.plot <- rbind(eval.auc.thick[["250"]][,1:8], eval.auc.thick[["500"]][,1:8], eval.auc.thick[["1000"]][,1:8])
auc.plot$species <- rep(species.plot, 3)
auc.plot$res[1:6] <- rep("250", 6)
auc.plot$res[7:12] <- rep("500", 6)
auc.plot$res[13:18] <- rep("1000", 6)
auc.plot$res <- factor(auc.plot$res, levels = c("250", "500", "1000"))

png(file="output/figs/auc.plot.thick.png", res = 200, width = 1800, height = 1100, bg = "white")

ggplot(auc.plot, aes(x= species, y= mean, colour= species)) +
  geom_point(data= auc.plot, aes(x= species, y= fold_1, shape= res), colour= "#CCCCCC", size= 2, position = position_dodge(width = 0.7)) +
  geom_point(data= auc.plot, aes(x= species, y= fold_2, shape= res), colour= "#CCCCCC", size= 2, position = position_dodge(width = 0.7)) +
  geom_point(data= auc.plot, aes(x= species, y= fold_3, shape= res), colour= "#CCCCCC", size= 2, position = position_dodge(width = 0.7)) +
  geom_point(data= auc.plot, aes(x= species, y= fold_4, shape= res), colour= "#CCCCCC", size= 2, position = position_dodge(width = 0.7)) +
  geom_point(data= auc.plot, aes(x= species, y= fold_5, shape= res), colour= "#CCCCCC", size= 2, position = position_dodge(width = 0.7)) +
  geom_point(data= auc.plot, aes(x= species, y= fold_6, shape= res), colour= "#CCCCCC", size= 2, position = position_dodge(width = 0.7)) +
  geom_point(data= auc.plot, aes(shape= res), size= 4, position = position_dodge(width = 0.7)) +
  labs(x= "species", y= "AUC") +
  scale_color_manual(values = c(cbPalette))+
  scale_shape_manual("grain", values= c(19, 15, 17)) +
  theme_classic(base_size = 14)
dev.off()

## To see the difference between 250 and 1000m AUCs
auc_1000_mean <- select(filter(auc.plot, res == "1000"), c(mean, species))
auc_1000_stack <- rbind(auc_1000_mean, auc_1000_mean, auc_1000_mean)

auc.plot %>% 
  mutate("difference" = mean - auc_1000_mean[,1])

#########################################
### spatial difference in predictions ###
#########################################

#plotting the ungridded difference layers
distribution.layers.thick <- readRDS("output/distribution.layers.thick.RData")

difference.thick <- list()
for (s in 1:6) {
  sp.resample <- resample(distribution.layers.thick[["1000"]][[s]], distribution.layers.thick[["250"]][[s]])
  difference.thick[[s]] <- (sp.resample- distribution.layers.thick[["250"]][[s]])
}
names(difference.thick) <- species.plot
diff.thick <- stack(difference.thick)
diff.stack.thick <- subset(diff.thick, order(species.plot))
saveRDS(diff.stack.thick, file= "output/1000-250.stack.thick.RData")

#adjuting the color scale
scale.range <- seq(from= -0.65, to= 0.65, by= 0.1)
col.ramp <- colorRampPalette(c("darkmagenta","white", "orangered"))

#plotting the figure
png(file="output/figs/SDM_1000-250_thick.png", res = 200, width = 1800, height = 900, bg = "white")

levelplot(diff.stack.thick, main= "1000m-250m thick", par.settings = rasterTheme(region = col.ramp(14)),  
          scales=list(y=list(draw=FALSE),x=list(draw=FALSE)), margin= FALSE, border= FALSE, at= scale.range, layout= c(3,2))
dev.off()

round(cellStats(diff.stack.thick, stat= "mean"), 3) #calculate the mean of each layer


################################
### Kulczynski's coefficient ###
################################
KUL.thick <- readRDS("output/KUL.thick.RData")

KUL.thick.all <- cbind(KUL.thick[["btp"]], KUL.thick[["egk"]], KUL.thick[["koa"]], KUL.thick[["rtp"]], KUL.thick[["bsw"]],KUL.thick[["wom"]]) #this reorders the species into the correct order
rownames(KUL.thick.all) <- c("1000", "500", "250")

png("output/figs/kulczynski.thick.png", width = 1200, height = 200, res = 200, pointsize = 6)
corrplot(KUL.thick.all, tl.col = "black", tl.srt = 45, na.label = " ", addCoef.col = "black", number.cex = 0.85, method = "shade", cl.pos = "n")
dev.off()


##################################
### Ranked variable importance ###
##################################

brt.list.whole.thick <- readRDS("output/brt.list.whole.thick.RData")

var.names <- brt.list.whole.thick[[1]][[1]]$var.names
rank.res <- NULL
rank.var <- NULL

for (s in 1:6) {
  for (r in 1:3) {  
    
    var.raw <- data.frame(brt.list.whole.thick[[r]][[s]]$contributions)
    
    rank <- match(var.names, var.raw[,"var"])
    
    rank.res <- cbind(rank.res, rank)
  }
  colnames(rank.res) <- res.names
  rank.df <- as.data.frame(rank.res)
  rank.df$variable <- var.names
  rank.df$species <- species.plot[s]
  rank.var <- rbind(rank.var, rank.df)
  rank.res <- NULL
}

rank.var %>% 
  gather("250", "500", "1000", key = "res", value = "rank" ) -> var.rank.long



var.rank.long$res <- factor(var.rank.long$res, levels = c("1000", "500", "250"))

png(file="output/figs/var.rank.png", res = 200, width = 1800, height = 1100, bg = "white")

ggplot(var.rank.long, aes(y= variable, x= res, fill=rank )) +
  geom_tile() +
  labs(x= "grain (per species)") +
  geom_text(aes(label= rank), colour= "white") +
  scale_fill_gradient(low ="#106003", high = "#83e664") +
  facet_grid(cols = vars(species)) 
dev.off()


#########################################
### Looking at the variable responses ###
#########################################
#Here I am lookin at the species response for the five most important variables across all species

brt.list.whole.thick <- readRDS("output/brt.list.whole.thick.RData")

# tree density
treedens.list <- list()
for (r in 1:3) {
  var.res <- NULL
  for (s in 1:6) {
    model <- brt.list.whole.thick[[r]][[s]]
    values <- plot.gbm(model, i.var="TREEDENS",return.grid=TRUE, type="response")
    colnames(values) <- c("x","y")
    values$name <- rep(paste(species.plot[s]), each=length(values[,2]))
    var.res <- rbind(var.res,values)
    
  }
  treedens.list[[r]] <- var.res
}  
names(treedens.list) <- res.names

png(file="output/figs/treedens.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(treedens.list[["1000"]], aes(x= x, y= y, colour= name)) +
  geom_line(aes(linetype= "1000")) +
  geom_line(data = treedens.list[["500"]], aes(x= x, y= y, colour= name, linetype= "500")) +
  geom_line(data = treedens.list[["250"]], aes(x= x, y= y, colour= name,linetype= "250")) +
  scale_color_manual("species", values = c(cbPalette))+
  scale_linetype_manual("resolution", values = c("1000" = 1, "500" = 2, "250" = 3)) +
  labs(x= "Percent Tree Cover", y= "Likelihood of Occurence") +
  facet_wrap(vars(name)) +
  theme_classic()
dev.off()

# Annual precipitation 
precan.list <- list()
for (r in 1:3) {
  var.res <- NULL
  for (s in 1:6) {
    model <- brt.list.whole.thick[[r]][[s]]
    values <- plot.gbm(model, i.var="PRECAN",return.grid=TRUE, type="response")
    colnames(values) <- c("x","y")
    values$name <- rep(paste(species.plot[s]), each=length(values[,2]))
    var.res <- rbind(var.res,values)
    
  }
  precan.list[[r]] <- var.res
}  
names(precan.list) <- res.names

png(file="output/figs/precan.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(precan.list[["1000"]], aes(x= x, y= y, colour= name)) +
  geom_line(aes(linetype= "1000")) +
  geom_line(data = precan.list[["500"]], aes(x= x, y= y, colour= name,  linetype= "500")) +
  geom_line(data = precan.list[["250"]], aes(x= x, y= y, colour= name,  linetype= "250")) +
  scale_color_manual("species", values = c(cbPalette))+
  scale_linetype_manual("resolution", values = c("1000" = 1, "500" = 2, "250" = 3)) +
  labs(x= "Annual Precipitation", y= "Likelihood of Occurence") +
  facet_wrap(vars(name)) +
  theme_classic()
dev.off()

#mean diurnal temp ranfe
tempmndirng.list <- list()
for (r in 1:3) {
  var.res <- NULL
  for (s in 1:6) {
    model <- brt.list.whole.thick[[r]][[s]]
    values <- plot.gbm(model, i.var="TEMPMNDIRNG",return.grid=TRUE, type="response")
    colnames(values) <- c("x","y")
    values$name <- rep(paste(species.plot[s]), each=length(values[,2]))
    var.res <- rbind(var.res,values)
    
  }
  tempmndirng.list[[r]] <- var.res
}  
names(tempmndirng.list) <- res.names

png(file="output/figs/tempmndirng.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(tempmndirng.list[["1000"]], aes(x= x, y= y, colour= name)) +
  geom_line(aes( linetype= "1000")) +
  geom_line(data = tempmndirng.list[["500"]], aes(x= x, y= y, colour= name,  linetype= "500")) +
  geom_line(data = tempmndirng.list[["250"]], aes(x= x, y= y, colour= name, linetype= "250")) +
  scale_color_manual("species", values = c(cbPalette))+
  scale_linetype_manual("resolution", values = c("1000" = 1, "500" = 2, "250" = 3)) +
  labs(x= "Mean Diurnal Temperature Range", y= "Likelihood of Occurence") +
  facet_wrap(vars(name)) +
  theme_classic()
dev.off()

#soil silt
soilsilt.list <- list()
for (r in 1:3) {
  var.res <- NULL
  for (s in 1:6) {
    model <- brt.list.whole.thick[[r]][[s]]
    values <- plot.gbm(model, i.var="SOILSILT",return.grid=TRUE, type="response")
    colnames(values) <- c("x","y")
    values$name <- rep(paste(species.plot[s]), each=length(values[,2]))
    var.res <- rbind(var.res,values)
    
  }
  soilsilt.list[[r]] <- var.res
}  
names(soilsilt.list) <- res.names

png(file="output/figs/soilsilt.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(soilsilt.list[["1000"]], aes(x= x, y= y, colour= name)) +
  geom_line(aes(linetype= "1000")) +
  geom_line(data = soilsilt.list[["500"]], aes(x= x, y= y, colour= name,  linetype= "500")) +
  geom_line(data = soilsilt.list[["250"]], aes(x= x, y= y, colour= name,  linetype= "250")) +
  scale_color_manual("species", values = c(cbPalette))+
  scale_linetype_manual("resolution", values = c("1000" = 1, "500" = 2, "250" = 3)) +
  labs(x= "Proportion of Silt Particles", y= "Likelihood of Occurence") +
  facet_wrap(vars(name)) +
  theme_classic()
dev.off()

#soil clay
soilclay.list <- list()
for (r in 1:3) {
  var.res <- NULL
  for (s in 1:6) {
    model <- brt.list.whole.thick[[r]][[s]]
    values <- plot.gbm(model, i.var="SOILCLAY",return.grid=TRUE, type="response")
    colnames(values) <- c("x","y")
    values$name <- rep(paste(species.plot[s]), each=length(values[,2]))
    var.res <- rbind(var.res,values)
    
  }
  soilclay.list[[r]] <- var.res
}  
names(soilclay.list) <- res.names

png(file="output/figs/soilclay.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(soilclay.list[["1000"]], aes(x= x, y= y, colour= name)) +
  geom_line(aes(linetype= "1000")) +
  geom_line(data = soilclay.list[["500"]], aes(x= x, y= y, colour= name,  linetype= "500")) +
  geom_line(data = soilclay.list[["250"]], aes(x= x, y= y, colour= name,  linetype= "250")) +
  scale_color_manual("species", values = c(cbPalette))+
  scale_linetype_manual("resolution", values = c("1000" = 1, "500" = 2, "250" = 3)) +
  labs(x= "Proportion of Clay Particles", y= "Likelihood of Occurence") +
  facet_wrap(vars(name)) +
  theme_classic()
dev.off()


##################################################
### Vistualision the species occurence records ###
##################################################
VIC_latlong <- st_read("data/VIC shape/VIC_GDA94LL_ADMIN_STATE.shp")
VIC_UTM <-  st_transform(VIC_latlong, CRS("+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs"))

species.plot <- c("Koala","Brushtail", "Kangaroo", "Ringtail", "Wallaby", "Wombat") 
species.clean <- read.csv("data/species_clean.csv")
scientific.names <- unique(species.clean$Scientific.Name)

species.occ <- NULL
for (s in 1:6) {
sp.sub <- subset(species.clean, species.clean$Scientific.Name == scientific.names[s])  
sp.sub$species <- species.plot[s]
species.occ <- rbind(species.occ, sp.sub)
}
 
png(file="output/figs/sp_records.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot() +
  geom_sf(data = VIC_UTM) +
  geom_point(data= species.occ, aes(x= x, y= y, colour= species), size= 0.7) +
  scale_color_manual(values = c(cbPalette))+
  facet_wrap(vars(species)) +
  theme_classic()
dev.off()


##########################################
### visalising the plain distributions ###
##########################################
distribution.layers.thick <- readRDS("output/distribution.layers.thick.RData") #reading in all of the gridded distributions

sp.distribution.thick_1000 <- stack(distribution.layers.thick[["1000"]]) #subsetting it to the 1000m resolution and renaming it
names(sp.distribution.thick_1000) <- species.plot
sp.distribution.thick_1000 <- subset(sp.distribution.thick_1000, order(species.plot))

terrain <- rasterTheme(region=terrain.colors(n= 20, rev= TRUE))

png(file="output/figs/SDM_distribution_1000_thick.png", res = 200, width = 1800, height = 900, bg = "white")
levelplot(sp.distribution.thick_1000, par.settings= terrain, scales=list(y=list(draw=FALSE), x=list(draw=FALSE)), 
          margin= FALSE, border= FALSE, layout= c(3,2))

dev.off()

###################################
### Thinned model visualisation ###
###################################

### COMPARING THE GRIDDED (THINNED) AND UNGRIDDED AUCS ###
eval.auc.thin <- readRDS("output/eval.auc.thin.RData") #both of these data files were created in the SDM_analysis script
eval.auc.thick <- readRDS("output/eval.auc.thick.RData")

res.names <- c("250", "500", "1000")
sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom")
species.plot <- c("Koala","Brushtail", "Kangaroo", "Ringtail", "Wallaby", "Wombat")

auc.sp.mean.thin <- as.data.frame(cbind(eval.auc.thin[["250"]][,8], eval.auc.thin[["500"]][,8], eval.auc.thin[["1000"]][,8]))
rownames(auc.sp.mean.thin) <- species.plot
names(auc.sp.mean.thin) <- c("250", "500", "1000")

auc.sp.mean.thick <- as.data.frame(cbind(eval.auc.thick[["250"]][,8], eval.auc.thick[["500"]][,8], eval.auc.thick[["1000"]][,8]))
rownames(auc.sp.mean.thick) <- species.plot
names(auc.sp.mean.thick) <- c("250", "500", "1000")

auc.whole <- rbind(auc.sp.mean.thin, auc.sp.mean.thick)
colnames(auc.whole) <- c("res_250", "res_500", "res_1000")
auc.whole$sp.name <- rep(species.plot, 2)
auc.whole$treatment <- c(rep("gridded", 6), rep("ungridded", 6))

png(file="output/figs/auc.plot.thin-thick_comp.png", res = 200, width = 1800, height = 1100, bg = "white")

ggplot(auc.whole, aes(x= sp.name, group = sp.name))+
  geom_point(data= auc.whole, aes(y= res_250, colour= treatment, shape= factor("250")), size= 3, position= position_dodge2(width = 0.9))+
  geom_point(data= auc.whole, aes(y= res_500, colour= treatment, shape= factor("500")), size= 3, position= position_dodge2(width = 0.9))+
  geom_point(data= auc.whole, aes(y= res_1000, colour= treatment, shape= factor("1000")), size= 3, position= position_dodge2(width = 0.5)) +
  scale_shape_manual("grain", values= c("250"= 19, "500"= 15, "1000"= 17)) +
  scale_color_manual(values= c("gridded" = "black", "ungridded" = "#CCCCCC")) +
  labs(x= "species", y= 'mean AUC')+
  theme_classic()
dev.off()

### KULCZYNSKI TABLE ###

KUL.thin <- readRDS("output/KUL.thin.RData")

KUL.thin.all <- cbind(KUL.thin[["btp"]], KUL.thin[["egk"]], KUL.thin[["koa"]], KUL.thin[["rtp"]], KUL.thin[["bsw"]],KUL.thin[["wom"]]) #this reorders the species into the correct order
rownames(KUL.thin.all) <- c("1000", "500", "250")

png("output/figs/kulczynski.thin.png", width = 1200, height = 200, res = 200, pointsize = 6)
corrplot(KUL.thin.all, tl.col = "black", tl.srt = 45, na.label = " ", addCoef.col = "black", number.cex = 0.85, method = "shade", cl.pos = "n")
dev.off()

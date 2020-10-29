library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(corrplot)

sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom")
species.plot <- c("Koala","Brushtail", "Kangaroo", "Ringtail", "Wallaby", "Wombat")
species <- c("Koala","Brushtail_Possum", "Kangaroo", "Ringtail_Possum", "Wallaby", "Wombat")
cbPalette <- c("#FFBD00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00")
res.names <- c("250", "500", "1000")

###########################
### Colllision AUC data ###
###########################
coll.auc_CS <- read.csv("output/collision/case_study/coll.auc_CS.csv") #this is data from the WVC analysis script

coll.glm.eval <- read.csv("output/collision/coll.glm.eval.csv")

coll.glm.eval %>%
  group_by(res) %>%
  group_by(species, add= TRUE) %>% 
  group_nest(.key = "AUC_list") %>% 
  mutate(mean_train_AUC= map_dbl(.$AUC_list, ~mean(.x$AUC_train))) %>% 
  mutate(mean_test_AUC= map_dbl(.$AUC_list, ~mean(.x$AUC_test)))  -> wvc.auc.plot

wvc.auc.plot$species <- factor(wvc.auc.plot$species, level= sp.order)
wvc.auc.plot$res <-  factor(wvc.auc.plot$res, levels = c("250", "500", "1000"))
coll.auc_CS$res <- factor(coll.auc_CS$res, levels = c("250", "500", "1000"))

png(file="output/figs/wvc.auc.plot.thick.png", res = 200, width = 1800, height = 1100, bg = "white")
  
ggplot(wvc.auc.plot, aes(x= species, shape= res)) +
  geom_point(data= wvc.auc.plot, aes(x= species, y= mean_train_AUC, colour= "VIC train"), size= 3, position = position_dodge(width = 0.6)) +
  geom_point(data= wvc.auc.plot, aes(x= species, y= mean_test_AUC, colour= "VIC test"), size= 3, position = position_dodge(width = 0.6)) +
  geom_point(data = coll.auc_CS, aes(x= species, y= AUC_CS, colour= "regional"), size= 3, position = position_dodge(width = 0.6)) +
  scale_shape_manual("grain", values = c("250" = 19, "500" = 15, "1000" = 17)) +
  scale_color_manual("model type", values=  c("VIC test" = "#3182bd", "VIC train" = "#99ccff", "regional"= "#e6550d")) +
  labs(y= "mean AUC") +
  theme_classic()
dev.off()  

#####################
### Density plots ###
#####################

diff.coll <- read.csv("output/collision/collision.density.csv")

png(file="output/figs/coll.density.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(diff.coll, aes(x= diff, colour= area)) +
  geom_density(aes(linetype= type)) +
  scale_color_manual("area", values = c("regional" = "#e6550d", "statewide" = "#3182bd")) + 
  scale_linetype_manual("data type", values =c("collision" = 1, "non-collision" = 5)) +
  facet_grid(vars(species), scales = "free_y") +
  labs(y= "relative density of difference", x= "absolute difference between the 250m and 1000m species occurences") +
  xlim(c(0,0.75)) +
  theme_classic()
dev.off()

##################################################
### Looking at the statewide variable response ###
##################################################
occurence.raw <- read.csv(file= "output/collision/occurence.master.csv")
tvol.raw <- read.csv(file= "output/collision/tvol.master.csv")
tspd.raw <- read.csv(file= "output/collision/tspd.master.csv")

occurence.master <- NULL
tvol.master <- NULL
tspd.master <- NULL

#I need to change the species column to a different name
for (s in 1:6) {
  occ.sub <-  filter(occurence.raw, species == sp.names[s])
  occ.sub$species <- species.plot[s]
  occurence.master <- rbind(occurence.master, occ.sub)
  
  tv.sub <-  filter(tvol.raw, species == sp.names[s])
  tv.sub$species <- species.plot[s]
  tvol.master <- rbind(tvol.master, tv.sub)
  
  ts.sub <-  filter(tspd.raw, species == sp.names[s])
  ts.sub$species <- species.plot[s]
  tspd.master <- rbind(tspd.master, ts.sub)
  }

occurence.master$resolution <- factor(occurence.master$resolution, levels = c("250", "500", "1000"))
tvol.master$resolution <- factor(tvol.master$resolution, levels = c("250", "500", "1000"))
tspd.master$resolution <- factor(tspd.master$resolution, levels = c("250", "500", "1000"))



#species occurence 
png(file="output/figs/species.occ.statewide.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(data = occurence.master, mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
  geom_line(aes(colour= species, linetype= resolution)) +
  geom_ribbon(alpha = 0.3) +
  ylab("") +
  xlab("Likelihood of Species Occurrence") +
  scale_color_manual("species", values = cbPalette) +
  facet_wrap(vars(species), scales = "free_y") +
  theme_classic() +
  theme(legend.key = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")) +
  theme(axis.title.x = element_text(margin = unit(c(0.3, 0, 0, 0), "cm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 0.3, 0 ,0), "cm"))) +
  theme(panel.grid.major = element_line(size = 0.1), panel.grid.minor = element_line(size = 0.1)) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0), lim = c(0, 1))
dev.off()


#Traffic Volume
png(file="output/figs/tvol.statewide.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(tvol.master, aes(x = x / 1000, y = y, ymin = ymin, ymax = ymax)) +
  geom_line(aes(colour= species, linetype= resolution)) +
  geom_ribbon(alpha = 0.3) +
  ylab("") +
  xlab("Traffic Volume (1000 vehicles/day)") +
  scale_color_manual("species", values = cbPalette) +
  facet_wrap(vars(species), scales = "free_y") +
  theme_classic() +
  theme(legend.key = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")) +
  theme(axis.title.x = element_text(margin = unit(c(0.3, 0, 0, 0), "cm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 0.3, 0 ,0), "cm"))) +
  theme(panel.grid.major = element_line(size = 0.1), panel.grid.minor = element_line(size = 0.1)) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0, 40, 5), expand = c(0, 0), lim = c(0, 40))
dev.off()

#Traffic speed
png(file="output/figs/tspd.statewide.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(tspd.master, aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
  geom_line(aes(colour= species, linetype= resolution)) +
  geom_ribbon(alpha = 0.3) +
  ylab("") +
  xlab("Traffic Speed (km/hour)") +
  scale_color_manual("species", values = cbPalette) +
  facet_grid(vars(species), scale= "free_y") +
  theme_classic() +
  theme(legend.key = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")) +
  theme(axis.title.x = element_text(margin = unit(c(0.3, 0, 0, 0), "cm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 0.3, 0 ,0), "cm"))) +
  theme(panel.grid.major = element_line(size = 0.1), panel.grid.minor = element_line(size = 0.1)) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0, 110, 10), expand = c(0, 0), lim = c(0, 110))
dev.off()

###################################
### Correlation of coefficients ###
###################################

coll.cor.raw<- read.csv(file= "output/collision/coll.cor.master.csv")

## Looking at the variable correlation per species and resolution
### I would expect the species-tvol and species=tspd to change, but not the tvol-tspd
#this is not neat, but it works to change the species names
coll.cor.raw$res <- factor(coll.cor.raw$res, levels = c("250", "500", "1000"))
coll.cor.raw %>% 
  gather(key = "comparison", "cor", sp_tvol, sp_tspd, tvol_tspd) -> coll.cor.raw

coll.cor.master <- NULL
for (s in 1:6) {
  sub <-  filter(coll.cor.raw, species == sp.names[s])
  sub$species <- species.plot[s]
  coll.cor.master <- rbind(coll.cor.master, sub)
}

png(file="output/figs/coll.cor.statewide.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(coll.cor.master, aes(x= comparison, y= as.numeric(cor), shape= res, colour= species)) +
  geom_point(position = position_dodge(width = 0.7), size= 3) +
  geom_hline(yintercept= 0, linetype="dashed", color = "grey") +
  labs(x= "parwise comparison of variables", y = "correlation coefficient") +
  scale_color_manual("species", values = cbPalette) +
  scale_shape_manual("resolution", values = c("250" = 17, "500" = 15, "1000" = 19)) +
  facet_wrap(vars(species)) +
  theme_classic() 
dev.off()

###################################
### Regional variable responses ###
###################################

occurence.cs <- read.csv("output/collision/case_study/occurence.master_CS.csv")
tvol.cs <- read.csv("output/collision/case_study/tvol.master_CS.csv")
tspd.cs<- read.csv("output/collision/case_study/tspd.master_CS.csv")

occurence.master_CS <- NULL
tvol.master_CS <- NULL
tspd.master_CS <- NULL

for (s in 1:6) {
  occ.sub <-  filter(occurence.cs, species == sp.names[s])
  occ.sub$species <- species.plot[s]
  occurence.master_CS <- rbind(occurence.master_CS, occ.sub)
  
  tv.sub <-  filter(tvol.cs, species == sp.names[s])
  tv.sub$species <- species.plot[s]
  tvol.master_CS <- rbind(tvol.master_CS, tv.sub)
  
  ts.sub <-  filter(tspd.cs, species == sp.names[s])
  ts.sub$species <- species.plot[s]
  tspd.master_CS <- rbind(tspd.master_CS, ts.sub)
}

occurence.master_CS$resolution <- factor(occurence.master_CS$resolution, levels = c("250", "500", "1000"))

tvol.master_CS$resolution <- factor(tvol.master_CS$resolution, levels = c("250", "500", "1000"))

tspd.master_CS$resolution <- factor(tspd.master_CS$resolution, levels = c("250", "500", "1000"))


#species occurence
png(file="output/figs/species.occ_CS.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(data = occurence.master_CS, mapping = aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
  geom_line(size = 0.2, aes(colour= species, linetype= resolution)) +
  geom_ribbon(alpha = 0.3) +
  ylab("") +
  xlab("Likelihood of Species Occurrence") +
  scale_color_manual("species", values = cbPalette) +
  facet_wrap(vars(species), scales = "free") +
  theme_classic() +
  theme(legend.key = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")) +
  theme(axis.title.x = element_text(margin = unit(c(0.3, 0, 0, 0), "cm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 0.3, 0 ,0), "cm"))) +
  theme(panel.grid.major = element_line(size = 0.1), panel.grid.minor = element_line(size = 0.1)) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0), lim = c(0, 1))
dev.off()

#Traffic Volume
png(file="output/figs/tvol_CS.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(tvol.master_CS, aes(x = x / 1000, y = y, ymin = ymin, ymax = ymax)) +
  geom_line(aes(colour= species, linetype= resolution)) +
  geom_ribbon(alpha = 0.3) +
  ylab("") +
  xlab("Traffic Volume (1000 vehicles/day)") +
  scale_color_manual("species", values = cbPalette) +
  facet_wrap(vars(species), scales = "free") +
  theme_classic() +
  theme(legend.key = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")) +
  theme(axis.title.x = element_text(margin = unit(c(0.3, 0, 0, 0), "cm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 0.3, 0 ,0), "cm"))) +
  theme(panel.grid.major = element_line(size = 0.1), panel.grid.minor = element_line(size = 0.1)) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0, 40, 5), expand = c(0, 0), lim = c(0, 40))
dev.off()

#Traffic speed
png(file="output/figs/tspd_CS.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot(tspd.master_CS, aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
  geom_line(aes(colour= species, linetype= resolution)) +
  geom_ribbon(alpha = 0.3) +
  ylab("") +
  xlab("Traffic Speed (km/hour)") +
  scale_color_manual("species", values = cbPalette) +
  facet_wrap(vars(species), scale= "free_y") +
  theme_classic() +
  theme(legend.key = element_blank()) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")) +
  theme(axis.title.x = element_text(margin = unit(c(0.3, 0, 0, 0), "cm"))) +
  theme(axis.title.y = element_text(margin = unit(c(0, 0.3, 0 ,0), "cm"))) +
  theme(panel.grid.major = element_line(size = 0.1), panel.grid.minor = element_line(size = 0.1)) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0, 110, 10), expand = c(0, 0), lim = c(0, 110))
dev.off()

############################################
### Plotting the distribution of records ###
############################################
### STATEWIDE ###
coll.dat <- readRDS("output/collision/collision.coords.RData") #transforming the collision coordinates list into a csv

coll.dat.sp <- NULL

for (s in 1:6) {
  sp.coll <- coll.dat[[s]]
  sp.coll$species <- species.plot[s]
  coll.dat.sp <- rbind(sp.coll, coll.dat.sp)
}

greater_daylesford <- st_read("data/VIC_shape/Daylesford_solid.shp")
VIC_latlong <- st_read("data/VIC shape/VIC_GDA94LL_ADMIN_STATE.shp")
VIC_UTM <-  st_transform(VIC_latlong, CRS("+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs"))

png(file="output/figs/coll.records.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot() +
  geom_sf(data= VIC_UTM) +
  geom_sf(data = greater_daylesford) +
  geom_point(data= coll.dat.sp, aes(x= x, y= y, colour= species), shape= 4) +
  scale_color_manual(values = c(cbPalette))+
  facet_wrap(vars(species)) +
  theme_classic()
dev.off()

### DAYLESFORD REGION ###

collision.coords_CS <- read.csv("output/collision/case_study/collision.coords_CS.csv")

png(file="output/figs/coll.records_CS.png", res = 200, width = 1800, height = 1100, bg = "white")
ggplot() +
  geom_sf(data = greater_daylesford) +
  geom_point(data= collision.coords_CS, aes(x= x, y= y, colour= species), shape= 4) +
  coord_sf(xlim = c(231982, 428083.7), ylim = c(5835056, 5968199)) +
  scale_color_manual(values = c(cbPalette))+
  facet_wrap(vars(species)) +
  theme_classic()
dev.off()


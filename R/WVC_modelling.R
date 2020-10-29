# This script runs through all of the WVC modelling

##################################
### Name vectors and functions ###
##################################
source("R/support_functions.R")

res.names <- c("250", "500", "1000")
spe.nam <- c("phas.cin","tri.vul","mac.gig","pse.per","wal.bic", "vom.urs")
sp.names <- c("koa", "btp", "egk", "rtp", "bsw", "wom") #three letter abbreviations of the spe.nam vector 
species <- c("Koala","Brushtail_Possum", "Kangaroo", "Ringtail_Possum", "Wallaby", "Wombat")

#############################################################
### Uploading the prediction layers to the spatial server ###
#############################################################

distribution.layers.thick <- readRDS("output/distribution.layers.thick.RData")


for (r in 1:3) {
  for (s in 1:6) {
    
    upload_brt_preds(distribution.layers.thick[[res.names[r]]][[spe.nam[s]]], paste0("vic_gda9455_grid_", sp.names[s], "_preds_brt_", res.names[r]))
  }
}

####################################
### Running the collision models ###
####################################

## This is a test to run one speccies model at one resolution
run_collision_model(species_name = "Kangaroo",
                    table_name = "vic_gda9455_grid_egk_preds_brt_500",
                    n_cores = 1,
                    n_rep= 1)

glm.test <- read.csv("output/collision/glm_metrics_Kangaroo_500.csv")

#This will run the collision models for all species/resolutions in parallel with (n_core * n_rep) repetitions 

for (r in 1:3) {
  for (s in 1:6) {
    
    run_collision_model(species_name = species[s],
                        table_name = paste0("vic_gda9455_grid_", sp.names[s], "_preds_brt_", res.names[r]), 
                        n_cores = 40,
                        n_rep = 25)
    
  }
}

###########################################
### whole data models for visualisation ###
###########################################

#These are all of the response variables
occurence.master <- NULL #species occurence
tvol.master <- NULL      #traffic volume
tspd.master <- NULL      #traffic speed
coll.cor.master <- NULL  # and this is to collect the correlation between covariates for each model

coll.data.res <- list()
coll.data.list <- list() #this is where the glm model is stored

drv <- DBI::dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- DBI::dbConnect(drv,
                      dbname = "qaeco_spatial",
                      user = "qaeco",
                      password = "Qpostgres15",
                      host = "boab.qaeco.com",
                      port = "5432")  #Connection to database server on Boab
on.exit(dbDisconnect(con))


for (r in 1:3) {
  for (s in 1:6) {
    
    
    
    table_name <- paste0("vic_gda9455_grid_", sp.names[s], "_preds_brt_", res.names[r])
    species_name <- species[s]
    
    
    roads <- data.table::as.data.table(dbGetQuery(con,paste0("
  SELECT a.uid, a.length, a.species, b.tvol, c.tspd
  FROM
  (SELECT r.uid AS uid, ST_Length(r.geom)/1000 AS length, sum((st_length(st_intersection(r.geom,g.geom))/st_length(r.geom)) * (g).val) AS species
    FROM gis_victoria.vic_gda9455_roads_state_orig_500 AS r,
    (SELECT (ST_PixelAsPolygons(rast)).val AS val, (ST_PixelAsPolygons(rast)).geom AS geom
    FROM gis_victoria.", table_name, ") AS g
    WHERE ST_Intersects(r.geom,g.geom)
    GROUP BY r.uid) AS a, 
    gis_victoria.vic_nogeom_roads_volpreds_500 AS b, 
    gis_victoria.vic_nogeom_roads_speedpreds_500 AS c
  WHERE
  a.uid = b.uid
  AND
  a.uid = c.uid
    ")))
    
    setkey(roads,uid)
    
    roads$coll <- as.integer(0)
    
    coll <- data.table::as.data.table(dbGetQuery(con, paste0("
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
    
    data <- copy(roads)
    data[coll, coll := i.coll]
    data <- na.omit(data)
    
    coll.data.res[[s]] <- data
    
    coll_glm <- glm(formula = coll ~ log(species) + log(tvol) + I(log(tvol)^2) + log(tspd), offset=log(length), family=binomial(link = "cloglog"), data = data)  #Fit regression model, offset accounts for road lengths
    
    occ_range <- seq(0, 1, 0.001)[-1]
    
    occ_fit <- predict.glm(coll_glm,
                           data.frame(species = occ_range,
                                      tvol = mean(data$tvol),
                                      tspd = mean(data$tspd),
                                      length = 0.5),
                           type = "response",
                           se.fit = TRUE)
    
    occ <- data.frame(x = occ_range,
                      y = occ_fit[["fit"]],
                      ymin = occ_fit[["fit"]] - 1.96 * occ_fit[["se.fit"]],
                      ymax = occ_fit[["fit"]] + 1.96 * occ_fit[["se.fit"]],
                      species = sp.names[s],
                      resolution = res.names[r])
    
    occurence.master <- rbind(occurence.master, occ)
    
    tvol_range <- seq(0, 40000, 100)[-1]
    
    tvol_fit <- predict.glm(coll_glm,
                            data.frame(species = mean(data$species),
                                       tvol = tvol_range,
                                       tspd = mean(data$tspd),
                                       length = 0.5),
                            type = "response",
                            se.fit = TRUE)
    
    tvol <- data.frame(x = tvol_range,
                       y = tvol_fit[["fit"]],
                       ymin = tvol_fit[["fit"]] - 1.96 * tvol_fit[["se.fit"]],
                       ymax = tvol_fit[["fit"]] + 1.96 * tvol_fit[["se.fit"]],
                       species = sp.names[s],
                       resolution = res.names[r])
    
    tvol.master <- rbind(tvol.master, tvol)
    
    tspd_range <- seq(0, 110, 1)[-1]
    
    tspd_fit <- predict.glm(coll_glm,
                            data.frame(species = mean(data$species),
                                       tvol = mean(data$tvol),
                                       tspd = tspd_range,
                                       length = 0.5),
                            type = "response",
                            se.fit = TRUE)
    
    tspd <- data.frame(x = tspd_range,
                       y = tspd_fit[["fit"]],
                       ymin = tspd_fit[["fit"]] - 1.96 * tspd_fit[["se.fit"]],
                       ymax = tspd_fit[["fit"]] + 1.96 * tspd_fit[["se.fit"]],
                       species = sp.names[s],
                       resolution = res.names[r])
    
    tspd.master <- rbind(tspd.master, tspd)
    
    saveRDS(coll_glm, paste0("output/collision/coll_glm_", sp.names[s], "_", res.names[r], ".RData"))
    
    cor <- cor(data[,3:5])
    
    cor.dat <- cor[lower.tri(cor)]
    
    cor.dat[4] <- sp.names[s]
    
    cor.dat[5] <- res.names[r]
    
    cor.dat <- t(cor.dat)
    
    colnames(cor.dat) <- c("sp_tvol", "sp_tspd", "tvol_tspd", "species", "res")
    
    coll.cor.master <- rbind(coll.cor.master, cor.dat)
    print(paste0("finished_", sp.names[s], "_", res.names[r]))
  }
  
  names(coll.data.res) <- sp.names
  coll.data.list[[r]] <- coll.data.res
}
names(coll.data.list) <- res.names

coll.cor.master <- as.data.frame(coll.cor.master)

write.csv(coll.cor.master, file= "output/collision/coll.cor.master.csv")
write.csv(occurence.master, file= "output/collision/occurence.master.csv")
write.csv(tvol.master, file= "output/collision/tvol.master.csv")
write.csv(tspd.master, file= "output/collision/tspd.master.csv")
saveRDS(coll.data.list, file = "output/collision/coll.data.list.RData")

#######################
### Regional Models ###
#######################
#These were only run once per species because there wasn't enough data for cross validation

#see above for explainations of the datafraes
occurence.master_CS <- NULL
tvol.master_CS <- NULL
tspd.master_CS <- NULL
coll.cor.master_CS <- NULL

coll.data.res <- list()
coll.data.list_CS <- list()


drv <- DBI::dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- DBI::dbConnect(drv,
                      dbname = "qaeco_spatial",
                      user = "qaeco",
                      password = "Qpostgres15",
                      host = "boab.qaeco.com",
                      port = "5432")  #Connection to database server on Boab
on.exit(dbDisconnect(con))

for (r in 1:3) {
  for (s in 1:6) {
    
    
    
    table_name <- paste0("vic_gda9455_grid_", sp.names[s], "_preds_brt_", res.names[r])
    species_name <- species[s]
    
    
    roads <- data.table::as.data.table(dbGetQuery(con,paste0("
     SELECT a.uid, a.length, a.species, b.tvol, c.tspd
    FROM
    (SELECT r.uid AS uid, ST_Length(r.geom)/1000 AS length, sum((st_length(st_intersection(r.geom,g.geom))/st_length(r.geom)) * (g).val) AS species
    FROM gis_victoria.vic_gda9455_roads_state_orig_500 AS r,
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
    (SELECT (ST_PixelAsPolygons(rast)).val AS val, (ST_PixelAsPolygons(rast)).geom AS geom
    FROM gis_victoria.", table_name, ") AS g
    WHERE ST_Intersects(r.geom,g.geom)
    AND
    ST_Within(r.geom, z.geom)
    GROUP BY r.uid) AS a, gis_victoria.vic_nogeom_roads_volpreds_500 AS b, gis_victoria.vic_nogeom_roads_speedpreds_500 AS c
  WHERE
  a.uid = b.uid
  AND
  a.uid = c.uid
  
    ")))
    
    
    
    setkey(roads,uid)
    
    roads$coll <- as.integer(0)
    
    coll <- data.table::as.data.table(dbGetQuery(con, paste0("
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
    
    
    data <- copy(roads)
    data[coll, coll := i.coll]
    data <- na.omit(data)
    
    coll.data.res[[s]] <- data
    
    #Running the Collision GLM
    coll_glm <- glm(formula = coll ~ log(species) + log(tvol) + I(log(tvol)^2) + log(tspd), offset=log(length), family=binomial(link = "cloglog"), data = data)  #Fit regression model, offset accounts for road lengths
    
    
    #Then we look at the variable responses
    
    #Collision risk response to species occurence
    occ_range <- seq(0, 1, 0.001)[-1]
    
    occ_fit <- predict.glm(coll_glm,
                           data.frame(species = occ_range,
                                      tvol = mean(data$tvol),
                                      tspd = mean(data$tspd),
                                      length = 0.5),
                           type = "response",
                           se.fit = TRUE)
    
    occ <- data.frame(x = occ_range,
                      y = occ_fit[["fit"]],
                      ymin = occ_fit[["fit"]] - 1.96 * occ_fit[["se.fit"]],
                      ymax = occ_fit[["fit"]] + 1.96 * occ_fit[["se.fit"]],
                      species = sp.names[s],
                      resolution = res.names[r])
    
    occurence.master_CS <- rbind(occurence.master_CS, occ)
    
    ## Traffic Volume
    tvol_range <- seq(0, 40000, 100)[-1]
    
    tvol_fit <- predict.glm(coll_glm,
                            data.frame(species = mean(data$species),
                                       tvol = tvol_range,
                                       tspd = mean(data$tspd),
                                       length = 0.5),
                            type = "response",
                            se.fit = TRUE)
    
    tvol <- data.frame(x = tvol_range,
                       y = tvol_fit[["fit"]],
                       ymin = tvol_fit[["fit"]] - 1.96 * tvol_fit[["se.fit"]],
                       ymax = tvol_fit[["fit"]] + 1.96 * tvol_fit[["se.fit"]],
                       species = sp.names[s],
                       resolution = res.names[r])
    
    tvol.master_CS <- rbind(tvol.master_CS, tvol)
    
    ##Traffic Speed
    tspd_range <- seq(0, 110, 1)[-1]
    
    tspd_fit <- predict.glm(coll_glm,
                            data.frame(species = mean(data$species),
                                       tvol = mean(data$tvol),
                                       tspd = tspd_range,
                                       length = 0.5),
                            type = "response",
                            se.fit = TRUE)
    
    tspd <- data.frame(x = tspd_range,
                       y = tspd_fit[["fit"]],
                       ymin = tspd_fit[["fit"]] - 1.96 * tspd_fit[["se.fit"]],
                       ymax = tspd_fit[["fit"]] + 1.96 * tspd_fit[["se.fit"]],
                       species = sp.names[s],
                       resolution = res.names[r])
    
    tspd.master_CS <- rbind(tspd.master_CS, tspd)
    
    
    saveRDS(coll_glm, paste0("output/collision/case_study/coll_glm_", sp.names[s], "_", res.names[r], "_CS.RData"))
    
    ## Now I want to look at how the correlation of the independant variables at the road segments
    cor <- cor(data[,3:5])
    
    cor.dat <- cor[lower.tri(cor)]
    
    cor.dat[4] <- sp.names[s]
    
    cor.dat[5] <- res.names[r]
    
    cor.dat <- t(cor.dat)
    
    colnames(cor.dat) <- c("sp_tvol", "sp_tspd", "tvol_tspd", "species", "res")
    
    coll.cor.master_CS <- rbind(coll.cor.master_CS, cor.dat)
    
  }
  
  names(coll.data.res) <- sp.names
  coll.data.list_CS[[r]] <- coll.data.res  
}

names(coll.data.list_CS) <- res.names

coll.cor.master_CS <- as.data.frame(coll.cor.master_CS)

write.csv(coll.cor.master_CS, file= "output/collision/case_study/coll.cor.master_CS.csv")
write.csv(occurence.master_CS, file= "output/collision/case_study/occurence.master_CS.csv")
write.csv(tvol.master_CS, file= "output/collision/case_study/tvol.master_CS.csv")
write.csv(tspd.master_CS, file= "output/collision/case_study/tspd.master_CS.csv")
saveRDS(coll.data.list_CS, file= "output/collision/case_study/coll.data.list_CS.RData")
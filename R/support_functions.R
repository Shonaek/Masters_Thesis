#These are the functions to run the collision model script (see WVC_modelling)

require(raster)
require(RPostgreSQL)
require(gbm)
require(dismo)
require(data.table)
require(ggplot2)
require(foreach)
require(doMC)
require(rpostgis)

roc <- function (obsdat, preddat){
  if (length(obsdat) != length(preddat)) 
    stop("obs and preds must be equal lengths")
  n.x <- length(obsdat[obsdat == 0])
  n.y <- length(obsdat[obsdat == 1])
  xy <- c(preddat[obsdat == 0], preddat[obsdat == 1])
  rnk <- rank(xy)
  roc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * n.y)
  return(roc)
}

invcloglog <- function (x) {1-exp(-exp(x))}


# Standard convention for naming tables on the spatial data server is:
# statecode_projection_geometrytype_speciescode_values_modeltype_resolution - e.g. vic_gda9455_grid_egk_preds_brt_500

upload_brt_preds <- function (raster_layer, table_name, write_local = FALSE, file = NULL) {
  drv <- DBI::dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- DBI::dbConnect(drv,
                        dbname = "qaeco_spatial",
                        user = "qaeco",
                        password = "Qpostgres15",
                        host = "boab.qaeco.com",
                        port = "5432")  #Connection to database server on Boab
  on.exit(dbDisconnect(con))
  
  
  raster::writeRaster(raster_layer, filename = "preds.tif")
  
  if(write_local) {
    if (is.null(file)) stop("You must supply a filename/path if writing raster predictions locally")
    raster::writeRaster(raster_layer, filename = file)
  }
  
  system(paste0("raster2pgsql -d -C -I -M -s 28355 -t auto ", getwd(), "/preds.tif gis_victoria.", table_name," | PGPASSWORD=Qpostgres15 psql -d qaeco_spatial -h boab.qaeco.com -p 5432 -U qaeco -w"))
  
  unlink("preds.tif")
  
}


# EXAMPLE USE
# upload_brt_preds(model = kang_brt,
#                  vars = env_vars,
#                  table_name = "vic_gda9455_grid_egk_preds_brt_500")

run_collision_model <- function (species_name, table_name, n_cores, n_rep) {
  drv <- DBI::dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con_coll <- DBI::dbConnect(drv,
                        dbname = "qaeco_spatial",
                        user = "qaeco",
                        password = "Qpostgres15",
                        host = "boab.qaeco.com",
                        port = "5432")  #Connection to database server on Boab
  
  
  
  roads <- data.table::as.data.table(dbGetQuery(con_coll,paste0("
  SELECT a.uid, a.length, a.species, b.tvol, c.tspd
  FROM
  (SELECT r.uid AS uid, ST_Length(r.geom)/1000 AS length, sum((st_length(st_intersection(r.geom,g.geom))/st_length(r.geom)) * (g).val) AS species
    FROM gis_victoria.vic_gda9455_roads_state_orig_500 AS r,
    (SELECT (ST_PixelAsPolygons(rast)).val AS val, (ST_PixelAsPolygons(rast)).geom AS geom
    FROM gis_victoria.", table_name, ") AS g
    WHERE ST_Intersects(r.geom,g.geom)
    GROUP BY r.uid) AS a, gis_victoria.vic_nogeom_roads_volpreds_500 AS b, gis_victoria.vic_nogeom_roads_speedpreds_500 AS c
  WHERE
  a.uid = b.uid
  AND
  a.uid = c.uid
    ")))
  setkey(roads,uid)
  
  roads$coll <- as.integer(0)
  
  coll <- data.table::as.data.table(dbGetQuery(con_coll, paste0("
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
  
  
  resolution <- dbGetQuery(con_coll, paste0("
                                     SELECT 
                                    ST_PixelHeight(rast)
                                     FROM
                                     gis_victoria.", table_name))[1,]
  
  
  n_training_years <- 7
  years <- unique(coll$year)
  
  
  output_sp <- foreach(j = seq_len(n_rep), .combine = rbind) %do%{
    
    registerDoMC(cores = n_cores)
    
    output_table <- foreach(i = seq_len(n_cores), .combine= rbind, .packages = "RPostgreSQL") %dopar% {
      
      drv <- DBI::dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
      con_par <- DBI::dbConnect(drv,
                                 dbname = "qaeco_spatial",
                                 user = "qaeco",
                                 password = "Qpostgres15",
                                 host = "boab.qaeco.com",
                                 port = "5432")  #Connection to database server on Boab
      
      
      set.seed(j+i^5)
      
      
      years.train <- sample(2010:2018, n_training_years)
      years.test <- years[!(years %in% years.train)]
      coll.train <- coll[year %in% years.train,]
      coll.test <- coll[year %in% years.test,]
      
      data <- copy(roads)
      data[coll.train, coll := i.coll]
      data <- na.omit(data)
      
      val_data <- copy(roads)
      val_data[coll.test, coll := i.coll]
      val_data <- na.omit(val_data)

      
      
      # print(paste0("Fitting GLM model to ",
      #              length(data$coll[data$coll == 1]),
      #              " collisions and ",
      #              length(data$coll[data$coll == 0]),
      #              " inferred non-collisions (background)..."))
      
      
      
      coll_glm <- glm(formula = coll ~ log(species) + log(tvol) + I(log(tvol)^2) + log(tspd), offset=log(length), family=binomial(link = "cloglog"), data = data)  #Fit regression model, offset accounts for road length and years of data
      
      saveRDS(coll_glm, file= paste0("output/collision/coll_glm_", species_name,"_", resolution, ".RData"))
      
      #print(summary(coll_glm))
      
      #print(paste0("% Deviance Explained: ", round(((coll_glm$null.deviance - coll_glm$deviance) / coll_glm$null.deviance) * 100, 2)))
      
      
      preds <- predict(coll_glm, data, type="response")
      
      val_coll <- data.table::as.data.table(dbGetQuery(con_par, paste0("
                                         SELECT DISTINCT ON (p.id)
                                         r.uid AS uid, CAST(1 AS INTEGER) AS coll
                                         FROM
                                         gis_victoria.vic_gda9455_roads_state_orig_500 as r,
                                         (SELECT DISTINCT ON (geom)
                                         id, geom
                                         FROM
                                         gis_victoria.vic_gda9455_fauna_wv_2014_2018_coll2
                                         WHERE
                                         species = '", species_name ,"') AS p
                                         WHERE ST_DWithin(p.geom,r.geom, 30)
                                         ORDER BY p.id, ST_Distance(p.geom,r.geom)
                                         ")))
      setkey(val_coll, uid)

      
      # print(paste0("Validating GLM model with ",
      #              length(val_data$coll[val_data$coll==1]),
      #              " independent collisions and ",
      #              length(val_data$coll[val_data$coll==0]),
      #              " inferred non-collisions (background)..."))
      
      
      # print(paste0("Validation data AUC: ", round(roc(val_data$coll, predict(coll_glm, val_data, type="response")), 2)))
      
      
      #print(paste0("finished_", species_name, "_iteration_", i))
      dbDisconnect(con_par) 
      
      data.frame("Iteration"= paste0(j, "-", i),
                "N_train_Col"= length(data$coll[data$coll == 1]),
                 "N_train_NoCol"= length(data$coll[data$coll == 0]),
                 "N_test_Col"= length(val_data$coll[val_data$coll == 1]),
                 "N_test_NoCol"= length(val_data$coll[val_data$coll == 0]),
                 "Reduced_Deviance"=  round(((coll_glm$null.deviance - coll_glm$deviance) / coll_glm$null.deviance) * 100, 2),
                 "AUC_train"= roc(data$coll, predict(coll_glm, data, type="response")),
                 "AUC_test" = roc(val_data$coll, predict(coll_glm, val_data, type="response")))
    
       
      }
    
    dbDisconnect(con_coll)
    
    output_table
    
  }

  
  
  write.csv(output_sp, file = paste0("output/collision/glm_metrics_", species_name,"_", resolution, ".csv"))
  
  #assign(paste0("test.glm.", species_name,"_", resolution), output_table)
  
  write.csv(signif(summary(coll.glm)$coefficients, digits=4),"output/vic_coll_coef_500.csv",row.names=FALSE)
  write.csv(formatC(anova(coll.glm)[2:5,2]/sum(anova(coll.glm)[2:5,2]), format='f',digits=4),"output/vic_coll_anova_500.csv",row.names=FALSE)

}
# EXAMPLE USE:
# run_collision_model(species_name = "Kangaroo",
#                     table_name = "vic_gda9455_grid_egk_preds_brt_500")

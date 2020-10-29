# Masters_Thesis explaination of R scripts

1- data_exploration_cleaning
-	Exploring and cleaning the raw species occurrence data
-	Preparing the environmental covariates for the SDMs (training and bias-correcting stacks)
-	Creating the SDM background points
-	Investigating the collision data
-	Visualising the study site 

2- SDM_data_lists
-	Creates the spatial blocking structure and extracts the covariate information into lists
-	Prepares covariate information into lists for the whole-data SDMs 
-	Repeated for ungridded and gridded data

3- BRT_loop
-	Runs the spatially-blocked BRT models for evaluation of performance 
-	Runs the whole-data BRT models spatial predictions
-	Repeated for gridded models

4- -SDM_analysis
-	AUC evaluation (repeated for gridded and ungridded models)
-	Kulczynski’s coefficient calculations (for gridded and ungridded models)
-	Difference between the 1000m and 250m spatial predictions (gridded model only)

5- SDM_visualisation
Code to reproduce all SDM figures in the main text and appendices 

6- support_functions
All functions necessary for running the WVC models

7- WVC_modelling
-	Statewide collision models with cross-validation
-	Whole-data statewide models for collision predictions
-	Greater Daylesford region collision models

8- WVC_analysis
-	AUC values
-	Extracting the density of collision and non-collision roads from the difference in species prediction rasters
-	Standardising and exporting the collision predictions (for visualisation in QGIS)
-	Extracting the model parameter information

9- WVC_visualisation
Code for making all figures (except the maps of collision predictions) in the main text and appendices

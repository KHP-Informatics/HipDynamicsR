

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   _______ __         _____                              __                 +
#  |   |   |__|.-----.|     \.--.--.-----.---.-.--------.|__|.----.-----.    +
#  |       |  ||  _  ||  --  |  |  |     |  _  |        ||  ||  __|__ --|    +
#  |___|___|__||   __||_____/|___  |__|__|___._|__|__|__||__||____|_____|    +
#              |__|          |_____|                                         +
#                                                                            +
# TITLE:   HipDynamics - An analysis to deduce cell population dynamics      +
# VERSION: 1.1                                                               +
# AUTHOR:  Maximilian Kerz (kerz.maximilian@gmail.com)            	     +
#                                                                            +
# ACKNOWLEDGEMENTS: Amos Folarin (amosfolarin@gmail.com)                     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                  
library("gplots")
library("ggplot2")
library("ggdendro")
library("reshape")
library("RMySQL")

# +++++++++++++++++++++++++++++++ Settings +++++++++++++++++++++++++++++++++++

# mySQL db admin
db_host <- "localhost"
db_name <- "hipsci_final"
db_usern <- "hipsci"
db_passwd <- "data"
table_per_img <- "ExpFeeder_1to16_Per_Image"
table_per_obj <- "ExpFeeder_1to16_Per_Object"

# PlateResult path and filename - plate results must be supplied and the path 
#             adjusted accordingly. 
path_PR <- "~/Documents/projects/HIPSCI/Plate_Results/" 
path_out <- paste(getwd(), "/output/", sep="")
path_code <- paste(getwd(), "/.hipdynamics/", sep="")

# FN concentrations
fn_concentrations <- c(0,1,5,25)
# thresholds
threshold <- c(50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,
               800,850,900,950,1000,1050,1100,1150,1200)

# ----------------------------------- OPTIONS

# method for inferential analysis
infer_from_raw_bins <- TRUE
infer_from_normalisation <- FALSE

# export options
export_raw_data <- FALSE
export_normalised_data <- FALSE

# apply threshold on raw data
threshold_raw <- 50
apply_threshold_on_raw <- FALSE

# ++++++++++++++++++++++++++++++ Statics +++++++++++++++++++++++++++++++++++++
source(paste(path_code, "hipDynamics_statics.R", sep=""))

# +++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++
source(paste(path_code,"hipDynamics_fun.R", sep=''))

# +++++++++++++++++++++++++++++++++ SCRIPT +++++++++++++++++++++++++++++++++++
print(paste(aPrompt, "LOADING DATA..."))
source(paste(path_code, "hipDynamics_staging.R", sep=""))

print(paste(aPrompt, "BEGINNING ANALYSIS..."))
source(paste(path_code, "hipDynamics_analysis.R", sep=""))

print(paste(aPrompt, "ANALYSIS COMPLETE."))
summary()

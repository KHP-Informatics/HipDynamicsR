# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   _______ __         _____                              __                 +
#  |   |   |__|.-----.|     \.--.--.-----.---.-.--------.|__|.----.-----.    +
#  |       |  ||  _  ||  --  |  |  |     |  _  |        ||  ||  __|__ --|    +
#  |___|___|__||   __||_____/|___  |__|__|___._|__|__|__||__||____|_____|    +
#              |__|          |_____|                                         +
#                                                                            +
# TITLE:   HipDynamics - Making high-content, live- and end-point imaging    +
#                        comparable.                                         +
# VERSION: 2.0                                                               +
# AUTHOR:  Maximilian Kerz (kerz.maximilian@gmail.com)                       +
#                                                                            +
# ACKNOWLEDGEMENTS: Amos Folarin (amosfolarin@gmail.com)                     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                  

# ----------------------------------- OPTIONS ----------------------------------

# analyse within experiment per line accoring to well or line?
analysis_selector <- "FN" # Select from "Line" or "Well"
bin_number <- 20

# set working directory
setwd("/Path/to/HipDynamics")

# mySQL db admin
db_host <- "localhost"
db_name <- "example_data"
db_usern <- "hipsci"
db_passwd <- "data"
db_image_table <- "Example_Img"
db_object_table <- "Example_Obj"
db_column_selector <- "CellsFinal"

# contrinuous Database queryies during analysis
query = TRUE
# normalised heatmaps
normalise = TRUE

# --------------------------------- OPTIONS END  -------------------------------

# PlateResult path and filename - plate results must be supplied and the path 
#             adjusted accordingly. 
path_PR <- paste(getwd(), "/plateResults/", sep="")
path_out <- paste(getwd(), "/output/", sep="")
path_code <- paste(getwd(), "/hipdynamics/", sep="")

# thresholds
threshold_raw <- 50

library("gplots")
library("ggplot2")
library("ggdendro")
library("reshape")
library("RMySQL")

# ++++++++++++++++++++++++++++++ Statics +++++++++++++++++++++++++++++++++++++
source(paste(path_code, "statics.R", sep=""))
# +++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++
source(paste(path_code,"fun.R", sep=''))
# +++++++++++++++++++++++++++++++++ SCRIPT +++++++++++++++++++++++++++++++++++
print(paste(aPrompt, "Staging ..."))
source(paste(path_code, "staging.R", sep=""))

print(paste(aPrompt, "Starting analysis ..."))
source(paste(path_code, "analysis.R", sep=""))

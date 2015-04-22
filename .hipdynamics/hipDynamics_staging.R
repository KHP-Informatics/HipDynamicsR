
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   _______ __         _____                              __                 +
#  |   |   |__|.-----.|     \.--.--.-----.---.-.--------.|__|.----.-----.    +
#  |       |  ||  _  ||  --  |  |  |     |  _  |        ||  ||  __|__ --|    +
#  |___|___|__||   __||_____/|___  |__|__|___._|__|__|__||__||____|_____|    +
#              |__|          |_____|                                         +
#                                                                            +
# TITLE:   HipDynamics - An analysis to deduce cell population dynamics    	 +
# VERSION: 1.0                                                               +
# AUTHOR:  Maximilian Kerz (kerz.maximilian@gmail.com)             			     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if(!exists("aPrompt"))
{
  print("ANALYSER: ERROR (Make sure you load hipDynamics.R first).")
  quit(save="no")
}

# +++++++++++++++++++++++++++++++ SCRIPT ++++++++++++++++++++++++++++++++++++++

# make sure the mysql.server is running
# establishing connection and extracting tables
drv <- dbDriver("MySQL")
con <- dbConnect(drv, db_name, user='hipsci', pass='data')

statement <- paste('Select', 
                   img_no_idx, ",", 
                   file_idx, ",", 
                   well_idx, 
                   "from", table_per_img)
per.img <- dbGetQuery(con, statement=statement)

statement <- paste('Select', 
                   obj_img_no_idx, ",", 
                   obj_area_idx, ",", 
                   obj_compact_idx, ",", 
                   obj_eccent_idx, ",", 
                   obj_euler_idx, ",", 
                   obj_formFactor_idx, ",", 
                   obj_majorAxisLength_idx, ",", 
                   obj_maxRad_idx, ",", 
                   obj_meanRad_idx, ",", 
                   obj_medRad_idx, ",", 
                   obj_minAxisLength_idx, ",", 
                   obj_perimeter_idx,
                   "from", table_per_obj)
per.obj <- dbGetQuery(con, statement=statement)

print(paste(aPrompt, "Connected to", db_name, 
            "and successfully loaded tables."))

# adding 3 columns to list
headers <- c("Order_ID", "FN_conc", "Cell_Line")
for(i in 1:length(headers)){
  per.img[[length(per.img)+1]] <- NA
  head_names <- names(per.img)
  head_names[length(per.img)] <- headers[i]
  colnames(per.img) <- head_names
}

# adding order_idx to [,242]
order_id <- 1
well_id_before <- ""
for(i in 1:length(per.img[,1])){
  well_id <- as.character(per.img[i,well_idx])
  # order_id is used in conjungtion with Hour where the first hour is
  # 0 in order_id
  if (well_id != well_id_before){
    well_id_before <- well_id
    order_id <- 1
  } else {
    order_id <- order_id + 1
  }
  per.img[i, order_idx] <- order_id
}

print(paste(aPrompt, "Creating experiments img- and obj-tables..."))


# experiment with cell line SGJ from hipsci_db.ExpSGJ_Per_Image.
if (table_per_img=="_ExpSGJ_Per_Image"){
  PRDictSGJ <- c("A4:1:SJG", "A7:5:SJG", "A10:25:SJG")
  AddPR(PRDictSGJ)
  expSGJ.img <- per.img
  expSGJ.obj <- per.obj
} else {
  # all experiments in hipsci_db.ExpFeeder_Per_Image of 
  #   which Plate_Results.txt exist.
  ImportAllExperiments()
}

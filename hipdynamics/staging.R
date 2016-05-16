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

# Prompt
aPrompt <- "[HipDynamics]"

if(!exists("aPrompt"))
{
    print("ANALYSER: ERROR (Make sure you load hipDynamics.R first).")
    quit(save="no")
}

# +++++++++++++++++++++++++++++++ SCRIPT ++++++++++++++++++++++++++++++++++++++

# make sure the mysql.server is running
# establishing connection and extracting filesNames and ImaegNumber
print(paste(aPrompt, "Connecting to", db_name, "..."))
drv <- dbDriver("MySQL")
con <- dbConnect(drv, db_name, user='hipsci', pass='data', host=db_host)
statement <- paste('Select', 
                   mysql.imageNumber, ",", 
                   mysql.fileName,
                   "from", db_image_table)
db.overview <- dbGetQuery(con, statement=statement)
# close all db connections
lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)

print(paste(aPrompt, "Creating Database look-up table ..."))
# create lookUp table from filenames and imagenumber
db.lookUp <- createLookUpTable(db.overview)

print(paste(aPrompt, "Creating Plate Results look-up table ..."))
# query plate results files
pr.files <- dir(path_PR)
pr.lookUp <- lapply(pr.files, createPlateResultLookUpTable)

print(paste(aPrompt, "Annotating Database look-up table ..."))
# checking for match between ExpDate and Date from db and pr look-up tables
expDate.unique <- unique(db.lookUp$ExpDate)
for(i in 1:length(expDate.unique)){
    for(j in 1:length(pr.lookUp)){
        
        pr.lookUp.date <- as.character(factor(pr.lookUp[[j]]$Date[1]))
        if(as.character(expDate.unique[i]) == pr.lookUp.date){
            print(paste(aPrompt, "Match for", pr.lookUp[[j]]$Experiment[1], "found."))
            db.lookUp <- annotateDBLookUpTable(db.lookUp, pr.lookUp[[j]])
        }
        
    }
}
# removing non-annotated rows
lookUp <- db.lookUp
if(length(which(db.lookUp$Experiment == "-")) > 0){
    lookUp <- db.lookUp[-which(db.lookUp$Experiment == "-"),]   
}
if(length(which(db.lookUp$Line == "-")) > 0){
    lookUp <- db.lookUp[-which(db.lookUp$Line == "-"),]   
}
# remove duplicate hours
for(i in 2:length(lookUp[,1])){
    dup.len <- length(lookUp[,1])
    if(i < dup.len){
        if(lookUp$Hour[(i-1)] == lookUp$Hour[i]){
            lookUp <- lookUp[-i,]
        }   
    }
}
#cleaning memory
#rm(list=c("pr.lookUp.date", "pr.lookUp", "db.lookUp", "db.overview", "pr.files"))

### Database query
# initiate database connection
drv <- dbDriver("MySQL")
con <- dbConnect(drv, db_name, user='hipsci', pass='data', host=db_host)
# query database names
statement <- paste("show columns from",
                   db_object_table, ";")
db.cols <- dbGetQuery(con, statement=statement)
# select columns names containing "CellsDoubleFinal2"
db.cols.select <- db.cols$Field[grepl(db_column_selector, db.cols$Field)]
# collapse vector into comma delimited string
db.cols.str <- paste(db.cols.select, collapse=", ")
# collapse vector into max() query
db.cols.max.query <- paste("select max(", 
                           paste(db.cols.select, collapse="), max("),
                           ") from ", db_object_table, ";")
# query max values of data set
print(paste(aPrompt, "Querying max object measures ..."))
db.cols.max <- dbGetQuery(con, statement=db.cols.max.query)
# collapse vector into min() query
db.cols.min.query <- paste("select min(", 
                           paste(db.cols.select, collapse="), min("),
                           ") from ", db_object_table, ";")
# query max values of data set
print(paste(aPrompt, "Querying min object measures ..."))
db.cols.min <- dbGetQuery(con, statement=db.cols.min.query)
# close all db connections
lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)


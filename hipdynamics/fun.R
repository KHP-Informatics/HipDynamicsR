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


# ++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++

# decomposes filename into a vector and adds the ImageNumber to it
#  20150630_B8_1_2015y07m01d_08h00m.tif -> "1" "20150630" "B8" "1" "2015y07m01d" "08"  
fileNameDecomposition <- function(file){
    idx <- file[1]
    str.v <- strsplit(as.character(file[2]),"_")
    # extract hour
    str.hr <- strsplit(as.character(str.v[[1]][length(str.v[[1]])]), "h")
    str.v[[1]][length(str.v[[1]])] <- str.hr[[1]][1]
    v <- c(idx, str.v[[1]])
    return(as.vector(unlist(v)))
}

# uses fileNameDecomposition function to create a lookup-table for all images
createLookUpTable <- function(overview){
    decomp <- fileNameDecomposition(overview[1,])
    decomp.ext <- c(decomp, "-", "-", "-")
    decomp.mtrx <- t(as.matrix(decomp.ext))
    colnames(decomp.mtrx) <- c("ImageNumber", "ExpDate", "Well", "Plate", "Date", "Hour", "Experiment", "Line", "FN")
    for(i in 2:length(overview[,1])){
        decomp.mtrx <- rbind(decomp.mtrx, c(fileNameDecomposition(overview[i,]), "-", "-", "-"))
    }
    # transforms mtrx into df
    decomp.df <- data.frame(decomp.mtrx)
    return(decomp.df)
}

# extract well, [FN] and cell line from plate results row
plateResultExtract <- function(row){
    rows <- c("A", "B", "C", "D", "E", "F", "G", "H")
    well_id <- paste(rows[as.integer(row$V1)], 
                     row$V2, sep="")
    fn_conc <- as.character(row$V60)
    cell_line <- as.character(factor(row$V61))
    # if less fields are present
    if(is.na(fn_conc)){
        fn_conc <- as.character(row$V11)
        cell_line <- as.character(factor(row$V10))
    }
    return(c(well_id, cell_line, fn_conc))
}

# based on the 'date' in filename of plate_result.txt 
#     corresponding experments are identified.
createPlateResultLookUpTable <- function(plate_results_txt){
    sub_name <- strsplit(plate_results_txt, "_")
    disjoint_date <- strsplit(sub_name[[1]][1], "-")
    date <- paste(disjoint_date[[1]][1], disjoint_date[[1]][2], 
                  disjoint_date[[1]][3], sep="")
    expId.idx <- which(regexpr("Exp", sub_name[[1]])>0)
    expId <- sub_name[[1]][expId.idx]
    
    # create plate results look up table
    file <- read.table(paste(path_PR, plate_results_txt, sep=""), 
                       header=FALSE, skip=9, sep="\t")  
    extract <- plateResultExtract(file[1,])
    extract.mtrx <- t(matrix(extract))
    colnames(extract.mtrx) <- c("Well", "Line", "FN")
    for(i in 2:length(file[,1])){
        extract.mtrx <- rbind(extract.mtrx, plateResultExtract(file[i,]))
    }
    
    file.len <- length(file[,1])
    Date <- rep(date, file.len)
    Experiment <- rep(expId, file.len)
    
    extract.mtrx <- cbind(Experiment, extract.mtrx)
    extract.mtrx <- cbind(Date, extract.mtrx)
    return(data.frame(extract.mtrx))
}

# annotating db look-up table with pr look up table information
annotateDBLookUpTable <- function(dblu, prlu){
    dblu$Experiment <- as.vector(dblu$Experiment)
    dblu$Line <- as.vector(dblu$Line)
    dblu$FN <- as.vector(dblu$FN)
    # Experiment
    prlu.date <- as.character(factor(prlu$Date[1]))
    expIdx <- which(dblu$ExpDate == prlu.date)
    experiment <- as.character(as.vector(prlu$Experiment[1]))
    dblu[expIdx, "Experiment"] <- experiment
    # Wells
    wells <- as.vector(unique(prlu$Well))
    for(k in 1:length(wells)){
        dbWellIdx <- which(dblu$Well[expIdx] == wells[k])
        prWellIdx <- which(prlu$Well == wells[k])
        dblu$Line[expIdx[dbWellIdx]] <- as.vector(prlu$Line[prWellIdx])
        dblu$FN[expIdx[dbWellIdx]] <- as.vector(prlu$FN[prWellIdx])
    }
    return(dblu)
}

ApplyThreshold <- function(exp_field,
                           bottom_thresh,
                           upper_thresh=0,
                           idx=FALSE){
    # BOOL 'idx' allows export of indexes rather than real values
    
    thresh_idx1 <- which(exp_field > bottom_thresh)
    thresh_final_idx <- thresh_idx1
    if(upper_thresh > bottom_thresh){
        thresh_idx2 <- which(exp_field[thresh_idx1] <= upper_thresh)
        thresh_final_idx <- thresh_idx1[thresh_idx2]
    }
    # BOOL 'idx' allows export of indexes rather than real values
    thresh_return <- exp_field[thresh_final_idx]
    if(idx==TRUE){
        thresh_return <- thresh_final_idx
    }
    return(thresh_return)
}

ApplyMultiThreshold <- function(exp_field,
                                thresh_list){
    
    # thresh_list convention as follows: thresh_list <- c(thresh1, thresh2, threshN)
    # 	where thresh2 acts as bottom threshold in 2nd cylce of For loop, and so on.
    bin_no <- length(thresh_list)-1
    bins <- list()
    for(i in 1:bin_no){
        bin <- ApplyThreshold(exp_field, thresh_list[i], thresh_list[i+1])
        if(length(bin)==0){
            bin <- NA
        }
        bins[[i]] <- bin
    }
    return(bins)
}

HistoOfBins <- function(bins){
    
    bin_no <- length(bins)
    histo <- c()
    for(i in 1:bin_no){
        histo[i] <- length(bins[[i]])
        if(is.na(bins[[i]][1])){
            histo[i] <- 0
        }
    }
    return(histo)
}

GetThresholdLabel <- function(threshold_in){
    thresh_label <- c()
    for(i in 1:length(threshold_in)-1){
        thresh_label[i] <- paste(threshold_in[i],"-",threshold_in[i+1])
    }
    return(thresh_label)
}

Normalise <- function(x){
    ret <- (x-min(x))/(max(x)-min(x))
    for(i in 1:length(ret)){
        if(is.nan(ret[i])){ret[i] <- 0}
    }
    return(ret)
}

PlotHeatHisto <- function(histo_input,
                         fn_conc,
                         experiment,
                         cell_line,
                         xlab="Time (in Hours)",
                         ylab="Cell Area",
                         main="",
                         overwrite_main=FALSE,
                         cluster=FALSE,
                         all_FNs=FALSE,
                         in_percent=FALSE){
    
    xlab <- "Time (in Hours)"
    ylab <- "Cell Area"
    if(in_percent == TRUE){
        ylab <- "Cell Area (in %)"
    }
    if(overwrite_main == FALSE){
        main <- paste(experiment, "-", cell_line,
                      ": Heatmap of cell area at",sep="")
        if(all_FNs == TRUE){
            main <- paste(main, "all [FN]'s")
        } else {
            main <- paste(main, "[FN] =", fn_conc)
        }
    }
    
    e <- melt(histo_input)
    e <- rename(e,c(X1="CellArea", X2="Hour"))
    e$CellArea <- with(e, reorder(CellArea))
    
    p <- ggplot(e, aes(Hour,CellArea))
    p <- p + geom_tile(aes(fill=value)) + scale_fill_gradient(low="white",
                                                              high="steelblue") + xlab(xlab) + ylab(ylab) + ggtitle(main)
    print(p)
}

buildHeatMatrix <- function(m1, m2){
    m1.ncol <- length(m1[1,])
    m2.ncol <- length(m2[1,])
    if(m1.ncol > m2.ncol){
        m1 <- m1[,1:m2.ncol]
        return(rbind(m1, m2))
    }
    if(m2.ncol > m1.ncol){
        m2 <- m2[,1:m1.ncol]
        return(rbind(m1, m2))
    }
    if(m2.ncol == m1.ncol){
        return(rbind(m1, m2))
    }
}

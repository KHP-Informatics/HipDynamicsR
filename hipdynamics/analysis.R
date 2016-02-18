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


if(!exists("aPrompt"))
{
    print("ANALYSER: ERROR (Make sure you load hipDynamics.R first).")
    quit(save="no")
}

# +++++++++++++++++++++++++++++++ SCRIPT ++++++++++++++++++++++++++++++++++++++

# transforming max and min factor vectors into regular vectors
db.max <- as.numeric(db.cols.max)
db.min <- as.numeric(db.cols.min)

# results table which will be exported as .csv
# creating a results table per experiment
results.table <- seq(1, (length(db.cols.select)*2+1))
results.table <- t(data.frame(results.table))
results.table.colnames <- c("Experiment", "Line", analysis_selector, 
                            paste(db.cols.select[2:length(db.cols.select)],"_gradient",sep=""),
                            paste(db.cols.select[2:length(db.cols.select)],"_yintercept", sep=""))
colnames(results.table) <- results.table.colnames

# Analyse per experiment
experiments <- unique(lookUp$Experiment)
for(experiment in experiments){
    print(paste(aPrompt, "Analysing", experiment, "..."))
    expIdx <- which(lookUp$Experiment == experiment)
   
    # Analyse per line
    lines <- unique(lookUp$Line[expIdx])
    
    # ++++++++ DEV +++++++++
    #   lines <- lines[1] #+
    # ++++++++++++++++++++++
    
    for(line in lines){
        print(paste(aPrompt, "  DB query for", line, "..."))
        lineIdx <- which(lookUp$Line[expIdx] == line)
        
        if(query == TRUE) res <- list()
       
        # Analyse per selector (see hipDynamicsPro.R file) 
        selected <- unique(lookUp[expIdx[lineIdx], analysis_selector])
        for(i in 1:length(selected)){  
            selIdx <- which(lookUp[expIdx[lineIdx], analysis_selector] == selected[i])
            
            if(query == TRUE){
                # ------ Selected DB Query -----------------------------------------
                # get image numbers for database query
                imageNumbers <- as.vector(lookUp$ImageNumber[expIdx[lineIdx[selIdx]]])
                # initiate database connection
                drv <- dbDriver("MySQL")
                con <- dbConnect(drv, db_name, user='hipsci', pass='data', host=db_host)
                
                statement <- paste("Select ImageNumber, ", db.cols.str, " from ", db_object_table, 
                                   " where ImageNumber >= ", imageNumbers[1],
                                   " and ImageNumber <= ", imageNumbers[length(imageNumbers)],
                                   " and ", db_column_selector ,"_AreaShape_Area > ",
                                   threshold_raw, ";", sep="")
                res[[i]] <- dbGetQuery(con, statement=statement)
                res[[i]] <- res[[i]][,-2]
                # close all db connections
                lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)
                # ------ Selected DB Query End--------------------------------------
            }
        }
        
        results <- list()
        # for selected (see hipDynamicsPro.R)
        for(s in 1:length(res)){
            print(paste(aPrompt, "    Histogramming", selected[s], "..."))
            
            selIdx <- which(lookUp[expIdx[lineIdx], analysis_selector] == selected[s])
            luIdx <- as.vector(expIdx[lineIdx[selIdx]])
            
            results[[s]] <- list()
            # for field
            for(field in 2:length(res[[s]][1,])){
                
                hours <- as.vector(unique(lookUp$Hour[luIdx]))
                results[[s]][[field]] <- list()
                
                # ------ Selected Objects according to hour --------------------
                # for hours 
                for(hour in hours){
                    imageNumbers <- as.numeric(as.vector(lookUp$Image[luIdx[which(lookUp$Hour[luIdx] == hour)]]))
                    # select image number from object table for hour
                    resIdx <- which(res[[s]]$ImageNumber == imageNumbers[1])
                    if(length(resIdx) > 1){
                        for(i in 2:length(imageNumbers)){
                            resIdx <- c(resIdx, which(res[[s]]$ImageNumber == imageNumbers[i]))
                        }  
                    }
                    results[[s]][[field]][[hour]] <- res[[s]][resIdx, field]
                }
                # ------ Selected Objects according to hour End ----------------
                
                # ------ Calculate thresholds ----------------------------------
                # binning can only be applied dynamically, since the range of the
                # values change, according to the measurements taken. The number
                # of bins can be specified and will be calculated according to 
                # the max and min of combined hours for each field.
                max <- db.max[field]
                min <- db.min[field]
                multiplier <- (max-min)/bin_number
                thresholds <- c(0:bin_number)*multiplier
                # ------ Calculate thresholds End ------------------------------
                
                # ------ Transform into time-dimensional matrix ----------------
                m.build <- c()
                for(hour in hours){
                    # ------ Bin and Histogram ---------------------------------
                    # bin each hour accoridng to threshholds and histogram for each hour
                    bins <- ApplyMultiThreshold(results[[s]][[field]][[hour]], thresholds)
                    m.build <- c(m.build, HistoOfBins(bins))
                    # ------ Bin and Histogram End -----------------------------
                }
                
                m.mtrx <- matrix(m.build, nrow=bin_number, byrow=FALSE)
                colnames(m.mtrx) <- hours
                rownames(m.mtrx) <- GetThresholdLabel(thresholds)
                
                results[[s]][[field]] <- m.mtrx
                # ------ Transform into time-dimensional matrix End ------------
                
                # ------ Normalise ---------------------------------------------
                if(normalise == TRUE){
                    # normalise per bin
                    for(bin in 1:bin_number){
                        results[[s]][[field]][bin,] <- Normalise(results[[s]][[field]][bin,])
                    }
                    # normalise per hour
                    for(hour in hours){
                        results[[s]][[field]][,hour] <- Normalise(results[[s]][[field]][,hour])
                    }
                }
                # ------ Normalise End -----------------------------------------
            }
        }
        
        print(paste(aPrompt, "  Thresholding & LnReg-ing..."))
        
        # initiating results rows
        expe.col <- matrix(rep(experiment, length(res)))
        line.col <- matrix(rep(line, length(res)))
        sele.col <- matrix(selected)
        start.col <- cbind(expe.col, line.col, sele.col)
        
        grad.col <- matrix(rep("-", length(res)))
        yint.col <- matrix(rep("-", length(res)))
        
        # initiating pdf outputs
        pdf(paste(path_out, "plots/", experiment, "_", 
                  line, "_", analysis_selector, "_bins", bin_number, ".pdf", sep=""),
            width=16, height=16, onefile=TRUE, family='Helvetica', paper='letter', 
            pointsize=12)
        
        for(field in 2:length(res[[s]][1,])){
            
            # ------ Heatmapping -----------------------------------------------
            # heatmap of all [FN]'s
            heat_mtrx <- buildHeatMatrix(results[[1]][[field]], results[[2]][[field]])
            if(length(results) > 2){
                for(s in 3:length(results)){
                    heat_mtrx <- buildHeatMatrix(heat_mtrx, results[[s]][[field]])
                }
            }
            # adding [FN] to rownames
            heat_mtrx_rows <- rownames(heat_mtrx)
            for(i in 1:length(selected)){
                idx_start <- bin_number*(i-1)+1
                idx_end <- bin_number*i
                heat_mtrx_rows[idx_start:idx_end] <- paste("[",selected[i],"] ",
                                                           heat_mtrx_rows[idx_start:idx_end],
                                                           sep="")
            }
            rownames(heat_mtrx) <- heat_mtrx_rows
            
            #PlotHeatHisto(heat_mtrx, 0, experiment, line, ylab=db.cols.select[field], all_FNs=TRUE)
            xlab <- "Time (in Hours)"
            ylab <- db.cols.select[field]
            main <- paste(experiment, "-", line,
                              ": Heatmap of ", ylab,sep="")
            
            # to maintain ordering
            h.cols <- c(1:length(heat_mtrx[1,]))
            colnames(heat_mtrx) <- h.cols[]
            e <- melt(heat_mtrx)
            e <- rename(e,c(X1="Field", X2="Hour"))
            e$Field <- with(e, reorder(Field))
            
            p <- ggplot(e, aes(Hour,Field))
            p <- p + geom_tile(aes(fill=value)) + scale_fill_gradient(low="white",
                                                                      high="steelblue") + xlab(xlab) + ylab(ylab) + ggtitle(main)
            print(p)
            #ggsave(paste(path_out, "plots/Heat_", experiment, "_", 
            #             line, "_", analysis_selector, "_bins", bin_number, "_field", field, ".pdf", sep=""),
            #       plot=p, width=16, height=16)
            
            # ------ Heatmapping End ------------------------------------------- 
            
            # initiating field column for line based results table
            yint.v <- c(1:length(res))
            grad.v <- c(1:length(res))
            
            # ------ Linear Regression -----------------------------------------
            par.len <- length(selected)
            if(par.len > 4) par.len <- 4
            par(mfrow=c(par.len,1))
            for(s in 1:length(results)){
                med <- c(1:length(results[[s]][[field]][1,]))
                
                for(i in 1:length(med)) med[i] <- median(log(results[[s]][[field]][,i]+1))      
                y <- med
                x <- c(1:length(med))
                fit <- lm(y ~ x)
                coefs <- coef(fit)
                y_intercept <- coefs[1]
                gradient <- coefs[2]
                field_title <- db.cols.select[field]
                # asigning values for field column
                yint.v[s] <- y_intercept
                grad.v[s] <- gradient
                
                boxplot(log(results[[s]][[field]]+1), 
                        main=paste(experiment, line, selected[s], "- Gradient:", round(gradient, 5),
                                   "Intercept:", round(y_intercept, 5)), 
                        ylab=paste("log", field_title), xlab="Time (in hours)")
                abline(fit ,col=2)
            }
            # ------ Linear Regression End -------------------------------------
            
            # adding to line based results.table for gradient and yintercept
            grad.col <- cbind(grad.col, matrix(grad.v))
            yint.col <- cbind(yint.col, matrix(yint.v))            
        }
        dev.off()
        
        # constructing line based results table
        results.table.line <- cbind(start.col,
                                    grad.col[,2:length(grad.col[1,])],
                                    yint.col[,2:length(grad.col[1,])])
        colnames(results.table.line) <- results.table.colnames
        # adding line based results table to full results table
        results.table <- rbind(results.table, results.table.line)
    }
}

# removing inititation row
results.table <- results.table[-1,]
expOut <- paste(experiments, collapse="_")
write.csv(results.table, file=paste(path_out, expOut, "_",analysis_selector, ".csv", sep=""))







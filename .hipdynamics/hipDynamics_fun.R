
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   _______ __         _____                              __                 +
#  |   |   |__|.-----.|     \.--.--.-----.---.-.--------.|__|.----.-----.    +
#  |       |  ||  _  ||  --  |  |  |     |  _  |        ||  ||  __|__ --|    +
#  |___|___|__||   __||_____/|___  |__|__|___._|__|__|__||__||____|_____|    +
#              |__|          |_____|                                         +
#                                                                            +
# TITLE:   HipDynamics - An analysis to deduce cell population dynamics    	 +
# VERSION: 1.1                                                               +
# AUTHOR:  Maximilian Kerz (kerz.maximilian@gmail.com)             			     +
#                                                                            +
# ACKNOWLEDGEMENTS: Amos Folarin (amosfolarin@gmail.com)                     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# +++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++


ImportPlateResults <- function(plate_results_txt){
  # reading PlateResults and returns a list in the following format:
  # ("A1:1", "A2:1", "A3:1", "A4:5", "A5:5", "A6:5", "A7:25", ..., "A8:25")
  plate_results = read.table(paste(path_PR,plate_results_txt, sep=""), 
                             sep=",")
  rows <- c("A", "B", "C", "D", "E", "F", "G", "H")
  PRDict <- c("")
  for(i in 10:length(plate_results[,1]))
  {
    str_to_split <- as.character(plate_results[i,1])
    Encoding(str_to_split) <- "latin1"
    str_split <- strsplit(str_to_split, "\t")
    # if line too early NA are introduced - 
    if(is.na(as.integer(str_split[[1]][1]))){
      print(paste(aPrompt, "SKIPPING row", i, "of", plate_results_txt))
    } else {
      well_id <- paste(rows[as.integer(str_split[[1]][1])], 
                       str_split[[1]][2], sep="")
      fn_conc <- str_split[[1]][60]
      cell_line <- str_split[[1]][61]
      # if less fields are present
      if(is.na(fn_conc)){
        fn_conc <- str_split[[1]][11]
        cell_line <- str_split[[1]][10]
      }
      PRDict[[length(PRDict)+1]] <- paste(well_id,":",fn_conc,":",
                                          cell_line, sep="")
    }
  }
  PRDict <- PRDict[2:length(PRDict)]
  return(PRDict)
}



AddPR <- function(PRDict, index_first=0, index_last=0){
  # input: ("A1:1", "A2:1", "A3:1", "A4:5", "A5:5", ..., "A8:25")
  # assign the former to appropriate fields in column fn_idx
  if (index_first == 0){
    index_first <- 0
    index_last <- length(per.img[,well_idx])
  } 
  exp_wells <- per.img[index_first:index_last,well_idx]
  
  for(i in 1:length(PRDict)){
    subDict <- PRDict[i]
    well_and_fn <- strsplit(subDict, ":")
    well_id <- well_and_fn[[1]][1]
    fn_conc <- as.integer(well_and_fn[[1]][2])
    cell_line <- well_and_fn[[1]][3]
    
    # add cell_liine to cell_lines if not present
    which_cl <- which(cell_lines == cell_line)
    if(length(which_cl)<1){
      cell_lines[length(cell_lines)+1] <<- cell_line
    }
    
    index_list <- which(exp_wells==well_id)
    index_start <- index_list[1] + index_first-1
    index_end <- index_list[length(index_list)] + index_first -1
    if(length(index_list) > 1){
      for(j in index_start:index_end)	{
        per.img[j,fn_idx] <<- fn_conc
        per.img[j,line_idx] <<- cell_line
      }
    }
  }
}



ImportExperiment <- function(plate_results_txt){
  
  # based on the 'date' in filename of plate_result.txt 
  # 	corresponding experments are identified.
  sub_name <- strsplit(plate_results_txt, "_")
  disjoint_date <- strsplit(sub_name[[1]][1], "-")
  date <- paste(disjoint_date[[1]][1], disjoint_date[[1]][2], 
                disjoint_date[[1]][3], sep="")
  
  index_list <- which(regexpr(date, per.img[,file_idx])>0)
  index_first <- index_list[1]
  index_last <- index_list[length(index_list)]
  
  if(!is.na(index_first)){
    # add experiment to global reference
    experiments[length(experiments)+1] <<- date
    print(paste(aPrompt,"FOUND experiment for", plate_results_txt,"- adding FN ..."))
    
    well_and_fn <- ImportPlateResults(plate_results_txt)
    AddPR(well_and_fn, index_first, index_last)
    
    # DEBUG: whether well, line and [FN arre assigned correctly
    # print(well_and_fn)
    
    print(paste(aPrompt,"DONE."))
    
  } else {
    # No matching experiment to date of plate_results.txt
    print(paste(aPrompt," ERROR (Experiment ", plate_results_txt, 
                " does not exist in hipsci_db.", table_per_img, ").", sep=""))
  }
}



ImportAllExperiments <- function(){
  # Resets experiments variable
  experiments <<- c()
  # cycles through all plate_results.txt's
  files_PR <- list.files(path_PR)
  for(file_PR in files_PR){
    ImportExperiment(file_PR)
  }
}



GetExperimentImg <- function(experiment_date){
  # extracts experiment table from per.img using
  # 	experiment's id (stored in experiments).
  index_list <- which(regexpr(experiment_date, per.img[,file_idx])>0)
  exp <- per.img[index_list,]
  return(exp)
}



GetExperimentObj <- function(exp_img){
  # extracts experiment table from per.obj using
  img_list <- exp_img[,img_no_idx]
  index_start <- which(per.obj[,obj_img_no_idx] == img_list[1])
  index_end <- which(per.obj[,obj_img_no_idx] == img_list[length(img_list)])
  index_list <- c(index_start[1]:index_end[length(index_end)])
  exp <- per.obj[index_list,]
  return(exp)
}



GetObj <- function(exp_obj, 
                   exp_img, 
                   img_idx){
  
  objs_idx <- which(exp_obj[,obj_img_no_idx] == exp_img[img_idx, img_no_idx])
	objs <- exp_obj[objs_idx,]
	return(objs)
}



GetField <- function(exp_obj, 
                     obj_idx){
  
	field <- exp_obj[,obj_idx]
	field_return <- field[!is.na(field)]	
	return(field_return)
}



GetFieldOfObj <- function(exp_obj, 
                          exp_img, 
                          img_idx,
                          obj_idx){
  
	objs <- GetObj(exp_obj, exp_img, img_idx)
	field <- GetField(objs, obj_idx)
	return(field)
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



Normalise <- function(x){
  ret <- (x-min(x))/(max(x)-min(x))
  for(i in 1:length(ret)){
    if(is.nan(ret[i])){ret[i] <- 0}
  }
  return(ret)
}



MeanStanDevOfBins <- function(bins){
  
	bin_no <- length(bins)
	mean <- c()
	sd <- c()
	
	for(i in 1:bin_no){
		mean <- c(mean, mean(bins[[i]]))
		sd <- c(sd, sd(bins[[i]]))
	}
	
	mean_sd <- c(mean, sd)
	mtrx_mean_sd <- matrix(mean_sd, nrow=bin_no)
	colnames(mtrx_mean_sd) <- c("Mean", "SD")
	return(mtrx_mean_sd)
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



GetImgsOfFN <- function(exp_img, 
                        fn_conc, 
                        as_type_matrix=TRUE, 
                        hour_suggest=c(24,23,20)){
  
	img_idx <- which(exp_img[,fn_idx] == fn_conc)
	if(as_type_matrix == TRUE){
    hours_test <- length(img_idx)
    for(h_sug in hour_suggest){
      if(hours_test%%h_sug == 0){
        img_hours <<- h_sug
        break
      }
    }
		img_idx <- matrix(img_idx, nrow=img_hours)
	}
	return(img_idx)
}



GetWellsOfFN <- function(exp_img, fn_conc){
	img_idx <- GetImgsOfFN(exp_img, fn_conc)
	wells <- exp_img[img_idx[1,], well_idx]
	return(wells)
}



GetThresholdLabel <- function(threshold_in){
	thresh_label <- c()
	for(i in 1:length(threshold_in)-1){
		thresh_label[i] <- paste(threshold_in[i],"-",threshold_in[i+1])
	}
	return(thresh_label)
}



GetHourLabel <- function(){
	hour_label <- c()
	for(i in 1:img_hours){
		hour_label[i] <- paste(i)
	}
	return(hour_label)
}



HClustFunc <- function(x, method = "complete", dmeth = "euclidean") {    
  hclust(dist(x, method = dmeth), method = method)
}



PlotExpHisto <- function(histo_input, 
                         thresh_list, 
                         fn_conc_list, 
                         experiment, 
                         cell_line, 
                         in_percent=FALSE){
  
	par(mfrow=c(1,3), oma=c(0,0,2,0))
	
  threshold_label <- GetThresholdLabel(thresh_list)
	bins_no <- length(thresh_list)-1
	legend_col <- c(2:bins_no+1)
	xlab <- "Time (in Hours)"
	ylab <- "Cell Area"
	if(in_percent == TRUE){
		ylab <- "Cell Count (in %)"
	}

	for(j in 1:length(histo_input)){
	
	main <- paste("Cell area over 24hrs at [FN] =", fn_conc_list[j])
	ymax_all <- c()
  for(i in 1:length(histo_input)){
    ymax_all <- max(histo_input[[i]])
  } 
  ymax <- max(ymax_all)
  if(in_percent == TRUE){
    ymax <- 100
  }
	
	plot(histo_input[[j]][1,], col=2, type='l', ylim=c(0,ymax), 
		ylab=ylab, xlab=xlab, main=main)
	legend("topright", threshold_label, col=legend_col,  lwd=c(1))
	for(i in 1:bins_no){
			lines(histo_input[[j]][i,], col=(i+1), type='l')
		}
	}
	title_out <- paste("Exp",experiment, "-", cell_line,sep="")
	title(title_out, outer=TRUE)
}



HeatExpHisto <- function(histo_input, 
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
	  main <- paste("Exp",experiment, "-", cell_line,
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



RawScatter <- function(raw_input, 
                       experiment, 
                       cell_line, 
                       fn_conc_input,
                       ymax, 
                       quantiles, 
                       index_name,
                       plot=TRUE, 
                       return_to_hour=img_hours){
  
  main <- paste("Exp", experiment,"-",cell_line, 
                ": Scatter plot of ", index_name," across time", sep="")
  ylab <- "Cell Area"
  xlab <- "Time (in Hours)"
  hours <- img_hours
  conc_no <- length(raw_input)/img_hours
  raw_col <- "lightsteelblue3"
  cols <- c("lightskyblue","royalblue","royalblue4")
  threshold_label <- c(paste(quantiles[1],"quantile"),paste(quantiles[2],"quantile"),
                       paste(quantiles[3],"quantile"))
  
  ret_row_label <- c()
  ret_quants <- list()
  r_quants1 <- c()
  r_quants2 <- c()
  r_quants3 <- c()
  
  for(conc in 1:conc_no){
    main <- paste("Exp", experiment,"-",cell_line, 
                  ": Scatter plot of ", index_name," at [FN]", 
                  fn_conc_input,"µg/ml",sep="")  
    quantile1 <- c()
    quantile2 <- c()
    quantile3 <- c()
    
    cell_no <- length(raw_input[[(conc-1)*hours + 1]])
    x_hour <- rep(1, cell_no)
    if(plot == TRUE){
      plot(y=raw_input[[(conc-1)*hours + 1]], x=x_hour, xlim=c(1,hours),
           ylim=c(0,ymax), xlab=xlab, ylab=ylab, main=main, col=raw_col)
      legend("topright", threshold_label, col=cols,  lwd=3)
    }
    quantile1 <- c(quantile1, quantile(raw_input[[(conc-1)*hours + 1]], quantiles[1]))
    quantile2 <- c(quantile2, quantile(raw_input[[(conc-1)*hours + 1]], quantiles[2]))
    quantile3 <- c(quantile3, quantile(raw_input[[(conc-1)*hours + 1]], quantiles[3]))
    
    for(hour in 2:hours){
      cell_no <- length(raw_input[[(conc-1)*hours + hour]])
      x_hour <- rep(hour, cell_no)
      if(plot == TRUE){
        points(y=raw_input[[(conc-1)*hours + hour]], x=x_hour, col=raw_col)
      }
      quantile1 <- c(quantile1, quantile(raw_input[[(conc-1)*hours+hour]], quantiles[1]))
      quantile2 <- c(quantile2, quantile(raw_input[[(conc-1)*hours+hour]], quantiles[2]))
      quantile3 <- c(quantile3, quantile(raw_input[[(conc-1)*hours+hour]], quantiles[3]))
    }
    if(plot == TRUE){
      lines(quantile1, lwd=3, col=cols[1])
      lines(quantile2, lwd=3, col=cols[2])
      lines(quantile3, lwd=3, col=cols[3])
    }
    
    # trim quantiles according to return_to_hour
    quantile1 <- quantile1[1:return_to_hour]
    quantile2 <- quantile2[1:return_to_hour]
    quantile3 <- quantile3[1:return_to_hour]
    
    # add quantiles to return vals
    ret_row_label <- c(ret_row_label, paste(experiment, "_", cell_line, "_FN",
                                            fn_concentrations[conc+1], sep=""))
    r_quants1 <- c(r_quants1, quantile1)
    r_quants2 <- c(r_quants2, quantile2)
    r_quants3 <- c(r_quants3, quantile3)
  }
  
  # create matrix with appropriate row names
  ret_quants[[1]] <- matrix(r_quants1, ncol=return_to_hour, byrow=TRUE)
  ret_quants[[2]] <- matrix(r_quants2, ncol=return_to_hour, byrow=TRUE)
  ret_quants[[3]] <- matrix(r_quants3, ncol=return_to_hour, byrow=TRUE)
  
  rownames(ret_quants[[1]]) <- ret_row_label
  rownames(ret_quants[[2]]) <- ret_row_label
  rownames(ret_quants[[3]]) <- ret_row_label
  colnames(ret_quants[[1]]) <- c(1:return_to_hour)
  colnames(ret_quants[[2]]) <- c(1:return_to_hour)
  colnames(ret_quants[[3]]) <- c(1:return_to_hour)
  
  return(ret_quants)
}



GetExpHistoMean <- function(fn_con_idx, 
                            thresh_list,
                            obj_idx_field,
                            get_SD=FALSE, 
                            get_raw=FALSE){
  
	col <- colnames(exp_obj)
  histo.mean <- list()
	histo.stde <- list()
	raw <- list()

	for(conc in 1:length(fn_con_idx)){

		print(paste(aPrompt,"  Getting",col[obj_idx_field],
                    "for [FN] =", fn_concentrations[conc+1]))

		histo_hour.mean <- c()
		histo_hour.stde <- c()

		for(hour in 1:img_hours){

			histo_well <- c()
			well_no <- length(fn_con_idx[[1]][hour,])
      raw_well <- c()
      
			for(well in 1:well_no){
				# gets cellarea per hour per [FN]
				c_idx <- fn_con_idx[[conc]][hour,well]
				a_idx <- obj_idx_field
				field <- GetFieldOfObj(exp_obj, exp_img, c_idx, a_idx)
        
        if(get_raw == FALSE){
				field.bins <- ApplyMultiThreshold(field, thresh_list)
				field.hist <- HistoOfBins(field.bins)
				histo_well <- c(histo_well, field.hist)
        } else {
          raw_well <- c(raw_well, field)
        }
			}
      
      if(get_raw == FALSE){
			bins_no <- length(thresh_list)-1
			histo_well <- matrix(histo_well, nrow=bins_no)
			histo_well_mean <- c()
			histo_well_stde <- c()
			for(bin in 1:bins_no){
				histo_well_mean <- c(histo_well_mean,mean(histo_well[bin,]))
				histo_well_stde <- c(histo_well_stde, sd(histo_well[bin,]))
			}

			histo_hour.mean <- c(histo_hour.mean, histo_well_mean)
			histo_hour.stde <- c(histo_hour.stde, histo_well_stde)
      } else {
        raw[[hour + (conc-1)*img_hours]] <- raw_well
      }
		}
		
    if(get_raw == FALSE){
		threshold_label <- GetThresholdLabel(thresh_list)
		hour_label <- GetHourLabel()
		
		mtrx_h.mean <- matrix(histo_hour.mean, nrow=bins_no)
		mtrx_h.stde <- matrix(histo_hour.stde, nrow=bins_no)
		rownames(mtrx_h.mean) <- threshold_label
		rownames(mtrx_h.stde) <- threshold_label
		colnames(mtrx_h.mean) <- hour_label
		colnames(mtrx_h.stde) <- hour_label
		histo.mean[[conc]] <- mtrx_h.mean
		histo.stde[[conc]] <- mtrx_h.stde
    }
	}

	ret <- histo.mean
	if(get_SD == TRUE){
		ret <- histo.stde
	}
  if(get_raw == TRUE){
    ret <- raw
  }
	return(ret)
}



HeatMapping <- function(){
  
  print(paste(aPrompt, 
              "  Histogramming [FN]s..."))
  
  exp_histo <- GetExpHistoMean(fn_exp_idx, threshold, get_raw=FALSE)
  
  # percent value histo
  exp_histo.perc <- exp_histo
  for(i in 1:length(exp_histo)){
    bins_no <- length(exp_histo[[i]][,1])  
    # calculate percentage (normalise) bin after prior normalisation
    for(j in 1:img_hours){
      exp_histo.perc[[i]][,j] <- Normalise(exp_histo.perc[[i]][,j], 100, 0)
    }
    # calculate percentage accross time
    for(j in 1:bins_no){
      exp_histo.perc[[i]][j,] <- Normalise(exp_histo.perc[[i]][j,], 100, 0)
    }
  }

  # single heatmaps of percent
  #print(paste(aPrompt, 
  #            "  Heatmapping of mean histogram of different bins per [FN] in percent."))
  #for(i in 1:length(exp_histo.perc)){
  #  HeatExpHisto(exp_histo.perc[[i]], fn_concentrations[i+1],
  #               experiment, line, FALSE, TRUE)
  #}
  
  # heatmap of all [FN]'s
  heat_mtrx <- rbind(exp_histo.perc[[1]],exp_histo.perc[[2]],exp_histo.perc[[3]])
  # adding [FN] to rownames
  heat_mtrx_rows <- rownames(heat_mtrx)
  for(i in 1:length(exp_histo.perc)){
    bins_no <- length(exp_histo.perc[[i]][,1])
    idx_start <- bins_no*(i-1)+1
    idx_end <- bins_no*i
    heat_mtrx_rows[idx_start:idx_end] <- paste("[",fn_concentrations[i+1],"] ",
                                               heat_mtrx_rows[idx_start:idx_end],
                                               sep="")
  }
  rownames(heat_mtrx) <- heat_mtrx_rows
  HeatExpHisto(heat_mtrx, fn_concentrations[i+1],
               experiment, line, TRUE, TRUE)
}

Proportionalise <- function(mtrx){
  array_mtrx <- mtrx
  for(i in 1:length(mtrx[,1])){
    hours_array <- array_mtrx[,i]
    max <- sum(hours_array)
    prop_array <- hours_array / max
    array_mtrx[,i] <- prop_array
  }
  return(array_mtrx)
}



DynamicCutoff <- function(input_mtrx, component, print_cutoffs=TRUE, display_graphs=FALSE){
  comp <- component
  mtrx <- input_mtrx
  
  title <- paste("Model of Noise of Component: ",comp)
  title2 <- paste("Graph of Significant Pix. Intensities: ",comp)
  
  yRange <- range(mtrx[,comp])
  dim <- dim(mtrx)
  den <- density(mtrx[,comp])
  intenseMax <- den$x[which.max(den$y)]
  densityMax <- den$y[which.max(den$y)]
  
  #calc for modelling
  denX <- c(1:512)
  denY <- c(1:512)
  cDenY <- c(1:512)
  
  # IMPORTANT: this model assumes the maximum density in lower half of density 
  #            plot, to avoid errors only the lower half is tested for maximum
  #            y-value and put forward for further computation.
  count_idx <- 256
  y_coord <- den$y
  count <- which.max(y_coord[1:count_idx])
  
  #noise curve values
  for(i in 1:count)cDenY[i] <- den$y[i]
  for(i in count:1)cDenY[(count-i)+count+1] <- den$y[i]
  for(i in (2*count):512)cDenY[i] <- 0.0
  
  #density curve values
  for(i in 1:512) denY[i] <- den$y[i]
  #x values
  for(i in 1:512) denX[i] <- den$x[i]
  
  
  #calc for subtraction
  mDenY <- c(1:512)
  aDenY <- c(1:512)
  #momentum curve
  for(i in 1:512)mDenY[i] <- denY[i]-cDenY[i]
  #additive curve
  incr <- 0
  for(i in 1:512){
    incr <- incr + mDenY[i]
    aDenY[i] <- incr
  }
  
  #calc for cutoff point
  
  #find intensity with maximum density increase
  momentMax <- denX[which.max(mDenY)]
  
  #find intensity of a density with a slope of less than 1e-4
  slopeX <- 1.0
  for(i in 512:which.max(mDenY)){
    slope <- (aDenY[i-20]-aDenY[i])/(denX[i-20]-denX[i])
    if(!is.na(slope) && slope<1e-4) slopeX <- denX[i-10]
  }
  
  #momentMax and slopeX are guidlines (blue) - its avg is used as guide to det. cutoff
  cutoffS <- (slopeX+momentMax)/2
  
  #cutoffM utilises momentMax and the maximum x value to det. the avg
  cutoffM <- (denX[512]*5/6+momentMax)/2
  
  #the cutoff is generated from the avg of cutoffS and cutoffM
  cutoff <- (cutoffS+cutoffM)/2
  
  cutoffs <- c(cutoff, momentMax, slopeX, cutoffS, cutoffM)
  
  if(print_cutoffs == TRUE){
    print("")
    print(paste("CUTOFF:      Hour", comp))
    print(paste("Cutoff:     ", round(cutoff, digit=2)))
    print(paste("MomentMax:  ", round(momentMax, digit=2)))
    print(paste("Slope(1e-4):", round(slopeX, digit=2)))
    print(paste("(Guides:    ", round(cutoffS, digit=2),
                " & ", round(cutoffM, digit=2),")"))
    print("")
  }
  
  if(display_graphs == TRUE){
    # no gradient
    cutoffST <- c(1:dim[1])
    for(i in 1:dim[1]) cutoffST[i] <- slopeX
    
    cutoffMT <- c(1:dim[1])
    for(i in 1:dim[1]) cutoffMT[i] <- momentMax
    
    cutoffT <- c(1:dim[1])
    for(i in 1:dim[1]) cutoffT[i] <- cutoff
    
    
    #GRAPHS
    par(mfrow=c(1,3))
    #GRAPH 1: density plots
    plot(den, type="l", main=title, ylab="Density",
         xlab="Pix.Intensity")
    #representing the lower values of the noise - curve will be modelled along it
    lines(denX,cDenY, type="l", col=3)    
    #vertical line showing the cutoff
    abline(v=intenseMax, col=2)
    text(x=intenseMax+yRange[2]/3, y=den$y[which.max(den$y)],labels=intenseMax, col=2)
    
    #GRAPH 2: significant values (cutoff)
    plot(denX, aDenY, type="l", col=1, main=title2, ylab="(Density-Modelled_Noise)", 
         xlab="Pix. Intensity")
    #momentum
    lines(denX, mDenY, col=3)
    #momentum max
    abline(v=momentMax, col=2)
    abline(v=slopeX, col=4)
    abline(v=cutoff, col=3)
    text(x=cutoff/3, y=aDenY[which.max(aDenY)]/10*9,labels=round(momentMax, digit=2), col=2)
    
    #GRAPH 3: adjusted curve
    plot(mtrx[,comp])
    text(x=dim[1]/4, y=cutoff/3*4,labels=paste("Cutoff :",
                                                round(momentMax, digit=2)), col=2)
    lines(cutoffST, col=4)
    lines(cutoffMT, col=2)
    lines(cutoffT, col=3)
  }
  
  return(cutoffs)
}



ComputeTrend <- function(input_mtrx, 
                         remove_debris=FALSE, 
                         print_result=FALSE, 
                         exp="00000000", 
                         line="llll", 
                         field="ffff",
                         fn_conc_input=0){
  
  mtrx <- input_mtrx
  scatter_array <- list()
  
  for(i in 1:length(mtrx[1,])){
    
    # iss: upper threshold?
    cutoffs <- DynamicCutoff(mtrx, i, print_cutoffs=FALSE, display_graphs=FALSE)
    maxMomentum <- cutoffs[2]
    
    threshed_idx <- ApplyThreshold(mtrx[,i], maxMomentum, idx=TRUE)
    if(length(threshed_idx)==0){
      scatter_array[i] <- list(c(0))
    } else {
      scatter_array[i] <- list(as.vector(threshed_idx))
    }
    
  }
  
  # removes lowest line (coordinate y=1, x=i), aka debris
  if(remove_debris==TRUE){
    for(i in 1:length(scatter_array)){
      threshed <- ApplyThreshold(scatter_array[[i]], 1)
      scatter_array[[i]] <- threshed 
    }
  }
  
  quants <- RawScatter(scatter_array, exp, line, fn_conc_input, length(mtrx[,1]), c(0.3,0.5,0.7), field)
  
  # linear regression and display of graph
  y <- as.vector(quants[[2]])
  x <- c(1:length(quants[[1]]))
  fit <- lm(y ~ x)
  coefs <- coef(fit)
  y_intercept <- coefs[1]
  gradient <- coefs[2]
  
  abline(fit ,col=2)
  
  if(print_result==TRUE){
    print(paste(aPrompt, "COMPUTE-TREND:", exp, line, field))
    print(paste(aPrompt, "Gradient:     ", gradient))
    print(paste(aPrompt, "Intercept:    ", y_intercept))
  }
  
  output = c(gradient, y_intercept)
  return(output)
}



saveCSV <- function(name, table, path_to_dir){
  
  if (!is.na(name)) {
    file_name = name
    write.csv(table, file = paste(path_to_dir, file_name, "_result.csv", sep=""), 
              row.names = TRUE)
  } else {
    print(paste(aPrompt, "The data.frame submitted is empty.No Data has been written."))
  }
}

summary <- function(){
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("+   _______ __         ______                             __               +")
  print("+  |   |   |__|.-----.|      |--.--.-----.---.-.--------.|__|.----.-----.  +")
  print("+  |       |  ||  _  ||  --  |  |  |     |  _  |        ||  ||  __|__ --|  +")
  print("+  |___|___|__||   __||_____/|___  |__|__|___._|__|__|__||__||____|_____|  +")
  print("+              |__|          |_____|                                       +")
  print("+                                                                          +")
  print("+ TITLE:   HipDynamics - An analysis to deduce cell population dynamics    +")
  print("+ VERSION: 1.1                                                             +")
  print("+ AUTHOR:  Maximilian Kerz (kerz.maximilian@gmail.com)                     +")
  print("+                                                                          +")
  print("+ ACKNOWLEDGEMENTS: Amos Folarin (amosfolarin@gmail.com)                   +")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("")
  print("---------------------------- SUMMARY ---------------------------------------")
  print("DATABASE Settings:")
  print(paste("Database: ",db_name))
  print(paste("User:     ",db_usern))
  print(paste("Img table:",table_per_img))
  print(paste("Obj table:",table_per_obj))
  print("")
  print("OPTION Settings:")
  print(paste("Raw-bin Inferential Analysis:   ",infer_from_raw_bins))
  print(paste("Normalised Inferential Analysis:",infer_from_normalisation))
  print(paste("Export raw data:                ",export_raw_data))
  print(paste("Export normalised data:         ",export_normalised_data))
  print("")
  print("ANALYSIS:")
  print(paste("No of Experiments:", length(unique(final_output_table[,1]))))
  print(paste("No of Lines:      ", length(unique(final_output_table[,2]))))
  print(paste("No of dif. [FN]s: ", length(unique(final_output_table[,3]))))
  print("")
  print("All information is held within the final_output_table. For graphs")
  print(paste("and exports please view the folder loacted at:", path_out))
  print("------------------------------- End ----------------------------------------")
  
}













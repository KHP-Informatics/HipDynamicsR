
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

if(!exists("aPrompt"))
{
  print("ANALYSER: ERROR (Make sure you load hipDynamics.R first).")
  quit(save="no")
}

# only one analysis method can be selcted: tested for
if(infer_from_raw_bins == infer_from_normalisation){
  print(paste(aPrompt, "Both analysis methods were selected. Please only select one."))
  stopifnot(infer_from_raw_bins != infer_from_normalisation)
}

# +++++++++++++++++++++++++++++++ SCRIPT ++++++++++++++++++++++++++++++++++++++

# instantiating new final output table
final_output_table <- data.frame()

for(experiment in experiments){
  
  print(paste(aPrompt, "For experiment:",experiment)) # PRINT
  exp_img_full <- GetExperimentImg(experiment)
  exp_obj_full <- GetExperimentObj(exp_img_full)
  
  # instantiating new experiment table
  experiment_output_table <- data.frame()
  
  for(idx_no in 1:fields_list){
    col <- colnames(exp_obj_full)
    index_name <- obj_index_list[idx_no]
    index_vars <- strsplit(index_name, "_")
    index_var <- index_vars[[1]][3]
    
    #dev.off() writes it to the files
    pdf(file=paste(path_out,index_var, "_", experiment,"_Thresh",threshold_raw,"_plots_&_heatmaps.pdf",
                   sep=""), height=9, width=7, onefile=TRUE, family='Helvetica', paper='letter', 
        pointsize=12)
    
    lines <- unique(exp_img_full[,line_idx])
    # as plate result do not contain all information for all rows NA are omitted
    lines <- lines[!is.na(lines)]
    
    for(line in lines){
      
      print(paste(aPrompt, "For line:",line)) # PRINT
      line_img_idx <- which(exp_img_full[,line_idx] == line)
      exp_img <- exp_img_full[line_img_idx,]
      exp_obj <- GetExperimentObj(exp_img)
      
      # get well index for [FN] = 1, 5, 25 (0)
      fn_exp_idx <- list()
      fn_exp_idx[[1]] <- GetImgsOfFN(exp_img, fn_concentrations[2], hour_suggest=img_hours_suggest)
      fn_exp_idx[[2]] <- GetImgsOfFN(exp_img, fn_concentrations[3], hour_suggest=img_hours_suggest)
      fn_exp_idx[[3]] <- GetImgsOfFN(exp_img, fn_concentrations[4], hour_suggest=img_hours_suggest)
      
      exp_histo <- GetExpHistoMean(fn_exp_idx, threshold, obj_index_list[idx_no], get_raw=FALSE)
      print(paste(aPrompt, 
                  "  Histogramming [FN]s..."))
      
      # percent value histo
      exp_histo.perc <- exp_histo
      for(i in 1:length(exp_histo)){
        bins_no <- length(exp_histo[[i]][,1])  
        # normalise per bin
        for(j in 1:bins_no){
          exp_histo.perc[[i]][j,] <- Normalise(exp_histo.perc[[i]][j,])
        }
        # normalise per hour
        for(j in 1:img_hours){
          exp_histo.perc[[i]][,j] <- Normalise(exp_histo.perc[[i]][,j])
        }
      }
      
      if(export_raw_data==TRUE){
        save(exp_histo, file=paste(path_out,experiment,line,"-exo_histo.RData", sep=""))
      }
      
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
      
      # save heat_mtrx dataframe
      if(export_normalised_data==TRUE){
        save(heat_mtrx, file=paste("~/Desktop/",experiment,line,"-heat_mtrx.RData", sep=""))
      }
      
      fn_con <- fn_concentrations[2:4]
      
      # SPECIFY IN SETTINGS
      if(infer_from_raw_bins==TRUE){       
        fn_output_table <- data.frame()
        par(mfrow=c(3,1))
        
        for(d in 1:3){
          med <- c(1:length(exp_histo[[d]][1,]))
          
          for(i in 1:length(med)) med[i] <- median(log(exp_histo[[d]][,i]+1))      
          y <- med
          x <- c(1:length(med))
          fit <- lm(y ~ x)
          coefs <- coef(fit)
          y_intercept <- coefs[1]
          gradient <- coefs[2]
          field <- index_var
          exp <- experiment
          
          boxplot(log(exp_histo[[d]]+1), 
                  main=paste(experiment, line, fn_con[d], "Gradient:", round(gradient, 5),
                             "Intercept", round(y_intercept, 5)))
          abline(fit ,col=2)
          
          output_table <- data.frame(exp, line, fn_con[d], gradient, y_intercept, field)
          
          if(length(fn_output_table)==0){
            fn_output_table <- output_table
          } else{
            fn_output_table <- rbind(fn_output_table, output_table)
          }
        }
        
        # OUTPUT
        print(paste(aPrompt, "Raw-Bin Inferential Analysis"))
        print(paste(aPrompt, fn_output_table[1,1],fn_output_table[1,2],fn_output_table[1,3], " ",
                    fn_output_table[1,4],fn_output_table[1,5],fn_output_table[1,6]))
        print(paste(aPrompt, fn_output_table[2,1],fn_output_table[2,2],fn_output_table[2,3], " ",
                    fn_output_table[2,4],fn_output_table[2,5],fn_output_table[2,6]))
        print(paste(aPrompt, fn_output_table[3,1],fn_output_table[3,2],fn_output_table[3,3], "",
                    fn_output_table[3,4],fn_output_table[3,5],fn_output_table[3,6]))
        
        if(length(experiment_output_table)==0){
          experiment_output_table <- fn_output_table
        } else {
          experiment_output_table <- rbind(experiment_output_table, fn_output_table)
        }
        if(length(final_output_table)==0){
          final_output_table <- output_table
        } else {
          final_output_table <- rbind(final_output_table, fn_output_table)
        }
      }
      
      # SPECIFY IN SETTINGS
      if(infer_from_normalisation==TRUE){
        
        fn1_mtrx <- heat_mtrx[1:23,]
        fn5_mtrx <- heat_mtrx[24:46,]
        fn25_mtrx <- heat_mtrx[47:69,]
        
        DynamicCutoff(fn1_mtrx, 20, print_cutoffs=FALSE, display_graphs=FALSE)
        
        par(mfrow=c(3,1))
        output1  <- ComputeTrend(fn1_mtrx, remove_debris=TRUE, exp=experiment, line=line, 
                                 field=index_var, fn_conc_input=1)
        output5  <- ComputeTrend(fn5_mtrx, remove_debris=TRUE, exp=experiment, line=line, 
                                 field=index_var, fn_conc_input=5)
        output25 <- ComputeTrend(fn25_mtrx, remove_debris=TRUE, exp=experiment, line=line, 
                                 field=index_var, fn_conc_input=25)
        print(aPrompt)
        
        exp <- c(experiment, experiment, experiment)
        line <- c(line, line, line)
        field <- c(index_var, index_var, index_var)
        gradient <- c(output1[[1]], output5[[1]], output25[[1]])
        y_intercept <- c(output1[[2]], output5[[2]], output25[[2]])
        output_table <- data.frame(exp, line, fn_con, gradient, y_intercept, field)
        
        # OUTPUT
        print(paste(aPrompt, "Normalised Inferential Analysis"))
        print(paste(aPrompt, output_table[1,1],output_table[1,2],output_table[1,3], " ",
                    output_table[1,4],output_table[1,5],output_table[1,6]))
        print(paste(aPrompt, output_table[2,1],output_table[2,2],output_table[2,3], " ",
                    output_table[2,4],output_table[2,5],output_table[2,6]))
        print(paste(aPrompt, output_table[3,1],output_table[3,2],output_table[3,3], "",
                    output_table[3,4],output_table[3,5],output_table[3,6]))
        
        if(length(experiment_output_table)==0){
          experiment_output_table <- output_table
        } else {
          experiment_output_table <- rbind(experiment_output_table, output_table)
        }
        if(length(final_output_table)==0){
          final_output_table <- output_table
        } else {
          final_output_table <- rbind(final_output_table, output_table)
        }
      }
    
      print(aPrompt)
    }
    print(aPrompt)
    dev.off()
    
  }
  saveCSV(experiment, experiment_output_table, path_out)
}
saveCSV("FINAL", final_output_table, path_out)



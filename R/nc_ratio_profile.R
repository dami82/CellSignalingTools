nc_ratio_profile <-
function(profile_tab, bckgr_subtr = "baseline", mode = "auc", condition = "", antibody = "Ab - I", show_labels = TRUE, noise = 0.05){
  #
  my_data <- profile_tab
  my_data <- my_data[!is.na(my_data[,"X1"]),]
  #
  for (i in 2:ncol(my_data)){
    #remove NAs, col by col
    my_data[,i][is.na(my_data[,i])] <- min(my_data[,i][!is.na(my_data[,i])])
  }
  res_tab <- sapply(1:((ncol(my_data)-1)/2), (function(jj){
    #
    bl_v <- paste("B", jj, sep = '')
    gr_v <- paste("G", jj, sep = '')
    #
    plot(1,1, 
         xlim = c(0,(max(my_data$X1) + max(my_data$X1)/9)),
         ylim = c(0,1.2),
         main = paste(condition, ": cell #", jj, sep =''), 
         xlab = "x-position", 
         ylab = "normalized signal", 
         type= "n")  
    #
    green_profile <- scale_01(my_data[,gr_v])
    blue_profile <- scale_01(my_data[,bl_v])
    my_range <- get_dapi_range(blue_profile)
    #
    green_matrix <- cbind(my_data[,"X1"],green_profile) 
    colnames(green_matrix) <- c("pos", "signal")
    #
    left_profile <- green_matrix[my_range[1]:(my_range[2]),] 
    cent_profile <- green_matrix[(my_range[2]+1):(my_range[3]-1),] 
    rite_profile <- green_matrix[my_range[3]:my_range[4],] 
    #
    left_ave_sig <- mean(left_profile[left_profile[,2] > noise ,2])
    cent_ave_sig <- mean(cent_profile[cent_profile[,2] > noise ,2])
    rite_ave_sig <- mean(rite_profile[rite_profile[,2] > noise ,2])
    #
    left_auc <- get_auc(left_profile)
    rite_auc <- get_auc(rite_profile)
    cyto_auc <- left_auc + rite_auc
    cent_auc <- get_auc(cent_profile)
    #
    baseline <- (min(rite_profile[,"pos"]) - max(left_profile[,"pos"])) * (mean(left_ave_sig, rite_ave_sig))
    base_corr_cent <- cent_auc - baseline 
    #
    if (bckgr_subtr == "confocal") {
      polygon(c(0,left_profile[,"signal"],0) ~ c(min(left_profile[,"pos"]),left_profile[,"pos"],max(left_profile[,"pos"])), col = "lightgreen")
      polygon(c(0,rite_profile[,"signal"],0) ~ c(min(rite_profile[,"pos"]),rite_profile[,"pos"],max(rite_profile[,"pos"])), col = "lightgreen")
      #
      polygon(c(0,cent_profile[,"signal"],0) ~ c(min(cent_profile[,"pos"]),cent_profile[,"pos"],max(cent_profile[,"pos"])), col = "limegreen")
    } else {
      polygon(c(0, left_profile[,"signal"],
                left_ave_sig, rite_ave_sig, 
                rite_profile[,"signal"], 0) ~ 
                c(min(left_profile[,"pos"]), left_profile[,"pos"],
                  max(left_profile[,"pos"]), min(rite_profile[,"pos"]),
                  rite_profile[,"pos"], max(rite_profile[,"pos"])), 
              col = "lightgreen")
      #
      polygon(c(left_ave_sig,cent_profile[,"signal"],rite_ave_sig) ~ c(min(cent_profile[,"pos"]),cent_profile[,"pos"],max(cent_profile[,"pos"])), col = "limegreen")
    }
    #
    lines(scale_01(my_data[,bl_v]) ~ my_data$X1, col = "darkblue", lwd = 3.75, lty = 1)
    #
    corr_x <- 0.015
    segments(corr_x,1.1,0.095,1.1, lwd = 15, col = "limegreen")
    segments(0.09,1.1,(0.15-corr_x),1.1, lwd = 15, col = "lightgreen")
    segments(0,1,0.15,1, lwd = 4, col = "blue")
    polygon(c(0,0,0.1,0.1),c(0.88,0.92,0.92,0.88), col = "limegreen")
    polygon(c(0,0,0.1,0.1),c(0.78,0.82,0.82,0.78), col = "lightgreen")
    text(0.15,1.1,labels = antibody, pos = 4, cex = 0.75, font = 2)
    text(0.15,1,labels = "DAPI", pos = 4, cex = 0.75, font = 2)
    text(0.1,0.9,labels = "nucl. signal", pos = 4, cex = 0.75, font = 2)
    text(0.1,0.8,labels = "cyto. signal", pos = 4, cex = 0.75, font = 2)
    #
    print_ratio <- 0
    #
    if (mode == "auc") {
      my_mode <- "area under curve"
      #return auc and ratio
      if (bckgr_subtr == "baseline") {
        print_ratio <- base_corr_cent/ (cyto_auc + baseline)
        vec_to_return <- c(base_corr_cent, cyto_auc + baseline, print_ratio)  
      } else {
        print_ratio <- cent_auc/cyto_auc
        vec_to_return <- c(cent_auc, cyto_auc, print_ratio)  
      }
    } else {
      my_mode <- "average signal"
      if (bckgr_subtr == "baseline") {
        #return averages and ratio
        print_ratio <- (cent_ave_sig-mean(rite_ave_sig, left_ave_sig)) /mean(rite_ave_sig, left_ave_sig)
        vec_to_return <- c(cent_ave_sig, mean(rite_ave_sig, left_ave_sig), print_ratio)
      } else {
        #return averages and ratio
        print_ratio <- cent_ave_sig/mean(rite_ave_sig, left_ave_sig)
        vec_to_return <- c(cent_ave_sig, mean(rite_ave_sig, left_ave_sig) , print_ratio)
      }
    }
    #
    if (bckgr_subtr == "baseline") {
      my_correction <- "avg. baseline"
    } else {
      my_correction <- "no correction"
    }
    #
    if(cent_ave_sig - mean(left_ave_sig, rite_ave_sig) > 0.2) {
      segments(max(my_data$X1) + max(my_data$X1)/27, cent_ave_sig, 
               max(my_data$X1) + max(my_data$X1)/9, cent_ave_sig, 
               lwd = 4, col = "red")
      segments(max(my_data$X1) + max(my_data$X1)/27, mean(left_ave_sig, rite_ave_sig), 
               max(my_data$X1) + max(my_data$X1)/9, mean(left_ave_sig, rite_ave_sig), 
               lwd = 4, col = "red")
      text(max(my_data$X1) + max(my_data$X1)/7.25, 
           mean(left_ave_sig, rite_ave_sig) + 0.05, 
           "cyto signal", pos = 2, cex = 0.75, font = 2)
      text(max(my_data$X1) + max(my_data$X1)/7.25, 
           cent_ave_sig +  0.05, 
           "nucl signal", pos = 2, cex = 0.75, font = 2)
      #
      if (show_labels == TRUE) {
        text(max(my_data$X1) + max(my_data$X1)/7.25, 
             1.15, paste("Nucl/Cyto Ratio:", round(print_ratio, digits = 2)), pos = 2, cex = 0.85, font = 2)
        text(max(my_data$X1) + max(my_data$X1)/7.25, 
             1.05, paste("Bckgrnd subtr.:",my_correction), pos = 2, cex = 0.685, font = 3)
        text(max(my_data$X1) + max(my_data$X1)/7.25, 
             0.982, paste("Values:", my_mode), pos = 2, cex = 0.685, font = 3)
      }
      #
    } else {
      if (show_labels == TRUE) {
        
        text(max(my_data$X1) + max(my_data$X1)/7.25, 1.05, 
             "no evident difference between nucl and cyto signal", pos = 2, cex = 0.75, font = 2)
      }
    }
    vec_to_return
  }))
  res_tab <- t(res_tab)
  colnames(res_tab) <- c("nucl", "cyto", "ratio")
  return(res_tab)
}

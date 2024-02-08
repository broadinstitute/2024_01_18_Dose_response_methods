# Process WTT Results from BMDExpress
# February 7, 2024
# Jessica Ewald

# The purpose of this script is to process and compile the WTT output from BMDExpress. The files 
# created by the software have added empty columns and character type, so here are some tricks 
# for processing them. This will have to be adapted for each creation of the WTT results: different 
# computing systems will be able to handle larger batches for the httr data. 

# httr 2 reps
httr.2.files <- list.files("../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/httr_2reps",
                           full.names = TRUE)

res <- data.frame()
for(i in c(1:length(httr.2.files))){

  col.nms <- unlist(read.table(httr.2.files[i], sep = "\t", fill = T, nrow = 3)[1, ])
  col.len <- length(col.nms)
  col.nms[(col.len-1):col.len] <- c("X1", "X2")


  temp <- read.table(httr.2.files[i], sep = "\t", fill = T, skip = 3, col.names = col.nms)

  temp <- temp[,c("Analysis", "Probe.ID", "Unadjusted.P.Value", "Adjusted.P.Value", "Max.Fold.Change")]
  temp[,3:5] <- apply(temp[,3:5], 2, as.numeric)

  res <- rbind(res, temp)

}
write.table(res, "../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/processed/httr_2reps_wtt.txt",
            row.names = F, sep = "\t", quote = FALSE)

# httr 3 reps
httr.3.files <- list.files("../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/httr_3reps",
                           full.names = TRUE)

res <- data.frame()
for(i in c(1:length(httr.3.files))){
  
  col.nms <- unlist(read.table(httr.3.files[i], sep = "\t", fill = T, nrow = 3)[1, ])
  col.len <- length(col.nms)
  col.nms[(col.len-1):col.len] <- c("X1", "X2")
  
  
  temp <- read.table(httr.3.files[i], sep = "\t", fill = T, skip = 3, col.names = col.nms)
  
  temp <- temp[,c("Analysis", "Probe.ID", "Unadjusted.P.Value", "Adjusted.P.Value", "Max.Fold.Change")]
  temp[,3:5] <- apply(temp[,3:5], 2, as.numeric)
  
  res <- rbind(res, temp)
  
}
write.table(res, "../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/processed/httr_3reps_wtt.txt",
            row.names = F, sep = "\t", quote = FALSE)

# htpp 2 reps
col.nms <- read.table("../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/htpp_2reps_wtt.txt",
                      sep = "\t", fill = T, nrows = 3)
col.nms <- unlist(col.nms[1,])
col.len <- length(col.nms)
col.nms[(col.len-1):col.len] <- c("X1", "X2")

res <- read.table("../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/htpp_2reps_wtt.txt",
                  sep = "\t", fill = T, skip = 3, col.names = col.nms)
res <- res[,c("Analysis", "Probe.ID", "Unadjusted.P.Value", "Adjusted.P.Value", "Max.Fold.Change")]
res[,3:5] <- apply(res[,3:5], 2, as.numeric)

write.table(res, "../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/processed/htpp_2reps_wtt.txt",
            row.names = F, sep = "\t", quote = FALSE)

# htpp 3 reps
col.nms <- read.table("../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/htpp_3reps_wtt.txt",
                      sep = "\t", fill = T, nrows = 3)
col.nms <- unlist(col.nms[1,])
col.len <- length(col.nms)
col.nms[(col.len-1):col.len] <- c("X1", "X2")

res <- read.table("../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/htpp_3reps_wtt.txt",
                  sep = "\t", fill = T, skip = 3, col.names = col.nms)
res <- res[,c("Analysis", "Probe.ID", "Unadjusted.P.Value", "Adjusted.P.Value", "Max.Fold.Change")]
res[,3:5] <- apply(res[,3:5], 2, as.numeric)

write.table(res, "../2024_01_18_Dose_response_methods-data/data/2c_WTT_results/processed/htpp_3reps_wtt.txt",
            row.names = F, sep = "\t", quote = FALSE)

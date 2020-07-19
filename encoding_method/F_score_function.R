F_score_fn <- function(dataframe){
  
  # change logic feature to numeric feature
  for (i in 1:(ncol(dataframe))){
    if (is.logical(dataframe[,i]) == TRUE){
      dataframe[,i] <- as.numeric(dataframe[,i])
    }
  }
  
  # variables
  half_index <- nrow(dataframe)/2
  pos_part <- dataframe[c(1:half_index),]
  neg_part <- dataframe[-c(1:half_index),]
  pos_mean <- colMeans(pos_part) 
  neg_mean <- colMeans(neg_part) 
  total_mean <- colMeans(dataframe) 
  df_pos_mean <- matrix(NA, ncol = ncol(pos_part), nrow = nrow(pos_part)) 
  df_neg_mean <- matrix(NA, ncol = ncol(pos_part), nrow = nrow(neg_part))
  for (i in 1:nrow(pos_part)){
    df_pos_mean[i,] <- pos_mean
    df_neg_mean[i,] <- neg_mean
  }
  
  # numeritor
  numerator <- (pos_mean - total_mean)^2 + (neg_mean - total_mean)^2 
  # denominator
  denominator <- 1/(half_index - 1)*((colSums(pos_part - df_pos_mean)^2) + (colSums(neg_part - df_neg_mean)^2))
  
  # F score 
  F_score <- as.data.frame(numerator / denominator)
  
  return(F_score)
}


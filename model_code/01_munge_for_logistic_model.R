
source("00_functions.R")

# Import essential meta features including 10 bp surrounding sequence
process_sample <- function(sample, editor){
  df_list <- subset_data_to_balance(importForLogistic(sample, editor, pad = 5, chrs = 1:20))
  saveRDS(df_list, file = paste0("../logistic_processed/", editor, "_", sample, "_dfs_for_logistic.rds"))
  paste0(sample, editor)
}

# Uniform processing on 7 May 2019
if(FALSE){
  # Process the ABE samples
  process_sample("243B", "ABE")
  process_sample("243C", "ABE")
  
  #process_sample("244B", "ABE")
  #process_sample("244C", "ABE")
  
  #process_sample("247B", "ABE")
  #process_sample("247C", "ABE")
  
  #process_sample("156B", "ABE")
  #process_sample("157B", "ABE")
  
  # Process the CBE samples - A3A
  process_sample("160F", "CBE")
  #process_sample("161F", "CBE")
  
  # Process CBE samples -- BE3
  process_sample("89B", "CBE")
  #process_sample("90B", "CBE")
}


library(dplyr)
source("00_functions.R")

io <- function(sample, editor){
  print(sample)
  x <- importEssential(sample, editor)
  saveRDS(x, file = paste0("processed_df/", editor, "-", sample, "-essentialDF.rds"))
  sample
}

io("160B", "CBE")
io("160F", "CBE")
io("89A", "CBE")
io("89B", "CBE")

io("243A", "ABE")
io("243B", "ABE")
io("243C", "ABE")
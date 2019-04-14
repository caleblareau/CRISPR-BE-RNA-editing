library(dplyr)
library(tools)

brc <- list.files("../../../", pattern = "all_positions.txt.gz", recursive = TRUE, full.names = TRUE)
abs_brc <- sapply(brc, file_path_as_absolute) %>% unname()
filepath_final <- abs_brc[grepl("exp", abs_brc) & !grepl("-nextseq-qc", abs_brc) & !grepl("exp1-novaseq", abs_brc) & !grepl("104", abs_brc)]
id <- gsub(".all_positions.txt.gz", "", gsub("EJUL", "", basename(filepath_final)))

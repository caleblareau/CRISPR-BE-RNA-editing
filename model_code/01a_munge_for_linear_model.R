
source("00_functions.R")

importForLinear <- function(sample, editor, chrs = 1:20){
  
  # Import sequence fasta files
  dir_seq <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-sequence")
  dir_struct <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-secondary")
  
  seq_fastas <- mixedsort(sort(list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)))
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas)][chrs]
  struct_files <- mixedsort(sort(list.files(dir_struct, pattern = paste0(sample), full.names = TRUE)))
  struct_files <- struct_files[grepl(".gz", struct_files)][chrs]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
  # Process structure data
  struct_input <- sapply(struct_files, function(faf) fread(paste0("zcat < ", faf), header = FALSE)[[1]][c(FALSE, TRUE)]) %>% unlist()
  paired <- substring(struct_input, 51, 51) != "."; rm(struct_input)
  
  # Get all sequence attributes
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  
  # Process meta to convert characters to numeric
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  meta <- data.frame(meta, sequence = unname(unlist(fasta_input)), paired)
  meta <- meta[meta$editRate > 0.02,]
  rownames(meta) <- NULL
  return(meta)
}


# Import essential meta features including 10 bp surrounding sequence
process_sample_linear <- function(sample, editor, chrs = 1:20){
  df <- importForLinear(sample, editor, chrs = chrs)
  saveRDS(df, file = paste0("../linear_processed/", editor, "_", sample, "_dfs_for_linear.rds"))

  paste0(sample, editor)
}

# Uniform processing on 7 May 2019
if(TRUE){
  # Process the ABE samples
  process_sample_linear("243C", "ABE")
  process_sample_linear("243B", "ABE")


  #process_sample("244B", "ABE")
  #process_sample("244C", "ABE")
  
  #process_sample("247B", "ABE")
  #process_sample("247C", "ABE")
  
  #process_sample("156B", "ABE")
  #process_sample("157B", "ABE")
  
  # Process the CBE samples - A3A
  process_sample_linear("160F", "CBE")
  #process_sample("161F", "CBE")
  
  # Process CBE samples -- BE3
  process_sample_linear("89B", "CBE")

  #process_sample("90B", "CBE")
}


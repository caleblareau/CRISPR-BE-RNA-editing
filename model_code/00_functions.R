library(seqinr)
library(stringr)
library(reshape2)
library(data.table)
library(dplyr)
library(gtools)

importForLogistic <- function(sample, editor, pad = 1, proximal_GC_radius = 10, chrs = 1:20){
  
  # Reorder chromosomes for fast subsetting
  order_chrs <- c(1, 12, 16:22, 2:11, 13:15, 23)
  
  # Import sequence fasta files
  dir_seq <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-sequence")
  dir_struct <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-secondary")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)[order_chrs]
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas)][chrs]
  struct_files <- list.files(dir_struct, pattern = paste0(sample), full.names = TRUE)[order_chrs]
  struct_files <- struct_files[grepl(".gz", struct_files)][chrs]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
  # Process structure data
  struct_input <- sapply(struct_files, function(faf) fread(paste0("zcat < ", faf), header = FALSE)[[1]][c(FALSE, TRUE)]) %>% unlist()
  paired <- substring(struct_input, 51, 51) != "."; rm(struct_input)
  
  # Get all sequence attributes
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  seq_proximal <- str_split_fixed(sapply(fasta_input, function(x) substring(x, first = 51 - pad, last = 51 + pad)) %>% unlist() %>% unname(), "", (2*pad + 1))
  
  # Pull GC content
  proximal_GC <- sapply(fasta_input, function(x) substring(x, first = 51 - proximal_GC_radius, last = 51 + proximal_GC_radius)) %>%
    unlist() %>% unname() %>% 
    DNAStringSet() %>% letterFrequency( "GC")
  
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  
  # Process meta to convert characters to numeric
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  meta <- data.frame(meta, seq_proximal, paired, proximal_GC = proximal_GC/sum(proximal_GC_radius*2))
  meta$isEdited <- meta$editRate > 0.05
  boo <- meta$editRate > 0.05 | meta$editRate == 0
  meta$chr <- stringr::str_split_fixed(meta$chr_pos, "-", 2)[,1]
  rownames(meta) <- NULL
  return(meta[boo,])
}


importForDeep_noStructure <- function(sample, editor){
  
  saveForDeepLift <- c("chr21", "chr22", "chrX")
  
  # Import sequence fasta files
  dir_seq <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-sequence")
  #dir_struct <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-structure")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas) & !grepl(paste(saveForDeepLift,collapse="|"), seq_fastas)]
  #struct_files <- list.files(dir_struct, pattern = paste0(sample), full.names = TRUE)
  #struct_files <- struct_files[grepl(".gz", struct_files)]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
  # Process structure data
  #struct_input <- sapply(struct_files, function(faf) fread(cmd = paste0("zcat < ", faf), header = FALSE)[[1]][c(FALSE, TRUE)]) %>% unlist() %>% unname()
  
  # Get all sequence attributes
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  meta$chr <- stringr::str_split_fixed(meta$chr_pos, "-", 2)[,1]
  seqs <- unlist(fasta_input) %>% unname()
  
  # Process meta to convert characters to numeric
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  rownames(meta) <- NULL
  meta$sequence <- seqs
  #meta$structure <- struct_input
  meta$isEdited <- meta$editRate > 0.05
  boo <- meta$editRate > 0.05 | meta$editRate == 0
  return(meta[boo,])
}

subset_data_to_balance <- function(all_data, seed = 5, n_excess = 5){
  
  stopifnot(all(c("isEdited") %in% colnames(all_data)))
  
  # Subset the data based on chromosome
  test_chromosomes <- c("chr16", "chr17", "chr18", "chr19", "chr20")
  set.seed(seed)
  put_in_train <- !grepl(paste(test_chromosomes,collapse="|"), all_data$chr)
  put_in_test <- grepl(paste(test_chromosomes,collapse="|"), all_data$chr)
  
  # Downsample
  n_pos_in_train <- sum(put_in_train & all_data[,"isEdited"])
  n_pos_in_test <- sum(put_in_test & all_data[,"isEdited"])
  n_neg_in_train <- sum(put_in_train & all_data[,"isEdited"])
  n_neg_in_test <- sum(put_in_test & all_data[,"isEdited"])
  
  # Extract indices
  set.seed(seed)
  n_excess <- n_excess
  idx_train <- sample(c(which(put_in_train & all_data[,"isEdited"]), 
                        sample(which(put_in_train & !all_data[,"isEdited"]), size = min(n_neg_in_train,n_excess*n_pos_in_train))))
  idx_test <- sample(c(which(put_in_test & all_data[,"isEdited"]), 
                       sample(which(put_in_test & !all_data[,"isEdited"]), size = min(n_neg_in_test,n_excess*n_pos_in_test))))
  return(list("train" =  all_data[idx_train,], 
              "test" =  all_data[idx_test,]))
}


make_one_hot_eight_channel <- function(x) {
  
  # Split into each individual matrix
  seq_mat <-str_split_fixed( x[,1],  "", 101)
  struct_mat <- str_split_fixed(x[,2],  "", 101)
  
  # Encode a fifth channel for sequence -- editing event
  seq_mat[,51] <- "Z"
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- rbind(reshape2::melt(seq_mat), reshape2::melt(struct_mat)); rm$new <- 1
  aa <- acast(rm, Var1  ~  Var2  ~  value, fill = 0, value.var = "new")  
  aa
}


make_one_hot_five_channel_string <- function(x) {
  
  # Split sequences into individual characters
  spm <- str_split_fixed(x, "", 101)
  
  # Encode a fifth channel
  spm[,51] <- "Z"
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- reshape2::melt(spm); rm$new <- 1
  aa <- acast(rm, Var1  ~  Var2  ~  value, fill = 0, value.var = "new")  
  aa
}


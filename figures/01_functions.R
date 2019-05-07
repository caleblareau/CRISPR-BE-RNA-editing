library(BuenColors)
library(seqinr)
library(stringr)
library(dplyr)
library(data.table)
library(gtools)
library(Biostrings)
options(datatable.fread.input.cmd.message=FALSE)


# Reorder chromosomes for fast subsetting
order_chrs <- c(1, 12, 16:22, 2:11, 13:15, 23)

importMetaOnly <- function(sample, editor, chrs = 1:23){

  # Pull meta data from 
  dir_struct <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-secondary")
  struct_files <- list.files(dir_struct, pattern = paste0(sample), full.names = TRUE)[order_chrs]
  struct_files <- struct_files[grepl(".gz", struct_files)][chrs]
  struct_input <- sapply(struct_files, function(faf) fread(paste0("zcat < ", faf), header = FALSE)[[1]][c(TRUE, FALSE)]) %>% unlist() %>% unname()
  
  meta <- data.frame(stringr::str_split_fixed(struct_input, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage", "annotation", "gene")
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  meta$chr <- stringr::str_split_fixed(meta$chr_pos, "-", 2)[,1]
  return(meta)
}

importEssential <- function(sample, editor, pad = 1, proximal_GC_radius = 10, chrs = 1:23){
  
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
  rownames(meta) <- NULL
  return(meta)
}

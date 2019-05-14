library(seqinr)
library(stringr)
library(reshape2)
library(data.table)
library(dplyr)
library(gtools)

"%ni%" <- Negate("%in%")

importForSecure <- function(sample, editor){
  
  saveForDeepLift <- c("chr21", "chr22", "chrX")
  
  # Import sequence fasta files
  dir_seq <- paste0("../../CRISPR-",editor,"-RNA-DATA/fastas/", editor, "-sequence")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas) & !grepl(paste(saveForDeepLift,collapse="|"), seq_fastas)]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
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
  return(meta)
}

import_sample_secure <- function(library){
  
  # Import control sample to filter out anything post-hoc
  a_sample <- paste0(library, "A")
  nCas9 <- importForSecure(a_sample, "ABE")
  rm_nozero_loci <- nCas9 %>% filter(editRate > 0) %>% pull(chr_pos) %>% as.character()
  
  c_sample <- paste0(library, "C")
  monoABE <- importForSecure(c_sample, "ABE")
  
  d_sample <- paste0(library, "D")
  monoABE_mut1 <- importForSecure(d_sample, "ABE")
  
  e_sample <- paste0(library, "E")
  monoABE_mut2 <- importForSecure(e_sample, "ABE")
  
  # Subset and remove stuff / merge
  monoABE_mut1 <- monoABE_mut1[,c("chr_pos", "editRate")]; colnames(monoABE_mut1) <- c("chr_pos", "mut_editRate")
  monoABE_mut2 <- monoABE_mut2[,c("chr_pos", "editRate")]; colnames(monoABE_mut2) <- c("chr_pos", "mut_editRate")
  mdf_mut1 <- inner_join(monoABE, monoABE_mut1, by = "chr_pos") %>% filter(editRate >0 | mut_editRate >0) %>% filter(chr_pos %ni% rm_nozero_loci)
  mdf_mut2 <- inner_join(monoABE, monoABE_mut2, by = "chr_pos")  %>% filter(editRate >0 | mut_editRate >0)%>% filter(chr_pos %ni% rm_nozero_loci)
  
  # Export
  mdf_mut1$Library <- paste0(library, "-K20AR21A")
  mdf_mut2$Library <- paste0(library, "-V82G")
  saveRDS(mdf_mut1, file = paste0("../secure_processed/", library, "-K20AR21A", ".rds"))
  saveRDS(mdf_mut2, file = paste0("../secure_processed/", library, "-V82G", ".rds"))
  
}
import_sample_secure("243")
# import_sample_secure("244")
import_sample_secure("247")

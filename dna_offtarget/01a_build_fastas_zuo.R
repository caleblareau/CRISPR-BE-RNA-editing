library(dplyr)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(reshape2)
library(stringr)
library(data.table)
library(diffloop)
library(seqinr)


# Crude pass that pulls from gDNA -- does not consider any potential splicing
getSequenceSimpleMM10 <- function(gr_query, reverse_complement = FALSE){
  
  gs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr_query)
  
  # Flip the sequence if the user requests it
  if(reverse_complement){
    gs_rc <- reverseComplement(gs)
    gs <- gs_rc
  }
  
  return(as.character(gs))
}

tab <- fread("data/zuo_BE3_offtargets.txt")
breaks <- which(is.na(tab[["V2"]]))

# Make altered df of essential attributes
idx <- 1:dim(tab)[1]
boo <- !(idx %in% breaks)
experiments <- tab$V1[breaks]
df <- data.frame(chr = tab$V1[boo],
                 pos = tab$V2[boo],
                 rate = round(tab$V10[boo], 3),
                 anno = gsub("_", "-", as.character(tab$V4)[boo]),
                 coverage = tab$V6[boo] + tab$V7[boo], 
                 strand = ifelse(grepl("G", as.character(tab$V3[boo])), "-", "+"),
                 experiment = c(rep("BE3_1", 279 - 1 - 1), rep("BE3_2", 418-280 -1 ),
                                rep("BE3_TyrC1", 746 - 419 -1 ), rep("BE3_TyrC2", 1104 - 747 -1),
                                rep("BE3_TyrD1", 1441 - 1105 -1), rep("BE3_TyrD2", max(idx) - 1442 )
                                )) %>% mutate(start = pos -50, end = pos + 50)
df_f <- df %>% filter(strand == "+")
df_r <- df %>% filter(strand == "-")

seqs_forward <- getSequenceSimpleMM10(makeGRangesFromDataFrame(df_f))
seqs_reverse <- getSequenceSimpleMM10(makeGRangesFromDataFrame(df_f))

# Establish headers and write out
header_pos <- paste0("JIN_", df_f$chr, "-", df_f$pos, "_F_", as.character(df_f$rate), "_", as.character(df_f$coverage), "_", df_f$anno, "_", df_f$experiment)
header_rev <- paste0("JIN_", df_r$chr, "-", df_r$pos, "_R_", as.character(df_r$rate), "_", as.character(df_r$coverage), "_", df_r$anno, "_", df_r$experiment)

fasta_file <- "fastas/zuo_all_BE.fasta"
write.fasta(as.list(seqs_forward, seqs_reverse), c(header_pos, header_rev), fasta_file, open = "w")

system(paste0("gzip ", fasta_file))


library(seqinr)
library(caret)
library(speedglm)
library(precrec)
library(BuenColors)
library(stringr)
library(gtools)

source("../figures/01_functions.R")

# Transform and stuff
transformInput <- function(df, library){
  df %>%
    mutate(transformed_outcome = logit(editRate)) %>%
    mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
    mutate(library = library)
}

importDNA <- function(sample,  pad = 1, proximal_GC_radius = 10){
  # Import sequence fasta files
  dir_seq <- paste0("fastas/")

  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas)]

  # Import data into a list
  fasta_input <-  unlist(read.fasta(seq_fastas, as.string = TRUE))
  
  # Process structure data
  paired <- TRUE # because its DNA
  
  # Get all sequence attributes
  seq_names <- names(fasta_input)
  seq_proximal <- str_split_fixed(substring(unname(fasta_input), first = 51 - pad, last = 51 + pad) %>% unlist() %>% unname(), "", (2*pad + 1))
  
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


# Do a simple linear model
CBE1 <- transformInput(importEssential("89B", "CBE", pad = 3), "Train")
linearMod <- speedlm(transformed_outcome ~ G.C + upstream + downstream + as.numeric(paired), data=CBE1 %>% filter(transformed_outcome <10 & transformed_outcome > -10))

Zuo_df <- transformInput(importDNA("Zuo_all_BE", pad = 3), "DNA-test") %>% filter(transformed_outcome <10 & transformed_outcome > -10)
Zuo_df$predicted <- predict(linearMod, Zuo_df)

cor(Zuo_df$transformed_outcome, Zuo_df$predicted)
p1 <- ggplot(shuf(Zuo_df), aes(x = transformed_outcome, y = predicted, color = gene)) +
  pretty_plot(fontsize = 9) + L_border() + scale_color_manual(values = c("firebrick", "green4", "purple3", "orange3", "black", "dodgerblue3")) +
  geom_point() + labs(x = "Observed edit rate (logit scale)", y = "Predicted edit rate (logit scale)", color = "Zuo et. al\nLibrary")
cowplot::ggsave(p1, file = "zuo_etal_predicted.pdf", width = 4, height = 3.5)






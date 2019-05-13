library(dplyr)
library(BuenColors)
library(ggseqlogo)
library(reshape2)
library(stringr)

make_one_hot_from_df <- function(x) {
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- melt(x, id.vars = "chr_pos"); rm$new <- 1
  aa <- acast(rm, chr_pos  ~  variable+value, fill = 0, value.var = "new")  
  aa[x$chr_pos,]
}


make_grid_motif_enrich <- function(editor, sample, name){
  
  # Import and make one hot encoding of proximal  nucleotides
  list_df <- readRDS(paste0("../linear_processed/",editor,"_",sample,"_dfs_for_linear.rds"))
  list_df <- list_df[list_df$editRate > 0.05 & list_df$editRate < 1,]
  list_df$transformedEditRate <- logit(list_df$editRate)
  seq_mat <- str_split_fixed(as.character(list_df$sequence), "", 101)
  essential_df <- data.frame(paired = list_df$paired, seq_mat, chr_pos = list_df$chr_pos)
  one_hot <- make_one_hot_from_df(essential_df)
  x3mers <- substring(as.character(list_df$sequence), 50, 52)
  print("making plots")
  
  lapply(unique(sort(x3mers)), function(x3mr){
    n = sum(x3mr == x3mers)
    
    data.frame(seq_mat[x3mr == x3mers,], x = 1) %>% reshape2::melt(id.vars = "x") %>%
      group_by(variable, value) %>% summarize(count = n()) %>%
      ungroup() %>% group_by(variable) %>%
      mutate(freq = count/sum(count)) %>%
      reshape2::dcast(variable~value, value.var = "freq") -> freq_df
    colnames(freq_df) <- toupper(colnames(freq_df))
    
    p1 <- ggseqlogo(t(freq_df[,c("A", "C", "G", "T")]))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) + L_border() +
      scale_y_continuous(expand = c(0,0)) +
      ggtitle(paste0(name, " cluster ", as.character(x3mr), " n=", as.character(n)))
    p1
  }) -> list_of_plots
  pdf(paste0("output_figures/cluster_motifs/", name, "_", editor, "_", sample, "_5p3pMotifs", ".pdf"), 10, 2)
  lapply(list_of_plots, print)
  dev.off()
}

make_grid_motif_enrich("CBE", "160F", "A3A")
make_grid_motif_enrich("CBE", "89B", "BE3")
make_grid_motif_enrich("ABE", "243B", "ABEmax")
make_grid_motif_enrich("ABE", "243C", "miniABEmax")

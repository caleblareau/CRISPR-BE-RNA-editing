library(uwot)
library(irlba)
library(reshape2)
library(Matrix)

make_one_hot_from_df <- function(x) {
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- melt(x, id.vars = "chr_pos"); rm$new <- 1
  aa <- acast(rm, chr_pos  ~  variable+value, fill = 0, value.var = "new")  
  aa[x$chr_pos,]
}


make_umap_grid_plot <- function(editor, sample, name){
  # Import and make one hot encoding of proximal  nucleotides
  list_df <- readRDS(paste0("../logistic_processed/",editor,"_",sample,"_dfs_for_logistic.rds"))
  df <- rbind(list_df[["train"]], list_df[["test"]])
  df <- df[df$editRate > 0 & df$editRate < 1,]
  df <- head(df, 5000)
  df$transformedEditRate <- logit(df$editRate)
  one_hot <- make_one_hot_from_df(df[,c(paste0("X", as.character(1:11)), "paired", "chr_pos")])
  
  
  # Run LSI
  nfreqs <- t(t(one_hot) / Matrix::colSums(one_hot))
  idf <- as(log(1 + ncol(one_hot) / Matrix::rowSums(one_hot)), "sparseVector")
  tf_idf_counts <- as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs
  SVD_go <-  irlba(tf_idf_counts, 20, 20)
  
  # Perform a supervised UMAP embedding
  unsupervised_umap <- umap(SVD_go$u)
  
  # Visualize
  colnames(unsupervised_umap) <- c("UMAP1", "UMAP2")
  plot_df <- data.frame(unsupervised_umap, SVD_go$u[,c(1,2, 3,4,5,6)], transformed_rate = df$transformedEditRate,
                        GC = df$G.C, paired = df$paired, two5p = df$X4, base5p = df$X5,base3p = df$X7, two3p = df$X8)
  
  dna_color_vec <- c("a" = "green3", "c" = "dodgerblue3", "g" = "orange3", "t" = "red")
  
  p1 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = base5p)) +
    geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
    scale_color_manual(values = dna_color_vec)
  
  p2 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = base3p)) +
    geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
    scale_color_manual(values = dna_color_vec)
  
  p3 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = paired)) +
    geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border()
  
  p4 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = two5p)) +
    geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
    scale_color_manual(values = dna_color_vec)
  
  p5 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = two3p)) +
    geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
    scale_color_manual(values = dna_color_vec)
  
  p6 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = transformed_rate)) +
    geom_point(size = 0.5) + scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
    pretty_plot(fontsize = 8) + L_border() + labs(color = "rate")
  
  cowplot::ggsave(cowplot::plot_grid(p1, p2, p3, p4, p5, p6), 
                  filename = paste0("output_figures/umaps/",name,"_embedding.pdf"), 
                  width = 10, height = 6)
}

make_umap_grid_plot("ABE", "243B", "miniABEmax")
make_umap_grid_plot("ABE", "243C", "ABEmax")
make_umap_grid_plot("CBE", "89B", "BE3")
make_umap_grid_plot("CBE", "160F", "A3A")



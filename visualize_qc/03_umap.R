library(uwot)
library(irlba)
library(reshape2)
library(Matrix)
library(stringr)
library(ggseqlogo)


make_one_hot_from_df <- function(x) {
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- melt(x, id.vars = "chr_pos"); rm$new <- 1
  aa <- acast(rm, chr_pos  ~  variable+value, fill = 0, value.var = "new")  
  aa[x$chr_pos,]
}


make_umap_grid_plot <- function(editor, sample, name){
  
  # Import and make one hot encoding of proximal  nucleotides
  list_df <- readRDS(paste0("../linear_processed/",editor,"_",sample,"_linearDF.rds"))
  list_df <- list_df[list_df$editRate > 0.05 & list_df$editRate < 1,]
  list_df$transformedEditRate <- logit(list_df$editRate)
  seq_mat <- str_split_fixed(as.character(list_df$sequence), "", 101)
  essential_df <- data.frame(paired = list_df$paired, seq_mat, chr_pos = list_df$chr_pos)
  one_hot <- make_one_hot_from_df(essential_df)
  
  # Run LSI
  nfreqs <- t(t(one_hot) / Matrix::colSums(one_hot))
  idf <- as(log(1 + ncol(one_hot) / Matrix::rowSums(one_hot)), "sparseVector")
  tf_idf_counts <- as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs
  SVD_go <-  irlba(tf_idf_counts, 50, 50)
  
  if(FALSE){
    df <- data.frame(SVD_go$v, what = colnames(tf_idf_counts), letter = str_split_fixed(colnames(tf_idf_counts), "_", 2)[,2])
    df <- df[df$letter %in% c("a", "c", "g", "t"),]
    p1 <- ggplot(df, aes(x = X1, y = X2, label = what, color = letter)) + 
      geom_text() + labs(x = "LSI1", y = "LSI2")
    p2 <- ggplot(df, aes(x = X3, y = X4, label = what, color = letter)) + 
      geom_text() + labs(x = "LSI3", y = "LSI4")
    p3 <- ggplot(df, aes(x = X5, y = X6, label = what, color = letter)) + 
      geom_text() + labs(x = "LSI5", y = "LSI6")
    p4 <- ggplot(df, aes(x = X7, y = X8, label = what, color = letter)) + 
      geom_text() + labs(x = "LSI7", y = "LSI8")
    
    cowplot::ggsave(cowplot::plot_grid(p1, p2, p3, p4, nrow = 2), 
                    file = "output_figures/umaps/LSI4.pdf", width = 10, height = 10)
  }
  
  # Perform a supervised UMAP embedding + community detection
  knn <- FNN::get.knn(SVD_go$u, algo="kd_tree", k = 10)[["nn.index"]]
  igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(knn)[,c("Var1", "value")]), directed=FALSE)), mode = "undirected")
  louvain_output <- igraph::cluster_louvain(igraphObj)
  per_edit_cluster <- igraph::membership(louvain_output)
  table(per_edit_cluster)
  
  
  lapply(sort(unique(per_edit_cluster)), function(cluster){
    n = sum(cluster == per_edit_cluster)
    
    data.frame(seq_mat[cluster == per_edit_cluster,], x = 1) %>% reshape2::melt(id.vars = "x") %>%
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
      ggtitle(paste0(name, " cluster ", as.character(cluster), " n=", as.character(n)))
    p1
  }) -> list_of_plots
  pdf(paste0("output_figures/cluster_motifs/", name, "_", editor, "_", sample, "_clusterMotifs", ".pdf"), 10, 2)
  list_of_plots
  dev.off()
  
  # Do umap things
  unsupervised_umap <- umap(SVD_go$u)
  if(typeof(unsupervised_umap) == "list") unsupervised_umap <- unsupervised_umap$layout
  
  # Visualize
  colnames(unsupervised_umap) <- c("UMAP1", "UMAP2")
  plot_df <- data.frame(unsupervised_umap, transformed_rate = df$transformedEditRate, cluster = as.character(per_edit_cluster), 
                        GC = df$G.C, paired = df$paired, two5p = df$X4, base5p = df$X5,base3p = df$X7, two3p = df$X8)
  
  dna_color_vec <- c("a" = "green3", "c" = "dodgerblue3", "g" = "orange3", "t" = "red")
  
  p0 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border()
  
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
  
  cowplot::ggsave(cowplot::plot_grid(p0, p1, p2, p6, p4, p5), 
                  filename = paste0("output_figures/umaps/",name,"_embedding.pdf"), 
                  width = 10, height = 6)
}

make_umap_grid_plot("ABE", "243B", "ABEmax")
make_umap_grid_plot("ABE", "243C", "miniABEmax")
make_umap_grid_plot("CBE", "89B", "BE3")
make_umap_grid_plot("CBE", "160F", "A3A")



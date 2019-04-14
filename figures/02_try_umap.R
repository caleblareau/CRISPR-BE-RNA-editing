library(uwot)
library(irlba)
library(reshape2)
library(Matrix)
library(ggrastr)

source("01_functions.R")

make_one_hot_from_df <- function(x) {
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- melt(x, id.vars = "seqid"); rm$new <- 1
  aa <- acast(rm, seqid  ~  variable+value, fill = 0, value.var = "new")  
  aa[x$seqid,]
}

# Import and make one hot encoding of proximal 20 nucleotides
df <- importEssential("146D", "CBE", pad = 20)
df_filt <- df %>% filter(editRate > 0 & editRate < 0.9) %>% mutate(transformed = logit(editRate), seqid = paste0("R", 1:n()))
one_hot <- make_one_hot_from_df(df_filt[,c(paste0("X", as.character(1:41)), "paired", "seqid")])

# Run LSI
nfreqs <- t(t(one_hot) / Matrix::colSums(one_hot))
idf <- as(log(1 + ncol(one_hot) / Matrix::rowSums(one_hot)), "sparseVector")
tf_idf_counts <- as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs
SVD_go <-  irlba(tf_idf_counts, 10, 10)

# Perform a supervised UMAP embedding
# nn <- FNN::get.knn(SVD_go$u, k = 5)
# names(nn) <- c("idx", "dist")
supervised_umap <- umap(SVD_go$u, y = df_filt$transformed,  n_epochs = 200, target_weight = 1)

# Visualize
colnames(supervised_umap) <- c("UMAP1", "UMAP2")
plot_df <- data.frame(supervised_umap, SVD_go$u[,c(1,2, 3,4,5,6)], rate = df_filt$transformed,
                      GC = df_filt$G.C, paired = df_filt$paired, X22 = df_filt$X22)
ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = X22)) +
  geom_point() 

ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = rate)) +
  geom_point() + scale_color_gradientn(colors = jdb_palette("brewer_spectra"))



svd_feature <- data.frame(SVD_go$v, names = colnames(tf_idf_counts))
ggplot(svd_feature, aes(x = X1, y = X5, label = names)) +
  geom_text()

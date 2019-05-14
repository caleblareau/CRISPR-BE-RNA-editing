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



# Import and make one hot encoding of proximal  nucleotides
list_df <- readRDS(paste0("../secure_processed/243-K20AR21A.rds"))
list_df <- list_df %>% filter(coverage > 100 & editRate >= 0.05) 
x3mers1 <- substring(as.character(list_df$sequence), 50, 52)
data.frame(three = x3mers1, WT = list_df$editRate, MUT = list_df$mut_editRate) %>%
  group_by(three) %>% summarize(mean_WT = round(mean(WT)*100,2), mean_K20R21A = round(mean(MUT)*100,2)) -> K20AR21A_mut

list_df2 <- readRDS(paste0("../secure_processed/243-V82G.rds"))
list_df2 <- list_df2 %>% filter(coverage > 100 & editRate >= 0.05) 

x3mers2 <- substring(as.character(list_df2$sequence), 50, 52)
data.frame(three = x3mers2, WT = list_df2$editRate, MUT = list_df2$mut_editRate) %>%
  group_by(three) %>% summarize(mean_WT = round(mean(WT)*100,2), mean_V82G = round(mean(MUT)*100,2)) -> V82G


final_df <- data.frame(K20AR21A_mut, mean_V82G = V82G$mean_V82G)


ggplot(final_df, aes(x = mean_K20R21A, y = mean_V82G, label = three)) +
  geom_text() +
  labs(x = "Mean K20AR21A editing rate", y = "mean V82G editing rate") +
  pretty_plot() + L_border() + scale_y_continuous(limits = c(0.5,2.2)) +
  scale_x_continuous(limits = c(0.5,2.2)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

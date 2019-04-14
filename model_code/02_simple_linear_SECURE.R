library(seqinr)
library(caret)
library(speedglm)
library(precrec)
library(BuenColors)
library(stringr)
library(gtools)

source("../figures/01_functions.R")

# Import essential meta features
x147A <- importEssential("147A", "CBE", pad = 3)
x147B <- importEssential("147B", "CBE", pad = 3)
x147C <- importEssential("147C", "CBE", pad = 3)
x147D <- importEssential("147D", "CBE", pad = 3)

x147A_clip <- x147A %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

x147B_clip <- x147B %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

x147C_clip <- x147C %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

x147D_clip <- x147D %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

compare_hist_df <- data.frame(
  value = c(x147A_clip$transformed_outcome, x147B_clip$transformed_outcome, x147C_clip$transformed_outcome, x147D_clip$transformed_outcome),
  library = c(rep("Control", dim(x147A_clip)[1]), rep("BE3", dim(x147B_clip)[1]), rep("R33A", dim(x147C_clip)[1]), rep("R33AK34A", dim(x147D_clip)[1]))
)

p1 <- ggplot(compare_hist_df, aes(x = value, color = library)) + 
  geom_density(adjust = 3) + labs(x = "logit-transformed editing rate", y = "Density") +
  scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 8) + L_border() + 
  scale_color_manual(values = c("firebrick", "black", "green4", "dodgerblue3")) +
  theme(legend.position = c(0.8, 0.6))
cowplot::ggsave(p1, file = "../figures/output_figures/03a_147_density.pdf", width = 2, height = 2)


summarize_tile <- function(df, library){
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  df %>% group_by(X3,X5) %>%
    summarize(edit = mean(editRate)) %>% 
    ungroup() %>%
    mutate(library=library, editNorm = range01(edit) )
}

tile_plot_df <- rbind(summarize_tile(x147A, "Control"),
                      summarize_tile(x147B, "BE3"),
                      summarize_tile(x147C, "R33A"),
                      summarize_tile(x147D, "R33AK34A"))

p2 <- ggplot(tile_plot_df, aes(x = X3, y = X5, fill = editNorm)) +
  geom_tile() + pretty_plot(fontsize = 6) + L_border() + labs(x = "5' base", y = "3' base") +
  scale_y_discrete(expand = c(0,0), labels = rev(c("A", "C", "G", "U"))) +
  scale_x_discrete(expand = c(0,0), labels = c("A", "C", "G", "U")) +
  scale_fill_gradientn(colors = jdb_palette("brewer_red")) +
  facet_wrap(~library) + 
  theme(legend.position = "none") 
cowplot::ggsave(p2, file = "../figures/output_figures/03b_147_tile4.pdf", width = 5, height = 5)




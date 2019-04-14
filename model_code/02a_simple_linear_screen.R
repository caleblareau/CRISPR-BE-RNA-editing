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

# Import essential meta features
x118A <- transformInput(importEssential("118A", "CBE", pad = 3), "Control")
x118B <- transformInput(importEssential("118B", "CBE", pad = 3), "BE3")
x118C <- transformInput(importEssential("118C", "CBE", pad = 3), "E63Q")
x119A <- transformInput(importEssential("119A", "CBE", pad = 3), "P29F")
x119B <- transformInput(importEssential("119B", "CBE", pad = 3), "P29T")
x119C <- transformInput(importEssential("119C", "CBE", pad = 3), "L182A")
x120A <- transformInput(importEssential("120A", "CBE", pad = 3), "R33A")
x120B <- transformInput(importEssential("120B", "CBE", pad = 3), "K34A")
x120C <- transformInput(importEssential("120C", "CBE", pad = 3), "R33AK34A")

compare_hist_df <- data.table::rbindlist(list(x118A, x118B, x118C, x119A, x119B, x119C, x120A, x120B, x120C)) %>%
  data.frame()

p1 <- ggplot(compare_hist_df, aes(x = transformed_outcome, color = library)) + 
  stat_ecdf(adjust = 3) + labs(x = "logit-transformed editing rate", y = "Density") +
  scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = c(0.8, 0.6))
cowplot::ggsave(p1, file = "../figures/output_figures/03c_screen_density.pdf", width = 5, height = 5)

# Do a simple linear model
CBE1 <- transformInput(importEssential("89B", "CBE", pad = 3), "Train")
linearMod <- speedlm(transformed_outcome ~ G.C + upstream + downstream + as.numeric(paired), data=CBE1 %>% filter(transformed_outcome <10 & transformed_outcome > -10))

get_cor <- function(df){
  df2 <- df %>% filter(transformed_outcome <10 & transformed_outcome > -10)
  cor(predict(linearMod, df2), df2$transformed_outcome)
}

levels <- c("nCas9", "BE3", "E63Q", "P29F", "P29T", "L182A", "R33A", "K34A", "R33AK34A")
tile2_df <- data.frame(
  Factor = factor(levels, level = levels),
  mean_edit_rate = round(sapply(list(x118A, x118B, x118C, x119A, x119B, x119C, x120A, x120B, x120C), function(x) mean(x$editRate)) * 100, 2),
  correlation = round(sapply(list(x118A, x118B, x118C, x119A, x119B, x119C, x120A, x120B, x120C), get_cor),2),
  n_hits = sapply(list(x118A, x118B, x118C, x119A, x119B, x119C, x120A, x120B, x120C), function(x) sum(x$editRate > 0.05))
)

tile2_df$mean_edit_rate_cap <- ifelse(tile2_df$mean_edit_rate > 1, 1, tile2_df$mean_edit_rate)
tile2_df$correlation_cap <- ifelse(tile2_df$correlation > 0.5, 0.5, tile2_df$correlation)
tile2_df$n_hits_cap <- ifelse(tile2_df$n_hits > 5000, 5000, tile2_df$n_hits)

pA <- ggplot(tile2_df, aes(x = Factor, y = 1, fill = mean_edit_rate_cap, label = mean_edit_rate)) +
  geom_tile(color = "black") +
  labs(x = "", fill = "") + geom_text(size = 2) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_red")) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pB <- ggplot(tile2_df, aes(x = Factor, y = 1, fill = n_hits_cap, label = n_hits)) +
  geom_tile(color = "black") +
  labs(x = "", fill = "") + geom_text(size = 2) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_blue"))+
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pC <- ggplot(tile2_df, aes(x = Factor, y = 1, fill = correlation_cap, label = correlation)) +
  geom_tile(color = "black") +
  labs(x = "", fill = "") + geom_text(size = 2) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_purple"), limits = c(0,0.5))+
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cowplot::ggsave(cowplot::plot_grid(pA, pB, pC,ncol = 1), file = "../figures/output_figures/03d_screen_tiles.pdf", width = 2.9, height = 4.8)


x160A <- transformInput(importEssential("nova160A", "CBE", pad = 3), "P2A")
x160B <- transformInput(importEssential("nova160B", "CBE", pad = 3), "nCas9-BE3")
x160C <- transformInput(importEssential("nova160C", "CBE", pad = 3), "nCas9-BE4")
x160D <- transformInput(importEssential("nova160D", "CBE", pad = 3), "BE3")
x160E <- transformInput(importEssential("nova160E", "CBE", pad = 3), "BE4max")
x160F <- transformInput(importEssential("nova160F", "CBE", pad = 3), "A3A")
x160G <- transformInput(importEssential("nova160G", "CBE", pad = 3), "eA3A")
x160H <- transformInput(importEssential("nova160H", "CBE", pad = 3), "AID")

levels <- c("P2A", "nCas9-BE3", "nCas9-BE4", "BE3", "BE4max", "A3A", "eA3A", "AID")
tile3_df <- data.frame(
  Factor = factor(levels, level = levels),
  mean_edit_rate = round(sapply(list(x160A, x160B, x160C, x160D, x160E, x160F, x160G, x160H), function(x) mean(x$editRate)) * 100, 2),
  correlation = round(sapply(list(x160A, x160B, x160C, x160D, x160E, x160F, x160G, x160H), get_cor),2),
  n_hits = sapply(list(x160A, x160B, x160C, x160D, x160E, x160F, x160G, x160H), function(x) sum(x$editRate > 0.05))
)

tile3_df$mean_edit_rate_cap <- ifelse(tile3_df$mean_edit_rate > 1, 1, tile3_df$mean_edit_rate)
tile3_df$correlation_cap <- ifelse(tile3_df$correlation > 0.5, 0.5, tile3_df$correlation)
tile3_df$n_hits_cap <- ifelse(tile3_df$n_hits > 5000, 5000, tile3_df$n_hits)

pA <- ggplot(tile3_df, aes(x = Factor, y = 1, fill = mean_edit_rate_cap, label = mean_edit_rate)) +
  geom_tile(color = "black") +
  labs(x = "", fill = "") + geom_text(size = 2) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_red")) +
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pB <- ggplot(tile3_df, aes(x = Factor, y = 1, fill = n_hits_cap, label = n_hits)) +
  geom_tile(color = "black") +
  labs(x = "", fill = "") + geom_text(size = 2) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_blue"))+
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pC <- ggplot(tile3_df, aes(x = Factor, y = 1, fill = correlation_cap, label = correlation)) +
  geom_tile(color = "black") +
  labs(x = "", fill = "") + geom_text(size = 2) + 
  scale_fill_gradientn(colors = jdb_palette("brewer_purple"), limits = c(0,0.5))+
  pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cowplot::ggsave(cowplot::plot_grid(pA, pB, pC,ncol = 1), file = "../figures/output_figures/03f_nova_tiles.pdf", width = 2.9, height = 4.8)




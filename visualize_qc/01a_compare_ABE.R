library(BuenColors)
library(lmtest)

source("01_functions.R")

# Import both but also filter down to only the treated sample
ABEmax <- importEssential("89B", "CBE", chrs = 1:2)

plot_df_tile <- ABEmax %>% group_by(X1, X3) %>% 
  summarize(count_5 = sum(editRate >= 0.05), n = n(), avg_edit = mean(editRate)) %>%
  mutate(prop_5 = count_5/n)

# Get attributes for plot together
plot_df_tile$X3 <- factor(as.character(plot_df_tile$X3), levels = c("t", "g", "c", "a"))
max_val <- round(max(plot_df_tile$prop_5) + 0.05, 1)
n <- dim(ABEmax)

p2 <- ggplot(plot_df_tile, aes(x = X1, y = X3, fill = prop_5)) +
  geom_tile() + pretty_plot(fontsize = 6) + L_border() + labs(x = "5' base", y = "3' base") +
  scale_y_discrete(expand = c(0,0), labels = rev(c("A", "C", "G", "U"))) +
  scale_x_discrete(expand = c(0,0), labels = c("A", "C", "G", "U")) +
  scale_fill_gradientn(colors = jdb_palette("brewer_red"), limits = c(0,max_val)) +
  theme(legend.position = "none") + ggtitle("A3A")
p2

cowplot::ggsave(p2, file= "output_figures/01b_tile-5percent_ABE.pdf", width = 1.1, height = 1.2)


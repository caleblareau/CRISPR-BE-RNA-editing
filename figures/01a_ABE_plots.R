library(BuenColors)
library(lmtest)

source("01_functions.R")

# Import both but also filter down to only the treated sample
both <- rbind(importEssential("156B", "ABE", chrs = 1:2), importEssential("156A", "ABE", chrs = 1:2))
meta <- both %>% filter(Library == "156B")

both %>% group_by(paired, Library) %>% summarize(count = n(), gZ = sum(editRate > 0)) %>%
  mutate(rate = gZ/count)

both$paired <- factor(as.character(both$paired))
both <- both %>% mutate(four = paste0(paired,"_",as.character(Library)))
p1 <- ggplot(both, aes(x = editRate + 0.001, color = four, linetype = paired) ) +
  stat_ecdf() + scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1.0)) +
  pretty_plot(fontsize = 7) + L_border() + labs(x = "A>I editing rate", y = "Empirical cumulative density") +
  scale_linetype_manual(values=c( "FALSE" = "solid", "TRUE" = "dashed")) +
  scale_color_manual(values = c("TRUE_156A" = "red", "FALSE_156A" = "darkgrey", "FALSE_156B" = "firebrick", "TRUE_156B" = "black")) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("ABE (n = 3.3M)")
cowplot::ggsave(p1, file= "output_figures/01_ecdf_ABE.pdf", width = 2.5, height = 1.6)

plot_df_tile <- meta %>% group_by(X1, X3) %>% 
  summarize(count_5 = sum(editRate >= 0.05), n = n(), avg_edit = mean(editRate)) %>%
  mutate(prop_5 = count_5/n)

# Get attributes for plot together
plot_df_tile$X3 <- factor(as.character(plot_df_tile$X3), levels = c("t", "g", "c", "a"))
max_val <- round(max(plot_df_tile$prop_5) + 0.05, 1)
n <- dim(meta)

p2 <- ggplot(plot_df_tile, aes(x = X1, y = X3, fill = prop_5)) +
  geom_tile() + pretty_plot(fontsize = 6) + L_border() + labs(x = "5' base", y = "3' base") +
  scale_y_discrete(expand = c(0,0), labels = rev(c("A", "C", "G", "U"))) +
  scale_x_discrete(expand = c(0,0), labels = c("A", "C", "G", "U")) +
  scale_fill_gradientn(colors = jdb_palette("brewer_red"), limits = c(0,max_val)) +
  theme(legend.position = "none") + ggtitle("ABE")
cowplot::ggsave(p2, file= "output_figures/01b_tile-5percent_ABE.pdf", width = 1.1, height = 1.2)



# Demonstrate non-linearity
plot_df_tile2 <- meta %>% group_by(X1, X3, paired) %>% 
  summarize( avg_edit = mean(editRate)) %>%
  ungroup() %>% reshape2::dcast(X1 + X2 ~ paired, value.var = "avg_edit") %>%
  mutate(fold_change = `FALSE`/`TRUE`)

max_val2 <- round(max(plot_df_tile2$fold_change) + 0.05, 1)
plot_df_tile$X3 <- factor(as.character(plot_df_tile$X3), levels = c("t", "g", "c", "a"))

p3 <- ggplot(plot_df_tile2, aes(x = X1, y = X3, fill = fold_change)) +
  geom_tile() + pretty_plot(fontsize=6) + L_border() + labs(x = "5' base", y = "3' base") +
  scale_y_discrete(expand = c(0,0), labels = rev(c("A", "C", "G", "U"))) +
  scale_x_discrete(expand = c(0,0), labels = c("A", "C", "G", "U")) +
  scale_fill_gradientn(colors = jdb_palette("brewer_blue"), limits = c(1,max_val2)) +
  theme(legend.position = "none") + ggtitle("ABE")
cowplot::ggsave(p3, file= "output_figures/01c_tile-FC_ABE.pdf", width = 1.1, height = 1.2)

meta$transformed_outcome <- logit(meta$editRate)
meta_clip <- meta %>% filter(transformed_outcome < 10 & transformed_outcome > -10)
lm1 <- lm(formula = transformed_outcome ~ X1 + X2, data = meta_clip)
lm2 <- lm(formula = transformed_outcome ~ X1 + X2 + X1*X2, data = meta_clip)
str(lrtest(lm1, lm2))


n
max_val
max_val2


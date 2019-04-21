library(BuenColors)
library(data.table)
library(dplyr)

dt_ABE <- fread("ABE-outfile.tsv", stringsAsFactors = TRUE, col.names = c("sequence", "base", "position", "value"))

ABE_plot <- ggplot(dt_ABE %>% filter(base != "e"), aes(x = (position- 51), y = value, color = base, group = (position- 51))) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  geom_hline(yintercept = 0) + 
  facet_wrap(~base, nrow = 4) + scale_y_continuous(limits = c(-1, 1)) +
  pretty_plot() + L_border() + 
  theme(legend.position = "none") + scale_color_manual(values = c("green3", "blue", "orange3", "firebrick")) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Position relative to edit", y = "Deeplift Importance")


dt_CBE <- fread("CBE-outfile.tsv", stringsAsFactors = TRUE, col.names = c("sequence", "base", "position", "value"))

CBE_plot <- ggplot(dt_CBE %>% filter(base != "e"), aes(x = (position- 51), y = value, color = base, group = (position- 51))) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  geom_hline(yintercept = 0) + 
  facet_wrap(~base, nrow = 4) + scale_y_continuous(limits = c(-1, 1)) +
  pretty_plot() + L_border() + 
  theme(legend.position = "none") + scale_color_manual(values = c("green3", "blue", "orange3", "firebrick")) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Position relative to edit", y = "Deeplift Importance")

ABE_mean <- dt_ABE %>% group_by(position -51) %>%
  summarize(ABE_importance = mean(abs(value)))

CBE_mean <- dt_CBE %>% group_by(position -51) %>%
  summarize(CBE_importance = mean(abs(value)))

plot_df <- data.frame(
  position = rep(c(-50:50),2),
  importance = c(ABE_mean$ABE_importance, CBE_mean$CBE_importance),
  editor= c(rep("ABE", 101), rep("CBE", 101))
)

cowplot::ggsave(CBE_plot, filename = "plots/CBE_distribution.pdf", width = 5, height = 5)

p1 <- ggplot(plot_df, aes(x = position, y = importance, color = editor)) +
  geom_line() + pretty_plot() + L_border() + 
  geom_vline(xintercept = 0, linetype = 2) + scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  labs(x = "Position relative to edit", y = "Deeplift Absolute Importance")

cowplot::ggsave(ABE_plot, filename = "plots/ABE_distribution.pdf", width = 6, height = 6)
cowplot::ggsave(CBE_plot, filename = "plots/CBE_distribution.pdf", width = 6, height = 6)
cowplot::ggsave(p1, filename = "plots/absolute_importance_compare.pdf", width = 6, height = 3)

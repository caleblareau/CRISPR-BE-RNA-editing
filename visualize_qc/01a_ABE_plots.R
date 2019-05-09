library(BuenColors)
library(lmtest)

source("00_functions.R")

# Import both but also filter down to only the treated sample
make_ecdf <- function(editor, treated_sample, control_sample){
  meta <- readRDS(paste0("processed_df/",editor,"-",treated_sample,"-essentialDF.rds")); meta$what <- "Treated"
  control <- readRDS(paste0("processed_df/",editor,"-",control_sample,"-essentialDF.rds")); control$what <- "Control"
  both <- rbind(meta, control)
  
  both %>% group_by(paired, Library) %>% summarize(count = n(), gZ = sum(editRate > 0)) %>%
    mutate(rate = gZ/count)
  
  both$paired <- factor(as.character(both$paired))
  both <- both %>% mutate(four = paste0(paired,"_",as.character(what)))
  p1 <- ggplot(both, aes(x = editRate + 0.001, color = four, linetype = paired) ) +
    stat_ecdf() + scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1.0)) +
    pretty_plot(fontsize = 6) + L_border() + labs(x = "X>X editing rate", y = "Cumulative Density") +
    scale_linetype_manual(values=c( "FALSE" = "solid", "TRUE" = "dashed")) +
    scale_color_manual(values = c("TRUE_Control" = "red", "FALSE_Control" = "darkgrey", "FALSE_Treated" = "firebrick", "TRUE_Treated" = "black")) +
    scale_y_continuous(expand = c(0,0)) +  theme(legend.position='none')
  cowplot::ggsave(p1, file= paste0("output_figures/01_ecdf_",editor, "_", treated_sample, "_",control_sample,".pdf"), width = 1.3, height = 1.3)
}
make_ecdf("CBE", "160F", "160B")
make_ecdf("CBE", "89B", "89A")
make_ecdf("ABE", "243B", "243A")
make_ecdf("ABE", "243C", "243A")

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
cowplot::ggsave(p2, file= paste0("output_figures/01_tile-5percent_ABE.pdf"), width = 1.1, height = 1.2)



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


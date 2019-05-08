library(speedglm)
library(precrec)
library(BuenColors)
library(dplyr)
library(data.table)

format_up_down <- function(df, pad){
  if(pad == 1){
    df2 <- df %>% mutate(upstream = paste0(X5), downstream = paste0(X7))
    return(df2)
  }
  if(pad == 3){
    df2 <- df %>% mutate(upstream = paste0(X3,X4,X5), downstream = paste0(X7,X8,X9))
    return(df2)
  }
  if(pad == 5){
    df2 <- df %>% mutate(upstream = paste0(X1,X2,X3,X4,X5), downstream = paste0(X7,X8,X9,X10,X11))
    return(df2)
  }
}


fit_logistic_model <- function(editor, sample, window, specific_editor){
  df_list <- readRDS(paste0("../logistic_processed/", editor, "_", sample, "_dfs_for_logistic.rds"))
  train_df <- format_up_down(df_list[["train"]], pad = window)
  test_df <- format_up_down(df_list[["test"]], pad = window)
  
  logisticMod <- speedglm(isEdited ~ G.C + upstream + downstream + as.numeric(paired), data=train_df, family=binomial('logit'))
  predicted_train <- predict(logisticMod, train_df)
  predicted_test <- predict(logisticMod, test_df)
  
  # Compute AUROC / AUPRC
  mm_data_obj_train <- mmdata(predicted_train, train_df$isEdited)
  eval_mod_curves_train <- evalmod(mm_data_obj_train)
  df_aucs_train <- auc(eval_mod_curves_train)
  df_aucs_train
  
  mm_data_obj_test <- mmdata(predicted_test, test_df$isEdited)
  eval_mod_curves_test <- evalmod(mm_data_obj_test)
  df_aucs_test <- auc(eval_mod_curves_test)
  
  
  output_df <- data.frame(
    sample,
    editor, 
    specific_editor,
    window,
    value = c(df_aucs_train$aucs, df_aucs_test$aucs),
    metric = c("AUROC", "AUPRC", "AUROC", "AUPRC"),
    dataset = c("Train", "Train", "Test", "Test")
  )
  return(output_df)
}

# ABE Max
rbindlist(lapply(c(1,3,5), function(i) {
  print(i)
  fit_logistic_model("ABE", "243B", i, "ABEmax")
})) %>% data.frame() -> ABEmaxdf
ABEmaxdf$metric <- factor(as.character(ABEmaxdf$metric), levels = c("AUROC", "AUPRC"))
ABEmaxdf


p1 <- ggplot(ABEmaxdf, aes(x = window, y = value, fill = rev(dataset))) +
  geom_bar(stat = "identity", color = "black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1)) +
  facet_wrap(~metric) + scale_fill_manual(values = c("dodgerblue3", "lightgrey")) +
  scale_x_continuous(breaks = c(1,3,5)) +
  pretty_plot() + L_border() + ggtitle("ABEmax") +
  labs(x = "Sequence modeled around edit (bp)", y = "value", fill = "")

# Mini ABE
rbindlist(lapply(c(1,3,5), function(i) {
  fit_logistic_model("ABE", "243C", i, "miniABEmax")
})) %>% data.frame() -> miniABEmaxdf
miniABEmaxdf$metric <- factor(as.character(miniABEmaxdf$metric), levels = c("AUROC", "AUPRC"))

p2 <- ggplot(miniABEmaxdf, aes(x = window, y = value, fill = rev(dataset))) +
  geom_bar(stat = "identity", color = "black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1)) +
  facet_wrap(~metric) + scale_fill_manual(values = c("dodgerblue3", "lightgrey")) +
  scale_x_continuous(breaks = c(1,3,5)) +
  pretty_plot() + L_border() + ggtitle("miniABEmax") +
  labs(x = "Sequence modeled around edit (bp)", y = "value", fill = "")

# BE3
rbindlist(lapply(c(1,3), function(i) {
  print(i)
  fit_logistic_model("CBE", "89B", i, "BE3")
})) %>% data.frame() -> BE3df
BE3df
BE3df$metric <- factor(as.character(BE3df$metric), levels = c("AUROC", "AUPRC"))

p3 <- ggplot(BE3df, aes(x = window, y = value, fill = rev(dataset))) +
  geom_bar(stat = "identity", color = "black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1)) +
  facet_wrap(~metric) + scale_fill_manual(values = c("dodgerblue3", "lightgrey")) +
  scale_x_continuous(breaks = c(1,3,5)) +
  pretty_plot() + L_border() + ggtitle("BE3") +
  labs(x = "Sequence modeled around edit (bp)", y = "value", fill = "")

# A3A
rbindlist(lapply(c(1,3, 5), function(i) {
  print(i)
  fit_logistic_model("CBE", "160F", i, "A3A")
})) %>% data.frame() -> A3Adf
A3Adf$metric <- factor(as.character(A3Adf$metric), levels = c("AUROC", "AUPRC"))

p4 <- ggplot(A3Adf, aes(x = window, y = value, fill = rev(dataset))) +
  geom_bar(stat = "identity", color = "black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1)) +
  facet_wrap(~metric) + scale_fill_manual(values = c("dodgerblue3", "lightgrey")) +
  scale_x_continuous(breaks = c(1,3,5)) +
  pretty_plot() + L_border() + ggtitle("A3A") +
  labs(x = "Sequence modeled around edit (bp)", y = "value", fill = "")

cowplot::ggsave(p1, file = "../output/plots/logistic_ABEmax_AUC.pdf", width = 5, height = 3)
cowplot::ggsave(p2, file = "../output/plots/logistic_miniABEmax_AUC.pdf", width = 5, height = 3)
cowplot::ggsave(p3, file = "../output/plots/logistic_BE3_AUC.pdf", width = 5, height = 3)
cowplot::ggsave(p4, file = "../output/plots/logistic_A3A_AUC.pdf", width = 5, height = 3)

write.table(ABEmaxdf, file = "../output/logistic_summary_statistics/logistic_ABEmax_AUC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(miniABEmaxdf, file = "../output/logistic_summary_statistics/logistic_miniABEmax_AUC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(BE3df, file = "../output/logistic_summary_statistics/logistic_BE3_AUC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(A3Adf, file = "../output/logistic_summary_statistics/logistic_A3A_AUC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

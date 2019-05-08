library(keras)
library(tensorflow)
library(seqinr)
library(stringr)
library(reshape2)
library(precrec)
library(data.table)
library(dplyr)
library(gtools)
library(BuenColors)
source("00_functions.R")

batch_size <- 256 # batch size for training 
n_epochs <- 10 # training epochs
model_dir <- "../model_outputs/ABE-ABEmax-27April-8Channel"

# Create output directory
dir.create(model_dir)

# import the data from fasta files
all_data <- importForDeep_wStructure("243B", "ABE")

# Split into training / test data
list_two <- subset_data_to_balance(all_data)

training_data <- list_two[["train"]][,c("sequence", "structure")] 
training_labels <- list_two[["train"]][,"isEdited"] %>% as.numeric()
validation_data <-list_two[["test"]][,c("sequence", "structure")]
validation_labels <- list_two[["test"]][,"isEdited"] %>% as.numeric()

# Setup one-hot encoding w/ standardization
training_array <- make_one_hot_eight_channel(training_data)
validation_array <- make_one_hot_eight_channel(validation_data)  

# Fit CNN via keras / tensorflow
model <-keras_model_sequential() %>%   
  layer_conv_1d(filters = 50, kernel_size = 11, activation = "relu", input_shape = c(101, 8), name = "conv1") %>% 
  layer_conv_1d(filters = 50, kernel_size = 11, activation = "relu", name = "conv2") %>% 
  layer_global_average_pooling_1d() %>%
  layer_dense(units = 50, activation = "relu") %>%
  layer_dropout(rate = 0.5, name = "dropout1") %>%
  layer_dense(units = 3, activation = "relu") %>%
  layer_dense(units = 1, activation = "linear") %>%
  layer_activation("sigmoid")

# Compile the model 
model %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_rmsprop(),
  metrics = c("accuracy")
)

# Clean up memory
rm(all_data)
rm(validation_data)
rm(training_data)
gc()

# Evaluate with traing and test data
sequence_fit <- model %>% fit(
  x=training_array, training_labels,
  batch_size = batch_size,
  epochs = n_epochs,
  validation_data = list(validation_array, validation_labels),
  shuffle = TRUE
)

# Compute explicit model predictions
training_prediction <- predict_proba(model, training_array)
validation_prediction <- predict_proba(model, validation_array)

# Determine AUROC / AUPRC
training_df <- mmdata(training_prediction, training_labels, modnames = "Training_data") %>% evalmod() %>% auc()
validation_df <- mmdata(validation_prediction, validation_labels, modnames = "Validation_data") %>% evalmod() %>% auc()

# Export model attributes
sum_stats <- rbind(training_df,validation_df)[,c(1,3,4)]
sum_stats
write.table(sum_stats, file = paste0(model_dir, "/", "model-statistics.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Save predictions for ROC / PRC plotting
out_df_pred <- data.frame(
  prediction = c(training_prediction, validation_prediction), 
  label = c(training_labels, validation_labels), 
  which = c(rep("training", length(training_labels)), rep("training", length(validation_labels)))
)

saveRDS(out_df_pred, paste0(model_dir, "/", "for_AUROC_AUPRC.rds"))

# Save Keras output
p <- plot(sequence_fit, method='ggplot2', smooth=TRUE) 
cowplot::ggsave(p, filename = paste0(model_dir, "/", "training_loss_perEpoch.pdf"), width=4, height=8, units='in')
model$save(paste0(model_dir, "/", 'overall_model.h5'))
model$save_weights(paste0(model_dir, "/", 'model_weights.h5'))
write.table(model$to_yaml(), paste0(model_dir, "/", 'model_attributes.yaml'), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)




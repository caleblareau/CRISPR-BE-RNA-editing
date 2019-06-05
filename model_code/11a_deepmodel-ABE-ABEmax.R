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
n_epochs <- 3 # training epochs
model_dir <- "../model_outputs/ABE-ABEmax-5June"

# Create output directory
dir.create(model_dir)

# import the data from fasta files
all_data <- importForDeep_noStructure("243B", "ABE")

# Split into training / test data
list_two <- subset_data_to_balance(all_data)

training_data <- list_two[["train"]][,"sequence"] %>% as.character()
training_labels <- list_two[["train"]][,"isEdited"] %>% as.numeric()
validation_data <-list_two[["test"]][,"sequence"] %>% as.character()
validation_labels <- list_two[["test"]][,"isEdited"] %>% as.numeric()

# Setup one-hot encoding w/ standardization
training_array <- make_one_hot_five_channel_string(training_data) 
validation_array <- make_one_hot_five_channel_string(validation_data)  

# Fit residualbind via keras / tensorflow
input_layer1 <-  layer_input(shape = c(101,5), name = 'input1')

model_i1 <- input_layer1 %>% 
  layer_conv_1d(filters = 96, kernel_size = 72, name = "conv_branch1") %>%
  layer_batch_normalization() %>% 
  layer_activation_relu() 

model_i2  <- model_i1 %>%
  layer_conv_1d(filters = 96, kernel_size = 1) %>% 
  layer_batch_normalization(name = "a") %>%
  layer_activation_relu() %>%
  layer_dropout(rate = 0.3, name = "dropout1") %>%
  layer_conv_1d(filters = 96, kernel_size = 1) %>% 
  layer_batch_normalization(name = "b") 

# Concatenate && flatten  
predictions <- layer_add(list(model_i1, model_i2)) %>%
  layer_activation_relu() %>%
  layer_average_pooling_1d(pool_size = 10) %>%
  layer_conv_1d(filters = 196, kernel_size = 3) %>% 
  layer_batch_normalization(name = "c") %>%
  layer_activation_relu() %>%
  layer_flatten() %>%
  
  layer_dense(units = 1, activation = "sigmoid")

model <- keras_model(inputs = input_layer1, outputs = predictions)
model


# Compile the model 
model %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adam(clipnorm = 1e-6),
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
training_prediction <- predict(model, training_array)
validation_prediction <- predict(model, validation_array)

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




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
model_dir <- "../model_outputs/ABE-ABEmax-linear-13may"

# Create output directory
dir.create(model_dir)

#Import and munge
df <- readRDS("../linear_processed/ABE_243B_linearDF.rds")
df_clip <- df %>%
  mutate(transformed_outcome = logit(editRate)) %>%
  filter(transformed_outcome < 10 & transformed_outcome > -10)

# Split into training / test data
train_chr <- paste0("chr", as.character(1:15))
test_chr <- paste0("chr", as.character(16:20))
df_train <- df_clip[df_clip$chr %in% train_chr,]
df_test <- df_clip[df_clip$chr %in% test_chr,]

training_data <- df_train[,"sequence"] %>% as.character()
training_labels <- df_train[,"editRate"] %>% as.numeric() *100
validation_data <- df_test[,"sequence"] %>% as.character()
validation_labels <- df_test[,"editRate"] %>% as.numeric() *100

# Setup one-hot encoding w/ standardization
training_array <- make_one_hot_five_channel_string(training_data) 
validation_array <- make_one_hot_five_channel_string(validation_data)  

# Fit CNN via keras / tensorflow
model <- keras_model_sequential() %>%   
  layer_conv_1d(filters = 50, kernel_size = 11, activation = "relu", input_shape = c(101, 5), name = "conv1") %>% 
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

compare_2 <- data.frame(predicted = validation_prediction, observed = df_test$editRate) # inv.logit

cor(compare_2)

smoothScatter(compare_2[,1]*100, compare_2[,2]*100,  xlim = c(0,35), ylim = c(0,35),
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "ABEmax - predicted", ylab = "ABEmax - observed")

# Save Keras output
p <- plot(sequence_fit, method='ggplot2', smooth=TRUE) 
cowplot::ggsave(p, filename = paste0(model_dir, "/", "training_loss_perEpoch.pdf"), width=4, height=8, units='in')
model$save(paste0(model_dir, "/", 'overall_model.h5'))
model$save_weights(paste0(model_dir, "/", 'model_weights.h5'))
write.table(model$to_yaml(), paste0(model_dir, "/", 'model_attributes.yaml'), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)




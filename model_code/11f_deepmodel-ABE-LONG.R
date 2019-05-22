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

batch_size <- 32 # batch size for training 
n_epochs <- 5 # training epochs
model_dir <- "../model_outputs/ABE-long-21May"

# Create output directory
dir.create(model_dir)

# Eval long 
make_one_hot_five_channel_string_long <- function(x) {
  
  # Split sequences into individual characters
  spm <- str_split_fixed(x, "", 1001)
  
  # Encode a fifth channel
  spm[,501] <- "Z"
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- reshape2::melt(spm); rm$new <- 1
  aa <- acast(rm, Var1  ~  Var2  ~  value, fill = 0, value.var = "new")  
  aa
}


importForDeep_noStructure_long <- function(sample){
  
  saveForDeepLift <- c("chr21", "chr22", "chrX")
  
  # Import sequence fasta files
  dir_seq <- paste0("../../crispr-be-rnaediting-data-long/ABE")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas) & !grepl(paste(saveForDeepLift,collapse="|"), seq_fastas)]
  seq_fastas <- seq_fastas[c(11,15, 7)]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
  # Get all sequence attributes
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  meta$chr <- stringr::str_split_fixed(meta$chr_pos, "-", 2)[,1]
  seqs <- unlist(fasta_input) %>% unname()
  
  # Process meta to convert characters to numeric
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  rownames(meta) <- NULL
  meta$sequence <- seqs
  meta$isEdited <- meta$editRate > 0.05
  boo <- meta$editRate > 0.05 | meta$editRate == 0
  return(meta[boo,])
}

# import the data from fasta files
all_data <- importForDeep_noStructure_long("long243B")

# Split into training / test data
list_two <- subset_data_to_balance(all_data, n_excess = 10)

training_data <- list_two[["train"]][,"sequence"] %>% as.character()
training_labels <- list_two[["train"]][,"isEdited"] %>% as.numeric()
validation_data <-list_two[["test"]][,"sequence"] %>% as.character()
validation_labels <- list_two[["test"]][,"isEdited"] %>% as.numeric()

# Setup one-hot encoding w/ standardization
training_array <- make_one_hot_five_channel_string_long(training_data) 
validation_array <- make_one_hot_five_channel_string_long(validation_data)  

# Fit CNN via keras / tensorflow
RB <- function(x, N = 32, W = 11){
  layer_batch_normalization(x) %>%
    layer_activation_relu() %>%
    layer_conv_1d(filters = N, kernel_size = W) %>% 
    layer_batch_normalization() %>%
    layer_activation_relu() %>%
    layer_conv_1d(filters = N, kernel_size = W)
}

# Do the SpliceAI model 
input_layer1 <-  layer_input(shape = c(1001,5), name = 'input1')
model_i <- input_layer1  %>%
  layer_conv_1d(filters = 32, kernel_size = 1, name = "conv1")
model_i1 <- model_i %>%  layer_conv_1d(filters = 32, kernel_size = 1, name = "conv_branch1")

model_i2  <- model_i %>%
  RB() %>%
  RB() %>%
  RB() %>%
  RB() %>%
  layer_conv_1d(filters = 32, kernel_size = 11)

# Concatenate && flatten  
predictions <- layer_concatenate(list(model_i1, model_i2), axis = 1) %>%
  layer_flatten(name = "flatten1") %>% 
  layer_dropout(rate = 0.3, name = "dropout2") %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dropout(rate = 0.3, name = "dropout3") %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dropout(rate = 0.3, name = "dropout4") %>%
  layer_dense(units = 1, activation = "linear") %>%
  layer_activation(activation = "sigmoid")

model <- keras_model(inputs = input_layer1, outputs = predictions)

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
  validation_data = list(training_array, training_labels),
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
  prediction = c(training_prediction, training_prediction), 
  label = c(training_labels, training_labels), 
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


# Try to look at what is still failing
plot_df <- data.frame(
  str_split_fixed(substr(list_two[["test"]][,c("sequence")], 500, 502), "", 3),
  predicted = validation_prediction,
  annotation = list_two[["test"]][,"annotation"],
  truth = as.factor(validation_labels)
)
colnames(plot_df) <- c("x5p", "e", "x3p", "predictedValue", "annotation","isEdited")

p1 <- ggplot(plot_df, aes(x = isEdited, y = predictedValue)) +
  facet_grid(~x5p) +
  geom_violin() + pretty_plot() + 
  L_border() + ggtitle("Stratified by 5' base")

cowplot::ggsave(p1, file = "../output/plots/violins_sample_deep.pdf", width = 5, height = 2)

plot_df %>% group_by(x5p, annotation, isEdited) %>%
  summarize(count = n(), percentile20score = quantile(predictedValue, 0.2), medianScore = median(predictedValue)) %>%
  ungroup() %>% group_by(annotation) %>% mutate(anno_prop = count / sum(count))

p5_df <- plot_df %>% filter(x5p == "g")
mmdata(p5_df$predictedValue, p5_df$isEdited, modnames = "Validation_data") %>% evalmod() %>% auc()




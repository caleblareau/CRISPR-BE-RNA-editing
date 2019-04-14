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

batch_size <- 256 # batch size for training 
n_epochs <- 6 # training epochs
model_dir <- "../model_outputs/ABE-19FEB"

# Create output directory
dir.create(model_dir)

importForDeep <- function(sample, editor){
  # Import sequence fasta files
  dir_seq <- paste0("../fastas/", editor, "-sequence")
  dir_struct <- paste0("../fastas/", editor, "-secondary")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas)]
  struct_files <- list.files(dir_struct, pattern = paste0(sample), full.names = TRUE)
  struct_files <- struct_files[grepl(".gz", struct_files)]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  
  # Process structure data
  struct_input <- sapply(struct_files, function(faf) fread(cmd = paste0("zcat < ", faf), header = FALSE)[[1]][c(FALSE, TRUE)]) %>% unlist() %>% unname()
  
  # Get all sequence attributes
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  seqs <- unlist(fasta_input) %>% unname()
  
  # Process meta to convert characters to numeric
  meta$editRate <- as.numeric(meta$editRate); meta$coverage <- as.numeric(meta$coverage)
  rownames(meta) <- NULL
  meta$sequence <- seqs; meta$structure <- struct_input
  meta_clip <- meta %>% 
    mutate(transformed_outcome = logit(editRate)) %>%
    filter(transformed_outcome < 10 & transformed_outcome > -10)
  return(meta_clip)
}

all_data <- importForDeep("156B", "ABE")

# Randomly subset the data
idx <- sample(1:(dim(all_data)[1]), round(dim(all_data)*0.8))

training_data <- all_data[idx,c('sequence', 'structure')]
training_labels <- all_data[idx,'transformed_outcome']
validation_data <- all_data[-idx,c('sequence', 'structure')]
validation_labels <- all_data[-idx,'transformed_outcome']

make_one_hot_eight_channel <- function(x) {
  
  # Split into each individual matrix
  seq_mat <-str_split_fixed( x[,1],  "", 101)
  struct_mat <- str_split_fixed(x[,2],  "", 101)
  
  # Encode a fifth channel for sequence -- editing event
  seq_mat[,51] <- "Z"
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- rbind(reshape2::melt(seq_mat), reshape2::melt(struct_mat)); rm$new <- 1
  aa <- acast(rm, Var1  ~  Var2  ~  value, fill = 0, value.var = "new")  
  aa
}


make_one_hot_five_channel <- function(x) {
  
  # Split sequences into individual characters
  spm <- str_split_fixed(x[,1], "", 101)
  
  # Encode a fifth channel
  spm[,51] <- "Z"
  
  # Rehsape into nseqs x sequence length x nucleotide
  rm <- reshape2::melt(spm); rm$new <- 1
  aa <- acast(rm, Var1  ~  Var2  ~  value, fill = 0, value.var = "new")  
  aa
}


# Setup one-hot encoding
# the minus helps standardize the mean for early epochs
training_array <- make_one_hot_five_channel(training_data) 
validation_array <- make_one_hot_five_channel(validation_data) 

ray_kernel <- c(8,5,5,3)
kkern <- ray_kernel

model <- keras_model_sequential() %>%   
  layer_conv_1d(filters = 1024, kernel_size = kkern[1], activation = "relu", input_shape = c(101, 5), name = "conv1") %>% 
  layer_max_pooling_1d(pool_size=1, name='downsample1') %>%
  layer_conv_1d(filters = 64, kernel_size = kkern[2], activation = "relu", name = "conv2") %>% 
  layer_max_pooling_1d(pool_size=1, name='downsample2') %>%
  layer_dropout(rate = 0.3, name = "dropout1") %>%
  layer_conv_1d(filters = 64, kernel_size = kkern[3], activation = "relu", name = "conv3") %>% 
  layer_conv_1d(filters = 64, kernel_size = kkern[4], activation = "relu", name = "conv4") %>% 
  layer_max_pooling_1d(pool_size=1, name='downsample3') %>%
  layer_flatten(name = "flatten1") %>% 
  layer_dropout(rate = 0.3, name = "dropout2") %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dropout(rate = 0.3, name = "dropout3") %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dropout(rate = 0.3, name = "dropout4") %>%
  layer_dense(units = 1, activation = "linear")

# Compile the model 
model %>% compile(
  loss = "mse",
  optimizer = optimizer_rmsprop(),
  metrics = list("mean_absolute_error")
)

# Evaluate with traing and test data
sequence_fit <- model %>% fit(
  x=training_array, training_labels,
  batch_size = batch_size,
  epochs = n_epochs,
  validation_data = list(validation_array, validation_labels),
  shuffle = FALSE
)

# Compute explicit model predictions
training_prediction <- predict(model, training_array)
validation_prediction <- predict(model, validation_array)

ABE <- data.frame(predicted = inv.logit(c(training_prediction, validation_prediction)),
                  observed = inv.logit(c(training_labels, validation_labels))) %>% 
  filter(predicted > 0.01 | observed > 0.01)

cor(ABE)

smoothScatter(ABE[,1]*100, ABE[,2]*100,  xlim = c(0,35), ylim = c(0,35),
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "Deep - ABE - predicted", ylab = "Deep - ABE - observed")


# Export model attributes
p <- plot(sequence_fit, method='ggplot2', smooth=TRUE) 
cowplot::ggsave(p, filename = paste0(model_dir, "/", "training_loss_perEpoch.pdf"), width=4, height=8, units='in')
model$save(paste0(model_dir, "/", 'overall_model.h5'))
model$save_weights(paste0(model_dir, "/", 'model_weights.h5'))
write.table(model$to_yaml(), paste0(model_dir, "/", 'model_attributes.yaml'), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)




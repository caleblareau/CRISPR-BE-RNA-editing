library(seqinr)
library(caret)
library(speedglm)
library(precrec)
library(BuenColors)
library(stringr)
library(gtools)
set.seed(100)  # for repeatability of samples

# Import fasta sequences
rep1 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/weighted_fastas/JUL89_weights.fasta.gz", as.string = TRUE)
rep2 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/weighted_fastas/JUL90_weights.fasta.gz", as.string = TRUE)

master_df <- data.frame(
  upstream = c(substr(unlist(rep1), 48, 50), substr(unlist(rep2), 48, 50)), 
  downstream = c(substr(unlist(rep1), 52, 53), substr(unlist(rep2), 52, 53)),
  GC = ((c(substr(unlist(rep1), 41, 61), substr(unlist(rep2), 41, 61)) %>% DNAStringSet() %>% letterFrequency( "GC"))/20)[,1]
)

# Append additional information in the read name
stat2 <- rbind(str_split_fixed(names(rep1), "_", 3), str_split_fixed(names(rep2), "_", 3))[,c(2,3)] 
master_df$edit_raw <- as.numeric(stat2[,1])
master_df$coverage <- as.numeric(stat2[,2])
master_df$edit <- logit(master_df$edit_raw)

# Filter extreme values
master_df <- master_df %>% filter(edit_raw > 0.005 & edit_raw < 0.995)
rm(rep1); rm(rep2)

# Split train /test
keep_train_boo <- 1:dim(master_df)[1] %in% sample(1:dim(master_df)[1], 0.7 * dim(master_df)[1])
train_data <- master_df[keep_train_boo,]
test_data <- master_df[!keep_train_boo,]

linearMod <- speedlm(edit ~ GC + upstream + downstream, data=train_data)
predicted_linear_model <- predict(linearMod, test_data)

if(FALSE){
  smoothScatter(x = inv.logit(predicted_linear_model), y = inv.logit(test_data$edit), xlim = c(0,1), ylim = c(0,1))
}



# For Apobec + the SECURE variants
if(FALSE){
  
  p1 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/variant_fasta/JUL146C_evaluate.fasta", as.string = TRUE)
  p2 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/variant_fasta/JUL147C_evaluate.fasta", as.string = TRUE)
  p3 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/variant_fasta/JUL148C_evaluate.fasta", as.string = TRUE)
  
  R33 <- data.frame(
    upstream = c(substr(unlist(p1), 48, 50), substr(unlist(p2), 48, 50),  substr(unlist(p3), 48, 50)), 
    downstream = c(substr(unlist(p1), 52, 53), substr(unlist(p2), 52, 53), substr(unlist(p3), 52, 53)),
    GC = ((c(substr(unlist(p1), 41, 61), substr(unlist(p2), 41, 61), substr(unlist(p3), 41, 61)) %>% DNAStringSet() %>% letterFrequency( "GC"))/20)[,1] 
  )
  stat2 <- rbind(str_split_fixed(names(p1), "_", 3), str_split_fixed(names(p2), "_", 3), str_split_fixed(names(p3), "_", 3))[,c(2,3)] 
  R33$edit_raw <- as.numeric(stat2[,1])
  R33$coverage <- as.numeric(stat2[,2])
  R33$edit <- logit(R33$edit_raw)
  predict_R33 <- predict(linearMod, R33)
  R33$score <- predict_R33
  R33$predicted <- inv.logit(predict_R33)
  R33$seq <-  c(substr(unlist(p1), 30, 53), substr(unlist(p2), 30, 53),  substr(unlist(p3), 30, 53))
  
  R33 %>% filter(downstream %in% c("aa", "ga")) %>% 
    arrange(desc(predicted)) %>% head(20)
  
  p1 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/variant_fasta/JUL146D_evaluate.fasta", as.string = TRUE)
  p2 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/variant_fasta/JUL147D_evaluate.fasta", as.string = TRUE)
  p3 <- read.fasta("/data/joung/caleb/base_editing/exp2-hiseq/testCNN/variant_fasta/JUL148D_evaluate.fasta", as.string = TRUE)
  
  R33K34 <- data.frame(
    upstream = c(substr(unlist(p1), 48, 50), substr(unlist(p2), 48, 50),  substr(unlist(p3), 48, 50)), 
    downstream = c(substr(unlist(p1), 52, 54), substr(unlist(p2), 52, 54), substr(unlist(p3), 52, 54)),
    GC = ((c(substr(unlist(p1), 41, 61), substr(unlist(p2), 41, 61), substr(unlist(p3), 41, 61)) %>% DNAStringSet() %>% letterFrequency( "GC"))/20)[,1] 
  )
  predict_R33K34 <- predict(linearMod, R33K34)
  
  baseline_df <- data.frame(
    predicted_value = predicted_linear_model, 
    Truth =  c("Test-data")
  )
  
  variant_df <- data.frame(
    predicted_value = c(predict_R33, predict_R33K34),
    Truth = c(rep("R33", length(predict_R33)), rep("R33K34", length(predict_R33K34)))
  )
  
  plot_df <- rbind(baseline_df, variant_df)
  
  pDensity <- ggplot(plot_df, aes(x = predicted_value, color = Truth)) +
    geom_density() + scale_color_manual(values = c("firebrick", "dodgerblue3", "green3", "purple3")) +
    pretty_plot() + L_border() + labs(x = "Predicted value", y = "density") +
    ggtitle("Cytosine base editor + Engineered Variants")
  
  cowplot::ggsave(pDensity, file = "../plots/cytosine-density-wSECURE-linear.pdf", width = 5, height = 3)
  
  mat <- data.frame(summary(linearMod)$coefficients)
  mat$what <- rownames(summary(linearMod)$coefficients)
  mat %>% arrange(desc(coef)) %>% head()
}

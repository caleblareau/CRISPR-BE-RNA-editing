library(seqinr)
library(caret)
library(speedglm)
library(BuenColors)
library(stringr)
library(gtools)


# Import essential meta features
df <- readRDS("../linear_processed/ABE_243B_linearDF.rds")

df_clip <- df %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = substring(sequence, 48,50), downstream = substring(sequence, 52,54)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

train_chr <- paste0("chr", as.character(1:15))
test_chr <- paste0("chr", as.character(16:20))

df_train <- df_clip[df_clip$chr %in% train_chr,]
df_test <- df_clip[df_clip$chr %in% test_chr,]

linearMod <- speedlm(editRate ~upstream + downstream + as.numeric(paired), data=df_train)
predicted_edit <- predict(linearMod, df_test)
compare_2 <- data.frame(predicted = (predicted_edit), observed = df_test$editRate) # inv.logit

cor(compare_2)
smoothScatter(compare_2[,1]*100, compare_2[,2]*100,  xlim = c(0,35), ylim = c(0,35),
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "ABEmax - predicted", ylab = "ABEmax - observed")

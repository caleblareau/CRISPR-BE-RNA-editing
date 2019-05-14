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

linearMod <- speedlm(transformed_outcome ~ G.C + upstream + downstream + as.numeric(paired), data=ABE1_clip)
predicted_ABE <- predict(linearMod, ABE2_clip)
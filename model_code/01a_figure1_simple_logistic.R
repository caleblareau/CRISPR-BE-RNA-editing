library(seqinr)
library(caret)
library(speedglm)
library(precrec)
library(BuenColors)
library(stringr)
library(gtools)
library(precrec)


source("../figures/01_functions.R")

# Import essential meta features
ABE1 <- importEssential("243C", "ABE", pad = 3, chrs = 1:10)
ABE2 <- importEssential("244C", "ABE", pad = 3, chrs = 11:20)

ABE1_clip <- ABE1 %>% 
  mutate(is_edited = editRate > 0) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(editRate > 0.05 | editRate == 0)

ABE2_clip <- ABE2 %>% 
  mutate(is_edited = editRate > 0) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(editRate > 0.05 | editRate == 0)

logisticMod <- speedglm(is_edited ~ G.C + upstream + downstream + as.numeric(paired), data=ABE1_clip,  family=binomial('logit'))
predicted_ABE <- predict(logisticMod, ABE2_clip)
mm_data_obj <- mmdata(predicted_ABE, ABE2_clip$is_edited)

eval_mod_curves <- evalmod(mm_data_obj)
df_aucs <- auc(eval_mod_curves)
df_aucs


#------------
# CBE predict
#------------

# Import essential meta features
CBE1 <- importEssential("160F", "CBE", pad = 3, chrs = 1:5)
CBE2 <- importEssential("161F", "CBE", pad = 3, chrs = 11:15)

CBE1_clip <- CBE1 %>% 
  mutate(is_edited = editRate > 0) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(editRate > 0.05 | editRate == 0)

CBE2_clip <- CBE2 %>% 
  mutate(is_edited = editRate > 0) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(editRate > 0.05 | editRate == 0)

logisticMod <- speedglm(is_edited ~ G.C + upstream + downstream + as.numeric(paired), data=CBE1_clip,  family=binomial('logit'))
predicted_CBE <- predict(logisticMod, CBE2_clip)
mm_data_obj <- mmdata(predicted_CBE, CBE2_clip$is_edited)

eval_mod_curves <- evalmod(mm_data_obj)
df_aucs <- auc(eval_mod_curves)
df_aucs



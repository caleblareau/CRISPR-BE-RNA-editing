library(seqinr)
library(caret)
library(speedglm)
library(precrec)
library(BuenColors)
library(stringr)
library(gtools)

source("../figures/01_functions.R")

# Import essential meta features
ABE1 <- importEssential("156B", "ABE", pad = 3)
ABE2 <- importEssential("157B", "ABE", pad = 3)

ABE1_clip <- ABE1 %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

ABE2_clip <- ABE2 %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

<
ABE <- data.frame(predicted = inv.logit(predicted_ABE), observed = ABE2_clip$editRate) %>% 
  filter(predicted > 0.01 | observed > 0.01)

cor(ABE[,1], ABE[,2])

pdf("../figures/output_figures/01e_ABE-Rep2-predicted.pdf", width = 3, height = 3)
smoothScatter(ABE[,1]*100, ABE[,2]*100,  xlim = c(0,35), ylim = c(0,35),
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "ABE - predicted", ylab = "ABE - observed")
dev.off()

#------------
# CBE predict
#------------

# Import essential meta features
CBE1 <- importEssential("89B", "CBE", pad = 3)
CBE2 <- importEssential("90B", "CBE", pad = 3)

CBE1_clip <- CBE1 %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

CBE2_clip <- CBE2 %>% 
  mutate(transformed_outcome = logit(editRate)) %>%
  mutate(upstream = paste0(X1,X2,X3), downstream = paste0(X5,X6,X7)) %>% 
  filter(transformed_outcome < 10 & transformed_outcome > -10)

linearMod <- speedlm(transformed_outcome ~ G.C + upstream + downstream + as.numeric(paired), data=CBE1_clip)
predicted_CBE <- predict(linearMod, CBE2_clip)
CBE <- data.frame(predicted = inv.logit(predicted_CBE), observed = CBE2_clip$editRate) %>% 
  filter(predicted > 0.01 | observed > 0.01)

cor(CBE[,1], CBE[,2])

pdf("../figures/output_figures/01e_CBE-Rep2-predicted.pdf", width = 3, height = 3)
smoothScatter(CBE[,1]*100, CBE[,2]*100,  xlim = c(0,70), ylim = c(0,70),
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "CBE - predicted", ylab = "CBE - observed")
dev.off()


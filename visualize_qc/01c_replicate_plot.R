library(BuenColors)

source("01_functions.R")

ABE1 <- importMetaOnly("156B", "ABE")
ABE2 <- importMetaOnly("157B", "ABE")

ABE <- inner_join(ABE1[,c("chr_pos", "editRate")], ABE2[,c("chr_pos", "editRate")], by = "chr_pos")

ABE %>% filter(editRate.x > 0.01 | editRate.y > 0.01) -> ABE
pdf("output_figures/01d_ABE-replicate.pdf", width = 3, height = 3)
smoothScatter(ABE[,2]*100, ABE[,3]*100, xlim = c(0,35), ylim = c(0,35),
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "ABE - Rep. 1 (A>I %)", ylab = "ABE - Rep. 2 (A>I %)")
dev.off()

cor(ABE[,2], ABE[,3])

CBE1 <- importMetaOnly("89B", "CBE")
CBE2 <- importMetaOnly("90B", "CBE")

CBE <- inner_join(CBE1[,c("chr_pos", "editRate")], CBE2[,c("chr_pos", "editRate")], by = "chr_pos")

CBE %>% filter(editRate.x > 0.01 | editRate.y > 0.01) -> CBE

pdf("output_figures/01d_CBE-replicate.pdf", width = 3, height = 3)
smoothScatter(CBE[,2]*100, CBE[,3]*100,  xlim = c(0,70), ylim = c(0,70),
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "CBE - Rep.1 (C>T %)", ylab = "CBE - Rep.2 (C>T %)")
dev.off()


cor(CBE[,2], CBE[,3])

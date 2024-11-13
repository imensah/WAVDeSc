library(SAVER)
saver <- read.csv("~/path/to/data/saver.csv")
denoise_saver <- saver(saver, ncores = 4)
saverdenoised <- denoise.saver$estimate
View(saverDenoised)

write.csv(saverDenoised,"~/path/to/save/denoised_saver5000.csv",
           row.names = FALSE)

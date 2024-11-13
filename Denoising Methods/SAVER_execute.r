library(SAVER)
saver <- read.csv("~/Desktop/Final_datasets_used/5000*3000_2pop/
                  nonUMI/saver/saver.csv")
denoise_saver <- saver(saver, ncores = 4)
saverdenoised <- denoise.saver$estimate
View(saverDenoised)

write.csv(saverDenoised,"~/Desktop/Final_datasets_used/5000*3000_2pop/
                    nonUMI/saver/denoised_saver5000.csv", row.names = FALSE)

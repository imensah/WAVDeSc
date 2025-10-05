# library(SAVER)
# saver <- read.csv("~/path/to/data/saver.csv")
# denoise_saver <- saver(saver, ncores = 4)
# saverdenoised <- denoise.saver$estimate
# View(saverDenoised)
# 
# write.csv(saverDenoised,"~/path/to/save/denoised_saver5000.csv",
#            row.names = FALSE)

library(SAVER)

repo_dir = "/home/jay/Desktop/WAVDeSc_repo" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/Data Simulations/Simulations.r"))
# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"

work_dir = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations" 
setwd(work_dir)


algorithm = "saver"
inputDir = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations/"
outputDir = paste0("denoised_",algorithm,"/")


# if(!dir.exists(outputDir)){dir.create(outputDir)}
# 
# # Get a list of all files in the input directory that contain "noisy" in their name
noisy_files <- list.files(inputDir, pattern = "noisy")
 
# write("experiment,method,cpuTime,elapsedTime,RAMpeak(MiB)", 
#       file=paste0(outputDir,"runtime.csv"), append=T)

for (files in noisy_files ){

  # Extract the experiment name by removing "noisy" and ".csv"
  experiment_name = gsub("noisy|\\.csv$", "", files)
  print(experiment_name)
  tmp = read.csv(file=paste0(inputDir, files))
  View(files)

  # denoise_saver <- saver(saver, ncores = 4)
  # saverdenoised <- denoise.saver$estimate
  # View(saverDenoised)
  # 
  # write.csv(saverDenoised,"~/path/to/save/denoised_saver5000.csv",
  #            row.names = FALSE)
  denoise_saver <-saver(tmp, ncores = 3)
  saverdenoised <- denoise.saver$estimate
  View(saverdenoised )
}





data <- read.csv("/home/jay/Desktop/WAVDeSc_repo/Data Simulations/Noisy_files/noisy counts UMI (100 x 50).csv")


data_matrix <- as.matrix(data)

attr(data_matrix,"class")<-"matrix"

denoise_saver <-saver(data_matrix, ncores = 3)
str(denoise_saver)
saverdenoised <- denoise_saver$estimate
View(saverdenoised )


# Ensure that data_filtered is numeric
data_filtered <- apply(data_filtered, 2, as.numeric)

# Check for any NA or NaN values and handle them
if(any(is.na(data_filtered))){
  cat("There are NA values in the data. Replacing them with 0.\n")
  data_filtered[is.na(data_filtered)] <- 0
}

# Confirm the structure and type of data_filtered
str(data_filtered)

# Run SAVER on the filtered data
denoise_saver <- saver(data_filtered, ncores = 3)







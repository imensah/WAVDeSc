library(SAVER)
library(Matrix)

repo_dir = "/home/jay/Desktop/WAVDeSc_repo" # Set this to the folder containing the cloned repository (i.e. local path to "review-scRNA-seq-DENOISING")
source(paste0(repo_dir, "/Data Simulations/Simulations.r"))
# Set this to the path of the folder containing the csv of datasets with noise
# e.g., to denoise synthetic datasets with the most variables genes, set this to "localPath/review-scRNA-seq-DENOISING/final_data/csv_mostvar"

work_dir = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations" 
setwd(work_dir)


algorithm = "saver"
inputDir = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations/"
outputDir = paste0("denoised_",algorithm,"/")


# Get a list of all files in the input directory that contain "noisy" in their name
noisy_files <- list.files(inputDir, pattern = "noisy")

for (files in noisy_files ){

  # Extract the experiment name by removing "noisy" and ".csv"
  experiment_name = gsub("noisy|\\.csv$", "", files)
  print(experiment_name)
  denoised_name <- paste0(outputDir, experiment_name, "_denoised_",algorithm,".csv")
  if(file.exists(denoised_name)){next}
   tmp = read.csv(file=paste0(inputDir, files))
#   tmp <- as.matrix(tmp)
#   View(tmp)

    # if (!any(grepl("matrix", class(x), ignore.case = TRUE))) {
    #     tmp = read.csv(file=paste0(inputDir, files))
    #     tmp <- as.matrix(tmp)
    # }



#  denoise_saver <- saver(tmp, ncores = 3)
#  saverdenoised <- denoise_saver$estimate
  denoise_saver <-saver(tmp, ncores = 3)
  saverdenoised <- denoise.saver$estimate
  View(saverdenoised )
}
class(tmp)

tmp_list<- lapply(tmp,function(x) as.numeric(x))
tmp_list
dim(tmp_list)
tmp_num <- apply(tmp,function(x) as.nummeric(x))
tmp_num <- apply(tmp,MARGIN = 1, FUN=sum)
dim(tmp_num)
class(tmp)


tmp_mat = as.matrix(tmp)
class(tmp_mat)

AA = SAVER::saver(as.numeric(tmp),ncores = 3)
str(tmp_mat)


tmp_data <- read.csv("/home/jay/Desktop/WAVDeSc_repo/Data Simulations/noisy counts UMI (100 x 50).csv")
# tmp_data_mat <- as.matrix(tmp_data)

tmp_data <- apply(tmp_data, 2, as.numeric)


# class(tmp_data_mat)
str(tmp_data)
View(tmp_data)

Ab <- SAVER::saver(tmp_data, ncores = 3)
str(tmp_data_mat)
type(tmp_data_mat)



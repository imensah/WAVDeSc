library(SAVER)
library(parallel)

cl <- makeCluster(3)  # Create a cluster with 3 cores
registerDoParallel(cl)  # Register the cluster for parallel operations

repo_dir = "/home/jay/Desktop/WAVDeSc_repo" # Set this to the folder containing the cloned repository (i.e. local path to "WAVDeSc_repo")
source(paste0(repo_dir, "/Data Simulations/Simulations.r"))
# Set this to the path of the folder containing the csv of datasets with noise
work_dir = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations" 
setwd(work_dir)


algorithm = "saver"
inputDir = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations/Noisy_files/"
outputDir = paste0("denoised_",algorithm,"/")


if(!dir.exists(outputDir)){dir.create(outputDir)} #Creates directory if it doesn't exist

noisy_files <- list.files(inputDir, pattern = "noisy")

write("experiment,method,cpuTime,elapsedTime,RAMpeak(MiB)", 
      file=paste0(outputDir,"runtime.csv"), append=T)

peakRAM_custom <- function(...){
  # Capture R expressions or function calls
  args <- c(as.list(match.call(expand.dots = FALSE)$`...`), NULL)
  Function <- sapply(args, function(e) paste0(deparse(e), collapse = ""))
  Function <- sapply(Function, function(e) gsub(" ", "", e))
  # Initialize containers for output
  numfunc <- length(args)
  totaltime <- vector("numeric", numfunc)
  usertime <- vector("numeric", numfunc)
  RAMused <- vector("numeric", numfunc)
  RAMpeak <- vector("numeric", numfunc)
  i <- 1
  for(arg in args){
    # Reset garbage collector and save baseline
    start <- gc(verbose = FALSE, reset = TRUE)
    # Evaluate regular and anonymous functions
    evalTime <- system.time(result <- eval.parent(arg))
    if(inherits(result, "function")){
      evalTime <- system.time(output <- result())
      rm(result)
    }else{
      output <- result
      rm(result)
    }
    # Call garbage collector and save post-eval
    end <- gc(verbose = FALSE, reset = FALSE)
    rm(output)
    # Calculate total and peak RAM used
    totaltime[i] <- as.numeric(evalTime["elapsed"])
    usertime[i] <- as.numeric(evalTime["user.self"])
    RAMused[i] <- end[2, 2] - start[2, 2]
    RAMpeak[i] <- end[2, 6] - start[2, 6]
    i <- i + 1
  }
  data.frame("Function_Call" = Function,
             "Elapsed_Time_sec" = totaltime,
             "user_Time_sec" = usertime,
             "Total_RAM_Used_MiB" = RAMused,
             "Peak_RAM_Used_MiB" = RAMpeak,
             row.names = 1:numfunc)
}

#Looping through files to be denoised
for (files in noisy_files ){
  # Extract the experiment name by removing "noisy counts" and ".csv"
  experiment_name = gsub("noisy counts|\\.csv$", "", files)
  print(experiment_name)
  # Define denoised file path
  denoised_path <- paste0(outputDir, experiment_name, "_denoised_",algorithm,".csv")
  
  if(file.exists(denoised_path)){next} # Skip if denoised file already exists
  
  # Load the Data
  data = read.csv(file=paste0(inputDir, files))
  data_matrix <- as.matrix(data)
  attr(data_matrix,"class")<-"matrix"
  
  # Run the saver algorithm
  print(paste0("Running ", algorithm))
  start = Sys.time() # record start time

  t <- try(tmpD <- saver(data_matrix, ncores = 3), silent = T)
  peakMemory = peakRAM_custom(t)
  if(class(t) != "try-error") {
    # Record the execution time and memory usage
    end = Sys.time()
    duration = end-start
    duration_sec = as.numeric(end) - as.numeric(start)
    print("Saving results")

    # Round and save denoised data
    tmpD <- round(tmpD$estimate, digits = 2)
    write.table(tmpD, file = denoised_path, sep = ",", col.names = T, row.names = T)
    write(paste(experiment_name,
                algorithm,
                peakMemory$user_Time_sec,
                peakMemory$Elapsed_Time_sec,
                peakMemory$Peak_RAM_Used_MiB,  sep=','),
          file=paste0(outputDir,"runtime.csv"), append=T)
  }else {
    print("ERROR!")
    write(paste(experiment_name, algorithm, "error", "error", "error", sep=','), file=paste0(outputDir,"runtime.csv"), append=T)
  }
}

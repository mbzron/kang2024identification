use_existing_files <- TRUE
data_subset_proportion <- 0.25

start <- proc.time()

source("kang2024_testing_1.r")

rm(list = ls())

source("kang2024_testing_2.r")

source("kang2024_testing_3.r")

rm(list = ls())

source("kang2024_testing_4.r")

end <- proc.time()

#display elapsed time
elapsed_time <- end - start
print(paste0("Total elapsed time: ", elapsed_time))
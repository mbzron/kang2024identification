source("logging.r")
library(copykat)
n_cores <- 4

sam_name <- "LUAD"
sample.name <- paste(sam_name, "_copykat_", sep = "")
luad_filename <- paste(sample.name, "CNA_raw_results.txt", sep = "")
log_message(paste("Loading existing", luad_filename, "..."))
RNA.copycat <- read.delim(luad_filename, sep = "\t", header = TRUE)

Aj <- copykat::convert.all.bins.hg20(
  DNA.mat = copykat::DNA.hg20,
  RNA.mat = RNA.copycat,
  n.cores = n_cores
)

saveRDS(Aj, file = "Aj.rds")
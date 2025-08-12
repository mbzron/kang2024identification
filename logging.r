# Define log file
log_file <- "my_log.txt"

# Function to write to log
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_message <- paste0("[", timestamp, "] ", message, "\n")
  print(full_message)
  cat(full_message, file = log_file, append = TRUE)
}
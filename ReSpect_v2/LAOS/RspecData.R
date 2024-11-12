# Function to extract strain amplitude and omega from filename
extract_parameters <- function(filename) {
  # Use regular expressions to extract the strain amplitude and omega (frequency)
  strain_amp <- as.numeric(sub(".*_(\\d+)%Om\\d+.*", "\\1", filename))  # Extract strain amplitude like '5%' or '100%'
  omega <- as.numeric(sub(".*Om(\\d+).*", "\\1", filename))  # Extract omega like 'Om1' or 'Om50'
  
  if (is.na(strain_amp) || is.na(omega)) {
    stop("Error in extracting parameters from filename: ", filename)
  }
  
  return(list(strain_amp = strain_amp, omega = omega))
}

# Function to process each file
process_file <- function(filepath) {
  # Extract strain amplitude and omega from filename
  filename <- basename(filepath)
  params <- extract_parameters(filename)
  omega <- params$omega
  
  # Read the data from the file
  data <- read.table(filepath, header = FALSE, sep = "\t")
  
  # Assign column names: time, total stress, strain rate, viscous stress, elastic stress
  colnames(data) <- c("time", "total_stress", "strain_rate", "viscous_stress", "elastic_stress")
  
  # Calculate G' and G'' (normalized by strain rate)
  G_prime <- data$elastic_stress / data$strain_rate
  G_double_prime <- data$viscous_stress / data$strain_rate
  
  # Return a data frame with omega, Gi(wi) = G_prime, and Gii(wi) = G_double_prime
  result <- data.frame(omega = omega, Gi_wi = mean(G_prime), Gii_wi = mean(G_double_prime))
  return(result)
}

# Main function to process multiple files in a directory
process_files_in_directory <- function(directory_path) {
  # List all files in the directory
  file_list <- list.files(directory_path, full.names = TRUE)
  
  # Process each file and combine results
  results <- do.call(rbind, lapply(file_list, process_file))
  
  # Return the results
  return(results)
}



# Usage
directory_path <- "C:/temp/Laos/KBKZ/WWF_62%_LAOS/"  # Set your directory path
results <- process_files_in_directory(directory_path)

# View results
print(results)

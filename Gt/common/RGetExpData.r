# Function: GetExpData
# Reads in the experimental data from the input file
# Input:  fname = name of file that contains G(t) in 2 columns [t Gt]
# Output: A n*1 vector "t", and a n*1 vector Gt

GetExpData <- function(fname) {
  
  # Load data from file
  data <- read.table(fname, header = FALSE)
  to <- data[[1]]   # Extract first column as 'to'
  Gto <- data[[2]]  # Extract second column as 'Gto'
  
  # Remove duplicates
  to <- unique(to)
  Gto <- Gto[match(to, data[[1]])]  # Match 'Gto' to unique 'to' values
  
  # Sanitize the input by spacing it out. Using linear interpolation
  t <- 10^seq(log10(min(to)), log10(max(to)), length.out = length(to))
  Gt <- approx(to, Gto, xout = t, method = "linear", rule = 2)$y
  
  return(list(t = t, Gt = Gt))
}

# Example usage:
# Replace 'your_data_file.txt' with the actual file name containing your data
result <- GetExpData('your_data_file.txt')
t <- result$t
Gt <- result$Gt

# Function: GetExpData
#
# Reads in the experimental data from the input file
#
# Input:  fname = name of file that contains G(t) in 2 columns [t Gt]
#
# Output: A n*1 vector "t", and a n*1 vector Gt

GetExpData <- function(fname) {
  data <- read.table(fname)  # Read data from file
  to <- data[, 1]            # Extract time vector
  Gto <- data[, 2]           # Extract Gt vector
  
  # Remove duplicate time values and interpolate to uniform t vector
  to_unique <- unique(to)
  Gto_unique <- Gto[match(to_unique, to)]
  
  # Interpolate using logspace to ensure spacing is uniform in log scale
  t <- logspace(log10(min(to_unique)), log10(max(to_unique)), length(to_unique))
  Gt <- approx(to_unique, Gto_unique, xout = t, method = "linear", rule = 2)$y  # Linear interpolation
  
  return(list(t = t, Gt = Gt))
}



# Function: GetExpData
# Purpose: Reads in the experimental data from the input file
# Input: fname = name of file that contains G*(w) in 3 columns [w Gp Gpp]
#        Space it evenly on a log scale
# Output: A n*1 vector "w", and a 2n*1 vector Gexp = [Gp; Gpp]
GetExpData <- function(fname) {
  data <- read.table(fname)
  wo <- data[, 1]
  Gpo <- data[, 2]
  Gppo <- data[, 3]
  
  # Remove duplicate frequency values
  unique_indices <- unique(wo, incomparables = FALSE, fromLast = FALSE)
  wo <- wo[unique_indices]
  Gpo <- Gpo[unique_indices]
  Gppo <- Gppo[unique_indices]
  
  # Interpolate on log-spaced frequencies
  w <- exp(seq(log(min(wo)), log(max(wo)), length.out = length(wo)))
  Gp <- approx(wo, Gpo, xout = w)$y
  Gpp <- approx(wo, Gppo, xout = w)$y
  
  # Use super-smoother to clean it up further (optional)
  Gp <- supsmu(w, Gp)$y
  Gpp <- supsmu(w, Gpp)$y
  
  Gexp <- c(Gp, Gpp)  # Combine Gp and Gpp into a single vector
  
  return(list(w = w, Gexp = Gexp))
}

# Example usage:
# Ensure the `supsmu` function is properly defined and loaded in your R environment.
# Replace `fname` with the actual file path containing your experimental data.
# Call GetExpData function to retrieve `w` (frequency vector) and `Gexp` (experimental data).
# result <- GetExpData(fname)
# w <- result$w
# Gexp <- result$Gexp


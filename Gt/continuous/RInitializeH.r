# Function: InitializeH
# Input:  Gexp = n*1 vector [Gt],
#         t = n*1 vector contains times,
#         s = relaxation modes,
# Output: H = guessed H
InitializeH <- function(Gexp, t, s) {
  
  # Guess initial spectrum
  H <- -5.0 + sin(pi * s)
  
  # Initial lambda value for regularization
  lambda <- 1e0
  
  # First optimization with large lambda
  Hlam <- LevenMarq(lambda, Gexp, H, t, s)
  
  # Successively improve the initial guess until you have a reasonably good guess for low lambda
  lambda <- 1e-8
  H <- LevenMarq(lambda, Gexp, Hlam, t, s)
  
  return(H)
}

# Example usage:
# Replace 'your_data_file.txt' with the actual file name containing your data
# Replace 't' and 's' with appropriate vectors for your data
# Replace 'Gexp' with the actual observed data vector G(t)
t <- seq(0, 10, by = 1)  # Example time vector
s <- seq(1, 5, by = 1)   # Example relaxation modes
Gexp <- rep(1, length(t))  # Example observed data vector

# Call the InitializeH function
H <- InitializeH(Gexp, t, s)

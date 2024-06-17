# Function: kernel
# Outputs the n*1 dimensional vector K(H)(t) which is comparable to Gexp = Gt
# Modifying kernel for unevenly spaced s_i
# Input: H = substituted CRS,
#        t = n*1 vector contains times
#        s = relaxation modes
kernel <- function(H, t, s) {
  
  ns <- length(s)
  hs <- numeric(ns)
  
  # Compute hs using trapezoidal rule
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  # Create meshgrid of s and t
  S <- matrix(rep(s, each = length(t)), nrow = length(t), byrow = TRUE)
  T <- matrix(rep(t, each = length(s)), nrow = length(s), byrow = FALSE)
  
  # Compute kernel
  kern <- exp(-T / S)
  
  # Clear unnecessary variables
  rm(S, T)
  
  # Calculate K
  K <- kern %*% (hs * exp(H))
  
  return(K)
}

# Example usage:
# Replace 'your_data_file.txt' with the actual file name containing your data
# Replace 'H_value' with the actual substituted CRS value
t <- seq(0, 10, by = 1)  # Example time vector
s <- seq(1, 5, by = 1)   # Example relaxation modes
H <- 0.5  # Example substituted CRS value

# Call the kernel function
K <- kernel(H, t, s)


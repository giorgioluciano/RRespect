# Function to calculate tau and g values
taugcal <- function(p, Gp1, r, taumin, N) {
  # Initialize vectors
  taui <- numeric(N)
  H <- numeric(N)
  gi <- numeric(N)
  
  # Calculate tau values
  taui[1] <- taumin
  for (i in 2:N) {
    taui[i] <- taui[i-1] * r
  }
  
  # Calculate G1 and H values
  G1 <- 2 * Gp1 * gamma(1 + p) / p / pi * sin(p * pi / 2)
  H <- G1 * taui^(-p) / gamma(p)
  
  # Calculate g values
  gi <- H * (r^(p/2) - r^(-p/2)) / p
  
  return(list(taui = taui, gi = gi))
}

# Function to calculate Eabz
Eabzcal <- function(a, b, z, k) {
  return(z^k / gamma(a * k + b))
}

# Function to calculate R2
CalcR2 <- function() {
  # Set file paths
  FilePath <- file.path(getwd(), "WWF_62%_LAOS")
  FilePath2 <- getwd()
  
  # Define omega and amplitude vectors
  omegaV <- c(1, 10, 50, 100)
  AmpV <- c(0.05, 0.25, 0.5, 1.0, 2.0)
  
  # Initialize error variables
  error <- 0.0
  ErrorAroundMean <- 0.0
  R2 <- 0.0
  
  # Loop through amplitudes and frequencies
  for (j in rev(seq_along(AmpV))) {
    for (i in seq_along(omegaV)) {
      omega <- omegaV[i]
      Amp <- AmpV[j]
      AmpP <- Amp * 100
      
      # Construct file names
      FileName <- sprintf("Exp_WGF62%%_%d%%Om%d.txt", AmpP, omega)
      FileName2 <- sprintf("Theo_WGF62%%_%d%%Om%d.txt", AmpP, omega)
      
      # Load experimental and theoretical data
      q <- read.table(file.path(FilePath, FileName))
      qtheo <- read.table(file.path(FilePath2, FileName2))
      
      # Interpolate theoretical values
      xinterp <- seq_len(nrow(q))
      theoElStrvals <- approx(0:(nrow(qtheo)-1), qtheo[,5] / max(qtheo[,2]), 
                              xout = (xinterp-1) * (nrow(qtheo)-1) / nrow(q))$y
      
      # Calculate errors for elastic stress
      error <- error + sum((theoElStrvals - q[,5] / max(q[,2]))^2) / nrow(q)
      MeanEl <- mean(q[,5] / max(q[,2]))
      ErrorAroundMean <- ErrorAroundMean + sum((MeanEl - q[,5] / max(q[,2]))^2) / nrow(q)
      R2 <- R2 + 0.5 * (1 - error / ErrorAroundMean)
      
      # Interpolate theoretical values for viscous stress
      theoVisStrvals <- approx(0:(nrow(qtheo)-1), qtheo[,4] / max(qtheo[,2]), 
                               xout = (xinterp-1) * (nrow(qtheo)-1) / nrow(q))$y
      
      # Calculate errors for viscous stress
      error <- error + sum((theoVisStrvals - q[,4] / max(q[,2]))^2) / nrow(q)
      MeanVis <- mean(q[,4] / max(q[,2]))
      ErrorAroundMean <- ErrorAroundMean + sum((MeanVis - q[,4] / max(q[,2]))^2) / nrow(q)
      R2 <- R2 + 0.5 * (1 - error / ErrorAroundMean)
    }
  }
  
  # Calculate final error and R2 values
  error <- error / length(AmpV) / length(omegaV) / 2
  R2 <- R2 / length(AmpV) / length(omegaV)
  
  return(R2)
}
CalcR2 <- function() {
  
  # Define file paths (adjust as per your directory structure)
  FilePath <- paste0(getwd(), '/WWF_62%_LAOS/')
  FilePath2 <- getwd()
  
  # Define vectors of omega and Amp values
  omegaV <- c(1, 10, 50, 100)
  AmpV <- c(0.05, 0.25, 0.5, 1.0, 2.0)
  
  # Initialize variables for error calculations
  error <- 0.0
  ErrorAroundMean <- 0.0
  R2 <- 0.0
  
  # Nested loops over AmpV and omegaV
  for (j in length(AmpV):1) {
    for (i in 1:length(omegaV)) {
      
      # Extract current omega and Amp values
      omega <- omegaV[i]
      Amp <- AmpV[j]
      AmpP <- Amp * 100
      
      # Construct file names based on current AmpP and omega
      FileName <- paste0('Exp_WGF62%_', AmpP, '%Om', omega, '.txt')
      FileName2 <- paste0('Theo_WGF62%_', AmpP, '%Om', omega, '.txt')
      
      # Load experimental and theoretical data
      q <- read.table(paste0(FilePath, FileName))
      qtheo <- read.table(paste0(FilePath2, '/', FileName2))
      
      # Generate interpolation points
      xinterp <- seq(1, nrow(q), by = 1)
      
      # Interpolate theoretical elastic strain values
      theoElStrvals <- (approx(x = seq(0, nrow(qtheo) - 1), y = qtheo[, 5] / max(qtheo[, 2]), 
                               xout = (xinterp - 1) * (nrow(qtheo) - 1) / nrow(q))$y)
      
      # Calculate error for elastic strain
      error <- error + sum((theoElStrvals - q[, 5] / max(q[, 2]))^2) / nrow(q)
      
      # Calculate mean elastic strain
      MeanEl <- mean(q[, 5] / max(q[, 2]))
      
      # Update ErrorAroundMean with elastic strain error
      ErrorAroundMean <- ErrorAroundMean + sum((MeanEl - q[, 5] / max(q[, 2]))^2) / nrow(q)
      
      # Update R2 with elastic strain error
      R2 <- R2 + 0.5 * (1 - error / ErrorAroundMean)
      
      # Interpolate theoretical visible strain values
      theoVisStrvals <- (approx(x = seq(0, nrow(qtheo) - 1), y = qtheo[, 4] / max(qtheo[, 2]), 
                                xout = (xinterp - 1) * (nrow(qtheo) - 1) / nrow(q))$y)
      
      # Calculate error for visible strain
      error <- error + sum((theoVisStrvals - q[, 4] / max(q[, 2]))^2) / nrow(q)
      
      # Calculate mean visible strain
      MeanVis <- mean(q[, 4] / max(q[, 2]))
      
      # Update ErrorAroundMean with visible strain error
      ErrorAroundMean <- ErrorAroundMean + sum((MeanVis - q[, 4] / max(q[, 2]))^2) / nrow(q)
      
      # Update R2 with visible strain error
      R2 <- R2 + 0.5 * (1 - error / ErrorAroundMean)
    }
  }
  
  # Average errors and R2 over all combinations of Amp and omega
  error <- error / length(AmpV) / length(omegaV) / 2
  R2 <- R2 / length(AmpV) / length(omegaV)
  
  # Return the calculated R2 value
  return(R2)
}

calculate_error <- function(FilePath, omegaV, AmpV, pattern_exp, FilePath2 = NULL, pattern_theo = NULL, show_axes = TRUE) {
  
  error <- 0.0
  
  # Set up A4-sized plot window (210 x 297 mm in inches)
  windows(width = 8.27, height = 11.69)  # A4 size in inches
  
  # Set up the layout for plotting
  par(mfrow = c(length(AmpV), length(omegaV)))
  
  for (j in rev(seq_along(AmpV))) {
    for (i in seq_along(omegaV)) {
      
      omega <- omegaV[i]
      Amp <- AmpV[j]
      AmpP <- round(Amp * 100)  # Convert to integer percentage format
      
      # Generate the file name for the experimental data
      FileName <- paste0(FilePath, sprintf(pattern_exp, AmpP, omega))
      
      # Print which experimental file is being loaded
      cat('Attempting to load experimental file:', FileName, '\n')
      
      # Load experimental data with error handling
      if (file.exists(FileName)) {
        q <- as.matrix(read.table(FileName))
        cat('Successfully loaded experimental file:', FileName, '\n')
      } else {
        warning(paste('Experimental file not found:', FileName))
        next  # Skip to the next iteration if the file is not found
      }
      
      # Plot experimental data
      plot(q[, 3] / max(q[, 3]), q[, 2] / max(q[, 2]), type = 'l', lwd = 0.8, 
           main = paste('Amp:', AmpP, 'Omega:', omega), xlab = '', ylab = '')
      lines(q[, 3] / max(q[, 3]), q[, 4] / max(q[, 2]), col = 'red', lty = 2, lwd = 0.8)
      
      # Hide axes if the option is set to FALSE
      if (!show_axes) {
        axis(1, labels = FALSE, ticks = FALSE)
        axis(2, labels = FALSE, ticks = FALSE)
      }
      
      # If theoretical file path and pattern are provided, plot theoretical data and calculate error
      if (!is.null(FilePath2) && !is.null(pattern_theo)) {
        
        # Correct pattern for theoretical file names
        FileName2 <- paste0(FilePath2, sprintf(pattern_theo, AmpP, omega))
        
        # Print which theoretical file is being loaded
        cat('Attempting to load theoretical file:', FileName2, '\n')
        
        # Load theoretical data with error handling
        if (file.exists(FileName2)) {
          qtheo <- as.matrix(read.table(FileName2))
          cat('Successfully loaded theoretical file:', FileName2, '\n')
          
          # Plot theoretical data
          lines(qtheo[, 3] / max(qtheo[, 3]), qtheo[, 2] / max(qtheo[, 2]), col = 'green', lwd = 0.8)
          lines(qtheo[, 3] / max(qtheo[, 3]), qtheo[, 4] / max(qtheo[, 2]), col = 'green', lty = 2, lwd = 0.8)
          
          # Interpolation and error calculation
          xinterp <- seq_len(nrow(q))
          theoElStrvals <- approx(seq(0, nrow(qtheo) - 1), qtheo[, 5] / max(qtheo[, 2]), 
                                  x = (xinterp - 1) * (nrow(qtheo) - 1) / nrow(q))$y
          error <- error + sum((theoElStrvals - q[, 5] / max(q[, 2]))^2) / nrow(q)
          
          theoVisStrvals <- approx(seq(0, nrow(qtheo) - 1), qtheo[, 4] / max(qtheo[, 2]), 
                                   x = (xinterp - 1) * (nrow(qtheo) - 1) / nrow(q))$y
          error <- error + sum((theoVisStrvals - q[, 4] / max(q[, 2]))^2) / nrow(q)
          
        } else {
          warning(paste('Theoretical file not found:', FileName2))
        }
      }
      
    }
  }
  
  # Final error output if theoretical data is provided
  if (!is.null(FilePath2) && !is.null(pattern_theo)) {
    error <- error / length(AmpV) / length(omegaV) / 2
    cat('Error =', sqrt(error), '\n')
  } else {
    cat('Only experimental data plotted. No error calculated.\n')
  }
}

# Example usage with experimental data and axis control
FilePath <- 'C:/temp/WWF_30_8_FMM/WWF_62%_LAOS/'
omegaV <- c(1, 10, 50, 100)
AmpV <- c(0.05, 0.25, 0.5, 1.0, 2.0)
pattern_exp <- 'Exp_WGF62%%_%d%%Om%d.txt'  # Adjusted pattern

# Run the function with hidden axes
calculate_error(FilePath, omegaV, AmpV, pattern_exp, show_axes = FALSE)




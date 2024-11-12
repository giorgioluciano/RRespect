plot_lissajous_curves <- function(file_path, omega_v, amp_v) {
  # Function to plot Lissajous curves for various amplitudes and frequencies
  
  # Set up the plotting layout
  layout_matrix <- matrix(1:(2 * length(amp_v) * length(omega_v)), 
                          nrow = 2 * length(amp_v), ncol = length(omega_v), byrow = TRUE)
  layout(layout_matrix[nrow(layout_matrix):1,])
  
  # Loop through amplitude and omega values
  for (j in length(amp_v):1) {
    for (i in 1:length(omega_v)) {
      # Extract current omega and amplitude values
      omega <- omega_v[i]
      amp <- amp_v[j]
      amp_p <- amp * 100
      
      # Generate file name
      file_name <- sprintf("Exp_WGF62%%_%g%%Om%g.txt", amp_p, omega)
      full_path <- file.path(file_path, file_name)
      
      # Load data
      data <- read.table(full_path, header = FALSE)
      
      # Normalize data for viscous plot (columns 3, 2, 4)
      x_norm_viscous <- data[,3] / max(data[,3])
      y1_norm_viscous <- data[,2] / max(data[,2])
      y2_norm_viscous <- data[,4] / max(data[,2])  # Corrected normalization
      
      # Normalize data for elastic plot (columns 1, 2, 5)
      x_norm_elastic <- data[,1] / max(data[,1])
      y1_norm_elastic <- data[,2] / max(data[,2])
      y2_norm_elastic <- data[,5] / max(data[,2])  # Corrected normalization
      
      # Create viscous plot
      plot(x_norm_viscous, y1_norm_viscous, type = "l", lwd = 0.8,
           main = sprintf("Viscous ($\\gamma$ = %g%%, $\\omega$ = %g rad/s)", amp_p, omega),
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", frame.plot = FALSE)
      lines(x_norm_viscous, y2_norm_viscous, lty = 2, col = "red", lwd = 0.8)
      
      # Create elastic plot
      plot(x_norm_elastic, y1_norm_elastic, type = "l", lwd = 0.8,
           main = sprintf("Elastic ($\\gamma$ = %g%%, $\\omega$ = %g rad/s)", amp_p, omega),
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", frame.plot = FALSE)
      lines(x_norm_elastic, y2_norm_elastic, lty = 2, col = "red", lwd = 0.8)
    }
  }
}

# Example usage:
file_path <- "C:/temp/WWF_30_8_FMM/WWF_62%_LAOS/"
omega_v <- c(1, 10, 50, 100)
amp_v <- c(0.05, 0.25, 0.5, 1.0, 2.0)

# Call the function to create both viscous and elastic plots
plot_lissajous_curves(file_path, omega_v, amp_v)

# Set LaTeX interpretation option if needed
options(tikzLatexPackages = c(getOption("tikzLatexPackages"), "\\usepackage{amsmath}"))
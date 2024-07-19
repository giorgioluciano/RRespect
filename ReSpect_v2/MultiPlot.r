# Function to plot nonlinear harmonics
plot_nonlinear_harmonics <- function(amps, omegas, dim3, dim4) {
  # Get dimensions of amplitude and omega vectors
  dim1 <- length(amps)
  dim2 <- length(omegas)
  
  # Initialize data matrix
  Alldata <- matrix(0, nrow = dim3, ncol = dim4)
  
  # Loop over amplitudes and omegas to import data
  for (i in 1:dim1) {
    for (j in 1:dim2) {
      AmpP <- amps[i] * 100
      omega <- omegas[j]
      FileName <- sprintf('NonlinearHarmonicsParameters%d%%Om%d.txt', AmpP, omega)
      Table <- read.table(FileName, header = FALSE)
      Alldata[(j - 1) * dim1 + i, ] <- as.numeric(Table)
    }
  }
  
  # Plot function
  plot_data <- function(col_index, plot_title) {
    xdata <- numeric(dim1)
    ydata <- numeric(dim1)
    
    plot.new()
    plot.window(xlim = range(xdata), ylim = range(ydata))
    title(main = plot_title)
    box()
    
    for (j in seq(1, dim2, by = 2)) {
      for (i in 1:dim1) {
        xdata[i] <- Alldata[(j - 1) * dim1 + i, 2] * 100
        ydata[i] <- Alldata[(j - 1) * dim1 + i, col_index]
      }
      xi <- seq(min(xdata), max(xdata), length.out = 150)
      yi <- spline(xdata, ydata, xout = xi)$y
      lines(xi, yi, lwd = 1.5)
    }
  }
  
  # Create a 2x2 plot layout
  par(mfrow = c(2, 2))
  
  # Call plot function for different columns and titles
  plot_data(4, "e3")
  plot_data(8, "v3")
  plot_data(11, "S data")
  plot_data(12, "T data")
}

# Function to create multi-panel plot for viscous and elastic data
create_multi_panel_plot <- function() {
  # Set file paths
  file_path <- "WWF_62%_LAOS/"
  file_path2 <- getwd()
  
  # Define omega and amplitude vectors
  omega_v <- c(1, 10, 50, 100)
  amp_v <- c(0.05, 0.25, 0.5, 1.0, 2.0)
  
  # Initialize error
  error <- 0.0
  
  # Create plots for viscous data
  plot_data <- function(is_viscous) {
    par(mfrow = c(length(amp_v), length(omega_v)))
    
    for (j in length(amp_v):1) {
      for (i in 1:length(omega_v)) {
        omega <- omega_v[i]
        amp <- amp_v[j]
        amp_p <- amp * 100
        
        file_name <- sprintf("Exp_WGF62%%_%d%%Om%d.txt", amp_p, omega)
        file_name2 <- sprintf("Theo_WGF62%%_%d%%Om%d.txt", amp_p, omega)
        
        q <- read.table(paste0(file_path, file_name))
        q_theo <- read.table(paste0(file_path2, file_name2))
        
        plot.new()
        plot.window(xlim = c(-1, 1), ylim = c(-1, 1))
        
        if (is_viscous) {
          plot(q[,3]/max(q[,3]), q[,2]/max(q[,2]), type = "l", lwd = 0.8)
          lines(q[,3]/max(q[,3]), q[,4]/max(q[,2]), lty = 2, col = "red", lwd = 0.8)
          lines(q_theo[,3]/max(q_theo[,3]), q_theo[,2]/max(q_theo[,2]), col = "green", lwd = 0.8)
          lines(q_theo[,3]/max(q_theo[,3]), q_theo[,4]/max(q_theo[,2]), lty = 2, col = "green", lwd = 0.8)
        } else {
          plot(q[,1]/max(q[,1]), q[,2]/max(q[,2]), type = "l", lwd = 0.8)
          lines(q[,1]/max(q[,1]), q[,5]/max(q[,2]), lty = 2, col = "red", lwd = 0.8)
          lines(q_theo[,1]/max(q_theo[,1]), q_theo[,2]/max(q_theo[,2]), col = "green", lwd = 0.8)
          lines(q_theo[,1]/max(q_theo[,1]), q_theo[,5]/max(q_theo[,2]), lty = 2, col = "green", lwd = 0.8)
        }
        
        axis(1, labels = FALSE, tick = FALSE)
        axis(2, labels = FALSE, tick = FALSE)
        box()
        
        # Calculate error
        if (is_viscous) {
          x_interp <- seq_len(nrow(q))
          theo_el_str_vals <- approx(seq_len(nrow(q_theo)) - 1, q_theo[,5] / max(q_theo[,2]), 
                                     xout = (x_interp - 1) * (nrow(q_theo) - 1) / nrow(q))$y
          error <- error + sum((theo_el_str_vals - q[,5] / max(q[,2]))^2) / nrow(q)
          
          theo_vis_str_vals <- approx(seq_len(nrow(q_theo)) - 1, q_theo[,4] / max(q_theo[,2]), 
                                      xout = (x_interp - 1) * (nrow(q_theo) - 1) / nrow(q))$y
          error <- error + sum((theo_vis_str_vals - q[,4] / max(q[,2]))^2) / nrow(q)
        }
      }
    }
  }
  
  # Plot viscous data
  plot_data(TRUE)
  
  # Plot elastic data
  plot_data(FALSE)
  
  # Calculate and return total error
  error <- error / (length(amp_v) * length(omega_v) * 2)
  return(sqrt(error))
}

# Main execution
amps <- c(0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.6, 1.0, 1.5, 2.0)
omegas <- c(1, 3, 10, 30, 100)

# Call the nonlinear harmonics plotting function
plot_nonlinear_harmonics(amps, omegas, 55, 12)

# Call the multi-panel plot function and print the error
error <- create_multi_panel_plot()
cat("Error:", error, "\n")
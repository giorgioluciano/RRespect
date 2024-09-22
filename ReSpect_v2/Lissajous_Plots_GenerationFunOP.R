# Function to plot pre-calculated files from a directory
Plot_Precalculated_Files <- function(directory, plot_elastic = TRUE, plot_viscous = TRUE) {
  # Get the list of pre-calculated files
  files <- list.files(directory, pattern = "^Theo_WGF62_.*\\.txt$", full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No pre-calculated files found in the directory.\n")
    return(NULL)
  }
  
  for (file in files) {
    cat("Processing file: ", file, "\n")
    
    # Read the data file
    data <- read.table(file, header = FALSE, sep = "\t")
    strain <- data[, 1]
    tau12t2 <- data[, 2]
    strainRt <- data[, 3]
    VisStr <- data[, 4]
    ElStr <- data[, 5]
    
    # Extract the parameters from the filename
    file_info <- gsub("^Theo_WGF62_|\\.txt$", "", basename(file))
    info_parts <- unlist(strsplit(file_info, "%Om"))
    AmpP <- as.numeric(info_parts[1])
    omega <- as.numeric(info_parts[2])
    
    if (plot_elastic) {
      # Plot elastic stress
      plot_name_elastic <- paste0("ElasticSt", AmpP, "%Om", omega, ".png")
      png(file = plot_name_elastic)
      plot(strain / max(strain), tau12t2 / max(tau12t2), type = "l", col = "black", lwd = 1, xlab = "", ylab = "")
      lines(strain / max(strain), ElStr / max(tau12t2), col = "red", lty = 2, lwd = 0.8)
      dev.off()
      cat("Elastic plot saved: ", plot_name_elastic, "\n")
    }
    
    if (plot_viscous) {
      # Plot viscous stress
      plot_name_viscous <- paste0("ViscousSt", AmpP, "%Om", omega, ".png")
      png(file = plot_name_viscous)
      plot(strainRt / max(strainRt), tau12t2 / max(tau12t2), type = "l", col = "black", lwd = 1, xlab = "", ylab = "")
      lines(strainRt / max(strainRt), VisStr / max(tau12t2), col = "red", lty = 2, lwd = 0.8)
      dev.off()
      cat("Viscous plot saved: ", plot_name_viscous, "\n")
    }
  }
  
  cat("Plotting completed.\n")
}


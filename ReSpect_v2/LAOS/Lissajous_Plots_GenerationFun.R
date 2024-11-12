Lissajous_Plots_GenerationFun <- function(aSosval, bSosval, omegaV, AmpV, save_files_only = FALSE) {
  dim <- length(omegaV) * length(AmpV)
  
  for (ij in 1:dim) {
    i <- floor((ij - 1) / length(AmpV)) + 1
    j <- ij - (i - 1) * length(AmpV)
    
    omega <- omegaV[i]
    Amp <- AmpV[j]
    AmpP <- Amp * 100
    
    # Call your KBKZ model function here (assuming it returns time2, strain, etc.)
    results <- KBKZmodelSimpsonSoskeyLissajous(omega, Amp, aSosval, bSosval)
    time2 <- results$time2
    strain <- results$strain
    strainRt <- results$strainRt
    tau12t2 <- results$tau12t2
    ElStr <- results$ElStr
    VisStr <- results$VisStr
    
    tau12t2 <- ElStr + VisStr
    
    # Save data to a file
    data <- cbind(strain, tau12t2, strainRt, VisStr, ElStr)
    file_name <- paste0('Theo_WGF62_', AmpP, '%Om', omega, '.txt')
    write.table(data, file = file_name, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Plot if save_files_only is FALSE
    if (!save_files_only) {
      # Elastic stress plot
      plot_name_elastic <- paste0('ElasticSt', AmpP, '%Om', omega, '.png')
      png(file = plot_name_elastic)
      plot(strain / max(strain), tau12t2 / max(tau12t2), type = "l", col = "black", lwd = 1, xlab = "", ylab = "")
      lines(strain / max(strain), ElStr / max(tau12t2), col = "red", lty = 2, lwd = 0.8)
      dev.off()
      
      # Viscous stress plot
      plot_name_viscous <- paste0('ViscousSt', AmpP, '%Om', omega, '.png')
      png(file = plot_name_viscous)
      plot(strainRt / max(strainRt), tau12t2 / max(tau12t2), type = "l", col = "black", lwd = 1, xlab = "", ylab = "")
      lines(strainRt / max(strainRt), VisStr / max(tau12t2), col = "red", lty = 2, lwd = 0.8)
      dev.off()
    }
    
    cat("Processed: ", ij, " of ", dim, "\n")
  }
}

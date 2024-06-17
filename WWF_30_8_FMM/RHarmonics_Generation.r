# Define omegaV and AmpV
omegaV <- c(1, 3, 10, 30, 100)
AmpV <- c(0.001, 0.005, 0.01, 0.03)

# Get dimensions of omegaV and AmpV
dim1 <- length(omegaV)
dim2 <- length(AmpV)

# Loop through combinations of omegaV and AmpV
for (i in 1:dim1) {
  for (j in 1:dim2) {
    
    # Extract current omega and Amp values
    omega <- omegaV[i]
    Amp <- AmpV[j]
    AmpP <- Amp * 100
    
    # Call KBKZmodelSimpsonLaunHarmonics function
    result <- KBKZmodelSimpsonLaunHarmonics(c(omega, Amp, 0.6658, 10.35, 278.9))
    
    # Extract values from result
    harmonicsI <- result$harmonicsI
    harmonicsR <- result$harmonicsR
    
    e1 <- harmonicsI[1]
    e3 <- -harmonicsI[2]
    GpM <- e1 - 3 * e3
    GpL <- e1 + e3
    
    v1 <- harmonicsR[1] / omegaV[i]
    v3 <- harmonicsR[2] / omegaV[i]
    etapM <- v1 - 3 * v3
    etapL <- v1 + v3
    
    S <- (GpL - GpM) / GpL
    T <- (etapL - etapM) / etapL
    
    # Construct file name
    FileName <- paste0('NonlinearHarmonicsParameters', AmpP, '%Om', omega, '.txt')
    
    # Save results to text file
    write.table(data.frame(omega, Amp, e1, e3, GpM, GpL, v1, v3, etapM, etapL, S, T),
                file = FileName, row.names = FALSE, col.names = TRUE)
    
    # Print progress (optional)
    cat("Processed:", i, j, "\n")
  }
}








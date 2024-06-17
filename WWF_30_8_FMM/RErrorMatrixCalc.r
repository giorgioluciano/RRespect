# Define initial values and intervals
aSoso <- 49.1698
aSosf <- 1400.0
bSoso <- 1.00261
bSosf <- 1.90
intervals <- 5

# Generate sequences of aSosVals and bSosVals
aSosVals <- seq(aSoso, aSosf, length.out = intervals + 1)
bSosVals <- seq(bSoso, bSosf, length.out = intervals + 1)

# Initialize ErrorMat with zeros
ErrorMat <- matrix(0, nrow = length(aSosVals), ncol = length(bSosVals))

# Loop through parameter combinations
for (i in 1:length(aSosVals)) {
  for (j in 1:length(bSosVals)) {
    # Simulate Lissajous plots and calculate error
    # Assuming Lissajous_Plots_GenerationFun(a, b) and CreateMultiPanelPlotViscousFun() are defined functions
    Lissajous_Plots_GenerationFun(c(aSosVals[i], bSosVals[j]))
    ErrorMat[i, j] <- CreateMultiPanelPlotViscousFun()
  }
}

# Display or further process ErrorMat as needed
print(ErrorMat)

# SetParameters.R
# Traduzione fedele di readInput / par di pyReSpect-freq

setParameters <- function(GstFile = "input.dat") {
  par <- list(
    GstFile    = GstFile,
    ns         = 100L,       # numero punti spettro s
    FreqEnd    = 2L,         # 1=exp(-pi/2), 2=1/w, 3=exp(+pi/2)
    lamC       = 0,          # 0 = auto L-curve
    lam_min    = 1e-10,
    lam_max    = 1e3,
    lamDensity = 3,          # punti per decade in lambda
    SmFacLam   = 0,          # shift lambda: 0=nessuno, >0=verso max, <0=verso min
    plateau    = FALSE,      # TRUE per gel/Laponite (G0)
    verbose    = TRUE,
    plotting   = FALSE
  )
  return(par)
}

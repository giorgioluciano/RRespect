
# SetParameters_time.R

SetParameters <- function(GtFile     = "Gt.dat",
                          ns         = 100,
                          lamC       = 0,
                          SmFacLam   = 0,
                          FreqEnd    = 1,
                          verbose    = 1,
                          plotting   = 1,
                          lam_min    = 1e-10,
                          lam_max    = 1e3,
                          lamDensity = 3,
                          rho_cutoff = 0,
                          MaxNumModes         = 0,
                          deltaBaseWeightDist = 0.2,
                          minTauSpacing       = 1.25,
                          condWt     = 0.5,
                          BaseDistWt = 0.5) {
  par <- list(
    GtFile     = GtFile,
    ns         = ns,
    lamC       = lamC,
    SmFacLam   = SmFacLam,
    FreqEnd    = FreqEnd,
    verbose    = verbose,
    plotting   = plotting,
    lam_min    = lam_min,
    lam_max    = lam_max,
    lamDensity = lamDensity,
    rho_cutoff = rho_cutoff,
    MaxNumModes         = MaxNumModes,
    deltaBaseWeightDist = deltaBaseWeightDist,
    minTauSpacing       = minTauSpacing,
    condWt     = condWt,
    BaseDistWt = BaseDistWt
  )
  if (par$BaseDistWt < 0 || par$BaseDistWt > 1) stop("BaseDistWt must be in [0,1]")
  if (par$condWt     < 0 || par$condWt     > 1) stop("condWt must be in [0,1]")
  if (par$SmFacLam   < -1 || par$SmFacLam  > 1) stop("SmFacLam must be in [-1,1]")
  if (par$lamDensity < 2)                        stop("lamDensity must be >= 2")
  return(par)
}


setwd("C:/test/Gt")

library(pracma)
library(nnls)
library(minpack.lm)

source("SetParameters.R")
source("common.R")
source("LevenMarq.R")
source("LCurve.R")
source("contSpec.R")
source("GetWeights.R")
source("GridDensity.R")
source("MaxwellModes.R")
source("discSpec.R")

base_dir <- "C:/test/Gt"

# prima controlla cosa c'è in tests/
cat("File in tests/:\n")
print(list.files("tests/"))

tests <- list.files("tests/", pattern = "\\.dat$", full.names = TRUE)

if (!dir.exists("output_R_time")) dir.create("output_R_time")

results <- data.frame()

for (f in tests) {
  cat("\n==============================\n")
  cat("Running R time-domain on:", f, "\n")
  cat("==============================\n")
  
  test_name <- gsub("\\.dat", "", basename(f))
  outdir    <- file.path(base_dir, "output_R_time", test_name)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # copia il .dat in outdir con percorso assoluto
  file.copy(f, file.path(outdir, basename(f)), overwrite = TRUE)
  
  par <- SetParameters(
    GtFile     = basename(f),
    ns         = 100,
    lamC       = 0,
    SmFacLam   = 0,
    FreqEnd    = 1,
    verbose    = 1,
    plotting   = 1,
    lam_min    = 1e-10,
    lam_max    = 1e3,
    lamDensity = 3,
    deltaBaseWeightDist = 0.2,
    minTauSpacing       = 1.25,
    MaxNumModes         = 0
  )
  
  old_wd <- getwd()
  setwd(outdir)
  crs <- contSpec(par)
  drs <- discSpec(par, crs)
  setwd(old_wd)   # torna SEMPRE a C:/test/Gt
  
  cat(sprintf("[R] %-15s : lamC=%.3e  Nopt=%d  err=%.4f  cond=%.3e\n",
              basename(f), crs$lamC, drs$Nopt, drs$error, drs$condKp))
  
  results <- rbind(results, data.frame(
    test      = basename(f),
    lamC      = crs$lamC,
    Nopt      = drs$Nopt,
    error     = drs$error,
    log10cond = log10(drs$condKp)
  ))
}

cat("\n===== RISULTATI BATCH =====\n")
print(results)
write.csv(results, "output_R_time/results_summary.csv", row.names = FALSE)
cat("(*) Salvato: output_R_time/results_summary.csv\n")

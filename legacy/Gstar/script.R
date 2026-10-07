
# script.R  —  ReSpectR main runner
# Usage: setwd() to the folder containing all .R files and Gst.dat, then source this file.

source("SetParameters.R")
source("common.R")
source("LevenMarq.R")
source("LCurve.R")
source("contSpec.R")
source("GetWeights.R")
source("GridDensity.R")
source("MaxwellModes.R")
source("discSpec.R")

# ── Set parameters (edit GstFile to point to your data) ──────────────────────
par <- SetParameters(
  GstFile    = "Gst.dat",
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

# ── Continuous Relaxation Spectrum ────────────────────────────────────────────
crs <- contSpec(par)
cat(sprintf("\n[contSpec] lamC = %.3e\n", crs$lamC))

# ── Discrete Relaxation Spectrum ──────────────────────────────────────────────
drs <- discSpec(par, crs)
cat(sprintf("\n[discSpec] Nopt = %d,  error = %.4f,  cond = %.3e\n",
            drs$Nopt, drs$error, drs$condKp))

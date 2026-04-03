# run_all.R

base_dir <- "C:/test/Gstar"
setwd(base_dir)

source(file.path(base_dir, "SetParameters.R"), encoding = "UTF-8")
source(file.path(base_dir, "common.R"),        encoding = "UTF-8")
source(file.path(base_dir, "contSpec.R"),      encoding = "UTF-8")
source(file.path(base_dir, "discSpec.R"),      encoding = "UTF-8")

test_dir  <- file.path(base_dir, "tests")
out_root  <- file.path(base_dir, "output_R_freq")
dir.create(out_root, showWarnings = FALSE)

all_dat <- sort(list.files(test_dir, pattern = "\\.dat$", full.names = FALSE))
cat("File trovati in tests/:\n"); print(all_dat); cat("\n")

plateau_tests <- c("test7.dat", "test1u.dat")

base_par <- list(
  ns                   = 100L,
  lamC                 = 0,
  SmFacLam             = 0,
  FreqEnd              = 1L,
  verbose              = TRUE,
  plotting             = FALSE,
  lam_min              = 1e-10,
  lam_max              = 1e3,
  lamDensity           = 3,
  plateau              = FALSE,
  deltaBaseWeightDist  = 0.2,
  minTauSpacing        = 1.25,
  MaxNumModes          = 0
)

results <- list()

for (fname in all_dat) {
  cat(strrep("=", 40), "\n")
  cat("Running on:", fname, "\n")
  cat(strrep("=", 40), "\n")

  tst_name <- sub("\\.dat$", "", fname)
  out_dir  <- file.path(out_root, tst_name)
  dir.create(file.path(out_dir, "output"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(test_dir, fname), file.path(out_dir, fname), overwrite = TRUE)

  par <- base_par
  par$GstFile <- file.path(test_dir, fname)
  par$plateau <- fname %in% plateau_tests

  old_wd <- getwd()
  setwd(out_dir)

  # ── contSpec ──
  res_cont <- tryCatch(
    contSpec(par),
    error = function(e) { cat("  ERRORE contSpec:", conditionMessage(e), "\n"); NULL }
  )

  # ── discSpec ──
  res_disc <- NULL
  if (!is.null(res_cont)) {
    res_disc <- tryCatch(
      getDiscSpecMagic(par),
      error = function(e) { cat("  ERRORE discSpec:", conditionMessage(e), "\n"); NULL }
    )
  }

  setwd(old_wd)

  if (!is.null(res_disc)) {
    cat(sprintf("  lamC=%.3e | Nopt=%d | error=%.4f\n\n",
                res_cont$lamC, res_disc$Nopt, res_disc$error))
  }

  results[[fname]] <- list(cont = res_cont, disc = res_disc)
}

cat("Batch completato.\n")

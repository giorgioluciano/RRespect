find_extdata_path_test <- function() {
  pkg_path <- system.file("extdata", package = "ReSpectR")
  if (nzchar(pkg_path) && dir.exists(pkg_path)) {
    return(pkg_path)
  }

  local_path <- normalizePath(file.path("..", "..", "inst", "extdata"), mustWork = FALSE)
  if (dir.exists(local_path)) {
    return(local_path)
  }

  stop("Could not locate extdata for tests.")
}

test_that("spectrum results have stable classes", {
  extdata <- find_extdata_path_test()
  time_file <- file.path(extdata, "time_tests", "test1.dat")

  par <- setParams(domain = "time", dataFile = time_file, verbose = FALSE)
  crs <- getContinuousSpectrum(par)
  drs <- getDiscreteSpectrum(par, crs = crs)

  expect_s3_class(crs, "respect_continuous_spectrum")
  expect_s3_class(drs, "respect_discrete_spectrum")
  expect_equal(crs$domain, "time")
  expect_equal(drs$domain, "time")
  expect_true(is.matrix(crs$Gfit))
  expect_true(is.list(drs$continuous))
})

test_that("plot methods run for time and frequency results", {
  extdata <- find_extdata_path_test()
  time_file <- file.path(extdata, "time_tests", "test1.dat")
  freq_file <- file.path(extdata, "freq_tests", "test1.dat")

  time_par <- setParams(domain = "time", dataFile = time_file, verbose = FALSE)
  freq_par <- setParams(domain = "frequency", dataFile = freq_file, verbose = FALSE)

  time_crs <- getContinuousSpectrum(time_par)
  time_drs <- getDiscreteSpectrum(time_par, crs = time_crs)
  freq_crs <- getContinuousSpectrum(freq_par)
  freq_drs <- getDiscreteSpectrum(freq_par, crs = freq_crs)

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)

  expect_invisible(plot(time_crs))
  expect_invisible(plot(time_drs))
  expect_invisible(plot(freq_crs))
  expect_invisible(plot(freq_drs))
})
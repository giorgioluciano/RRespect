test_that("getDiscreteSpectrum no longer requires forced write output", {
  root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
  par <- setParams(domain = "frequency", dataFile = "missing-input.dat")
  err <- tryCatch(
    {
      getDiscreteSpectrum(par, writeOutput = FALSE, projectDir = root)
      NULL
    },
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_false(grepl("requires writeOutput", conditionMessage(err), fixed = TRUE))
})

test_that("compute wrappers validate malformed par", {
  expect_error(getContinuousSpectrum(list()), "setParams")
  expect_error(getDiscreteSpectrum(list()), "setParams")
})

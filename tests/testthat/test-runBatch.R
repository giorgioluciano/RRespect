test_that("runBatch validates input source arguments", {
  expect_error(runBatch(domain = "time", dataFiles = NULL, dataDir = NULL), "Provide either")
})

test_that("runBatch returns row with error status when file fails", {
  fake <- tempfile(fileext = ".dat")
  root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
  res <- runBatch(
    domain = "frequency",
    dataFiles = fake,
    writeOutput = FALSE,
    progress = FALSE,
    projectDir = root
  )
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 1)
  expect_equal(res$status[[1]], "error")
})

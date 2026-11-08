test_that("API aliases are callable", {
  expect_equal(getParams(domain = "time")$domain, "time")
  expect_error(contSpectrum(list()), "setParams")
  expect_error(discSpectrum(list()), "setParams")
})

test_that("writeResults forwards correctly", {
  tmp <- tempfile("respect-out-")
  res <- writeResults(list(s = c(1, 2), H = c(0.1, 0.2)), outputDir = tmp)
  expect_true(is.character(res))
  expect_true(any(grepl("H\\.dat$", res)))
})

test_that("domain batch aliases set expected domain", {
  root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
  fake <- tempfile(fileext = ".dat")

  a <- runTimeBatch(dataFiles = fake, writeOutput = FALSE, progress = FALSE, projectDir = root)
  b <- runFrequencyBatch(dataFiles = fake, writeOutput = FALSE, progress = FALSE, projectDir = root)

  expect_equal(nrow(a), 1)
  expect_equal(nrow(b), 1)
  expect_equal(a$status[[1]], "error")
  expect_equal(b$status[[1]], "error")
})

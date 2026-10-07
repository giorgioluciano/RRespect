test_that("setParams returns normalized lists", {
  p_time <- setParams(domain = "time")
  expect_equal(p_time$domain, "time")
  expect_equal(p_time$dataFile, "Gt.dat")
  expect_true(is.numeric(p_time$deltaBaseWeightDist))

  p_freq <- setParams(domain = "frequency")
  expect_equal(p_freq$domain, "frequency")
  expect_equal(p_freq$dataFile, "Gst.dat")
  expect_true(is.logical(p_freq$plateau))
})

test_that("setParams validates constrained values", {
  expect_error(setParams(lamDensity = 1), "lamDensity")
  expect_error(setParams(smFacLam = 2), "smFacLam")
  expect_error(setParams(condWt = 2), "condWt")
})

test_that("compatibility wrappers map to setParams", {
  p1 <- setReSpectParameters(domain = "time")
  p2 <- getReSpectParameters(domain = "time")
  p3 <- setRespectParameters(domain = "time")
  p4 <- getParams(domain = "time")
  expect_equal(p1$domain, "time")
  expect_equal(p2$domain, "time")
  expect_equal(p3$domain, "time")
  expect_equal(p4$domain, "time")
})

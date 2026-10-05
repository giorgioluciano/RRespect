test_that("frequency-domain discrete spectrum reconstructs a Maxwell system", {
  omega <- 10^seq(-4, 4, length.out = 161)
  
  truth_g <- c(1500, 400, 80)
  truth_tau <- c(0.005, 0.1, 2)
  
  dat <- make_maxwell_data(
    omega = omega,
    G = truth_g,
    tau = truth_tau
  )
  
  input_file <- tempfile(
    pattern = "maxwell-frequency-",
    fileext = ".dat"
  )
  
  utils::write.table(
    dat,
    file = input_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  par <- setParams(
    domain = "frequency",
    dataFile = input_file,
    ns = 100L,
    lamC = 0,
    smFacLam = 0,
    freqEnd = 1L,
    verbose = FALSE,
    plotting = FALSE,
    lamMin = 1e-10,
    lamMax = 1e3,
    lamDensity = 3L,
    plateau = FALSE
  )
  
  fit <- getDiscreteSpectrum(
    par = par,
    writeOutput = FALSE,
    projectDir = tempdir()
  )
  
  expect_s3_class(fit, "respect_discrete_spectrum")
  expect_equal(fit$domain, "frequency")
  expect_equal(fit$result_type, "discrete")
  
  expect_true(all(is.finite(fit$g)))
  expect_true(all(is.finite(fit$tau)))
  expect_true(all(fit$g > 0))
  expect_true(all(fit$tau > 0))
  
  wt <- outer(fit$continuous$w, fit$tau, "*")
  
  fitted_gp <- rowSums(
    sweep(
      wt^2 / (1 + wt^2),
      2,
      fit$g,
      "*"
    )
  )
  
  fitted_gpp <- rowSums(
    sweep(
      wt / (1 + wt^2),
      2,
      fit$g,
      "*"
    )
  )
  
  n <- length(fit$continuous$w)
  observed_gp <- fit$continuous$Gexp[seq_len(n)]
  observed_gpp <- fit$continuous$Gexp[n + seq_len(n)]
  
  rel_rmse <- function(observed, fitted) {
    sqrt(mean((observed - fitted)^2)) / max(observed)
  }
  
  expect_lt(rel_rmse(observed_gp, fitted_gp), 1e-3)
  expect_lt(rel_rmse(observed_gpp, fitted_gpp), 1e-3)
  
  expect_equal(
    sum(fit$g),
    sum(truth_g),
    tolerance = 0.1
  )
  for (i in seq_along(truth_tau)) {
    j <- which.min(abs(log(fit$tau) - log(truth_tau[i])))
    
    expect_equal(
      log(fit$tau[j]),
      log(truth_tau[i]),
      tolerance = 0.08
    )
    
    expect_equal(
      fit$g[j],
      truth_g[i],
      tolerance = 0.05 * truth_g[i]
    )
  }
})
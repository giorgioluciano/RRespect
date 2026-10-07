test_that("frequency kernel agrees with independent trapezoidal integration", {
  s <- 10^seq(-6, 6, length.out = 101)
  w <- 10^seq(-4, 4, length.out = 81)
  x <- log(s)

  density <- 1000 * exp(
    -0.5 * ((x - log(0.1)) / 0.8)^2
  )

  H <- log(density)

  trapezoid <- function(y) {
    sum(
      diff(x) *
        (head(y, -1L) + tail(y, -1L)) / 2
    )
  }

  reference_gp <- vapply(
    w,
    function(omega) {
      z <- omega * s
      trapezoid(density * z^2 / (1 + z^2))
    },
    numeric(1)
  )

  reference_gpp <- vapply(
    w,
    function(omega) {
      z <- omega * s
      trapezoid(density * z / (1 + z^2))
    },
    numeric(1)
  )

  kernel <- ReSpectR:::.gstarGetKernMat(s, w)

  predicted <- ReSpectR:::.gstarKernelPrestore(
    H,
    kernel
  )

  n <- length(w)

  expect_true(all(is.finite(predicted)))

  expect_lt(
    max(abs(
      predicted[seq_len(n)] / reference_gp - 1
    )),
    1e-10
  )

  expect_lt(
    max(abs(
      predicted[n + seq_len(n)] / reference_gpp - 1
    )),
    1e-10
  )

  plateau <- 50

  predicted_plateau <- ReSpectR:::.gstarKernelPrestore(
    H,
    kernel,
    G0 = plateau
  )

  reference_plateau <- c(
    reference_gp + plateau,
    reference_gpp
  )

  expect_lt(
    max(abs(
      predicted_plateau / reference_plateau - 1
    )),
    1e-10
  )
})

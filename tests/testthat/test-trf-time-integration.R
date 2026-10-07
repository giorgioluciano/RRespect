
test_that("experimental time API uses TRF for initialization and final fit", {
  local_mocked_bindings(
    .gtLevenMarq = function(...) {
      stop("Experimental workflow called legacy.")
    }
  )

  local({
    file <- tempfile(fileext = ".dat")
    on.exit(unlink(file), add = TRUE)

    t <- exp(seq(log(0.01), log(10), length.out = 30))

    for (plateau in c(FALSE, TRUE)) {
      G <- exp(-t / 0.5) + 0.4 * exp(-t / 10) +
        if (plateau) 0.2 else 0

      write.table(
        cbind(t, G), file,
        row.names = FALSE, col.names = FALSE
      )

      par <- setParams(
        domain = "time",
        dataFile = file,
        ns = 20L,
        lamC = 0.01,
        plateau = plateau,
        solver = "experimental",
        verbose = FALSE,
        plotting = FALSE
      )

      out <- getContinuousSpectrum(par, writeOutput = FALSE)

      expect_s3_class(out, "respect_continuous_spectrum")
      expect_identical(out$solver, "experimental")
      expect_true(out$diagnostics$initial$success)
      expect_true(out$diagnostics$final$success)
      expect_null(out$diagnostics$scan)
      expect_true(all(is.finite(out$H)))
      expect_true(all(is.finite(out$Gfit)))

      Hraw <- -5 * rep(1, length(out$s)) + sin(pi * out$s)
      initial <- .gtGetHExperimental(
        1, out$Gexp, out$wexp, Hraw, out$kernMat,
        G0 = if (plateau) min(out$Gexp) else NULL
      )
      final <- .gtGetHExperimental(
        out$lamC, out$Gexp, out$wexp,
        initial$H, out$kernMat, initial$G0
      )

      expect_equal(out$H, final$H, tolerance = 1e-10)
      expect_equal(
        out$G0,
        if (plateau) final$G0 else 0,
        tolerance = 1e-10
      )
      expect_identical(out$lam_C, out$lamC)
      expect_equal(out$G_fit, as.numeric(out$Gfit[, 2]))
      expect_null(out$dH)
      expect_null(out$H_lam)
      expect_equal(
        out$diagnostics$final$cost,
        final$diagnostics$cost,
        tolerance = 1e-10
      )

      r <- .gtResidualLM(
        c(out$H, out$G0),
        out$lamC, out$Gexp, out$wexp, out$kernMat
      )
      expect_equal(
        0.5 * sum(r^2),
        out$diagnostics$final$cost,
        tolerance = 1e-10
      )

      expect_equal(
        as.numeric(out$Gfit[, 2]),
        final$Gfit,
        tolerance = 1e-10
      )
    }
  })
})

test_that("experimental automatic lambda uses TRF throughout the scan", {
  local_mocked_bindings(
    .gtLevenMarq = function(...) {
      stop("Experimental lambda scan called legacy.")
    }
  )

  local({
    file <- tempfile(fileext = ".dat")
    on.exit(unlink(file), add = TRUE)

    t <- exp(seq(log(0.01), log(10), length.out = 30))

    for (plateau in c(FALSE, TRUE)) {
      G <- exp(-t / 0.5) + 0.4 * exp(-t / 10) +
        if (plateau) 0.2 else 0

      write.table(
        cbind(t, G), file,
        row.names = FALSE, col.names = FALSE
      )

      par <- setParams(
        domain = "time",
        dataFile = file,
        ns = 20L,
        lamC = 0,
        lamMin = 1e-4,
        lamMax = 10,
        lamDensity = 2L,
        plateau = plateau,
        solver = "experimental",
        verbose = FALSE,
        plotting = FALSE
      )

      out <- getContinuousSpectrum(par, writeOutput = FALSE)
      scan <- out$diagnostics$scan

      expect_identical(out$solver, "experimental")
      expect_true(is.finite(out$lamC))
      expect_gte(out$lamC, par$lamMin * (1 - 1e-12))
      expect_lte(out$lamC, par$lamMax * (1 + 1e-12))
      expect_true(length(scan) > 0)
      expect_equal(length(out$lam), length(scan))
      expect_equal(dim(out$H_lam), c(length(out$s), length(out$lam)))
      expect_equal(length(out$dH), length(out$s))
      expect_true(all(is.finite(out$dH)))
      expect_true(all(out$dH >= 0))
      expect_true(all(is.finite(out$log_P)))
      expect_equal(
        out$dH,
        .respectErrorBand(
          out$H_lam,
          exp(out$log_P) / sum(exp(out$log_P))
        ),
        tolerance = 1e-12
      )
      expect_true(all(vapply(
        scan, function(d) isTRUE(d$success), logical(1)
      )))
      expect_true(all(vapply(
        scan,
        function(d) all(is.finite(d$x)) && is.finite(d$cost),
        logical(1)
      )))
      expect_true(all(vapply(
        scan,
        function(d) all(diff(d$cost_history) <= 0),
        logical(1)
      )))
      expect_true(out$diagnostics$final$success)
      expect_true(all(is.finite(out$Gfit)))
    }
  })
})

test_that("experimental dispatcher rejects an unsuccessful TRF fit", {
  local_mocked_bindings(
    .gtGetHExperimental = function(...) {
      list(
        diagnostics = list(
          success = FALSE,
          status = 0L,
          nfev = 1L
        )
      )
    }
  )

  expect_error(
    .gtFitH(
      lam = 0.01,
      Gexp = 1,
      wexp = 1,
      H = rep(-5, 3),
      kernMat = matrix(1, 1, 3),
      solver = "experimental"
    ),
    "did not converge"
  )
})

test_that("experimental API validation accepts both domains", {
  for (domain in c("time", "frequency")) {
    par <- setParams(
      domain = domain,
      solver = "experimental",
      verbose = FALSE
    )
    expect_silent(.validatePar(par, operation = "continuous"))
    expect_silent(.validatePar(par, operation = "discrete"))
  }
})

test_that("updated TRF workflow is the public default", {
  for (domain in c("time", "frequency")) {
    par <- setParams(domain = domain, verbose = FALSE)
    expect_identical(par$solver, "experimental")
    expect_identical(par$lamDensity, 2L)
    expect_silent(.validatePar(par, operation = "continuous"))
    expect_silent(.validatePar(par, operation = "discrete"))
  }
})

test_that("explicit legacy selection uses its dispatcher", {
  par <- setParams(domain = "time", solver = "legacy")
  expect_identical(par$solver, "legacy")
  expect_identical(par$lamDensity, 3L)

  local_mocked_bindings(
    .gtLevenMarq = function(lam, Gexp, wexp, H, kernMat, G0 = NULL) {
      if (is.null(G0)) H else list(H = H, G0 = G0)
    },
    .gtGetHExperimental = function(...) {
      stop("Legacy workflow called experimental.")
    }
  )

  for (G0 in list(NULL, 0.2)) {
    out <- .gtFitH(
      lam = 0.01,
      Gexp = 1,
      wexp = 1,
      H = rep(-5, 3),
      kernMat = matrix(1, 1, 3),
      G0 = G0,
      solver = par$solver
    )

    expect_identical(out$solver, "legacy")
    expect_equal(out$H, rep(-5, 3))
    expect_equal(out$G0, G0)
    expect_null(out$diagnostics)
  }
})


test_that("experimental frequency continuous API bypasses legacy LM", {
  local_mocked_bindings(
    .gstarGetH = function(...) {
      stop("Experimental frequency workflow called legacy LM.")
    }
  )

  local({
    file <- tempfile(fileext = ".dat")
    on.exit(unlink(file), add = TRUE)

    w <- exp(seq(log(0.01), log(100), length.out = 30))
    tau <- c(0.5, 10)
    g <- c(1, 0.4)
    wt <- outer(w, tau)
    base_Gp <- as.vector((wt^2 / (1 + wt^2)) %*% g)
    Gpp <- as.vector((wt / (1 + wt^2)) %*% g)

    for (plateau in c(FALSE, TRUE)) {
      Gp <- base_Gp + if (plateau) 0.2 else 0

      write.table(
        cbind(w, Gp, Gpp, rep(1, length(w)), rep(1, length(w))),
        file,
        row.names = FALSE,
        col.names = FALSE
      )

      for (lambda in c(0.01, 0)) {
        par <- setParams(
          domain = "frequency",
          dataFile = file,
          ns = 20L,
          lamC = lambda,
          lamMin = 1e-4,
          lamMax = 10,
          lamDensity = 2L,
          plateau = plateau,
          solver = "experimental",
          verbose = FALSE,
          plotting = FALSE
        )

        out <- getContinuousSpectrum(par, writeOutput = FALSE)

        expect_identical(out$solver, "experimental")
        expect_identical(out$domain, "frequency")
        expect_true(out$diagnostics$initial$success)
        expect_true(out$diagnostics$final$success)
        expect_true(all(is.finite(out$H)))
        expect_true(all(is.finite(out$G_fit)))
        expect_equal(length(out$G_fit), 2 * length(out$w))
        expect_equal(
          out$G_fit,
          c(out$Gfit[, 2], out$Gfit[, 3]),
          tolerance = 1e-12
        )
        expect_identical(out$lam_C, out$lamC)

        if (!plateau) expect_equal(out$G0, 0)

        x <- if (plateau) c(out$H, out$G0) else out$H
        r <- .gstarResidualLM(
          x, out$lamC, out$Gexp, out$wexp, out$kernMat
        )

        expect_equal(
          0.5 * sum(r^2),
          out$diagnostics$final$cost,
          tolerance = 1e-10
        )

        if (lambda == 0) {
          expect_true(is.finite(out$lamC))
          expect_true(length(out$lam) > 0)
          expect_equal(
            dim(out$H_lam),
            c(length(out$s), length(out$lam))
          )
          expect_true(all(is.finite(out$dH)))
          expect_true(all(is.finite(out$log_P)))
          expect_true(all(vapply(
            out$diagnostics$scan,
            function(d) isTRUE(d$success),
            logical(1)
          )))
        } else {
          expect_null(out$lam)
          expect_null(out$dH)
          expect_null(out$H_lam)
          expect_null(out$diagnostics$scan)
        }
      }
    }
  })
})

test_that("frequency plateau Jacobian affects storage rows only", {
  n <- 5L
  ns <- 4L
  K <- matrix(seq_len(2 * n * ns) / 100, 2 * n, ns)
  G <- rep(2, 2 * n)
  weights <- seq(1, 2, length.out = 2 * n)
  x <- c(rep(-2, ns), 0.3)

  J <- .gstarJacobianLM(x, 0.01, G, weights, K)

  expect_equal(J[seq_len(n), ns + 1L], -weights[seq_len(n)] / G[seq_len(n)])
  expect_equal(J[n + seq_len(n), ns + 1L], rep(0, n))
  expect_equal(J[2 * n + seq_len(ns - 2L), ns + 1L], rep(0, ns - 2L))
})

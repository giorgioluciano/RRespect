
test_that("experimental discrete API uses its own pipeline in both domains", {
  local_mocked_bindings(
    .gtDiscSpec = function(...) stop("Called legacy time discrete."),
    .gstarGetDiscSpecMagic = function(...) stop("Called legacy frequency discrete."),
    .gtFineTuneSolution = function(...) stop("Called legacy time refinement."),
    .gstarFineTuneSolution = function(...) stop("Called legacy frequency refinement.")
  )

  local({
    file <- tempfile(fileext = ".dat")
    on.exit(unlink(file), add = TRUE)

    for (domain in c("time", "frequency")) {
      axis <- exp(seq(log(0.01), log(10), length.out = 30))
      K <- .experimentalMaxwellKernel(c(0.5, 3), axis, domain)
      base <- as.vector(K %*% c(1, 0.4))

      for (plateau in c(FALSE, TRUE)) {
        G <- base

        if (plateau) {
          if (domain == "time") {
            G <- G + 0.2
          } else {
            G[seq_along(axis)] <- G[seq_along(axis)] + 0.2
          }
        }

        input <- if (domain == "time") {
          cbind(axis, G, rep(1, length(axis)))
        } else {
          n <- length(axis)
          cbind(
            axis, G[seq_len(n)], G[n + seq_len(n)],
            rep(1, n), rep(1, n)
          )
        }

        write.table(
          input, file,
          row.names = FALSE, col.names = FALSE
        )

        par <- setParams(
          domain = domain,
          dataFile = file,
          ns = 20L,
          lamC = 0.01,
          plateau = plateau,
          maxNumModes = 3L,
          solver = "experimental",
          verbose = FALSE,
          plotting = FALSE
        )

        crs <- getContinuousSpectrum(par)
        expect_silent(out <- getDiscreteSpectrum(par, crs = crs))

        expect_s3_class(out, "respect_discrete_spectrum")
        expect_identical(out$solver, "experimental")
        expect_equal(out$N, length(out$g))
        expect_equal(length(out$g), length(out$tau))
        expect_equal(length(out$dtau), length(out$tau))
        expect_true(all(is.finite(out$g)))
        expect_true(all(out$g >= 0))
        expect_true(all(is.finite(out$tau)))
        expect_true(all(diff(out$tau) >= 0))
        expect_true(all(is.finite(out$G_fit)))

        bounds <- .experimentalTauBounds(axis, domain)
        expect_true(all(out$tau >= bounds[1]))
        expect_true(all(out$tau <= bounds[2]))

        expect_equal(
          out$error,
          sum((crs$wexp * (out$G_fit / crs$Gexp - 1))^2),
          tolerance = 1e-10
        )

        expect_equal(length(out$wt_base), length(out$AIC_bst))
        expect_equal(length(out$N_bst), length(out$nz_N_bst))
        expect_true(all(out$nz_N_bst <= out$N_bst))
        if (!plateau) expect_equal(out$G0, 0)
      }
    }
  })
})

#' Plot ReSpectR Spectrum Results
#'
#' Native plotting methods for ReSpectR result objects. Continuous-spectrum
#' results show the recovered spectrum together with the fitted data curve.
#' Discrete-spectrum results show Maxwell modes overlaid on the continuous
#' envelope, plus the fitted response implied by the discrete modes.
#'
#' @param x A result object with class `"respect_continuous_spectrum"` or
#'   `"respect_discrete_spectrum"`.
#' @param y Unused.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @name plot_respect_spectrum
#' @examples
#' \dontrun{
#' time_file <- system.file("extdata", "time_tests", "test1.dat", package = "ReSpectR")
#' par <- setParams(domain = "time", dataFile = time_file, verbose = FALSE)
#' crs <- getContinuousSpectrum(par)
#' plot(crs)
#' drs <- getDiscreteSpectrum(par, crs = crs)
#' plot(drs)
#' }
NULL

.respectSpectrumAmplitude <- function(x) {
  exp(x$H)
}

.respectDiscreteFitTime <- function(x) {
  crs <- x$continuous
  fit <- as.vector(exp(-outer(crs$t, 1 / x$tau)) %*% x$g)
  cbind(crs$t, fit)
}

.respectDiscreteFitFrequency <- function(x) {
  crs <- x$continuous
  n   <- length(crs$w)
  ws  <- outer(crs$w, x$tau)
  ws2 <- ws^2
  gp  <- as.vector((ws2 / (1 + ws2)) %*% x$g)
  gpp <- as.vector((ws  / (1 + ws2)) %*% x$g)
  if (!is.null(crs$G0) && isTRUE(crs$G0 != 0)) {
    gp <- gp + crs$G0
  }
  cbind(crs$w, gp, gpp)
}

.respectModeOverlayPlot <- function(crs, drs, main) {
  amp  <- .respectSpectrumAmplitude(crs)
  ymin <- min(c(amp[amp > 0], drs$g[drs$g > 0]), na.rm = TRUE)

  graphics::plot(crs$s, amp,
                 type = "l", log = "xy", lwd = 2, col = "#1b6ca8",
                 xlab = expression(tau), ylab = "spectrum amplitude",
                 main = main)
  graphics::segments(drs$tau, rep(ymin, length(drs$tau)), drs$tau, drs$g,
                     col = "#c44e52", lwd = 2)
  graphics::points(drs$tau, drs$g, pch = 19, col = "#c44e52")
  graphics::legend("topright",
                   legend = c("continuous envelope", "discrete modes"),
                   col = c("#1b6ca8", "#c44e52"), lty = c(1, 1), pch = c(NA, 19),
                   bty = "n")
}

#' @rdname plot_respect_spectrum
#' @export
plot.respect_continuous_spectrum <- function(x, y = NULL, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = c(1, 2))

  amp <- .respectSpectrumAmplitude(x)
  graphics::plot(x$s, amp,
                 type = "l", log = "xy", lwd = 2, col = "#1b6ca8",
                 xlab = expression(tau), ylab = "spectrum amplitude",
                 main = sprintf("Continuous Spectrum (%s)", x$domain), ...)

  if (identical(x$domain, "time")) {
    graphics::plot(x$t, x$Gexp,
                   log = "xy", pch = 16, cex = 0.6, col = "#4c4c4c",
                   xlab = "t", ylab = "G(t)",
                   main = "Observed vs fitted")
    graphics::lines(x$Gfit[, 1], x$Gfit[, 2], lwd = 2, col = "#1b9e77")
    graphics::legend("topright", legend = c("observed", "fit"),
                     col = c("#4c4c4c", "#1b9e77"), pch = c(16, NA), lty = c(NA, 1), bty = "n")
  } else {
    n <- length(x$w)
    graphics::plot(x$w, x$Gexp[seq_len(n)],
                   log = "xy", pch = 16, cex = 0.6, col = "#4c4c4c",
                   xlab = expression(omega), ylab = expression(G^"*"(omega)),
                   main = "Observed vs fitted")
    graphics::points(x$w, x$Gexp[n + seq_len(n)], pch = 1, cex = 0.6, col = "#7f7f7f")
    graphics::lines(x$Gfit[, 1], x$Gfit[, 2], lwd = 2, col = "#1b9e77")
    graphics::lines(x$Gfit[, 1], x$Gfit[, 3], lwd = 2, col = "#d95f02")
    graphics::legend("topright",
                     legend = c("G' observed", "G'' observed", "G' fit", "G'' fit"),
                     col = c("#4c4c4c", "#7f7f7f", "#1b9e77", "#d95f02"),
                     pch = c(16, 1, NA, NA), lty = c(NA, NA, 1, 1), bty = "n")
  }

  invisible(x)
}

#' @rdname plot_respect_spectrum
#' @export
plot.respect_discrete_spectrum <- function(x, y = NULL, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = c(1, 2))

  crs <- x$continuous
  .respectModeOverlayPlot(crs, x, sprintf("Discrete Modes (%s)", x$domain))

  if (identical(x$domain, "time")) {
    fit <- .respectDiscreteFitTime(x)
    graphics::plot(crs$t, crs$Gexp,
                   log = "xy", pch = 16, cex = 0.6, col = "#4c4c4c",
                   xlab = "t", ylab = "G(t)",
                   main = "Discrete fit vs observed", ...)
    graphics::lines(fit[, 1], fit[, 2], lwd = 2, col = "#c44e52")
    graphics::legend("topright", legend = c("observed", "discrete fit"),
                     col = c("#4c4c4c", "#c44e52"), pch = c(16, NA), lty = c(NA, 1), bty = "n")
  } else {
    fit <- .respectDiscreteFitFrequency(x)
    n <- length(crs$w)
    graphics::plot(crs$w, crs$Gexp[seq_len(n)],
                   log = "xy", pch = 16, cex = 0.6, col = "#4c4c4c",
                   xlab = expression(omega), ylab = expression(G^"*"(omega)),
                   main = "Discrete fit vs observed", ...)
    graphics::points(crs$w, crs$Gexp[n + seq_len(n)], pch = 1, cex = 0.6, col = "#7f7f7f")
    graphics::lines(fit[, 1], fit[, 2], lwd = 2, col = "#c44e52")
    graphics::lines(fit[, 1], fit[, 3], lwd = 2, col = "#8172b3")
    graphics::legend("topright",
                     legend = c("G' observed", "G'' observed", "G' discrete fit", "G'' discrete fit"),
                     col = c("#4c4c4c", "#7f7f7f", "#c44e52", "#8172b3"),
                     pch = c(16, 1, NA, NA), lty = c(NA, NA, 1, 1), bty = "n")
  }

  invisible(x)
}
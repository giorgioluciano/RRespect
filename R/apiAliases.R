#' Continuous Spectrum Alias
#'
#' Alias for [getContinuousSpectrum()].
#'
#' @param ... Arguments passed to [getContinuousSpectrum()].
#' @return A list with continuous-spectrum results.
#' @export
contSpectrum <- function(...) {
  getContinuousSpectrum(...)
}

#' Discrete Spectrum Alias
#'
#' Alias for [getDiscreteSpectrum()].
#'
#' @param ... Arguments passed to [getDiscreteSpectrum()].
#' @return A list with discrete-spectrum results.
#' @export
discSpectrum <- function(...) {
  getDiscreteSpectrum(...)
}

#' Write Results Alias
#'
#' Alias for [writeSpectrumResults()].
#'
#' @param ... Arguments passed to [writeSpectrumResults()].
#' @return Invisibly returns written file paths.
#' @export
writeResults <- function(...) {
  writeSpectrumResults(...)
}

#' Time Batch Alias
#'
#' Convenience alias for running [runBatch()] in time domain.
#'
#' @param ... Arguments passed to [runBatch()].
#' @return Data frame with batch status and diagnostics.
#' @export
runTimeBatch <- function(...) {
  runBatch(domain = "time", ...)
}

#' Frequency Batch Alias
#'
#' Convenience alias for running [runBatch()] in frequency domain.
#'
#' @param ... Arguments passed to [runBatch()].
#' @return Data frame with batch status and diagnostics.
#' @export
runFrequencyBatch <- function(...) {
  runBatch(domain = "frequency", ...)
}

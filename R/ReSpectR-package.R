#' ReSpectR: Relaxation Spectrum Extraction in Time and Frequency Domains
#'
#' Tools for continuous and discrete relaxation spectrum extraction from
#' rheological data in both time-domain G(t) and frequency-domain G*(ω)
#' formulations.
#'
#' The algorithm applies Tikhonov regularisation with Bayesian lambda selection
#' via the L-curve method (continuous spectrum), followed by an AIC-guided
#' non-negative least-squares search for discrete Maxwell modes.
#'
#' ## Main workflow
#'
#' ```
#' par <- setParams(domain = "time", dataFile = "Gt.dat")
#' crs <- getContinuousSpectrum(par)
#' drs <- getDiscreteSpectrum(par, crs = crs)
#' writeSpectrumResults(drs, outputDir = "output")
#' ```
#'
#' ## Key functions
#'
#' | Function | Purpose |
#' |---|---|
#' | [setParams()] | Build a validated parameter list |
#' | [getContinuousSpectrum()] | Extract continuous H(τ) spectrum |
#' | [getDiscreteSpectrum()] | Extract discrete Maxwell modes |
#' | `plot()` | Plot classed spectrum result objects |
#' | [writeSpectrumResults()] | Write results to disk |
#' | [runBatch()] | Process multiple data files |
#'
#' Short-form aliases: [contSpectrum()], [discSpectrum()],
#' [writeResults()], [runTimeBatch()], [runFrequencyBatch()].
#'
#' Returned objects carry stable classes
#' `"respect_continuous_spectrum"` and `"respect_discrete_spectrum"`, so they
#' can be passed directly to `plot()` for built-in visualisation.
#'
#' @seealso
#' Package repository: \url{https://github.com/}
#'
#' @references
#' Shanbhag, S. (2019). pyReSpect: A computer program to extract discrete and
#' continuous spectra from stress relaxation experiments. *Macromolecular
#' Theory and Simulations*, 28(6), 1900005.
#'
#' @keywords internal
"_PACKAGE"

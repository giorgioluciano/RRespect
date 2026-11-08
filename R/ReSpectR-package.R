#' ReSpectR: Relaxation Spectrum Extraction in Time and Frequency Domains
#'
#' Tools for continuous and discrete relaxation spectrum extraction from
#' rheological data in the time domain and the frequency domain.
#'
#' Continuous spectra use Tikhonov regularisation with Bayesian selection of
#' the regularisation parameter, or a user-specified fixed value. Bayesian
#' selection should not be confused with selection from an L-curve corner.
#' The updated default workflow uses trust-region reflective least squares.
#' Discrete Maxwell spectra use an AIC-guided non-negative least-squares search,
#' nonlinear refinement, pruning and mode merging.
#'
#' The retained selector `solver = "experimental"` identifies the updated
#' default workflow. Historical solvers remain available with
#' `solver = "legacy"`; they do not require the removed prototype directory.
#' An optional equilibrium modulus can be fitted in either domain.
#'
#' ## Main workflow
#'
#' ```
#' par <- setParams(domain = "time", dataFile = "Gt.dat")
#' crs <- getContinuousSpectrum(par, writeOutput = FALSE)
#' drs <- getDiscreteSpectrum(par, crs = crs, writeOutput = FALSE)
#' writeSpectrumResults(drs, outputDir = "output")
#' ```
#'
#' ## Main functions
#'
#' - [setParams()]: construct a parameter list.
#' - [getContinuousSpectrum()]: estimate the continuous log-spectrum.
#' - [getDiscreteSpectrum()]: estimate discrete Maxwell modes.
#' - [writeSpectrumResults()]: explicitly write results.
#' - [runBatch()]: process multiple data files.
#'
#' Result objects support `plot()`. Short-form aliases include
#' [contSpectrum()], [discSpectrum()], [writeResults()], [runTimeBatch()]
#' and [runFrequencyBatch()].
#'
#' ## Documentation and interpretation
#'
#' The `ReSpectR-intro` vignette introduces the workflows.
#' The `ReSpectR-python-comparison` vignette illustrates comparison against
#' bundled historical Python output; it does not reproduce the separate
#' 23-configuration validation audit. Use `browseVignettes("ReSpectR")`
#' when vignettes are installed.
#'
#' The nominal AIC-selected mode count and the final surviving mode count
#' need not be equal. Local uncertainty estimates are not calibrated
#' confidence intervals. Agreement between implementations does not establish
#' uniqueness of the recovered spectrum or goodness of fit to observations.
#'
#' @seealso
#' Package repository: \url{https://github.com/giorgioluciano/RRespect}
#'
#' @references
#' Shanbhag, S. (2019). pyReSpect: A computer program to extract discrete and
#' continuous spectra from stress relaxation experiments. *Macromolecular
#' Theory and Simulations*, 28(6), 1900005.
#'
#' @keywords internal
"_PACKAGE"

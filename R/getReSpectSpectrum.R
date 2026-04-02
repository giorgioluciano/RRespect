#' Compute Continuous Spectrum
#'
#' Runs the legacy continuous-spectrum solver through a package API.
#' During migration this function bridges to domain-specific legacy code.
#'
#' @param par Parameter list from [setParams()].
#' @param writeOutput Logical. If `TRUE`, write solver outputs.
#' @param outputDir Output directory used when `writeOutput = TRUE`.
#' @param projectDir ReSpectR project root. Defaults to current working directory.
#'
#' @return A list with continuous-spectrum results.
#' @export
#' @examples
#' \dontrun{
#' par <- setParams(domain = "time", dataFile = "Gt.dat")
#' crs <- getContinuousSpectrum(par)
#' # crs$s  — relaxation-time grid
#' # crs$H  — log-spectrum values
#' }
getContinuousSpectrum <- function(
  par,
  writeOutput = FALSE,
  outputDir = "output",
  projectDir = NULL
) {
  .validatePar(par)
  domain    <- par$domain
  legacyPar <- .toLegacyParameters(par)

  if (domain == "time") {
    out      <- .gtContSpec(legacyPar, writeOutput = isTRUE(writeOutput),
                            outputDir = outputDir)
    out$Gfit <- cbind(out$t, .gtKernelPrestore(out$H, out$kernMat))
    return(.asRespectContinuousSpectrum(out, par))
  }

  out      <- .gstarContSpec(legacyPar, writeOutput = isTRUE(writeOutput),
                              outputDir = outputDir)
  Kfit <- .gstarKernelPrestore(out$H, out$kernMat, out$G0)
  n    <- length(out$w)
  out$Gfit <- cbind(out$w, Kfit[seq_len(n)], Kfit[n + seq_len(n)])
  .asRespectContinuousSpectrum(out, par)
}

#' Compute Discrete Spectrum
#'
#' Runs the legacy discrete-spectrum solver through a package API.
#' During migration this function bridges to domain-specific legacy code.
#'
#' @param par Parameter list from [setParams()].
#' @param crs Optional continuous-spectrum result object.
#' @param writeOutput Logical. If `TRUE`, write solver outputs.
#' @param outputDir Output directory used when `writeOutput = TRUE`.
#' @param projectDir ReSpectR project root. Defaults to current working directory.
#'
#' @return A list with discrete-spectrum results.
#' @export
#' @examples
#' \dontrun{
#' par <- setParams(domain = "time", dataFile = "Gt.dat")
#' crs <- getContinuousSpectrum(par)
#' drs <- getDiscreteSpectrum(par, crs = crs)
#' # drs$g    — mode weights
#' # drs$tau  — relaxation times
#' }
getDiscreteSpectrum <- function(
  par,
  crs = NULL,
  writeOutput = FALSE,
  outputDir = "output",
  projectDir = NULL
) {
  .validatePar(par)
  domain    <- par$domain
  legacyPar <- .toLegacyParameters(par)

  if (domain == "time") {
    if (is.null(crs)) {
      crs <- getContinuousSpectrum(
        par,
        writeOutput = FALSE,
        outputDir = outputDir,
        projectDir = projectDir
      )
    }
    out        <- .gtDiscSpec(legacyPar, crs,
                              writeOutput = isTRUE(writeOutput),
                              outputDir = outputDir)
    out$dmodes <- cbind(out$g, out$tau)
    return(.asRespectDiscreteSpectrum(out, par, crs))
  }

  if (is.null(crs)) {
    crs <- getContinuousSpectrum(
      par,
      writeOutput = isTRUE(writeOutput),
      outputDir = outputDir,
      projectDir = projectDir
    )
  }

  out        <- .gstarGetDiscSpecMagic(legacyPar,
                                       writeOutput = isTRUE(writeOutput),
                                       outputDir   = outputDir,
                                       contResult  = crs)
  out$dmodes <- cbind(out$g, out$tau)
  .asRespectDiscreteSpectrum(out, par, crs)
}

.asRespectContinuousSpectrum <- function(x, par) {
  x$domain      <- par$domain
  x$result_type <- "continuous"
  x$dataFile    <- par$dataFile
  class(x) <- c("respect_continuous_spectrum", "respect_spectrum_result", class(x))
  x
}

.asRespectDiscreteSpectrum <- function(x, par, crs) {
  x$domain      <- par$domain
  x$result_type <- "discrete"
  x$dataFile    <- par$dataFile
  x$continuous  <- crs
  class(x) <- c("respect_discrete_spectrum", "respect_spectrum_result", class(x))
  x
}

.validatePar <- function(par) {
  if (!is.list(par) || is.null(par$domain)) {
    stop("par must be a parameter list produced by setParams()")
  }
  if (!par$domain %in% c("time", "frequency")) {
    stop("par$domain must be either 'time' or 'frequency'")
  }
}

.toLegacyParameters <- function(par) {
  if (par$domain == "time") {
    return(list(
      GtFile = par$dataFile,
      ns = par$ns,
      lamC = par$lamC,
      SmFacLam = par$smFacLam,
      FreqEnd = par$freqEnd,
      verbose = par$verbose,
      plotting = par$plotting,
      lam_min = par$lamMin,
      lam_max = par$lamMax,
      lamDensity = par$lamDensity,
      rho_cutoff = 0,
      MaxNumModes = if (!is.null(par$maxNumModes)) par$maxNumModes else 0,
      deltaBaseWeightDist = if (!is.null(par$deltaBaseWeightDist)) par$deltaBaseWeightDist else 0.2,
      minTauSpacing = if (!is.null(par$minTauSpacing)) par$minTauSpacing else 1.25,
      condWt = if (!is.null(par$condWt)) par$condWt else 0.5,
      BaseDistWt = if (!is.null(par$baseDistWt)) par$baseDistWt else 0.5
    ))
  }

  list(
    GstFile = par$dataFile,
    ns = par$ns,
    lamC = par$lamC,
    SmFacLam = par$smFacLam,
    FreqEnd = par$freqEnd,
    verbose = par$verbose,
    plotting = par$plotting,
    lam_min = par$lamMin,
    lam_max = par$lamMax,
    lamDensity = par$lamDensity,
    plateau = if (!is.null(par$plateau)) par$plateau else FALSE,
    MaxNumModes = if (!is.null(par$maxNumModes)) par$maxNumModes else 0,
    deltaBaseWeightDist = if (!is.null(par$deltaBaseWeightDist)) par$deltaBaseWeightDist else 0.2,
    minTauSpacing = if (!is.null(par$minTauSpacing)) par$minTauSpacing else 1.25
  )
}

# Legacy bridge and project-root resolver removed: all domain code is now
# package-native. The projectDir parameter is kept for API compatibility but
# is no longer used internally.

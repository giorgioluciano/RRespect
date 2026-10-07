#' Compute Continuous Spectrum
#'
#' Runs a continuous-spectrum solver through the package API.
#' The experimental TRF option supports both continuous-spectrum domains.
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
  .validatePar(par, operation = "continuous")

  domain <- par$domain
  legacyPar <- .toLegacyParameters(par)

  if (domain == "time") {
    out <- .gtContSpec(
      legacyPar,
      writeOutput = isTRUE(writeOutput),
      outputDir = outputDir
    )

    out$Gfit <- cbind(
      out$t,
      .gtKernelPrestore(out$H, out$kernMat, out$G0)
    )

    return(.asRespectContinuousSpectrum(out, par))
  }

  out <- if (identical(legacyPar$solver, "experimental")) {
    .gstarContSpecExperimental(
      legacyPar,
      writeOutput = isTRUE(writeOutput),
      outputDir = outputDir
    )
  } else {
    .gstarContSpec(
      legacyPar,
      writeOutput = isTRUE(writeOutput),
      outputDir = outputDir
    )
  }

  Kfit <- .gstarKernelPrestore(
    out$H, out$kernMat, out$G0
  )

  n <- length(out$w)

  out$Gfit <- cbind(
    out$w,
    Kfit[seq_len(n)],
    Kfit[n + seq_len(n)]
  )

  .asRespectContinuousSpectrum(out, par)
}



#' Compute Discrete Spectrum
#'
#' Runs a discrete-spectrum solver through the package API.
#' The experimental option uses the updated Python-style workflow.
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

  if (identical(legacyPar$solver, "experimental")) {
    if (is.null(crs)) {
      crs <- getContinuousSpectrum(
        par,
        writeOutput = FALSE,
        outputDir = outputDir,
        projectDir = projectDir
      )
    }

    out <- .experimentalDiscSpec(
      legacyPar, crs, domain,
      writeOutput = isTRUE(writeOutput),
      outputDir = outputDir
    )
    out$dmodes <- cbind(out$g, out$tau, out$dtau)
    return(.asRespectDiscreteSpectrum(out, par, crs))
  }

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
  if (is.null(x$solver)) x$solver <- "legacy"

  if (
    identical(x$solver, "experimental") &&
    identical(par$domain, "time")
  ) {
    scan <- x$scan_result
    x$lam_C <- x$lamC
    x$G_fit <- as.numeric(x$Gfit[, 2])
    if (is.null(x$G0)) x$G0 <- 0.0

    x["lam"] <- list(if (is.null(scan)) NULL else scan$lam)
    x["rho"] <- list(if (is.null(scan)) NULL else scan$rho)
    x["eta"] <- list(if (is.null(scan)) NULL else scan$eta)
    x["log_P"] <- list(if (is.null(scan)) NULL else scan$logP)
    x["H_lam"] <- list(if (is.null(scan)) NULL else scan$Hlambda)
    x["dH"] <- list(if (is.null(scan)) NULL else scan$dH)

    x$reference <- list(
      repository = "shane5ul/pyReSpect-time",
      version = "2.1.0",
      commit = "cc2545c461deda614d0344c0731d9d5e55c04f43"
    )
  }

  x$result_type <- "continuous"
  x$dataFile    <- par$dataFile
  class(x) <- c("respect_continuous_spectrum", "respect_spectrum_result", class(x))
  x
}

.asRespectDiscreteSpectrum <- function(x, par, crs) {
  x$domain      <- par$domain
  x$solver <- if (is.null(par$solver)) "legacy" else par$solver
  x$result_type <- "discrete"
  x$dataFile    <- par$dataFile
  x$continuous  <- crs
  class(x) <- c("respect_discrete_spectrum", "respect_spectrum_result", class(x))
  x
}

.validatePar <- function(par, operation = "discrete") {
  if (!is.list(par) || is.null(par$domain)) {
    stop("par must be a parameter list produced by setParams()")
  }

  if (
    !is.character(par$domain) ||
    length(par$domain) != 1L ||
    is.na(par$domain) ||
    !par$domain %in% c("time", "frequency")
  ) {
    stop("par$domain must be either 'time' or 'frequency'")
  }

  solver <- if (is.null(par$solver)) "legacy" else par$solver

  if (
    !is.character(solver) ||
    length(solver) != 1L ||
    is.na(solver) ||
    !solver %in% c("legacy", "experimental")
  ) {
    stop("par$solver must be either 'legacy' or 'experimental'")
  }


}

.toLegacyParameters <- function(par) {
  if (par$domain == "time") {
    return(list(
      GtFile = par$dataFile,
	  solver = if (is.null(par$solver)) "legacy" else par$solver,
      ns = par$ns,
      lamC = par$lamC,
      SmFacLam = par$smFacLam,
      FreqEnd = par$freqEnd,
      verbose = par$verbose,
      plotting = par$plotting,
      lam_min = par$lamMin,
      lam_max = par$lamMax,
      lamDensity = par$lamDensity,
	  plateau = isTRUE(par$plateau),
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
	solver = if (is.null(par$solver)) "legacy" else par$solver,
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

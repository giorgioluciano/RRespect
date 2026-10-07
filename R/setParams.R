#' Build ReSpectR Parameters
#'
#' Creates a normalized parameter list for either time-domain or
#' frequency-domain spectrum workflows.
#'
#' @param domain Character scalar. Either "time" or "frequency".
#' @param dataFile Character scalar. Input data filename.
#' @param ns Integer. Number of spectrum grid points.
#' @param lamC Numeric. Regularization value. Use 0 for automatic selection.
#' @param smFacLam Numeric in \[-1, 1\]. Lambda smoothing shift.
#' @param freqEnd Integer in \{1, 2, 3\}. End-condition mode.
#' @param verbose Logical. Print progress messages.
#' @param plotting Logical. Enable plotting.
#' @param lamMin Numeric. Minimum lambda for scan.
#' @param lamMax Numeric. Maximum lambda for scan.
#' @param lamDensity Integer >= 2. Points per decade. If NULL, uses 3 for legacy and 2 for experimental.
#' @param plateau Logical. Infer a plateau modulus in either domain.
#' @param maxNumModes Integer >= 0. Hard cap on mode count.
#' @param deltaBaseWeightDist Numeric in (0, 1). AIC base-weight spacing.
#' @param minTauSpacing Numeric > 1. Minimum adjacent tau spacing.
#' @param condWt Numeric in \[0, 1\]. Condition-number weighting.
#' @param baseDistWt Numeric in \[0, 1\]. Base-distance weighting.
#' @param solver Solver selection: "experimental" (updated TRF workflow, default) or "legacy".
#'   Experimental workflows support continuous and discrete spectra in both
#'   domains. Their numerical compatibility is tested against Python 2.1.
#'
#' @return Named list of validated parameters.
#' @export
#' @examples
#' p <- setParams(domain = "time", dataFile = "Gt.dat")
#' p$ns        # 100
#' p$domain    # "time"
setParams <- function(
  domain = c("time", "frequency"),
  dataFile = NULL,
  ns = 100L,
  lamC = 0,
  smFacLam = 0,
  freqEnd = 1L,
  verbose = TRUE,
  plotting = FALSE,
  lamMin = 1e-10,
  lamMax = 1e3,
  lamDensity = NULL,
  plateau = FALSE,
  maxNumModes = 0L,
  deltaBaseWeightDist = 0.2,
  minTauSpacing = 1.25,
  condWt = 0.5,
  baseDistWt = 0.5,
  solver = c("experimental", "legacy")
) {
  domain <- match.arg(domain)
  solver <- match.arg(solver)

  if (is.null(lamDensity)) {
    lamDensity <- if (identical(solver, "experimental")) 2L else 3L
  }


  if (is.null(dataFile) || !nzchar(dataFile)) {
    dataFile <- if (domain == "time") "Gt.dat" else "Gst.dat"
  }

  if (lamDensity < 2) stop("lamDensity must be >= 2")
  if (smFacLam < -1 || smFacLam > 1) stop("smFacLam must be in [-1, 1]")
  if (condWt < 0 || condWt > 1) stop("condWt must be in [0, 1]")
  if (baseDistWt < 0 || baseDistWt > 1) stop("baseDistWt must be in [0, 1]")

  out <- list(
    domain = domain,
    solver = solver,
	dataFile = dataFile,
	ns = as.integer(ns),
    lamC = lamC,
    smFacLam = smFacLam,
    freqEnd = as.integer(freqEnd),
    verbose = isTRUE(verbose),
    plotting = isTRUE(plotting),
    lamMin = lamMin,
    lamMax = lamMax,
    lamDensity = as.integer(lamDensity),
    plateau = isTRUE(plateau),
    maxNumModes = as.integer(maxNumModes),
    deltaBaseWeightDist = deltaBaseWeightDist,
    minTauSpacing = minTauSpacing
  )

  if (domain == "time") {
    out$condWt <- condWt
    out$baseDistWt <- baseDistWt
  }

  out
}

#' Get ReSpectR Parameters
#'
#' Alias for [setParams()] for users who prefer a "get"-style constructor.
#'
#' @param ... Arguments passed to [setParams()].
#' @return Named list of validated parameters.
#' @export
#' @examples
#' p <- getParams(domain = "frequency", dataFile = "Gst.dat")
#' p$plateau   # FALSE
getParams <- function(...) {
  setParams(...)
}

# Compatibility wrappers for migration code.
setReSpectParameters <- function(...) {
  setParams(...)
}

getReSpectParameters <- function(...) {
  setParams(...)
}

setRespectParameters <- function(...) {
  setParams(...)
}

#' Run ReSpectR Batch Processing
#'
#' Executes continuous and discrete spectrum extraction for a set of input files.
#'
#' @param domain Character scalar. Either "time" or "frequency".
#' @param dataFiles Optional character vector of data file paths.
#' @param dataDir Optional directory to scan for data files when `dataFiles` is NULL.
#' @param pattern File pattern used with `dataDir`.
#' @param outputRoot Optional output root directory. If `NULL`, no files are written.
#' @param writeOutput Logical. Explicit control of writing outputs.
#' @param progress Logical. Show progress bar when available.
#' @param verbose Logical. Forwarded to solver parameters.
#' @param projectDir Deprecated; kept for backward compatibility. No longer used.
#' @param ... Extra parameters forwarded to [setParams()].
#'
#' @return Data frame summarizing per-file run status and diagnostics.
#' @export
#' @examples
#' \dontrun{
#' # Process all .dat files in a directory
#' results <- runBatch(
#'   domain    = "time",
#'   dataDir   = "path/to/data",
#'   outputRoot = "path/to/output"
#' )
#' results[results$status == "error", ]
#' }
runBatch <- function(
  domain = c("time", "frequency"),
  dataFiles = NULL,
  dataDir = NULL,
  pattern = "\\.dat$",
  outputRoot = NULL,
  writeOutput = !is.null(outputRoot),
  progress = TRUE,
  verbose = FALSE,
  projectDir = NULL,
  ...
) {
  domain <- match.arg(domain)

  files <- .resolveBatchFiles(dataFiles = dataFiles, dataDir = dataDir, pattern = pattern)
  if (length(files) == 0) stop("No data files found for batch processing")

  pb <- NULL
  if (isTRUE(progress) && requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      total = length(files),
      clear = FALSE,
      show_after = 0,
      format = "[runBatch] :bar :percent :current/:total :message"
    )
  }

  out <- vector("list", length(files))

  for (i in seq_along(files)) {
    filePath <- normalizePath(files[i], mustWork = FALSE)
    tag <- tools::file_path_sans_ext(basename(filePath))
    if (!is.null(pb)) pb$tick(tokens = list(message = tag))

    result <- tryCatch(
      {
        par <- setParams(
          domain = domain,
          dataFile = filePath,
          verbose = isTRUE(verbose),
          plotting = FALSE,
          ...
        )

        runOutDir <- if (isTRUE(writeOutput)) {
          if (is.null(outputRoot) || !nzchar(outputRoot)) {
            stop("outputRoot must be provided when writeOutput = TRUE")
          }
          dir <- file.path(outputRoot, tag, "output")
          dir.create(dir, recursive = TRUE, showWarnings = FALSE)
          dir
        } else {
          "output"
        }

        crs <- getContinuousSpectrum(
          par = par,
          writeOutput = isTRUE(writeOutput),
          outputDir = runOutDir,
          projectDir = projectDir
        )

        drs <- getDiscreteSpectrum(
          par = par,
          crs = crs,
          writeOutput = isTRUE(writeOutput),
          outputDir = runOutDir,
          projectDir = projectDir
        )

        data.frame(
          file = filePath,
          status = "ok",
          lamC = if (!is.null(crs$lamC)) crs$lamC else NA_real_,
          Nopt = if (!is.null(drs$Nopt)) drs$Nopt else if (!is.null(drs$g)) length(drs$g) else NA_real_,
          error = if (!is.null(drs$error)) drs$error else NA_real_,
          message = "",
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        data.frame(
          file = filePath,
          status = "error",
          lamC = NA_real_,
          Nopt = NA_real_,
          error = NA_real_,
          message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )

    out[[i]] <- result
  }

  do.call(rbind, out)
}

.resolveBatchFiles <- function(dataFiles = NULL, dataDir = NULL, pattern = "\\.dat$") {
  if (!is.null(dataFiles)) {
    return(as.character(dataFiles))
  }

  if (is.null(dataDir) || !nzchar(dataDir)) {
    stop("Provide either dataFiles or dataDir")
  }

  list.files(dataDir, pattern = pattern, full.names = TRUE)
}

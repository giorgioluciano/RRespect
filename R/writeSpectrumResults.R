#' Write Spectrum Results to Disk
#'
#' Writes a result object to disk. This is the intended explicit file-output
#' path for ReSpectR so compute functions can stay side-effect free.
#'
#' @param results List-like result object.
#' @param outputDir Character scalar output directory.
#' @param overwrite Logical. Overwrite existing files.
#' @return Invisibly returns a character vector of written file paths.
#' @export
#' @examples
#' res   <- list(s = c(0.01, 0.1, 1, 10), H = c(-3, 0, 0, -3))
#' paths <- writeSpectrumResults(res, outputDir = tempdir())
#' basename(paths)   # "H.dat"
writeSpectrumResults <- function(results, outputDir = "output", overwrite = TRUE) {
  if (!is.list(results)) stop("results must be a list")
  if (!nzchar(outputDir)) stop("outputDir must be a non-empty path")

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  written <- character(0)

  write_table_if_present <- function(obj, fileName) {
    fullPath <- file.path(outputDir, fileName)
    if (file.exists(fullPath) && !isTRUE(overwrite)) {
      stop("Refusing to overwrite existing file: ", fullPath)
    }
    utils::write.table(obj, fullPath, row.names = FALSE, col.names = FALSE)
    written <<- c(written, fullPath)
  }

  if (!is.null(results$H) && !is.null(results$s)) {
    write_table_if_present(cbind(results$s, results$H), "H.dat")
  }

  if (!is.null(results$Gfit)) {
    write_table_if_present(results$Gfit, "Gfit.dat")
  }

  if (!is.null(results$dmodes)) {
    write_table_if_present(results$dmodes, "dmodes.dat")
  }

  if (!is.null(results$aic)) {
    write_table_if_present(results$aic, "aic.dat")
  }

  invisible(written)
}

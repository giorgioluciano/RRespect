.createProgressBar <- function(total, format = "[:bar] :percent :current/:total") {
  if (!requireNamespace("progress", quietly = TRUE)) {
    return(NULL)
  }
  progress::progress_bar$new(
    total = total,
    clear = FALSE,
    show_after = 0,
    format = format
  )
}

.tickProgressBar <- function(pb) {
  if (!is.null(pb)) pb$tick()
  invisible(NULL)
}

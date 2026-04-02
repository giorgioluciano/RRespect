# Load package R files for repository-local test execution (dev mode only).
r_dir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
if (dir.exists(r_dir)) {
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) {
    source(f, local = FALSE)
  }
}

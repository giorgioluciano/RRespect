# Global Imports
library(ggplot2)      # For plotting
library(dplyr)        # For data manipulation
library(tidyr)        # For data tidying
library(pracma)       # For numerical methods, if needed


# Impostazioni di visualizzazione grafica
## Functions common to both discrete and continuous spectra

readInput <- function(fname = 'inp.dat') {
  # Reads data from the input file and populates the parameter list
  par <- list()
  
  # Read the input file
  lines <- readLines(fname)
  
  for (line in lines) {
    li <- trimws(line)
    
    # Skip comments and empty lines
    if (nchar(li) > 0 && substr(li, 1, 1) != "#") {
      parts <- strsplit(li, ":")[[1]]
      key <- trimws(parts[1])
      
      # Replace Python-style booleans with R-style booleans
      val <- trimws(parts[2])
      val <- gsub("False", "FALSE", val)
      val <- gsub("True", "TRUE", val)
      
      # Evaluate the parameter and add it to the list
      par[[key]] <- eval(parse(text = val))
    }
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists("output")) {
    dir.create("output")
  }
  
  return(par)
}

GetExpData <- function(fname) {
  
  data <- tryCatch({
    as.matrix(read.table(fname, header = FALSE))
  }, error = function(e) {
    stop("*Error*: G(t) data file is either not in the correct path or incorrectly formatted.")
  })
  cols <- ncol(data)
  if (cols == 2) {
    to <- data[, 1]    
    Gto <- data[, 2]   
    wGo <- rep(1, length(Gto))  
  } else if (cols == 3) {
    to <- data[, 1]    
    Gto <- data[, 2]   
    wGo <- data[, 3]   
  } else {
    stop("*Error*: G(t) data file is either not in the correct path, or incorrectly formatted")
  }
  
  unique_indices <- !duplicated(to)
  to <- to[unique_indices]
  Gto <- Gto[unique_indices]
  wGo <- wGo[unique_indices]
  
  return(list(t = to, Gt = Gto, wG = wGo))  
}

getKernMat <- function(s, t) {
  # Function: getKernMat(s, t)
  # Furnish kerMat() which helps faster kernel evaluation
  # Given s, t generates hs * exp(-T/S) [n * ns matrix], where hs = wi = weights
  # for trapezoidal rule integration.
  
  ns <- length(s)
  hsv <- numeric(ns)
  hsv[1] <- 0.5 * log(s[2] / s[1])
  hsv[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  
  if (ns > 2) {
    hsv[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  }
  
  S <- outer(s, rep(1, length(t)))  # Create a matrix of s values
  T <- outer(rep(1, length(s)), t)  # Create a matrix of t values
  res <- as.matrix(exp(-T / S) * hsv)
  return(res)  # Element-wise multiplication
  
}

kernel_prestore <- function(H, kernMat, G0 = 0) {
  # Function: kernel_prestore(H, kernMat, G0)
  # Turbocharging kernel function evaluation by prestoring kernel matrix.
  # Returns K*h, where h = exp(H).
  print("Entering kernel_prestore function")
  print(paste("Length of H:", length(H)))
  print(paste("Dimensions of kernMat:", dim(kernMat)))
  print(paste("G0:", G0))
  
  result <- kernMat %*% exp(H) + G0
  
  print("Exiting kernel_prestore function")
  print(paste("Length of result:", length(result)))
  return(result)
}


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


  if (cols == 3) {
    wo <- data[, 1]
    Gpo <- data[, 2]
    Gppo <- data[, 3]
    WG1 <- WG2 <- NULL  # Placeholder per mantenere la consistenza
  } else if (cols > 3) {
    wo <- data[, 1]
    Gpo <- data[, 2]
    Gppo <- data[, 3]
    WG1 <- data[, 4]
    WG2 <- data[, 5]
  } else {
    stop("*Error*: G(t) data file is either not in the correct path, or incorrectly formatted")
  }
  
 
  unique_indices <- !duplicated(wo)
  wo <- wo[unique_indices]
  Gpo <- Gpo[unique_indices]
  Gppo <- Gppo[unique_indices]
  if (!is.null(WG1)) {
    WG1 <- WG1[unique_indices]
    WG2 <- WG2[unique_indices]
  }
  

  if (cols == 3) {
    fp <- approxfun(wo, Gpo, rule = 2)
    fpp <- approxfun(wo, Gppo, rule = 2)
    w <- exp(seq(log(min(wo)), log(max(wo)), length.out = 100))
    Gp <- fp(w)
    Gpp <- fpp(w)
    Gst <- c(Gp, Gpp)  
    wexp <- rep(1, length(Gst))  
    
    return(list(w = w, Gst = Gst, wexp = wexp))
  } else {
    Gst <- c(Gpo, Gppo)
    wexp <- c(WG1, WG2)
    
    return(list(w = wo, Gst = Gst, wexp = wexp))
  }
}



getKernMat <- function(s, w) {
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
  
  S <- outer(s, rep(1, length(w)))  # Create a matrix of s values
  W <- outer(rep(1, length(s)), w)  # Create a matrix of w values
  S <- t(S)
  W <- t(W)
  ws <- S * W
  ws2 <- ws^2
  res <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))
  
  # Assumendo che hsv sia un vettore e vada replicato per righe, fai:
  hsv_matrix <- matrix(hsv, nrow = nrow(res), ncol = ncol(res), byrow = TRUE)
  
  # Ora moltiplica elemento per elemento
  res <- res * hsv_matrix
  return(res)
  
  
}

kernel_prestore <- function(H, kernMat, G0 = 0) {
 
  
  H=as.matrix(H)
  G0 <- as.matrix(G0)  #
  cat("Class of H:", class(H), "\n")
  cat("Type of H:", typeof(H), "\n")
  print(H)
  
  n <- nrow(kernMat) / 2
  
  
  if (!is.null(G0)) {
    G0v <- numeric(2 * n)  
    G0v[1:n] <- G0         
  } else {
    G0v <- 0  
  }
  
  # Calcolo del risultato: kernMat * exp(H) + G0v
  result <- kernMat %*% exp(H) + G0v
  
  return(result)
}
# 7/2023: allowing an optional weight column in the input data file
#         improving encapsulation of functions

# Help to find continuous spectrum
# March 2019 major update:
# (i)   added plateau modulus G0 (also in pyReSpect-time) calculation
# (ii)  following Hansen Bayesian interpretation of Tikhonov to extract p(lambda)
# (iii) simplifying lcurve (starting from high lambda to low)
# (iv)  changing definition of rho2 and eta2 (no longer dividing by 1/n and 1/nl)


InitializeH <- function(Gexp, wexp, s, kernMat, ...) {
  # Use message() for logging instead of print(), which is more flexible for debugging
  message("Debug: InitializeH called with Gexp dimensions: ", paste(dim(Gexp), collapse = "x"))
  message("Debug: wexp dimensions: ", paste(dim(wexp), collapse = "x"))
  message("Debug: s length: ", length(s))
  message("Debug: kernMat dimensions: ", paste(dim(kernMat), collapse = "x"))
  
  # Initial guess for H
  H <- -5.0 * rep(1, length(s)) + sin(pi * s)
  lam <- 1.0
  
  args <- list(...)
  
  # Error handling using tryCatch
  tryCatch({
    if (length(args) > 0) {
      G0 <- args[[1]]
      message("Debug: Calling getH with G0 = ", G0)
      
      # Assuming getH returns a list with Hlam and G0
      result <- getH(lam, Gexp, wexp, H, kernMat, G0)
      Hlam <- result$Hlam
      G0 <- result$G0
      
      # Log the returned shapes (or lengths) for clarity
      message("Debug: getH returned with Hlam length: ", length(Hlam))
      message("Debug: G0: ", G0)
      
      return(list(Hlam = Hlam, G0 = G0))
    } else {
      message("Debug: Calling getH without G0")
      
      # Assuming getH returns Hlam when G0 is not passed
      Hlam <- getH(lam, Gexp, wexp, H, kernMat)
      
      message("Debug: getH returned with Hlam length: ", length(Hlam))
      
      return(Hlam)
    }
  }, error = function(e) {
    # Catch and display errors
    message("Error in InitializeH: ", e$message)
    message("Debug: Error type: ", class(e)[1])
    message("Debug: Error details: ", capture.output(str(e)))
    stop(e)  # Re-raise the error after logging
  })
}


# getAmatrix <- function(ns) 
#             {
#             L <- matrix(0, ns, ns)
#             diag(L) <- -2
#             if (ns > 1) {
#             L[cbind(1:(ns-1), 2:ns)] <- 1 }
#             if (ns > 1) {
#             L[cbind(2:ns, 1:(ns-1))] <- 1 }
#             # Return the matrix A = t(L) %*% L
#             return(t(L) %*% L)}


getAmatrix <- function(ns) {
  # Check for valid input size
  if (ns < 2) {
    stop("Error: 'ns' must be greater than or equal to 2")
  }
  
  # Debug message
  message("Debug: getAmatrix called with ns = ", ns)
  
  # Create the tridiagonal matrix L
  L <- diag(-2, ns, ns)          # Main diagonal with -2
  if (ns > 1) {
    L[row(L) == col(L) - 1] <- 1 # Upper diagonal with 1
    L[row(L) == col(L) + 1] <- 1 # Lower diagonal with 1
  }
  
  # Submatrix L[2:(ns-1), ] as in Python version
  nl <- ns - 2
  L <- L[2:(nl+1), ]
  
  # Return the matrix A = t(L) %*% L
  tryCatch({
    A <- t(L) %*% L
    message("Debug: Matrix A successfully computed with dimensions: ", paste(dim(A), collapse = "x"))
    return(A)
  }, error = function(e) {
    message("Error in getAmatrix: ", e$message)
    stop(e)  # Propagate the error
  })
}

# 
# getBmatrix <- function(H, kernMat, Gexp, wexp, ...) {
#   # Get dimensions of kernMat
#   n  <- dim(kernMat)[1]
#   ns <- dim(kernMat)[2]
#   nl <- ns - 2
#   
#   # Initialize residual vector r of size n
#   r <- rep(0, n)
#   
#   # Create Kmatrix based on Gexp and wexp
#   Kmatrix <- (wexp / Gexp) %*% matrix(1, nrow = 1, ncol = ns)
#   
#   # Jacobian residual matrix
#   Jr <- -kernelD(H, kernMat) * Kmatrix
#   
#   # Check if optional argument G0 is provided
#   args <- list(...)
#   if (length(args) > 0) {
#     G0 <- args[[1]]
#     Kh <- kernel_prestore(H, kernMat, G0)
#     pritn("Kh")
#     print(Kh)
#     r <- wexp * (1 - Kh / Gexp)
# 
#   } else {
#     r <- wexp * (1 - kernel_prestore(H, kernMat) / Gexp)
#     print("r")
#     print(r)
#   }
#   
#   # Compute matrix B
#   B <- t(Jr) %*% Jr + diag(as.vector(t(r) %*% Jr))
#   
#   return(B)
# }

getBmatrix <- function(H, kernMat, Gexp, wexp, ...) {
  # Check input dimensions for consistency
  if (length(H) != dim(kernMat)[2] || length(Gexp) != dim(kernMat)[1] || length(wexp) != length(Gexp)) {
    stop("Error: Dimensions of H, kernMat, Gexp, or wexp do not match.")
  }
  
  if (any(Gexp == 0)) {
    stop("Error: Gexp contains zero, leading to division by zero.")
  }
  
  # Debugging information
  message("Debug: getBmatrix called with H length: ", length(H))
  message("Debug: kernMat dimensions: ", paste(dim(kernMat), collapse = "x"))
  message("Debug: Gexp length: ", length(Gexp), ", wexp length: ", length(wexp))
  
  # Get dimensions of kernMat
  n  <- dim(kernMat)[1]
  ns <- dim(kernMat)[2]
  nl <- ns - 2
  
  # Initialize residual vector r
  r <- rep(0, n)
  
  # Create Kmatrix based on Gexp and wexp
  Kmatrix <- (wexp / Gexp) %*% matrix(1, nrow = 1, ncol = ns)
  
  # Compute Jacobian residual matrix Jr using kernelD function
  Jr <- -kernelD(H, kernMat) * Kmatrix
  
  # Check if optional argument G0 is provided
  args <- list(...)
  tryCatch({
    if (length(args) > 0) {
      G0 <- args[[1]]
      message("Debug: G0 provided")
      Kh <- kernel_prestore(H, kernMat, G0)
      message("Debug: Kh calculated")
      r <- wexp * (1 - Kh / Gexp)
    } else {
      message("Debug: No G0 provided, calculating without it")
      r <- wexp * (1 - kernel_prestore(H, kernMat) / Gexp)
    }
    
    # Compute matrix B
    B <- t(Jr) %*% Jr + diag(as.vector(t(r) %*% Jr))
    
    message("Debug: B matrix successfully computed with dimensions: ", paste(dim(B), collapse = "x"))
    return(B)
  }, error = function(e) {
    message("Error in getBmatrix: ", e$message)
    stop(e)  # Propagate the error
  })
}

 
# 
# 
# lcurve <- function(Gexp, wexp, Hgs, kernMat, par, ...) {
#   # Optional argument handling
#   args <- list(...)
#   
#   # If plateau is set, extract G0
#   if (par$plateau) {
#     G0 <- args[[1]]
#   }
#   
#   # Number of lambda points
#   npoints <- as.integer(par$lamDensity * (log10(par$lam_max) - log10(par$lam_min)))
#   
#   # Generate lambda values
#   hlam <- (par$lam_max / par$lam_min)^(1 / (npoints - 1))
#   lam <- par$lam_min * hlam^(0:(npoints-1))
#   
#   # Initialize vectors
#   eta <- rep(0, npoints)
#   rho <- rep(0, npoints)
#   logP <- rep(0, npoints)
#   H <- Hgs
#   n <- length(Gexp)
#   ns <- length(H)
#   nl <- ns - 2
#   logPmax <- -Inf  # Set to negative infinity for comparison
#   Hlambda <- matrix(0, nrow = ns, ncol = npoints)
#   
#   # Error Analysis: Create Amatrix
#   Amat <- getAmatrix(length(H))
#   LogDetN <- determinant(Amat, logarithm = TRUE)$modulus
#   
#   # Iterate over lambda values (from largest to smallest)
#   for (i in seq(npoints, 1)) {
#     lamb <- lam[i]
#     
#       H <- result[[1]]
#       G0 <- result[[2]]
#       rho[i] <- norm(wexp * (1 - kernel_prestore(H, kernMat, G0) / Gexp), type = "2")
#       Bmat <- getBmatrix(H, kernMat, Gexp, wexp, G0)
#     } else {
#       H <- getH(lamb, Gexp, wexp, H, kernMat)
#       rho[i] <- norm(wexp * (1 - kernel_prestore(H, kernMat) / Gexp), type = "2")
#       Bmat <- getBmatrix(H, kernMat, Gexp, wexp)
#     }
#     
#     eta[i] <- norm(diff(H, differences = 2), type = "2")
#     Hlambda[, i] <- H
#     
#     # Compute the log determinant of lamb*Amat + Bmat
#     LogDetC <- determinant(lamb * Amat + Bmat, logarithm = TRUE)$modulus
#     V <- rho[i]^2 + lamb * eta[i]^2
#     
#     # Calculate logP (assuming a prior exp(-lam))
#     logP[i] <- -V + 0.5 * (LogDetN + ns * log(lamb) - LogDetC) - lamb
#     
#     if (logP[i] > logPmax) {
#       logPmax <- logP[i]
#     } else if (logP[i] < logPmax - 18) {
#       break
#     }
#   }
#   
#   # Truncate results to significant lambda values
#   lam <- lam[i:npoints]
#   logP <- logP[i:npoints]
#   eta <- eta[i:npoints]
#   rho <- rho[i:npoints]
#   logP <- logP - max(logP)
#   Hlambda <- Hlambda[, i:npoints]
#   
#   # Compute optimal lambda (lamM) using a probability-weighted mean
#   plam <- exp(logP)
#   plam <- plam / sum(plam)
#   lamM <- exp(sum(plam * log(lam)))
#   
#   # Smoothness factor adjustments
#   if (par$SmFacLam > 0) {
#     lamM <- exp(log(lamM) + par$SmFacLam * (max(log(lam)) - log(lamM)))
#   } else if (par$SmFacLam < 0) {
#     lamM <- exp(log(lamM) + par$SmFacLam * (log(lamM) - min(log(lam))))
#   }
#   
#   # Optional plotting
#   pdf(file = "output/logP.pdf")
#   if (par$plotting) {
#     plot(logP ~ lam, log = "x", type = "o", ylim = c(-20, 1), xlab = expression(lambda), ylab = expression(log(p(lambda))))
#     abline(v = lamM, col = "gray", lty = 2)
#     legend("topleft", legend = expression(lambda[M]), col = "gray", lty = 2)
#     #dev.copy2pdf(file = "output/logP.pdf")
#   }
#   dev.off()
#   return(list(lamM = lamM, lam = lam, rho = rho, eta = eta, logP = logP, Hlambda = Hlambda))
# }

      lcurve <- function(Gexp, wexp, Hgs, kernMat, par, ...) {
        # Debug information
        message("Debug: lcurve called with Gexp length: ", length(Gexp), ", Hgs length: ", length(Hgs))
        
        # Check if plateau parameter is TRUE and set G0 if available
        args <- list(...)
        if (par$plateau && length(args) > 0) {
          G0 <- args[[1]]
        }
        
        # Setup for lambda range
        npoints <- as.integer(par$lamDensity * (log10(par$lam_max) - log10(par$lam_min)))
        hlam <- (par$lam_max / par$lam_min)^(1.0 / (npoints - 1))
        lam <- par$lam_min * hlam^(0:(npoints - 1))
        
        # Initialize variables for eta, rho, logP
        eta <- numeric(npoints)
        rho <- numeric(npoints)
        logP <- numeric(npoints)
        H <- Hgs
        n <- length(Gexp)
        ns <- length(H)
        nl <- ns - 2
        logPmax <- -Inf
        Hlambda <- matrix(0, nrow = ns, ncol = npoints)
        
        # Furnish A_matrix for error analysis
        Amat <- getAmatrix(length(H))
        logDetN <- determinant(Amat)$modulus
        
        # Main loop over lambda values
        for (i in rev(1:length(lam))) {
          lamb <- lam[i]
          message("Debug: lcurve iteration ", i, ", lambda: ", lamb)
          
          # Compute H, rho, and Bmat based on plateau condition
          if (par$plateau && length(args) > 0) {
            H <- getH(lamb, Gexp, wexp, H, kernMat, G0)
            rho[i] <- sqrt(sum((wexp * (1 - kernel_prestore(H, kernMat, G0) / Gexp))^2))
            Bmat <- getBmatrix(H, kernMat, Gexp, wexp, G0)
          } else {
            H <- getH(lamb, Gexp, wexp, H, kernMat)
            rho[i] <- sqrt(sum((wexp * (1 - kernel_prestore(H, kernMat) / Gexp))^2))
            Bmat <- getBmatrix(H, kernMat, Gexp, wexp)
          }
          
          # Compute eta and Hlambda
          eta[i] <- sqrt(sum(diff(H, differences = 2)^2))
          Hlambda[, i] <- H
          
          # Compute determinant for regularization term
          logDetC <- determinant(lamb * Amat + Bmat)$modulus
          V <- rho[i]^2 + lamb * eta[i]^2
          
          # Compute log posterior probability
          logP[i] <- -V + 0.5 * (logDetN + ns * log(lamb) - logDetC) - lamb
          
          # Update logPmax and check for cutoff
          if (logP[i] > logPmax) {
            logPmax <- logP[i]
          } else if (logP[i] < logPmax - 18) {
            break
          }
        }
        
        # Truncate values to significant lambda range
        lam <- lam[i:length(lam)]
        logP <- logP[i:length(logP)]
        eta <- eta[i:length(eta)]
        rho <- rho[i:length(rho)]
        logP <- logP - max(logP)
        Hlambda <- Hlambda[, i:npoints]
        
        # Compute optimal lambda using weighted log(P)
        plam <- exp(logP)
        plam <- plam / sum(plam)
        lamM <- exp(sum(plam * log(lam)))
        
        # Adjust Smoothness Factor
        if (par$SmFacLam > 0) {
          lamM <- exp(log(lamM) + par$SmFacLam * (max(log(lam)) - log(lamM)))
        } else if (par$SmFacLam < 0) {
          lamM <- exp(log(lamM) + par$SmFacLam * (log(lamM) - min(log(lam))))
        }
        
        # Plotting if enabled
        if (par$plotting) {
          plot(log10(lam), logP, type = "o", col = "blue", xlab = expression(log(lambda)), ylab = expression(log(p(lambda))), main = "L-Curve")
          abline(v = log10(lamM), col = "gray", lwd = 2)
        }
        
        return(list(lamM = lamM, lam = lam, rho = rho, eta = eta, logP = logP, Hlambda = Hlambda))
      }
      
# getH <- function(lam, Gexp, wexp, H, kernMat, ...) {
#   print("Entering getH function")
#   print(paste("lam:", lam))
#   print(paste("Length of Gexp:", length(Gexp)))
#   print(paste("Length of wexp:", length(wexp)))
#   print(paste("Length of H:", length(H)))
#   print(paste("Dimensions of kernMat:", dim(kernMat)))
#   
#   library(minpack.lm)
#   
#   # Handle optional arguments
#   args <- list(...)
#   
#   # Try optimization with error handling
#   result <- tryCatch({
#     if (length(args) > 0) {
#       G0 <- args[[1]]  # Get G0 from the arguments
#       Hplus <- c(H, G0)  # Append G0 to H
#       
#       # Optimization with G0
#       res_lm <- nls.lm(
#         par = Hplus,
#         fn = residualLM,  # Residual function to minimize
#         lam = lam, Gexp = Gexp, wexp = wexp, kernMat = kernMat
#       )
#       print("Debug: Optimization completed with G0")
#       
#       # Return optimized H and G0
#       return(list(H_lam = res_lm$par[1:(length(Hplus) - 1)], G0 = res_lm$par[length(Hplus)]))
#       
#     } else {
#       # Optimization without G0
#       res_lm <- nls.lm(
#         par = H,
#         fn = residualLM,  # Residual function to minimize
#         lam = lam, Gexp = Gexp, wexp = wexp, kernMat = kernMat
#       )
#       print("Exiting getH function")
#       return(res_lm$par)
#     }
#   }, error = function(e) {
#     print(paste("Error in getH:", e$message))
#     NULL  # Return NULL if an error occurs
#   })
#   
#   return(result)
# }
      
      
      getH <- function(lam, Gexp, wexp, H, kernMat, ...) {
        message("Entering getH function")
        message(paste("lam:", lam))
        message(paste("Length of Gexp:", length(Gexp)))
        message(paste("Length of wexp:", length(wexp)))
        message(paste("Length of H:", length(H)))
        message(paste("Dimensions of kernMat:", dim(kernMat)))
        
        library(minpack.lm)
        
        # Handle optional arguments
        args <- list(...)
        
        # Try optimization with error handling
        result <- tryCatch({
          if (length(args) > 0) {
            G0 <- args[[1]]  # Get G0 from the arguments
            Hplus <- c(H, G0)  # Append G0 to H
            
            # Optimization with G0
            res_lm <- nls.lm(
              par = Hplus,
              fn = residualLM,
              # Residual function to minimize
              lam = lam,
              Gexp = Gexp,
              wexp = wexp,
              kernMat = kernMat
            )
            message("Debug: Optimization completed with G0")
            
            # Return optimized H and G0
            return(list(H_lam = res_lm$par[1:(length(Hplus) - 1)], G0 = res_lm$par[length(Hplus)]))
            
          } else {
            # Optimization without G0
            res_lm <- nls.lm(
              par = H,
              fn = residualLM,
              # Residual function to minimize
              lam = lam,
              Gexp = Gexp,
              wexp = wexp,
              kernMat = kernMat
            )
            message("Exiting getH function")
            return(res_lm$par)
          }
        }, error = function(e) {
          message(paste("Error in getH:", e$message))
          NULL  # Return NULL if an error occurs
        })
        
        return(result)
      }
      
# 
# residualLM <- function(H, lam, Gexp, wexp, kernMat) {
#   #    HELPER FUNCTION: Gets Residuals r
#   #    Input  : H       = guessed H,
#   #             lambda  = regularization parameter ,
#   #             Gexp    = experimental data,
#   #             wexp    = weighting factors,
#   #             kernMat = matrix for faster kernel evaluation
#   #             G0      = plateau (if provided)
#   #   
#   #    Output : a set of n+nl residuals,
#   #             the first n correspond to the kernel
#   #             the last  nl correspond to the smoothness criterion
#   print("Entering residualLM function")
#   print(paste("Length of H:", length(H)))
#   print(paste("lam:", lam))
#   print(paste("Length of Gexp:", length(Gexp)))
#   print(paste("Length of wexp:", length(wexp)))
#   print(paste("Dimensions of kernMat:", dim(kernMat)))
#   n <- nrow(kernMat)
#   ns <- ncol(kernMat)
#   nl <- ns - 2
#   
#   r <- rep(0, n + nl)  
#   
#   if (length(H) > ns) {
#     G0 <- tail(H, 1)  
#     H <- head(H, -1)  
#     Kh <- kernel_prestore(H, kernMat, G0)
#     
#     r[1:n] <- wexp * (Kh - Gexp)
#   } else {
#     Kh <- kernel_prestore(H, kernMat) 
#     
#     r[1:n] <- wexp * (Kh - Gexp)
#   }
#   
#   r[(n + 1):(n + nl)] <- sqrt(lam) * diff(H, differences = 2)
#   
#   print("r")
#   print(r)
#   
#   # Debug
#   
#   print("Exiting residualLM function")
#   print(paste("Length of H:", length(H)))
#   print(paste("Length of Gexp:", length(Gexp)))
#   print(paste("Dimensions of kernMat:", dim(kernMat)))
#   print(paste("Length of r:", length(r)))
#   return(r)
#   return(r)  
# }
      
      residualLM <- function(H, lam, Gexp, wexp, kernMat) {
        #    HELPER FUNCTION: Gets Residuals r
        #    Input  : H       = guessed H,
        #             lambda  = regularization parameter ,
        #             Gexp    = experimental data,
        #             wexp    = weighting factors,
        #             kernMat = matrix for faster kernel evaluation
        #             G0      = plateau (if provided)
        #   
        #    Output : a set of n+nl residuals,
        #             the first n correspond to the kernel
        #             the last  nl correspond to the smoothness criterion
        message("Entering residualLM function")
        message(paste("Length of H:", length(H)))
        message(paste("lam:", lam))
        message(paste("Length of Gexp:", length(Gexp)))
        message(paste("Length of wexp:", length(wexp)))
        message(paste("Dimensions of kernMat:", dim(kernMat)))
        
        n <- nrow(kernMat)
        ns <- ncol(kernMat)
        nl <- ns - 2
        
        r <- rep(0, n + nl)  
        
        if (length(H) > ns) {
          G0 <- tail(H, 1)  
          H <- head(H, -1)  
          Kh <- kernel_prestore(H, kernMat, G0)
          
          r[1:n] <- wexp * (Kh - Gexp)
        } else {
          Kh <- kernel_prestore(H, kernMat) 
          
          r[1:n] <- wexp * (Kh - Gexp)
        }
        
        r[(n + 1):(n + nl)] <- sqrt(lam) * diff(H, differences = 2)
        
        message("Residual vector r: ", paste(r, collapse = ", "))
        
        message("Exiting residualLM function")
        message(paste("Length of H:", length(H)))
        message(paste("Length of Gexp:", length(Gexp)))
        message(paste("Dimensions of kernMat:", dim(kernMat)))
        message(paste("Length of r:", length(r)))
        
        return(r)
      }     
 
      
# 
# jacobianLM <- function(H, lam, Gexp, wexp, kernMat) {
#   n <- nrow(kernMat)
#   ns <- ncol(kernMat)
#   nl <- ns - 2
#   
#   # Crea la matrice tridiagonale L
#   L <- matrix(0, nrow = ns, ncol = ns)
#   diag(L) <- -2
#   diag(L[,-1]) <- 1
#   diag(L[,-ncol(L)]) <- 1
#   L <- L[2:(nl+1), ]  
#   
#   Kmatrix <- (wexp / Gexp) %*% matrix(1, 1, ns)
#   
#   if (length(H) > ns) {
#     G0 <- tail(H, 1)  
#     H <- head(H, -1)  
#     
#     Jr <- matrix(0, n + nl, ns + 1)
#     
#     Jr[1:n, 1:ns] <- -kernelD(H, kernMat) * Kmatrix
#     Jr[1:n, ns + 1] <- -(wexp / Gexp)  
#     
#     
#     Jr[(n + 1):(n + nl), 1:ns] <- sqrt(lam) * L
#     Jr[(n + 1):(n + nl), ns + 1] <- 0  
#     
#   } else {
#     
#     Jr <- matrix(0, n + nl, ns)
#     
#     Jr[1:n, 1:ns] <- -kernelD(H, kernMat) * Kmatrix
#     Jr[(n + 1):(n + nl), 1:ns] <- sqrt(lam) * L
#   }
#   
#   
#   gradient <- colSums(Jr)
#   
#   return(gradient)
# }

# 
# kernelD <- function(H, kernMat) {
#   n <- nrow(kernMat)
#   ns <- ncol(kernMat)
#   
#   # Debug: Log the input shapes
#   cat("Debug: Entering kernelD function\n")
#   cat(paste("Length of H:", length(H), "\n"))
#   cat(paste("Dimensions of kernMat:", dim(kernMat), "\n"))
#   
#   Hsuper <- matrix(rep(exp(H), each = n), nrow = n, byrow = TRUE)
#   
#   DK <- kernMat * Hsuper
#   
#   # Debug: Log the output shape
#   cat(paste("Debug: Output DK dimensions:", dim(DK), "\n"))
#   
#   return(DK)
# }


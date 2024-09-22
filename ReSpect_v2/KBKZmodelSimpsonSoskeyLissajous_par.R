library(parallel)

KBKZmodelSimpsonSoskeyLissajous <- function(omega, Amp, aSos, bSos, gi, taui, verbose = FALSE) {
  PointsPerPeriod <- 125
  NoPeriods <- 200
  CutPeriods <- 100
  etai <- taui * gi
  
  mint <- 0
  maxt <- 2 * pi / omega * NoPeriods
  tspan <- c(mint, maxt)
  tstep <- 2 * pi / omega / PointsPerPeriod
  
  time <- seq(mint, maxt, by = tstep)
  dim <- length(time)
  tau12 <- numeric(dim)
  
  PreSum <- function(tp, t) {
    Amp * (sin(omega * t) - sin(omega * tp)) * 1.0 / (1.0 + aSos * (abs(Amp * (sin(omega * t) - sin(omega * tp))))^bSos)
  }
  
  Integrand <- function(tp, t, PreSumValue) {
    PreSumValue * sum(gi / taui * exp(-(t - tp) / taui))
  }
  
  SimpsonOdd <- function(integrand, h) {
    n <- length(integrand)
    
    if (n < 3) {
      return(sum(integrand) * h)
    }
    
    if (n == 3) {
      return((integrand[1] + 4*integrand[2] + integrand[3]) * h / 3)
    }
    
    return(
      (integrand[1] + 
         4 * sum(integrand[seq(2, n - 1, by = 2)]) + 
         2 * sum(integrand[seq(3, n - 2, by = 2)]) +
         integrand[n]) * h / 3
    )
  }
  
  SimpsonThreeEight <- function(integrand, h) {
    3 * h / 8 * (integrand[1] + 3 * integrand[2] + 3 * integrand[3] + integrand[4])
  }
  
  # Initialize parallel cluster
  cores <- detectCores() - 1  # Use all cores except 1
  cl <- makeCluster(cores)
  
  # Progress bar
  if (verbose) pb <- txtProgressBar(min = 0, max = dim, style = 3)
  
  for (i in 1:dim) {
    PreSumn <- sapply(1:i, function(j) PreSum(time[j], time[i]))
    integrand2 <- parLapply(cl, 1:i, function(j) Integrand(time[j], time[i], PreSumn[j]))
    
    if (i == 1) {
      tau12[i] <- 0
    } else if (i == 2) {
      tau12[i] <- (integrand2[[1]] + integrand2[[2]]) / 2 * tstep
    } else if (i %% 2 == 1) {
      tau12[i] <- SimpsonOdd(unlist(integrand2), tstep)
    } else {
      if (length(integrand2) >= 4) {
        tau12[i] <- SimpsonOdd(unlist(integrand2)[1:(length(integrand2) - 3)], tstep) +
          SimpsonThreeEight(unlist(integrand2)[(length(integrand2) - 3):length(integrand2)], tstep)
      } else {
        tau12[i] <- SimpsonOdd(unlist(integrand2), tstep)
      }
    }
    
    # Update progress bar
    if (verbose) {
      setTxtProgressBar(pb, i)
      if (i %% 100 == 0) {
        cat(sprintf("Processed step %d of %d\n", i, dim))
      }
    }
  }
  
  # Close the progress bar and cluster
  if (verbose) close(pb)
  stopCluster(cl)
  
  # Final calculations
  etas <- 0.01
  tau12s <- etas * Amp * omega * cos(omega * time)
  tau12t <- tau12 + tau12s
  
  time2 <- time[(PointsPerPeriod * CutPeriods + 1):length(time)]
  tau12t2 <- tau12t[(PointsPerPeriod * CutPeriods + 1):length(tau12t)]
  
  strain <- Amp * sin(omega * time2)
  strainRt <- Amp * omega * cos(omega * time2)
  
  Fs <- length(time2) / (2 * pi * (NoPeriods - CutPeriods))
  T <- 1 / Fs
  L <- length(time2)
  t <- (0:(L - 1)) * T
  S <- tau12t2
  X <- S
  
  Y <- fft(X)
  
  P2 <- abs(Y / L)
  P1 <- P2[1:(L %/% 2 + 1)]
  P1[2:(length(P1) - 1)] <- 2 * P1[2:(length(P1) - 1)]
  
  aSos <- 2 * pi * Fs * (0:(L %/% 2)) / L
  harmInd <- seq(1, 61, by = 2)
  harmonics <- approx(aSos, P1, xout = harmInd)$y
  
  harmonics <- harmonics / Amp
  
  return(list(time2 = time2, strain = strain, strainRt = strainRt, tau12t2 = tau12t2, harmonics = harmonics))
}

# Example usage with verbose output to track progress

gi <- c(177464, 54980, 35272.4, 22493.7, 16039.9, 13572.9, 9004.74, 6926.17, 527.529, 2.7341)
taui <- c(0.00020788, 0.00288704, 0.0234895, 0.163525, 1.03726, 5.24526, 15.6349, 35.5925, 167.211, 4810.48)

omegaV = 1
AmpV = 0.05
aSos = 818.19
bSos = 1.21

result <- KBKZmodelSimpsonSoskeyLissajous(omegaV, AmpV, aSos, bSos, gi, taui, verbose = TRUE)

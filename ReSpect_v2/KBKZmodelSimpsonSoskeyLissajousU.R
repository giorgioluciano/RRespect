library(parallel)

KBKZmodelSimpsonSoskeyLissajous <- function(omega, Amp, aSos, bSos, gi, taui, verbose=FALSE) {
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
  
  Integrand <- function(tp, t, PreSumn) {
    PreSumn * sum(gi / taui * exp(-(t - tp) / taui))
  }
  
  SimpsonOdd <- function(integrand, h) {
    n <- length(integrand)
    if (n < 3) {
      return(sum(integrand) * h)  # Trapezoidal rule for small n
    }
    if (n == 3) {
      return((integrand[1] + 4 * integrand[2] + integrand[3]) * h / 3)
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
  
  # Rilevamento del sistema operativo
  is_windows <- .Platform$OS.type == "windows"
  num_cores <- detectCores() - 1  # Numero di core disponibili
  
  # Configurazione per il parallelismo su Windows
  if (is_windows) {
    cl <- makeCluster(num_cores)
    clusterExport(cl, c("Integrand", "PreSum", "gi", "taui", "Amp", "omega", "aSos", "bSos", "tstep", "PreSum"))
  }
  
  # Progress bar
  if (verbose) pb <- txtProgressBar(min = 0, max = dim, style = 3)
  
  for (i in 1:dim) {
    if (is_windows) {
      # Parallelismo con parLapply per Windows
      integrand2 <- parLapply(cl, 1:i, function(j) Integrand(time[j], time[i], PreSum(time[j], time[i])))
    } else {
      # Parallelismo con mclapply per Unix-like
      integrand2 <- mclapply(1:i, function(j) Integrand(time[j], time[i], PreSum(time[j], time[i])), mc.cores = num_cores)
    }
    
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
    
    if (verbose && i %% 100 == 0) {
      cat(sprintf("Processed step %d of %d\n", i, dim))
    }
  }
  
  if (verbose) close(pb)
  
  # Chiudere il cluster su Windows
  if (is_windows) stopCluster(cl)
  
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
  
  P2r <- Re(Y / L)
  P1r <- P2r[1:(L %/% 2 + 1)]
  P1r[2:(length(P1r) - 1)] <- 2 * P1r[2:(length(P1r) - 1)]
  
  harmonics <- approx(aSos, P1r, xout = harmInd)$y
  harmonicsR <- harmonics / Amp
  
  P2i <- Im(Y / L)
  P1i <- P2i[1:(L %/% 2 + 1)]
  P1i[2:(length(P1i) - 1)] <- 2 * P1i[2:(length(P1i) - 1)]
  
  harmonics <- approx(aSos, P1i, xout = harmInd)$y
  harmonicsI <- -harmonics / Amp
  
  delta <- atan(harmonicsR[1] / harmonicsI[1])
  
  ElStr <- numeric(length(harmonicsI))
  VisStr <- numeric(length(harmonicsR))
  for (i in 1:length(harmonicsI)) {
    ElStr <- ElStr + Amp * (harmonicsI[i] * sin((2 * (i - 1) + 1) * omega * time2))
    VisStr <- VisStr + Amp * (harmonicsR[i] * cos((2 * (i - 1) + 1) * omega * time2))
  }
  
  return(list(time2 = time2, strain = strain, strainRt = strainRt, 
              tau12t2 = tau12t2, ElStr = ElStr, VisStr = VisStr))
}

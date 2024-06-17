KBKZmodelSimpsonSoskeyLissajous <- function(para) {
  # KBKZ model from Song 2020
  
  # gi and taui values
  gi <- c(177464, 54980, 35272.4, 22493.7, 16039.9, 13572.9, 9004.74, 
          6926.17, 527.529, 2.7341)
  taui <- c(0.00020788, 0.00288704, 0.0234895, 0.163525, 1.03726,
            5.24526, 15.6349, 35.5925, 167.211, 4810.48)
  
  # Constants and settings
  PointsPerPeriod <- 125
  NoPeriods <- 400
  CutPeriods <- 300
  
  etai <- taui * gi
  
  # Extract parameters from input vector para
  omega <- para[1]
  Amp <- para[2]
  aSos <- para[3]
  bSos <- para[4]
  
  mint <- 0
  maxt <- 2 * pi / omega * NoPeriods
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
    (integrand[1] + 
     4 * sum(integrand[2:seq(length(integrand) - 1, by = 2)]) + 
     2 * sum(integrand[3:seq(length(integrand) - 2, by = 2)]) +
     integrand[length(integrand)]) * h / 3
  }
  
  SimpsonThreeEight <- function(integrand, h) {
    3 * h / 8 * (integrand[1] +
                 3 * integrand[2] +
                 3 * integrand[3] +
                 integrand[4])
  }
  
  # Perform integration for each time point
  for (i in seq_along(tau12)) {
    integrand2 <- numeric(i)
    time2 <- seq(mint, time[i], by = tstep)
    PreSumn <- PreSum(time2, time[i])
    for (j in seq_along(integrand2)) {
      integrand2[j] <- Integrand(time[j], time[i], PreSumn[j])
    }
    
    if (i == 1) {
      tau12[i] <- 0
    } else if (i == 2) {
      tau12[i] <- (integrand2[1] + integrand2[2]) / 2 * tstep
    } else if (i %% 2 == 1) {
      tau12[i] <- SimpsonOdd(integrand2, tstep)
    } else {
      tau12[i] <- SimpsonOdd(integrand2[1:(length(integrand2) - 3)], tstep) +
                  SimpsonThreeEight(integrand2[(length(integrand2) - 3):length(integrand2)], tstep)
    }
  }
  
  # Adjust tau12 and time for CutPeriods
  time2 <- time[(PointsPerPeriod * CutPeriods + 1):length(time)]
  tau12s <- 0.01 * Amp * omega * cos(omega * time)
  tau12t <- tau12 + tau12s
  tau12t2 <- tau12t[(PointsPerPeriod * CutPeriods + 1):length(tau12t)]
  
  # Compute strain and strainRt
  strain <- Amp * sin(omega * time2)
  strainRt <- Amp * omega * cos(omega * time2)
  
  Fs <- length(time2) / (2 * pi * (NoPeriods - CutPeriods))
  L <- length(time2)
  T <- 1 / Fs
  
  # Fourier transform
  S <- tau12t2
  X <- S
  Y <- fft(X)
  
  P2 <- abs(Y / L)
  P1 <- P2[1:(L/2 + 1)]
  P1[2:(length(P1) - 1)] <- 2 * P1[2:(length(P1) - 1)]
  
  aSos <- 2 * pi * Fs * (0:(L/2)) / L
  
  # Calculate harmonics
  harmInd <- seq(1, 61, by = 2)
  harmonics <- interp1(aSos, P1, harmInd) / Amp
  
  P2r <- Re(Y) / L
  P1r <- P2r[1:(L/2 + 1)]
  P1r[2:(length(P1r) - 1)] <- 2 * P1r[2:(length(P1r) - 1)]
  
  harmonicsR <- interp1(aSos, P1r, harmInd) / Amp
  
  P2i <- Im(Y) / L
  P1i <- P2i[1:(L/2 + 1)]
  P1i[2:(length(P1i) - 1)] <- 2 * P1i[2:(length(P1i) - 1)]
  
  harmonicsI <- -interp1(aSos, P1i, harmInd) / Amp
  
  delta <- atan(harmonicsR[1] / harmonicsI[1])
  
  # Adjust time2 and tau12t2
  time2 <- time[(PointsPerPeriod * (NoPeriods - 1) + 1):length(time)]
  tau12t2 <- tau12t[(PointsPerPeriod * NoPeriods + 1):length(tau12t)]
  
  # Compute ElStr and VisStr
  ElStr <- sum(Amp * harmonicsI * sin((2 * (seq_along(harmonicsI) - 1) + 1) * omega * time2))
  VisStr <- sum(Amp * harmonicsR * cos((2 * (seq_along(harmonicsR) - 1) + 1) * omega * time2))
  
  return(list(time2 = time2, strain = strain, strainRt = strainRt, tau12t2 = tau12t2, ElStr = ElStr, VisStr = VisStr))
}

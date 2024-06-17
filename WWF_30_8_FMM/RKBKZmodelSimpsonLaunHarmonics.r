KBKZmodelSimpsonLaunHarmonics <- function(para) {
  # KBKZ model from Song 2020
  
  # Call taugcal function to get taui and gi
  taug_result <- taugcal(0.256, 9963/7537*4775.94, 3.162, 0.001, 15)
  taui <- taug_result$taui
  gi <- taug_result$gi
  
  # Constants and settings
  PointsPerPeriod <- 125
  NoPeriods <- 200
  CutPeriods <- 100
  
  etai <- taui * gi
  
  # Extract parameters from input vector para
  omega <- para[1]
  Amp <- para[2]
  f <- para[3]
  n1 <- para[4]
  n2 <- para[5]
  
  mint <- 0
  maxt <- 2 * pi / omega * NoPeriods
  tspan <- c(mint, maxt)
  tstep <- 2 * pi / omega / PointsPerPeriod
  
  time <- seq(mint, maxt, by = tstep)
  dim <- length(time)
  tau12 <- numeric(dim)
  
  PreSum <- function(tp, t) {
    Amp * (sin(omega * t) - sin(omega * tp)) *
      (f * exp(-n1 * abs(Amp * (sin(omega * t) - sin(omega * tp)))) +
       (1 - f) * exp(-n2 * abs(Amp * (sin(omega * t) - sin(omega * tp)))))
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
  tau12t <- tau12[(PointsPerPeriod * CutPeriods + 1):length(tau12)]
  
  # Additional calculations
  etas <- 0.01
  tau12s <- etas * Amp * omega * cos(omega * time)
  tau12t2 <- tau12t + tau12s
  
  # Calculate sampling frequency Fs and length L
  Fs <- length(time2) / (2 * pi * (NoPeriods - CutPeriods))
  L <- length(time2)
  T <- 1 / Fs
  t <- (0:(L - 1)) * T
  
  # Fourier transform
  S <- tau12t2
  Y <- fft(S)
  
  P2 <- abs(Y / L)
  P1 <- P2[1:(L/2 + 1)]
  P1[2:(length(P1) - 1)] <- 2 * P1[2:(length(P1) - 1)]
  
  f <- 2 * pi * Fs * (0:(L/2)) / L
  
  # Calculate harmonicsI and harmonicsR
  harmInd <- seq(1, 61, by = 2)
  harmonicsI <- -interp1(f, P1, harmInd) / Amp
  harmonicsR <- interp1(f, Re(Y) / L, harmInd) / Amp
  
  # Return harmonicsI and harmonicsR
  return(list(harmonicsI = harmonicsI, harmonicsR = harmonicsR))
}

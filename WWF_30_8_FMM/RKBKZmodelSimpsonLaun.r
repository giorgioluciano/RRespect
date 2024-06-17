KBKZmodelSimpsonLaun <- function(para) {
  # Extract parameters from input vector para
  omega <- para[1]
  Amp <- para[2]
  f <- para[3]
  n1 <- para[4]
  n2 <- para[5]
  
  # Call taugcal function to get taui and gi
  taug_result <- taugcal(0.256, 9963/7537*4775.94, 3.162, 0.001, 15)
  taui <- taug_result$taui
  gi <- taug_result$gi
  
  # Constants and settings
  PointsPerPeriod <- 40
  NoPeriods <- 200
  CutPeriods <- 100
  etas <- 0.01
  
  etai <- taui * gi
  tau12 <- numeric(length = PointsPerPeriod * NoPeriods)
  
  # Functions for integration using Simpson's rule
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
    time <- seq(0, 2 * pi / omega * NoPeriods, length.out = PointsPerPeriod * NoPeriods)
    PreSumn <- PreSum(time[1:i], time[i])
    for (j in seq_along(integrand2)) {
      integrand2[j] <- Integrand(time[j], time[i], PreSumn[j])
    }
    
    if (i == 1) {
      tau12[i] <- 0
    } else if (i == 2) {
      tau12[i] <- (integrand2[1] + integrand2[2]) / 2 * (time[2] - time[1])
    } else if (i %% 2 == 1) {
      tau12[i] <- SimpsonOdd(integrand2, time[2] - time[1])
    } else {
      tau12[i] <- SimpsonOdd(integrand2[1:(length(integrand2) - 3)], time[2] - time[1]) +
                  SimpsonThreeEight(integrand2[(length(integrand2) - 3):length(integrand2)], time[2] - time[1])
    }
  }
  
  # Adjust tau12 and time for CutPeriods
  time2 <- time[(PointsPerPeriod * CutPeriods + 1):length(time)]
  tau12t <- tau12[(PointsPerPeriod * CutPeriods + 1):length(tau12)]
  
  # Calculate sampling frequency Fs and length L
  Fs <- length(time2) / (2 * pi * (NoPeriods - CutPeriods))
  L <- length(time2)
  T <- 1 / Fs
  t <- (0:(L - 1)) * T
  
  # Calculate Fourier transform
  S <- tau12t
  Y <- fft(S)
  
  P2 <- abs(Y / L)
  P1 <- P2[1:(L/2 + 1)]
  P1[2:(length(P1) - 1)] <- 2 * P1[2:(length(P1) - 1)]
  
  f <- Fs * (0:(L/2)) / L
  
  # Calculate harmonicsR and harmonicsI
  harmonicsR <- c(interp1(f, P1, 1.0), interp1(f, P1, 3.0), interp1(f, P1, 5.0)) / Amp
  harmonicsI <- c(interp1(f, Re(Y) / L, 1.0), interp1(f, Re(Y) / L, 3.0), interp1(f, Re(Y) / L, 5.0)) / Amp
  
  # Return harmonicsR and harmonicsI
  return(list(harmonicsR = harmonicsR, harmonicsI = harmonicsI))
}

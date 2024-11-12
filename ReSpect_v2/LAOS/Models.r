# Main function to choose and run the KBKZ model
run_KBKZ_model <- function(model_type, para) {
  # Common parameters
  PointsPerPeriod <- 125
  NoPeriods <- 200
  CutPeriods <- 100
  etas <- 0.01

  # Choose model based on input
  if (model_type == "Laun") {
    result <- KBKZ_Laun(para, PointsPerPeriod, NoPeriods, CutPeriods, etas)
  } else if (model_type == "Lissajous") {
    result <- KBKZ_Lissajous(para, PointsPerPeriod, NoPeriods, CutPeriods, etas)
  } else if (model_type == "SoskeyLissajous") {
    result <- KBKZ_SoskeyLissajous(para, PointsPerPeriod, NoPeriods, CutPeriods, etas)
  } else {
    stop("Invalid model type")
  }
  
  return(result)
}

# Helper functions

taugcal <- function(p, Gp1, r, taumin, N) {
  taui <- numeric(N)
  H <- numeric(N)
  gi <- numeric(N)
  
  taui[1] <- taumin
  for (i in 1:(N-1)) {
    taui[i+1] <- taui[i] * r
  }
  
  G1 <- 2 * Gp1 * gamma(1 + p) / p / pi * sin(p * pi / 2)
  for (i in 1:N) {
    H[i] <- G1 * taui[i]^(-p) / gamma(p)
    gi[i] <- H[i] * (r^(p/2) - r^(-p/2)) / p
  }
  
  return(list(taui = taui, gi = gi))
}

calculate_stress <- function(time, PreSum, Integrand, gi, taui, tstep) {
  dim <- length(time)
  tau12 <- numeric(dim)
  
  for (i in 1:dim) {
    time2 <- seq(time[1], time[i], by = tstep)
    PreSumn <- PreSum(time2, time[i])
    integrand2 <- sapply(time2, function(tp) Integrand(tp, time[i], PreSumn[which(time2 == tp)], gi, taui))
    
    if (i == 1) {
      tau12[i] <- 0
    } else if (i == 2) {
      tau12[i] <- (integrand2[1] + integrand2[2])/2 * tstep
    } else if (i %% 2 == 1) {
      tau12[i] <- SimpsonOdd(integrand2, tstep)
    } else {
      tau12[i] <- SimpsonOdd(integrand2[1:(length(integrand2)-3)], tstep) + 
        SimpsonThreeEight(integrand2[(length(integrand2)-3):length(integrand2)], tstep)
    }
  }
  
  return(tau12)
}

calculate_harmonics <- function(tau12t2, Fs, Amp) {
  L <- length(tau12t2)
  Y <- fft(tau12t2)
  
  P2 <- abs(Y/L)
  P1 <- P2[1:(L/2+1)]
  P1[2:(length(P1)-1)] <- 2 * P1[2:(length(P1)-1)]
  
  f <- 2*pi*Fs*(0:(L/2))/L
  harmInd <- seq(1, 61, by=2)
  
  harmonics <- approx(f, P1, xout=harmInd)$y / Amp
  
  P2r <- Re(Y/L)
  P1r <- P2r[1:(L/2+1)]
  P1r[2:(length(P1r)-1)] <- 2 * P1r[2:(length(P1r)-1)]
  harmonicsR <- approx(f, P1r, xout=harmInd)$y / Amp
  
  P2i <- Im(Y/L)
  P1i <- P2i[1:(L/2+1)]
  P1i[2:(length(P1i)-1)] <- 2 * P1i[2:(length(P1i)-1)]
  harmonicsI <- -approx(f, P1i, xout=harmInd)$y / Amp
  
  return(list(harmonicsR = harmonicsR, harmonicsI = harmonicsI))
}

# Model-specific functions

KBKZ_Laun <- function(para, PointsPerPeriod, NoPeriods, CutPeriods, etas) {
  result <- taugcal(0.256, 9963/7537*4775.94, 3.162, 0.001, 15)
  taui <- result$taui
  gi <- result$gi
  
  omega <- para[1]
  Amp <- para[2]
  f <- para[3]
  n1 <- para[4]
  n2 <- para[5]
  
  mint <- 0
  maxt <- 2*pi/omega*NoPeriods
  tstep <- 2*pi/omega/PointsPerPeriod
  time <- seq(mint, maxt, by=tstep)
  
  PreSum <- function(tp, t) {
    Amp * (sin(omega*t) - sin(omega*tp)) * 
      (f*exp(-n1*abs(Amp*(sin(omega*t) - sin(omega*tp)))) + 
       (1-f)*exp(-n2*abs(Amp*(sin(omega*t) - sin(omega*tp)))))
  }
  
  Integrand <- function(tp, t, PreSumn, gi, taui) {
    PreSumn * sum(gi/taui * exp(-(t-tp)/taui))
  }
  
  tau12 <- calculate_stress(time, PreSum, Integrand, gi, taui, tstep)
  
  tau12s <- etas*Amp*omega*cos(omega*time)
  tau12t <- tau12 + tau12s
  
  time2 <- time[(PointsPerPeriod*CutPeriods+1):length(time)]
  tau12t2 <- tau12t[(PointsPerPeriod*CutPeriods+1):length(tau12t)]
  
  Fs <- length(time2) / (2*pi*(NoPeriods-CutPeriods))
  
  harmonics <- calculate_harmonics(tau12t2, Fs, Amp)
  
  return(harmonics)
}

KBKZ_Lissajous <- function(para, PointsPerPeriod, NoPeriods, CutPeriods, etas) {
  # Similar to KBKZ_Laun, but with different PreSum function
  # ...
}

KBKZ_SoskeyLissajous <- function(para, PointsPerPeriod, NoPeriods, CutPeriods, etas) {
  # Similar to KBKZ_Laun, but with different gi, taui, and PreSum function
  # ...
}

# Simpson's rule integration methods
SimpsonOdd <- function(integrand, h) {
  (integrand[1] + 4*sum(integrand[seq(2, length(integrand)-1, by=2)]) + 
   2*sum(integrand[seq(3, length(integrand)-2, by=2)]) + 
   integrand[length(integrand)]) * h/3
}

SimpsonThreeEight <- function(integrand, h) {
  3*h/8 * (integrand[1] + 3*integrand[2] + 3*integrand[3] + integrand[4])
}

# Usage example
para <- c(1, 0.1, 0.6658, 10.35, 278.9)
result <- run_KBKZ_model("Laun", para)
print(result)
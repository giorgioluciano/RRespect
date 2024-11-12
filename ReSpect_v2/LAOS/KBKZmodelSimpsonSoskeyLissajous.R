KBKZmodelSimpsonSoskeyLissajous <- function(omega, Amp, aSos, bSos, gi, taui) {
  # KBKZ model from Song 2020 
  
 
  PointsPerPeriod <- 125
  NoPeriods <- 200
  CutPeriods = 100;
  etai <- taui * gi
  
  mint <- 0
  maxt <- 2 * pi / omega * NoPeriods
  tspan <- c(mint, maxt)
  tstep <- 2 * pi / omega / PointsPerPeriod
  
  time <- seq(mint, maxt, by = tstep)
  
  dim <- length(time)
  
  tau12 <- numeric(dim)
  
  PreSum <- function(tp, t) {
    Amp * (sin(omega * t) - sin(omega * tp)) *
      1.0 / (1.0 + aSos * (abs(Amp * (sin(omega * t) - sin(omega * tp))))^bSos)
  }
  
  Integrand <- function(tp, t, PreSumn) {
    PreSumn * sum(gi / taui * exp(-(t - tp) / taui))
  }
  
  SimpsonOdd <- function(integrand, h) {
    n <- length(integrand)
    if (n < 3) {
      return(sum(integrand) * h)  # Use trapezoidal rule for small n
    }
    if (n == 3) {
      return((integrand[1] + 4*integrand[2] + integrand[3]) * h / 3)
    }
    (integrand[1] + 
        4 * sum(integrand[seq(2, n - 1, by = 2)]) + 
        2 * sum(integrand[seq(3, n - 2, by = 2)]) +
        integrand[n]) * h / 3
  }
  
  SimpsonThreeEight <- function(integrand, h) {
    3 * h / 8 * (integrand[1] + 3 * integrand[2] + 3 * integrand[3] + integrand[4])
  }
  
  for (i in 1:dim) {
    integrand2 <- numeric(i)
    time2 <- seq(mint, time[i], by = tstep)
    PreSumn <- PreSum(time2, time[i])
    for (j in 1:i) {
      integrand2[j] <- Integrand(time[j], time[i], PreSumn[j])
    }
    
    if (i == 1) {
      tau12[i] <- 0
    } else if (i == 2) {
      tau12[i] <- (integrand2[1] + integrand2[2]) / 2 * tstep
    } else if (i %% 2 == 1) {
      tau12[i] <- SimpsonOdd(integrand2, tstep)
    } else {
      if (length(integrand2) >= 4) {
        tau12[i] <- SimpsonOdd(integrand2[1:(length(integrand2) - 3)], tstep) +
          SimpsonThreeEight(integrand2[(length(integrand2) - 3):length(integrand2)], tstep)
      } else {
        tau12[i] <- SimpsonOdd(integrand2, tstep)
      }
    }
  }
  
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
  time2 <- time[(PointsPerPeriod * (NoPeriods - 1) + 1):length(time)]
  tau12t2 <- tau12t[(PointsPerPeriod * NoPeriods + 1):length(tau12t)]
  
  strain <- Amp * sin(omega * time2)
  strainRt <- Amp * omega * cos(omega * time2)
  
  ElStr <- 0.0
  VisStr <- 0.0
  dimHarm <- length(harmonicsI)
  for (i in 1:dimHarm) {
    ElStr <- ElStr + Amp * (harmonicsI[i] * sin((2 * (i - 1) + 1) * omega * time2))
    VisStr <- VisStr + Amp * (harmonicsR[i] * cos((2 * (i - 1) + 1) * omega * time2))
  }
  
  return(list(time2 = time2, strain = strain, strainRt = strainRt,
              tau12t2 = tau12t2, ElStr = ElStr, VisStr = VisStr))
}

# Example of how to use the function:
gi <- c(177464, 54980, 35272.4, 22493.7, 16039.9, 13572.9, 9004.74, 6926.17, 527.529, 2.7341)
taui <- c(0.00020788, 0.00288704, 0.0234895, 0.163525, 1.03726, 5.24526, 15.6349, 35.5925, 167.211, 4810.48)

# Example values for the parameters, replace with actual values
#omega = c(1, 10,50, 100)
#Amp= c(0.05, 0.25, 0.5, 1.0, 2.0)


omegaV = 1
AmpV= 0.05

aSos = 818.19
bSos = 1.21

result <- KBKZmodelSimpsonSoskeyLissajous(omega, Amp, aSos, bSos, gi, taui)

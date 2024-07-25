library(minpack.lm)

# Dati sperimentali (esempio)
shear_rate <- c(0.1, 0.5, 1, 5, 10, 50, 100)
stress <- c(10, 30, 50, 100, 150, 300, 400)

# Funzione del modello Winter-Soskey
soskey_model <- function(shear_rate, eta0, alpha, n) {
  eta0 * shear_rate / (1 + (alpha * abs(shear_rate))^n)
}

# Fit del modello
fit <- nlsLM(stress ~ soskey_model(shear_rate, eta0, alpha, n),
             start = list(eta0 = 100, alpha = 0.1, n = 0.5),
             data = data.frame(shear_rate, stress))

# Risultati
summary(fit)
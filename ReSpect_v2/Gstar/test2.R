library(minpack.lm)

# Funzione del modello
soskey_model <- function(gamma, a, b) {
  1 / (1 + (a * gamma)^b)
}


gamma <- dataset1$V1
h_experimental <- dataset1$V2

# Creazione del dataframe
data <- data.frame(gamma = gamma, h = h_experimental)

# Fit del modello
fit <- nlsLM(h ~ soskey_model(gamma, a, b),
             data = data,
             start = list(a = 0.01, b = 2.5))

# Riepilogo dei risultati
summary(fit)

# Estrazione dei parametri stimati
a_fitted <- coef(fit)["a"]
b_fitted <- coef(fit)["b"]

# Creazione di un grafico per visualizzare il fit
gamma_range <- seq(min(gamma), max(gamma), length.out = 100)
h_fitted <- soskey_model(gamma_range, a_fitted, b_fitted)

plot(gamma, h_experimental, log = "x", pch = 16, 
     xlab = "gamma", ylab = "h(gamma)", 
     main = "Fit del modello di Soskey")
lines(gamma_range, h_fitted, col = "red")
legend("topright", legend = c("Dati sperimentali", "Modello fittato"), 
       pch = c(16, NA), lty = c(NA, 1), col = c("black", "red"))

# Stampa dei parametri fittati
cat("Parametri fittati:\n")
cat("a =", a_fitted, "\n")
cat("b =", b_fitted, "\n")
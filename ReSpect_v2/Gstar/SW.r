library(minpack.lm)



soskey_model <- function(w, a, b) {
  1 / (1 + (a * abs(w)^b))
}

data <- data.frame(w = w, h = h)

fit <- nlsLM(h ~ soskey_model(w, a, b),
             data = data,
             start = list(a = 0.9, b = 0.5))

summary(fit)

a_fitted <- coef(fit)["a"]
b_fitted <- coef(fit)["b"]


w_range <- seq(min(w), max(w), length.out = 100)
h_fitted <- soskey_model(w_range, a_fitted, b_fitted)

plot(w, h, log = "x", pch = 16, 
     xlab = "w", ylab = "h(w)", 
     main = "Fit del modello di Soskey")
lines(w_range, h_fitted, col = "red")
legend("topright", legend = c("Dati sperimentali", "Modello fittato"), 
       pch = c(16, NA), lty = c(NA, 1), col = c("black", "red"))

# Stampa dei parametri fittati
cat("Parametri fittati:\n")
cat("a =", a_fitted, "\n")
cat("b =", b_fitted, "\n")
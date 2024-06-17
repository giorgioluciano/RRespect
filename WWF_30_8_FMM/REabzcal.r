Eabzcal <- function(a, b, z, k) {
  Eabz <- z^k / gamma(a * k + b)
  return(Eabz)
}
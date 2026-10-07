make_maxwell_data <- function(
    omega,
    G = c(1500, 400, 80),
    tau = c(0.005, 0.1, 2),
    noise_relative = 0,
    seed = 123
) {
  stopifnot(
    is.numeric(omega),
    all(is.finite(omega)),
    all(omega > 0),
    is.numeric(G),
    is.numeric(tau),
    length(G) == length(tau),
    all(is.finite(G)),
    all(is.finite(tau)),
    all(G > 0),
    all(tau > 0),
    length(noise_relative) == 1L,
    is.finite(noise_relative),
    noise_relative >= 0
  )

  wt <- outer(omega, tau, "*")

  Gp <- rowSums(
    sweep(
      wt^2 / (1 + wt^2),
      2,
      G,
      "*"
    )
  )

  Gpp <- rowSums(
    sweep(
      wt / (1 + wt^2),
      2,
      G,
      "*"
    )
  )

  if (noise_relative > 0) {
    set.seed(seed)

    Gp <- Gp * (
      1 + rnorm(
        length(Gp),
        mean = 0,
        sd = noise_relative
      )
    )

    Gpp <- Gpp * (
      1 + rnorm(
        length(Gpp),
        mean = 0,
        sd = noise_relative
      )
    )
  }

  data.frame(
    omega = omega,
    Gp = Gp,
    Gpp = Gpp
  )
}


evaluate_noisy_maxwell <- function(seed, noise_relative = 0.01) {
  omega <- 10^seq(-4, 4, length.out = 161)
  truth_g <- c(1500, 400, 80)
  truth_tau <- c(0.005, 0.1, 2)

  truth <- make_maxwell_data(
    omega = omega,
    G = truth_g,
    tau = truth_tau
  )

  dat <- make_maxwell_data(
    omega = omega,
    G = truth_g,
    tau = truth_tau,
    noise_relative = noise_relative,
    seed = seed
  )

  input_file <- tempfile(fileext = ".dat")
  on.exit(unlink(input_file), add = TRUE)

  utils::write.table(
    dat,
    file = input_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  par <- setParams(
    domain = "frequency",
    dataFile = input_file,
    ns = 100L,
    lamC = 0,
    smFacLam = 0,
    freqEnd = 1L,
    verbose = FALSE,
    plotting = FALSE,
    lamMin = 1e-10,
    lamMax = 1e3,
    lamDensity = 3L,
    plateau = FALSE
  )

  fit <- getDiscreteSpectrum(
    par = par,
    writeOutput = FALSE,
    projectDir = tempdir()
  )

  stopifnot(
    length(fit$g) > 0L,
    length(fit$g) == length(fit$tau),
    all(is.finite(fit$g)),
    all(is.finite(fit$tau)),
    all(fit$g > 0),
    all(fit$tau > 0)
  )

  wt <- outer(omega, fit$tau, "*")

  fitted_gp <- as.vector(
    (wt^2 / (1 + wt^2)) %*% fit$g
  )

  fitted_gpp <- as.vector(
    (wt / (1 + wt^2)) %*% fit$g
  )

  matched <- vapply(
    truth_tau,
    function(tau_true) {
      which.min(abs(log(fit$tau / tau_true)))
    },
    integer(1)
  )

  normalized_rmse <- function(observed, fitted) {
    sqrt(mean((observed - fitted)^2)) / max(observed)
  }

  list(
    distinct_matches = length(unique(matched)),
    gp_error_to_truth = normalized_rmse(
      truth$Gp,
      fitted_gp
    ),
    gpp_error_to_truth = normalized_rmse(
      truth$Gpp,
      fitted_gpp
    ),
    total_modulus_relative_error = abs(
      sum(fit$g) - sum(truth_g)
    ) / sum(truth_g),
    max_log_tau_error = max(
      abs(log(fit$tau[matched] / truth_tau))
    ),
    max_mode_relative_error = max(
      abs(fit$g[matched] - truth_g) / truth_g
    )
  )
}


evaluate_smooth_continuous <- function() {
  omega <- 10^seq(-4, 4, length.out = 161)

  total_g <- 1000
  center_tau <- 0.1
  sigma_log_tau <- 0.8

  density_at <- function(tau) {
    total_g * dnorm(
      log(tau),
      mean = log(center_tau),
      sd = sigma_log_tau
    )
  }

  trapezoid <- function(y, x) {
    sum(
      diff(x) *
        (head(y, -1L) + tail(y, -1L)) / 2
    )
  }

  truth_tau <- 10^seq(-6, 6, length.out = 4001)
  truth_x <- log(truth_tau)
  truth_density <- density_at(truth_tau)

  gp <- vapply(
    omega,
    function(w) {
      z <- w * truth_tau

      trapezoid(
        truth_density * z^2 / (1 + z^2),
        truth_x
      )
    },
    numeric(1)
  )

  gpp <- vapply(
    omega,
    function(w) {
      z <- w * truth_tau

      trapezoid(
        truth_density * z / (1 + z^2),
        truth_x
      )
    },
    numeric(1)
  )

  input_file <- tempfile(fileext = ".dat")
  on.exit(unlink(input_file), add = TRUE)

  utils::write.table(
    data.frame(
      omega = omega,
      Gp = gp,
      Gpp = gpp
    ),
    file = input_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  par <- setParams(
    domain = "frequency",
    dataFile = input_file,
    ns = 100L,
    lamC = 0,
    smFacLam = 0,
    freqEnd = 1L,
    verbose = FALSE,
    plotting = FALSE,
    lamMin = 1e-10,
    lamMax = 1e3,
    lamDensity = 3L,
    plateau = FALSE
  )

  fit <- getDiscreteSpectrum(
    par = par,
    writeOutput = FALSE,
    projectDir = tempdir()
  )

  cont <- fit$continuous

  stopifnot(
    all(c("s", "H", "w", "Gexp") %in% names(cont)),
    length(cont$s) == length(cont$H),
    all(is.finite(cont$s)),
    all(cont$s > 0),
    all(diff(cont$s) > 0),
    all(is.finite(cont$H))
  )

  x <- log(cont$s)
  estimated_density <- exp(cont$H)
  reference_density <- density_at(cont$s)

  stopifnot(
    all(is.finite(estimated_density)),
    all(estimated_density >= 0)
  )

  fitted_gp <- vapply(
    cont$w,
    function(w) {
      z <- w * cont$s

      trapezoid(
        estimated_density * z^2 / (1 + z^2),
        x
      )
    },
    numeric(1)
  )

  fitted_gpp <- vapply(
    cont$w,
    function(w) {
      z <- w * cont$s

      trapezoid(
        estimated_density * z / (1 + z^2),
        x
      )
    },
    numeric(1)
  )

  n <- length(cont$w)
  observed_gp <- cont$Gexp[seq_len(n)]
  observed_gpp <- cont$Gexp[n + seq_len(n)]

  estimated_total <- trapezoid(
    estimated_density,
    x
  )

  reference_total <- trapezoid(
    reference_density,
    x
  )

  stopifnot(
    is.finite(estimated_total),
    estimated_total > 0,
    is.finite(reference_total),
    reference_total > 0
  )

  estimated_center <- trapezoid(
    x * estimated_density,
    x
  ) / estimated_total

  reference_center <- trapezoid(
    x * reference_density,
    x
  ) / reference_total

  estimated_width <- sqrt(
    trapezoid(
      (x - estimated_center)^2 * estimated_density,
      x
    ) / estimated_total
  )

  normalized_rmse <- function(observed, fitted) {
    sqrt(mean((observed - fitted)^2)) / max(observed)
  }

  data.frame(
    gp_normalized_rmse = normalized_rmse(
      observed_gp,
      fitted_gp
    ),
    gpp_normalized_rmse = normalized_rmse(
      observed_gpp,
      fitted_gpp
    ),
    gp_relative_rms = sqrt(mean(
      (fitted_gp / observed_gp - 1)^2
    )),
    gpp_relative_rms = sqrt(mean(
      (fitted_gpp / observed_gpp - 1)^2
    )),
    density_relative_l1 = trapezoid(
      abs(estimated_density - reference_density),
      x
    ) / reference_total,
    total_modulus_relative_error = abs(
      estimated_total - total_g
    ) / total_g,
    center_log_error = abs(
      estimated_center - reference_center
    ),
    estimated_log_width = estimated_width,
    true_log_width = sigma_log_tau,
    truth_mass_on_fit_grid = reference_total / total_g
  )
}

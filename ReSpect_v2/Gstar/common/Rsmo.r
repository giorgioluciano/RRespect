supsmooth <- function(x, y, ...) {
  # supsmooth: Smoothing of scatterplots using Friedman's supersmoother algorithm.
  #
  # Syntax:
  #   Y_SMOOTH = supsmooth(X, Y, 'PropertyName', PropertyValue, ...)
  #
  # Inputs:
  #   X, Y are same-length vectors.
  #
  # Output:
  #   Y_SMOOTH is a smoothed version of Y.
  #
  # Example,
  #   x <- seq(0, 1, length.out = 201)
  #   y <- sin(2.5 * x) + 0.05 * rnorm(201)
  #   smo <- supsmooth(x, y)
  #   plot(x, y, col = "red", pch = 16)
  #   lines(x, smo, col = "blue")
  #
  # The supersmoother algorithm computes three separate smooth curves from
  # the input data with symmetric spans of 0.05*n, 0.2*n and 0.5*n, where n
  # is the number of data points. The best of the three smooth curves is
  # chosen for each predicted point using leave-one-out cross validation. The
  # best spans are then smoothed by a fixed-span smoother (span = 0.2*n) and
  # the prediction is computed by linearly interpolating between the three
  # smooth curves. This final smooth curve is then smoothed again with a
  # fixed-span smoother (span = 0.05*n).
  #
  # According to comments by Friedman, "For small samples (n < 40) or if
  # there are substantial serial correlations between observations close in
  # x-value, then a prespecified fixed span smoother (span > 0) should be
  # used. Reasonable span values are 0.2 to 0.4."
  #
  #
  # The following optional property/value pairs can be specified as arguments
  # to control the indicated behavior:
  #
  #   Property    Value
  #   ----------  ----------------------------------------------------------
  #   Weights     Vector of relative weights of each data point. Default is
  #               for all points to be weighted equally.
  #
  #   Span        Sets the width of a fixed-width smoothing operation
  #               relative to the number of data points, 0 < SPAN < 1.
  #               Setting this to be non-zero disables the supersmoother
  #               algorithm. Default is 0 (use supersmoother).
  #
  #   Period      Sets the period of periodic data. Default is Inf
  #               (infinity) which implies that the data is not periodic.
  #               Can also be set to zero for the same effect.
  #
  #   Alpha       Sets a small-span penalty to produce a greater smoothing
  #               effect. 0 < Alpha < 10, where 0 does nothing and 10
  #               produces the maximum effect. Default = 0.
  #
  #   Unsorted    If the data points are not already sorted in order of the X
  #               values then setting this to true will sort them for you.
  #               Default = false.
  #
  # All properties names are case-insensitive and need only be unambiguous.
  # For example,
  #
  #   Y_SMOOTH = supsmooth(X, Y, 'weights', rep(1, n), 'per', 2 * pi)
  #
  # is valid usage.

  # Friedman, J. H. (1984). A Variable Span Smoother. Tech. Rep. No. 5,
  # Laboratory for Computational Statistics, Dept. of Statistics, Stanford
  # Univ., California.

  # Version: 1.0, 12 December 2007
  # Author:  Douglas M. Schwarz
  # Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
  # Real_email = regexprep(Email, {'=', '*'}, {'@', '.'})

  # Load necessary library
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required.")
  }

  # Input checks.
  if (length(list(...)) %% 2 != 0) {
    stop("Properties must be specified by property/value pairs.")
  }

  # x and y must be vectors with same number of points (at least 5).
  if (!is.vector(x) || !is.vector(y) || length(x) != length(y) || length(y) < 5) {
    stop("X and Y must be equal-length vectors of at least 5 points.")
  }

  # Default property values
  prop <- list(
    weights = NULL,
    span = 0,
    period = Inf,
    alpha = 0,
    unsorted = FALSE
  )

  # Process inputs and set prop fields
  args <- list(...)
  properties <- names(prop)
  arg_index <- 1

  while (arg_index <= length(args)) {
    arg <- args[[arg_index]]
    if (is.character(arg)) {
      prop_index <- which(tolower(properties) == tolower(arg))
      if (length(prop_index) == 1) {
        prop[[properties[prop_index]]] <- args[[arg_index + 1]]
      } else {
        stop(sprintf("Property '%s' does not exist or is ambiguous.", arg))
      }
      arg_index <- arg_index + 2
    } else if (is.list(arg)) {
      for (i in seq_along(arg)) {
        prop_index <- which(tolower(properties) == tolower(names(arg)[i]))
        if (length(prop_index) == 1) {
          prop[[properties[prop_index]]] <- arg[[i]]
        } else {
          stop(sprintf("Property '%s' does not exist or is ambiguous.", names(arg)[i]))
        }
      }
      arg_index <- arg_index + 1
    } else {
      stop("Properties must be specified by property/value pairs or structures.")
    }
  }

  # Validate properties
  if (!is.null(prop$weights) && length(prop$weights) != length(x)) {
    stop("Weights property must be a vector of the same length as X and Y.")
  }

  if (!is.numeric(prop$span) || prop$span < 0 || prop$span >= 1) {
    stop("Span property must be a numeric scalar, 0 <= span < 1.")
  }

  if (!is.numeric(prop$period) || prop$period < 0) {
    stop("Period property must be a numeric scalar >= 0.")
  }

  if (!is.numeric(prop$alpha) || prop$alpha < 0 || prop$alpha > 10) {
    stop("Alpha property must be a numeric scalar between 0 and 10.")
  }

  if (!is.logical(prop$unsorted)) {
    stop("Unsorted property must be a logical scalar.")
  }

  # Function to smooth data
  smooth <- function(x, y, weights = NULL, span, period) {
    if (!is.null(weights)) {
      fit <- loess(y ~ x, weights = weights, span = span)
    } else {
      fit <- loess(y ~ x, span = span)
    }
    return(predict(fit, x))
  }

  # Make x and y into column vectors and sort if necessary
  if (prop$unsorted) {
    order <- order(x)
    x <- x[order]
    y <- y[order]
    if (!is.null(prop$weights)) {
      prop$weights <- prop$weights[order]
    }
  }

  # If prop$span > 0 then we have a fixed span smooth
  if (prop$span > 0) {
    smo <- smooth(x, y, prop$weights, prop$span, prop$period)
    return(smo)
  }

  spans <- c(0.05, 0.2, 0.5)
  nspans <- length(spans)
  n <- length(y)

  # Compute three smooth curves
  smo_n <- matrix(0, n, nspans)
  acvr_smo <- matrix(0, n, nspans)
  for (i in seq_along(spans)) {
    smo_n[, i] <- smooth(x, y, prop$weights, spans[i], prop$period)
    abs_cv_res <- abs(y - smo_n[, i])
    acvr_smo[, i] <- smooth(x, abs_cv_res, prop$weights, spans[2], prop$period)
  }

  # Select which smooth curve has smallest error using cross validation
  resmin <- apply(acvr_smo, 1, min)
  index <- apply(acvr_smo, 1, which.min)
  span_cv <- spans[index]

  # Apply alpha
  if (prop$alpha != 0) {
    small <- 1e-7
    tf <- resmin < acvr_smo[, 3] & resmin > 0
    span_cv[tf] <- span_cv[tf] + (spans[3] - span_cv[tf]) *
      (1 - exp(-prop$alpha * (log(acvr_smo[tf, 3] + small) - log(resmin[tf] + small))))
    span_cv[span_cv > spans[3]] <- spans[3]
  }

  # Smooth best spans
  smooth_spans <- smooth(x, span_cv, span = spans[2], period = prop$period)

  # Linearly interpolate three smooth curves
  smo_final <- numeric(n)
  for (i in 2:nspans) {
    tf <- smooth_spans >= spans[i - 1] & smooth_spans <= spans[i]
    weight <- (smooth_spans[tf] - spans[i - 1]) / (spans[i] - spans[i - 1])
    smo_final[tf] <- (1 - weight) * smo_n[tf, i - 1] + weight * smo_n[tf, i]
  }

  # Smooth result with fixed span smoother
  smo_final <- smooth(x, smo_final, span = spans[1], period = prop$period)
  return(smo_final)
}

# Example usage
x <- seq(0, 1, length.out = 201)
y <- sin(2.5 * x) + 0.05 * rnorm(201)
smo <- supsmooth(x, y)
plot(x, y, col = "red", pch = 16)
lines(x, smo, col = "blue")

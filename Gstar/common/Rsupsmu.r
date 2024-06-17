supsmu <- function(x, y, ...) {
  # Input checks
  if (length(x) != length(y) || length(x) < 5) {
    stop("X and Y must be vectors of equal length with at least 5 points.")
  }
  
  # Define default properties
  prop <- list(
    weights = NULL,
    span = 0,
    period = Inf,
    alpha = 0,
    unsorted = FALSE
  )
  
  # Process optional arguments
  args <- list(...)
  valid_props <- c("weights", "span", "period", "alpha", "unsorted")
  for (arg_name in names(args)) {
    prop_index <- match(tolower(arg_name), tolower(valid_props))
    if (is.na(prop_index)) {
      stop(sprintf("Property '%s' does not exist or is ambiguous.", arg_name))
    }
    prop[[valid_props[prop_index]]] <- args[[arg_name]]
  }
  
  # Validate weights property
  if (!is.null(prop$weights)) {
    if (length(prop$weights) != length(y)) {
      stop("Weights must be a vector of the same length as X and Y.")
    }
  }
  
  # Validate span property
  if (!is.numeric(prop$span) || length(prop$span) != 1 || prop$span < 0 || prop$span >= 1) {
    stop("Span must be a numeric scalar with 0 <= span < 1.")
  }
  
  # Validate period property
  if (!is.numeric(prop$period) || length(prop$period) != 1 || prop$period < 0) {
    stop("Period must be a numeric scalar >= 0.")
  }
  if (is.infinite(prop$period)) {
    prop$period <- 0
  }
  
  # Validate alpha property
  if (!is.numeric(prop$alpha) || length(prop$alpha) != 1 || prop$alpha < 0 || prop$alpha > 10) {
    stop("Alpha must be a numeric scalar with 0 <= alpha <= 10.")
  }
  
  # Validate unsorted property
  if (!is.logical(prop$unsorted) || length(prop$unsorted) != 1) {
    stop("Unsorted must be a logical scalar.")
  }
  
  # Select smoothing function based on properties
  if (is.null(prop$weights) && prop$period != 0) {
    smooth_fcn <- smooth_aper
  } else if (!is.null(prop$weights) && prop$period != 0) {
    smooth_fcn <- smooth_wt_aper
  } else if (is.null(prop$weights) && prop$period == 0) {
    smooth_fcn <- smooth_per
  } else {
    smooth_fcn <- smooth_wt_per
  }
  
  # Sort x and y if unsorted is true
  if (prop$unsorted) {
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    if (!is.null(prop$weights)) {
      prop$weights <- prop$weights[ord]
    }
  }
  
  # If span > 0, apply fixed span smoothing
  if (prop$span > 0) {
    smo <- smooth_fcn(x, y, prop$weights, prop$span, prop$period)
    return(smo)
  }
  
  # Define spans for supersmoother
  spans <- c(0.05, 0.2, 0.5)
  nspans <- length(spans)
  
  # Compute three smooth curves
  smo_n <- matrix(0, nrow = length(x), ncol = nspans)
  acvr_smo <- matrix(0, nrow = length(x), ncol = nspans)
  for (i in 1:nspans) {
    res <- smooth_fcn(x, y, prop$weights, spans[i], prop$period)
    smo_n[, i] <- res[[1]]
    acvr_smo[, i] <- res[[2]]
  }
  
  # Select best smooth curve using cross-validation
  resmin <- apply(acvr_smo, 1, min)
  index <- apply(acvr_smo, 1, which.min)
  span_cv <- spans[index]
  
  # Apply alpha
  if (prop$alpha != 0) {
    small <- 1e-7
    tf <- resmin < acvr_smo[, 3] & resmin > 0
    span_cv[tf] <- span_cv[tf] + (spans[3] - span_cv[tf]) *
      max(small, resmin[tf] / acvr_smo[tf, 3]) ^ (10 - prop$alpha)
  }
  
  # Smooth span_cv and clip at spans(1) and spans(end)
  smo_span <- smooth_fcn(x, span_cv, prop$weights, spans[2], prop$period)[[1]]
  smo_span <- pmin(pmax(smo_span, spans[1]), spans[length(spans)])
  
  # Interpolate each point
  smo_raw <- matrix(0, nrow = length(x), ncol = 1)
  for (i in 1:length(x)) {
    smo_raw[i] <- approx(spans, smo_n[i,], smo_span[i])$y
  }
  
  # Apply final smooth
  smo <- smooth_fcn(x, smo_raw, prop$weights, spans[1], prop$period)[[1]]
  
  # Reshape smo to original dimensions if unsorted
  if (prop$unsorted) {
    smo <- smo[order]
  }
  
  return(smo)
}

# Subfunctions
smooth_wt_aper <- function(x, y, w, span, period) {
  # Function logic for weighted aperiodic smoothing
  # This function should return a list containing smoothed values and acvr
  # Implement according to original MATLAB smooth_wt_aper function
}

smooth_wt_per <- function(x, y, w, span, period) {
  # Function logic for weighted periodic smoothing
  # This function should return a list containing smoothed values and acvr
  # Implement according to original MATLAB smooth_wt_per function
}

smooth_aper <- function(x, y, w, span, period) {
  # Function logic for unweighted aperiodic smoothing
  # This function should return a list containing smoothed values and acvr
  # Implement according to original MATLAB smooth_aper function
}

smooth_per <- function(x, y, w, span, period) {
  # Function logic for unweighted periodic smoothing
  # This function should return a list containing smoothed values and acvr
  # Implement according to original MATLAB smooth_per function
}


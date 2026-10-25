#' Integrated kernel plug-in estimator of AROCS
#'
#' Estimate the surface area by integrating three kernel density estimates,
#' rather than averaging estimated density ratios over observed pairs.
#' @param X1,X2,X3 Nonempty finite numeric vectors in ordinal class order.
#' @param kernel_type Either \code{"gauss"} or \code{"epan"}.
#' @param bwd One positive bandwidth, or three bandwidths in class order.
#' @param na.rm Remove missing observations within each class.
#' @param n_grid Initial number of equally spaced grid points (at least 3).
#' @param max_grid Maximum grid size. Refinement doubles the number of intervals.
#' @param rel.tol,abs.tol Relative and absolute quadrature tolerances.
#' @param tail_prob Upper bound on the omitted mass of each Gaussian KDE.
#'   Ignored for Epanechnikov, whose full support is used.
#' @param details Return numerical diagnostics as a list instead of a scalar.
#' @details The estimator is
#' \deqn{\int_{x<y}\sqrt{
#' [\hat f_1(x)\hat f_3(y)]^2+
#' [\hat f_2(x)\hat f_3(y)]^2+
#' [\hat f_1(x)\hat f_2(y)]^2}\,dx\,dy.}
#' This follows by multiplying the density ratios in the population expectation
#' representation by its joint density. It is a different estimator from the
#' empirical pair average in equation (3.1), implemented by \code{arocs_kernel}.
#' The ratio-free expression remains defined when kernel densities are zero.
#'
#' Densities are evaluated by \code{pdf_kernel}. Tensor trapezoidal quadrature
#' assigns half weight to the diagonal. Grids are refined until successive
#' estimates differ by at most \code{max(abs.tol, rel.tol * abs(estimate))}.
#' This difference is a convergence diagnostic, not a rigorous error bound;
#' very narrow kernels may require a larger initial grid. Failure to converge
#' raises a warning. The Gaussian truncation error is bounded by
#' \code{6 * tail_prob}, using the sum of the three density products as an
#' envelope. Numerical error is distinct from sampling error and smoothing bias.
#' No statistical confidence interval or bias correction is implied.
#'
#' Use the same bandwidth selection procedure in every bootstrap sample when
#' bandwidths are data dependent. Check smoothing sensitivity and coverage;
#' integration alone does not guarantee nominal confidence interval coverage.
#' @return A numeric scalar, or a list with \code{estimate},
#'   \code{grid_difference}, \code{n_grid}, \code{bounds},
#'   \code{tail_bound}, \code{converged}, and \code{bwd}.
#' @seealso \code{arocs_kernel}, \code{pdf_kernel}
#' @examples
#' set.seed(1)
#' x <- list(rnorm(50), rnorm(50, 1), rnorm(50, 2))
#' h <- vapply(x, bw.nrd0, numeric(1))
#' arocs_kernel_integral(x[[1]], x[[2]], x[[3]], bwd = h)
#' @export
arocs_kernel_integral <- function(X1, X2, X3,
    kernel_type = c("gauss", "epan"), bwd, na.rm = FALSE,
    n_grid = 257L, max_grid = 4097L, rel.tol = 1e-4, abs.tol = 1e-6,
    tail_prob = 1e-10, details = FALSE) {
  kernel_type <- match.arg(kernel_type)
  flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!flag(na.rm) || !flag(details)) stop("'na.rm' and 'details' must be logical scalars.")
  if (!is.numeric(bwd) || !is.null(dim(bwd)) ||
      !(length(bwd) %in% c(1L, 3L)) || any(!is.finite(bwd)) || any(bwd <= 0))
    stop("'bwd' must contain one or three positive finite bandwidths.")
  bwd <- rep(bwd, length.out = 3L)
  if (!scalar(n_grid) || n_grid < 3 || n_grid != floor(n_grid) ||
      !scalar(max_grid) || max_grid < 2*n_grid-1 || max_grid != floor(max_grid))
    stop("Grid sizes must be integers with n_grid >= 3 and max_grid >= 2*n_grid-1.")
  if (!scalar(rel.tol) || rel.tol < 0 || !scalar(abs.tol) || abs.tol < 0 ||
      max(rel.tol, abs.tol) == 0) stop("Invalid quadrature tolerances.")
  if (!scalar(tail_prob) || tail_prob <= 0 || tail_prob >= 1)
    stop("'tail_prob' must be between zero and one.")
  samples <- list(X1, X2, X3)
  for (k in seq_along(samples)) {
    x <- samples[[k]]
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || any(is.infinite(x)))
      stop("Samples must be nonempty numeric vectors without infinite values.")
    if (na.rm) x <- x[!is.na(x)]
    if (!length(x)) stop("A sample is empty after removing missing values.")
    samples[[k]] <- x
  }
  if (any(vapply(samples, anyNA, logical(1)))) return(NA_real_)
  # Standardize all groups together to keep density products on a safe scale.
  center <- median(unlist(samples, use.names = FALSE))
  unit <- max(c(vapply(samples, function(x) max(abs(x-center)), numeric(1)), bwd))
  if (!is.finite(unit)) stop("Data range is too large; rescale the observations.")
  samples <- lapply(samples, function(x) (x-center)/unit)
  h <- bwd/unit
  if (any(!is.finite(h)) || any(h <= 0))
    stop("Bandwidths are too small relative to the data range; rescale inputs.")
  radius <- if (kernel_type == "gauss") qnorm(tail_prob/2, lower.tail = FALSE) else 1
  lo <- min(vapply(samples, min, numeric(1)) - radius*h)
  hi <- max(vapply(samples, max, numeric(1)) + radius*h)
  evaluate <- function(n) {
    grid <- seq(lo, hi, length.out = n)
    f <- lapply(seq_len(3), function(k) pdf_kernel(grid, samples[[k]], kernel_type, h[k]))
    if (any(!is.finite(unlist(f)))) stop("Nonfinite KDE values; check bandwidths.")
    arocs_integral_grid_C(f[[1]], f[[2]], f[[3]], (hi-lo)/(n-1))
  }
  estimate <- evaluate(n_grid)
  difference <- Inf
  converged <- FALSE
  while (2*n_grid-1 <= max_grid) {
    n_grid <- 2*n_grid-1
    updated <- evaluate(n_grid)
    difference <- abs(updated-estimate)
    estimate <- updated
    if (difference <= max(abs.tol, rel.tol*abs(estimate))) {
      converged <- TRUE
      break
    }
  }
  if (!converged) warning("Quadrature did not converge; increase 'max_grid' and inspect 'details'.")
  if (!details) return(estimate)
  list(estimate = estimate, grid_difference = difference, n_grid = n_grid,
       bounds = c(lo, hi)*unit + center,
       tail_bound = if (kernel_type == "gauss") 6*tail_prob else 0,
       converged = converged, bwd = bwd)
}

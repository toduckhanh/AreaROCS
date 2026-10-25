log_sum_exp <- function(x) {
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

tanh_sinh_unit <- function(h = 0.2, K = 25L) {
  u <- (-K:K) * h
  a <- (pi / 2) * sinh(u)
  q <- (1 + tanh(a)) / 2
  w <- suppressWarnings(h * pi / 4 * cosh(u) / cosh(a)^2)
  keep <- is.finite(q) & is.finite(w) & q > 0 & q < 1 & w > 0
  list(nodes = q[keep], weights = w[keep])
}

TS <- tanh_sinh_unit()

#' Area of the ROC surface for three normal distributions
#'
#' Compute the area of the ROC surface (AROCS) for three normal class
#' distributions by numerical integration. Supplying estimated normal
#' parameters gives a normal-parametric plug-in estimate of AROCS.
#'
#' @param mu Numeric vector of length three containing the finite means of
#'   classes 1, 2 and 3, in that order.
#' @param sigma Numeric vector of length three containing the strictly positive,
#'   finite standard deviations of classes 1, 2 and 3, in that order.
#'   These are standard deviations, not variances.
#'
#' @details
#' Let \eqn{F_k} and \eqn{f_k} be the distribution and density functions of
#' \eqn{N(\mu_k, \sigma_k^2)}, respectively. The function evaluates
#' \deqn{
#' A = \int_0^1 \int_0^{1-F_3(F_1^{-1}(p_1))}
#' \sqrt{1 +
#' \left\{\frac{f_2(F_1^{-1}(p_1))}{f_1(F_1^{-1}(p_1))}\right\}^2 +
#' \left\{\frac{f_2(F_3^{-1}(1-p_3))}{f_3(F_3^{-1}(1-p_3))}\right\}^2}
#' \,dp_3\,dp_1.
#' }
#' This is the surface area, not the volume under the ROC surface (VUS).
#' Classes must be supplied in the intended ordinal order, with larger test
#' values interpreted as indicating higher classes. Parameters are not sorted
#' automatically. All parameter entries must be supplied without missing values;
#' the current implementation does not explicitly validate these requirements.
#'
#' Integration uses a fixed tanh-sinh quadrature rule on the probability scale,
#' with internal step size \code{h = 0.2} and truncation index \code{K = 25}.
#' Density ratios and quadrature contributions are combined on the log scale
#' to reduce numerical overflow and underflow. The returned value is a numerical
#' approximation: the function does not provide an integration error estimate
#' or an adaptive convergence check. Accuracy should be checked independently
#' for extreme separation or markedly different class standard deviations.
#'
#' Known population parameters give a numerical population AROCS value. For
#' normal-parametric estimation, replace them with estimates from the three
#' samples. The normal maximum likelihood standard deviation is
#' \code{sqrt(mean((x - mean(x))^2))}; \code{sd(x)} instead uses the divisor
#' \code{length(x) - 1}. The function does not fit a model or compute a
#' standard error or confidence interval.
#'
#' @return A numeric scalar containing the numerically approximated AROCS.
#' @references Nguyen, C.-K.-A., Phung, V.-T., To, D.-K. and Adimari, G.
#'   \emph{Area of the ROC surface}, equations (2.1) and (3.3), manuscript.
#' @seealso \code{\link{arocs_kernel}}, \code{\link{arocs_kernel_integral}},
#'   \code{\link[stats]{Normal}}
#' @examples
#' # Population AROCS: sigma contains standard deviations.
#' arocs_normal(mu = c(0, 1, 2), sigma = sqrt(c(1, 1.5, 1)))
#'
#' # Identical distributions: area is sqrt(3) / 2.
#' arocs_normal(mu = rep(0, 3), sigma = rep(1, 3))
#' sqrt(3) / 2
#'
#' # Normal-parametric plug-in estimate using maximum likelihood parameters.
#' set.seed(2510)
#' samples <- list(rnorm(50, 0, 1), rnorm(50, 1, sqrt(1.5)), rnorm(50, 2, 1))
#' mu_hat <- vapply(samples, mean, numeric(1))
#' sigma_hat <- vapply(samples, function(x) sqrt(mean((x - mean(x))^2)),
#'                     numeric(1))
#' arocs_normal(mu = mu_hat, sigma = sigma_hat)
#' @export
arocs_normal <- function(mu, sigma) {
  # Integrate first over p1 in (0,1), then over p3 in
  # (0, 1-F3(F1^{-1}(p1))). Tanh-sinh nodes resolve the steep boundary faces
  # considerably better than a rectangular probability grid.
  pieces <- vector("list", length(TS$nodes))
  for (i in seq_along(TS$nodes)) {
    p1 <- TS$nodes[i]
    x <- qnorm(p1, mu[1], sigma[1])
    upper <- pnorm(x, mu[3], sigma[3], lower.tail = FALSE)
    if (!is.finite(upper) || upper <= 0) next
    p3 <- upper * TS$nodes
    keep <- is.finite(p3) & p3 > 0
    if (!any(keep)) next
    p3 <- p3[keep]
    y <- qnorm(p3, mu[3], sigma[3], lower.tail = FALSE)
    log_r1 <- dnorm(x, mu[2], sigma[2], log = TRUE) -
      dnorm(x, mu[1], sigma[1], log = TRUE)
    log_r3 <- dnorm(y, mu[2], sigma[2], log = TRUE) -
      dnorm(y, mu[3], sigma[3], log = TRUE)
    m <- pmax(0, 2 * log_r1, 2 * log_r3)
    log_element <- 0.5 * (
      m + log(exp(-m) + exp(2 * log_r1 - m) + exp(2 * log_r3 - m))
    )
    pieces[[i]] <- log(TS$weights[i]) + log(upper) +
      log(TS$weights[keep]) + log_element
  }
  log_values <- unlist(pieces, use.names = FALSE)
  return(exp(log_sum_exp(log_values[is.finite(log_values)])))
}

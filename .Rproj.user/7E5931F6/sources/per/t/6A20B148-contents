#' Kernel density estimator of the area of the ROC surface
#'
#' Compute the empirical pair average in equation (3.1) of
#' \emph{Area of the ROC surface}, using \code{pdf_kernel} for the density
#' ratios. Groups must be supplied in their intended ordinal order.
#'
#' @param X1,X2,X3 Nonempty finite numeric vectors for classes 1, 2 and 3.
#' @param kernel_type Kernel used for all three density estimates: Gaussian
#'   (\code{"gauss"}) or Epanechnikov (\code{"epan"}).
#' @param bwd Positive finite bandwidth, either a scalar shared by all classes
#'   or a vector of length three in class order. Bandwidths are required;
#'   they can be selected with the package bandwidth selectors.
#' @param na.rm Logical; remove missing observations separately within each
#'   class. Otherwise, any missing observation gives \code{NA_real_}.
#'
#' @details
#' With kernel density estimates \eqn{\hat f_k}, the estimator is
#' \deqn{\frac{1}{n_1 n_3}\sum_{i=1}^{n_1}\sum_{j=1}^{n_3}
#' I(X_{1i}<X_{3j})\sqrt{1+
#' (\hat f_2(X_{1i})/\hat f_1(X_{1i}))^2+
#' (\hat f_2(X_{3j})/\hat f_3(X_{3j}))^2}.}
#' Each density uses the entire corresponding sample, including the evaluation
#' observation for the denominators (no leave-one-out correction). Thus the
#' denominators are positive for both supported kernels in ordinary floating
#' point arithmetic. Nonfinite densities or ratios and zero denominators cause
#' an error rather than silently changing the estimator with a density floor.
#'
#' Ties contribute zero. Normalization uses all \eqn{n_1 n_3} pairs, not just
#' ordered pairs. No clipping to population bounds is applied. In particular,
#' finite-sample estimates need not lie in the population range stated in the
#' manuscript. Pair contributions are computed one row at a time to avoid
#' allocating an \eqn{n_1} by \eqn{n_3} matrix.
#'
#' @return A numeric scalar containing the kernel estimate of AROCS.
#' @references Nguyen, C.-K.-A., Phung, V.-T., To, D.-K. and Adimari, G.
#'   \emph{Area of the ROC surface}, equation (3.1), manuscript.
#' @seealso \code{pdf_kernel}, \code{bwd_plugin}, \code{bwd_cv}
#' @examples
#' X1 <- c(-1.2, -0.5, 0, 0.4)
#' X2 <- c(-0.2, 0.6, 1, 1.5)
#' X3 <- c(0.5, 1.2, 1.8, 2.4)
#' arocs_kernel(X1, X2, X3, bwd = 0.5)
#' h <- vapply(list(X1, X2, X3), bwd_plugin, numeric(1), bwd_method = "RT")
#' arocs_kernel(X1, X2, X3, kernel_type = "epan", bwd = h)
#' @export
arocs_kernel <- function(X1, X2, X3, kernel_type = c("gauss", "epan"),
                         bwd, na.rm = FALSE) {
  kernel_type <- match.arg(kernel_type)
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm))
    stop("'na.rm' must be TRUE or FALSE.")
  if (!is.numeric(bwd) || !is.null(dim(bwd)) ||
      !(length(bwd) %in% c(1L, 3L)) ||
      any(!is.finite(bwd)) || any(bwd <= 0))
    stop("'bwd' must contain one or three positive finite bandwidths.")
  bwd <- rep(bwd, length.out = 3L)
  samples <- list(X1, X2, X3)
  for (k in seq_along(samples)) {
    x <- samples[[k]]
    if (!is.numeric(x) || !is.null(dim(x)) || length(x) == 0L)
      stop("'X", k, "' must be a nonempty numeric vector.")
    if (any(is.infinite(x))) stop("'X", k, 
                                  "' must not contain infinite values.")
    if (na.rm) x <- x[!is.na(x)]
    if (length(x) == 0L)
      stop("'X", k, "' has no observations after removing missing values.")
    samples[[k]] <- x
  }
  if (any(vapply(samples, anyNA, logical(1)))) return(NA_real_)
  X1 <- samples[[1L]]
  X2 <- samples[[2L]]
  X3 <- samples[[3L]]
  f11 <- pdf_kernel(X1, X1, kernel_type, bwd[1L])
  f21 <- pdf_kernel(X1, X2, kernel_type, bwd[2L])
  f23 <- pdf_kernel(X3, X2, kernel_type, bwd[2L])
  f33 <- pdf_kernel(X3, X3, kernel_type, bwd[3L])
  if (any(!is.finite(c(f11, f21, f23, f33))) ||
      any(f11 <= 0) || any(f33 <= 0))
    stop("Kernel densities are numerically invalid; rescale data or change bandwidths.")
  r21 <- f21 / f11
  r23 <- f23 / f33
  if (any(!is.finite(c(r21, r23))))
    stop("Density ratios overflowed; rescale data or change bandwidths.")
  row_means <- vapply(seq_along(X1), function(i) {
    j <- which(X1[i] < X3)
    if (!length(j)) return(0)
    # Scaled Euclidean norm avoids overflow when squaring large ratios.
    scale <- pmax(1, r21[i], r23[j])
    norm <- sqrt((1 / scale)^2 + (r21[i] / scale)^2 + (r23[j] / scale)^2)
    sum((scale / length(X3)) * norm)
  }, numeric(1))
  mean(row_means)
}

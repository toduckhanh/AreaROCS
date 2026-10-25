library(AreaROCS)

# Independent direct KDE and double-sum reference, including ties and
# unequal sample sizes. This also checks bandwidth-to-class mapping.
kde <- function(x, sample, h, type) {
  vapply(x, function(z) {
    u <- (z - sample) / h
    mean(if (type == "gauss") dnorm(u) else
      0.75 * pmax(1 - u^2, 0)) / h
  }, numeric(1))
}
reference <- function(a, b, c, h, type) {
  r <- kde(a, b, h[2], type) / kde(a, a, h[1], type)
  s <- kde(c, b, h[2], type) / kde(c, c, h[3], type)
  out <- 0
  for (i in seq_along(a)) for (j in seq_along(c)) {
    if (a[i] < c[j]) out <- out + sqrt(1 + r[i]^2 + s[j]^2)
  }
  out / (length(a) * length(c))
}
equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-12)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr, error = identity), "error"))
a <- c(-2, 0, 1)
b <- c(-1, 0.4, 0.9, 2)
c <- c(0, 1, 1.4, 3, 4)
h <- c(0.4, 0.7, 1.1)
for (type in c("gauss", "epan")) {
  equal(arocs_kernel(a, b, c, type, h), reference(a, b, c, h, type))
  equal(arocs_kernel(a, b, c, type, 0.8),
        reference(a, b, c, rep(0.8, 3), type))
  equal(arocs_kernel(rev(a), rev(b), rev(c), type, h),
        arocs_kernel(a, b, c, type, h))
  equal(arocs_kernel(2*a + 5, 2*b + 5, 2*c + 5, type, 2*h),
        arocs_kernel(a, b, c, type, h))
  # Identical KDEs give ratios one; strict ties still contribute zero.
  x <- c(-1, 0, 1, 2)
  equal(arocs_kernel(x, x, x, type, 1), sqrt(3) * 3/8)
  equal(arocs_kernel(1, 1, 1, type, 1), 0)
  equal(arocs_kernel(c(3, 4), b, c(0, 1), type, h), 0)
  equal(arocs_kernel(c(a, NA), c(NA, b), c(c, NaN), type, h, na.rm = TRUE),
        arocs_kernel(a, b, c, type, h))
}
# Compact disjoint supports give zero numerator densities, hence value one
# when every pair is ordered. Do not replace this estimator by an integral
# or force it to the population upper bound of three.
equal(arocs_kernel(0, 10, 20, "epan", 1), 1)
stopifnot(identical(arocs_kernel(c(a, NA), b, c, bwd = h), NA_real_))
for (bad in list(0, -1, Inf, NA_real_, c(1, 2), c(1, 2, 0), "a"))
  fails(arocs_kernel(a, b, c, bwd = bad))
for (bad in list(numeric(), "a", Inf, matrix(1:4, 2)))
  fails(arocs_kernel(bad, b, c, bwd = h))
fails(arocs_kernel(NA_real_, b, c, bwd = h, na.rm = TRUE))
fails(arocs_kernel(a, b, c, bwd = h, na.rm = NA))
fails(arocs_kernel(a, b, c, "unknown", h))

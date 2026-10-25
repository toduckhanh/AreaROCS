library(AreaROCS)
near <- function(x, y, tol = 3e-5) stopifnot(abs(x-y) < tol)
for (kernel in c("gauss", "epan")) {
  x <- c(-1, -0.3, 0.5, 2)
  near(arocs_kernel_integral(x, x, x, kernel, 0.7, rel.tol = 1e-6,
                            max_grid = 8193), sqrt(3)/2)
  a <- c(-1, 0, 1); b <- c(-0.2, 0.8); c <- c(0, 1, 2, 3)
  h <- c(0.6, 0.8, 0.9)
  v <- arocs_kernel_integral(a,b,c,kernel,h)
  near(v, arocs_kernel_integral(3*a+10,3*b+10,3*c+10,kernel,3*h))
  near(v, arocs_kernel_integral(rev(a),rev(b),rev(c),kernel,h))
}
# Single-observation Gaussian KDEs are exactly normal densities. Check against
# independent adaptive integration rather than the package's grid integrator.
f <- function(x, y) sqrt((dnorm(x,0,1)*dnorm(y,2,1))^2 +
                         (dnorm(x,1,sqrt(1.5))*dnorm(y,2,1))^2 +
                         (dnorm(x,0,1)*dnorm(y,1,sqrt(1.5)))^2)
oracle <- integrate(function(x) vapply(x, function(z)
  integrate(function(y) f(z,y), z, 14, rel.tol=1e-8)$value,
  numeric(1)), -12, 14, rel.tol=1e-8)$value
near(arocs_kernel_integral(0,1,2,bwd=sqrt(c(1,1.5,1)),rel.tol=1e-6,
                          max_grid=8193), oracle, 3e-6)
# Three ordered, disjoint compact supports: all three surface faces count.
near(arocs_kernel_integral(0,4,8,"epan",1,rel.tol=1e-6,max_grid=8193),3)
d <- arocs_kernel_integral(0,1,2,bwd=1,details=TRUE)
stopifnot(d$converged, d$grid_difference <= 1e-4*d$estimate,
          d$n_grid >= 513, d$tail_bound == 6e-10)
stopifnot(is.na(arocs_kernel_integral(NA_real_,1,2,bwd=1)))
near(arocs_kernel_integral(c(0,NA),1,2,bwd=1,na.rm=TRUE),
     arocs_kernel_integral(0,1,2,bwd=1))
fails <- function(expr) stopifnot(inherits(tryCatch(expr,error=identity),"error"))
fails(arocs_kernel_integral(numeric(),1,2,bwd=1))
fails(arocs_kernel_integral(0,1,2,bwd=0))
fails(arocs_kernel_integral(0,1,2,bwd=1,max_grid=257))
fails(arocs_kernel_integral(0,1,2,bwd=1,tail_prob=0))

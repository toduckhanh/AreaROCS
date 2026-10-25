### Name: arocs_normal
### Title: Area of the ROC surface for three normal distributions
### Aliases: arocs_normal

### ** Examples

# Population AROCS: sigma contains standard deviations.
arocs_normal(mu = c(0, 1, 2), sigma = sqrt(c(1, 1.5, 1)))

# Identical distributions: area is sqrt(3) / 2.
arocs_normal(mu = rep(0, 3), sigma = rep(1, 3))
sqrt(3) / 2

# Normal-parametric plug-in estimate using maximum likelihood parameters.
set.seed(2510)
samples <- list(rnorm(50, 0, 1), rnorm(50, 1, sqrt(1.5)), rnorm(50, 2, 1))
mu_hat <- vapply(samples, mean, numeric(1))
sigma_hat <- vapply(samples, function(x) sqrt(mean((x - mean(x))^2)),
                    numeric(1))
arocs_normal(mu = mu_hat, sigma = sigma_hat)




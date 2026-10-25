# Integrated KDE estimator of AROCS

The existing `arocs_kernel()` remains the empirical-pair estimator (3.1).
The new `arocs_kernel_integral()` is the full KDE plug-in functional:

\[
\widehat A_{\mathrm{int}}=\int_{x<y}
\sqrt{[\widehat f_1(x)\widehat f_3(y)]^2+
       [\widehat f_2(x)\widehat f_3(y)]^2+
       [\widehat f_1(x)\widehat f_2(y)]^2}\,dx\,dy.
\]

Multiplying the population ratio integrand by f1(x)f3(y) gives this expression.
Replacing all densities with KDEs gives a different finite-sample estimator
from (3.1). It uses `pdf_kernel()` without dividing estimated densities.
For compact kernels this expression also handles zero densities directly.

The old denominator includes its own observation, contributing K(0)/(n*h).
Also, empirical pair averaging samples only the observed class-1/class-3
locations; large ratios in sparsely sampled tails can have substantial weight
in the target integral. Both mechanisms can contribute to finite-sample bias.
The simulation below establishes the combined improvement, not a separate
causal decomposition of those two mechanisms. Applying exp() to log-ratios
does not fix either statistical issue. Leave-one-out alone can introduce zero
denominators for compact kernels or extremely small ones for Gaussian kernels.

## Validation run

Seed 25102026; 1,000 independent samples; n1=n2=n3=50;
mu=(0,1,2); sigma=sqrt(1,1.5,1); Gaussian kernels;
bandwidths selected by bw.nrd0 separately in each sample and group.
True area 1.5378919 (rounded), independently checked by adaptive integration
on the density-product expression and by refining the normal tanh-sinh rule.

| Estimator | Mean | Bias | SD | RMSE | Oracle-SD diagnostic |
|---|---:|---:|---:|---:|---:|
| Normal parametric | 1.545461 | 0.007569 | 0.080955 | 0.081268 | 0.948 |
| Kernel empirical pairs (3.1) | 1.440215 | -0.097676 | 0.068582 | 0.119330 | 0.683 |
| Kernel integrated plug-in | 1.513568 | -0.024324 | 0.076454 | 0.080194 | 0.934 |

The last column reproduces the vignette's diagnostic using the SD over all
Monte Carlo replications. It is not coverage of an implementable confidence
interval based on one data set. Do not describe 0.934 as verified bootstrap
coverage. No bootstrap coverage study was run for this report.

The new estimator reduces absolute bias by about 75% in this scenario, but
does not eliminate smoothing bias. No bandwidth was tuned against the truth.
Normal parametric estimation is appropriate when the normal model is correct;
performance in this scenario does not establish performance for other models.

## Usage and numerical checks

```r
h <- vapply(list(X1, X2, X3), bw.nrd0, numeric(1))
arocs_kernel_integral(X1, X2, X3, bwd=h)
arocs_kernel_integral(X1, X2, X3, bwd=h, details=TRUE)
```

The numerical method refines a common grid, computes densities once per grid,
and sums the triangular integral in C++. No quadratic-sized matrix is stored.
The diagonal has half weight as part of quadrature, not as an empirical tie
correction. Gaussian tails are truncated with omitted-area bound 6*tail_prob;
Epanechnikov integration includes the full support. Grid differences are
convergence diagnostics, not rigorous error bounds. Check a finer starting
grid when bandwidths are very narrow or very different across groups.

Tests cover the identical-density area sqrt(3)/2, three ordered disjoint
Epanechnikov supports (area 3), agreement with adaptive normal integration,
common affine transformations, permutations, missing data and invalid inputs.
The original empirical-pair regression tests also pass.

## Confidence intervals

`compare_arocs_estimators.R` contains an optional per-sample bootstrap study
with independent within-group resampling and bandwidth reselection. It reports
percentile and basic interval coverage, Monte Carlo standard errors, and mean
widths. Use at least 1,000 outer replications and 999 bootstrap replications
for a more informative study; this is computationally more expensive.

Integration fixes how the fitted KDE functional is evaluated, not its
smoothing bias or bootstrap validity. Undersmoothing and explicit bias
correction require separate evaluation for this functional; a result for KDE
confidence bands is not automatically a theorem for AROCS. For background on
why smoothing bias complicates bootstrap inference, see Hall and Horowitz,
[A simple bootstrap method for constructing nonparametric confidence bands
for functions](https://arxiv.org/abs/1309.4864).

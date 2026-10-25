# Run from the package root after installing the updated package.
# Optional environment variables: AROCS_MC_REPS, AROCS_BOOT_MC, AROCS_BOOT_REPS.
library(AreaROCS)
set.seed(25102026)
mu <- c(0,1,2)
sigma <- sqrt(c(1,1.5,1))
n <- 50L
env_int <- function(key, default) as.integer(Sys.getenv(key, default))
M <- env_int("AROCS_MC_REPS", 1000)
MC_boot <- env_int("AROCS_BOOT_MC", 0)
B <- env_int("AROCS_BOOT_REPS", 199)
out_dir <- "tmp/arocs-comparison"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Refine the existing normal quadrature locally without changing package state.
normal_refined <- arocs_normal
normal_env <- new.env(parent=environment(arocs_normal))
normal_env$TS <- AreaROCS:::tanh_sinh_unit(h=0.05,K=100L)
environment(normal_refined) <- normal_env
truth <- normal_refined(mu,sigma)
cat("Normal default:",arocs_normal(mu,sigma),"refined:",truth,"\n")

estimate <- function(x) {
  h <- vapply(x,bw.nrd0,numeric(1))
  c(normal=normal_refined(vapply(x,mean,numeric(1)),vapply(x,sd,numeric(1))),
    pair=arocs_kernel(x[[1]],x[[2]],x[[3]],bwd=h),
    integral=arocs_kernel_integral(x[[1]],x[[2]],x[[3]],bwd=h))
}
draw <- function() lapply(seq_len(3),function(k) rnorm(n,mu[k],sigma[k]))
values <- replicate(M,estimate(draw()))
mc_sd <- apply(values,1,sd)
summary <- data.frame(method=rownames(values), truth=truth,
  mean=rowMeans(values), bias=rowMeans(values)-truth, sd=mc_sd,
  rmse=sqrt(rowMeans((values-truth)^2)),
  oracle_sd_coverage=rowMeans(abs(values-truth) <= qnorm(.975)*mc_sd))
print(summary)
write.csv(summary,file.path(out_dir,"mc-summary.csv"),row.names=FALSE)
saveRDS(values,file.path(out_dir,"mc-estimates.rds"))

# Genuine per-sample nonparametric bootstrap intervals. Resample all three
# groups independently and reselect bandwidths within every bootstrap sample.
# Basic and percentile intervals are both diagnostic candidates, not a
# guarantee of 95% coverage. MC_boot=200 gives coverage MCSE about 0.015 at .95.
if (MC_boot > 0) {
  cover <- width <- array(NA_real_,c(3,2,MC_boot),
    dimnames=list(names(estimate(draw())),c("percentile","basic"),NULL))
  for (r in seq_len(MC_boot)) {
    x <- draw()
    point <- estimate(x)
    boots <- replicate(B,estimate(lapply(x,function(z) sample(z,length(z),replace=TRUE))))
    for (k in seq_len(3)) {
      q <- quantile(boots[k,],c(.025,.975),names=FALSE)
      intervals <- cbind(percentile=q,basic=2*point[k]-rev(q))
      cover[k,,r] <- intervals[1,] <= truth & truth <= intervals[2,]
      width[k,,r] <- intervals[2,]-intervals[1,]
    }
    if (r %% 10 == 0) cat("Bootstrap MC",r,"/",MC_boot,"\n")
  }
  coverage <- apply(cover,c(1,2),mean)
  result <- list(coverage=coverage,mcse=sqrt(coverage*(1-coverage)/MC_boot),
                 mean_width=apply(width,c(1,2),mean),M=MC_boot,B=B,truth=truth)
  print(result)
  saveRDS(result,file.path(out_dir,"bootstrap-summary.rds"))
}

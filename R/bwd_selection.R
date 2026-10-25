# plugin method ----
#' @export
bwd_plugin <- function(X, bwd_method = c("RT", "AZZ"), hc = NULL){
  bwd.method <- match.arg(bwd_method)
  if(bwd_method %in% c("RT", "AZZ") && is.null(hc)) {
    hc <- switch (bwd_method,
      RT = 0.9,
      AZZ = 1.3
    )
  }
  bwd <- switch(bwd_method,
                RT = hc*min(c(sd(X), IQR(X)/1.34))*length(X)^(-1/5),
                AZZ = hc*sd(X)*length(X)^(-1/3)
  )
  return(bwd)
}

# leave-one-out cross validation method ----
#' @export
bwd_cv <- function(X, kernel_type = c("gauss", "epan"),
                   n_bwd = 101, n_x = 151){
  if(any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  x_grid <- seq(min(X), max(X), length.out = n_x)
  hx <- diff(x_grid)[1]
  range_X <- max(X) - min(X)
  bwd_seq <- seq(range_X/200, range_X/2, length.out = n_bwd)
  res <- cv_bwd_C(x = x_grid, X = X, Ktype = kernel_type_value, hx = hx,
                  bwd_seq = bwd_seq)
  return(res)
}

# PCO method ----
#' @export
bwd_pco <- function(X, kernel_type = c("gauss", "epan"), n_bwd = 101,
                    h_min = NULL){
  if(any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  n_X <- length(X)
  range_X <- max(X) - min(X)
  if(is.null(h_min)) h_min <- 1/(sqrt(2*pi)*n_X)
  bwd_seq <- seq(h_min, range_X/2, length.out = n_bwd)
  res <- PCO_bwd_C(X = X, Ktype = kernel_type_value, bwd_seq = bwd_seq)
  return(res)
}

# k-fold cross validation method for Projection ----
#' @export
kcv_proj <- function(X, k_fold = 5, M = NULL, n_x = 151){
  if(any(is.na(X))) return(NA)
  ab <- range(X)
  n_X <- length(X)
  if(is.null(M)) M <- c(1:(min(100, n_X)/2))
  nk_fold <- ceiling(n_X/k_fold)
  id_cv <- lapply(1:(k_fold - 1), function(i){
    (1 + (i - 1)*nk_fold) : (i*nk_fold)
  })
  id_cv[[k_fold]] <- (1 + (k_fold - 1)*nk_fold):n_X
  x_grid <- seq(min(X), max(X), length.out = n_x)
  hx <- diff(x_grid)[1]
  cv_out <- sapply(M, function(y){
    temp <- sapply(1:k_fold, function(i) {
      temp1 <- sapply(x_grid, function(u) mean(u - X[id_cv[[i]]] >= 0))
      temp2 <- cdf_proj(x = x_grid, X = X[-id_cv[[i]]], M = y, ab = ab)
      temp3 <- (temp1 - temp2)^2
      return(simpson(temp3, n_x - 1, hx))
    })
    cv_proj <- mean(temp)
    return(cv_proj)
  })
  out <- list()
  out$M = M;
  out$cv_out = cv_out;
  out$m_opt = M[which.min(cv_out)];
  return(out)
}

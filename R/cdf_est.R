#' @importFrom Rcpp evalCpp
#' @useDynLib AreaROCS, .registration = TRUE

#' @export
cdf_kernel <- function(x, X, kernel_type = c("gauss", "epan"), bwd){
  if(any(is.na(x)) | any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  res <- cdf_kernel_C(x, X, kernel_type_value, bwd)
  return(res)
}

#' @export
cdf_proj <- function(x, X, M, ab = NULL){
  if(any(is.na(x)) | any(is.na(X))) return(NA)
  if(is.null(ab)) ab <- range(X)
  res <- cdf_proj_C(x, X, ab[1], ab[2], M)
  return(res)
}

#' @export
pdf_kernel <- function(x, X, kernel_type = c("gauss", "epan"), bwd){
  if(any(is.na(x)) | any(is.na(X))) return(NA)
  kernel_type <- match.arg(kernel_type)
  kernel_type_value <- switch(kernel_type,
                              "gauss" = 0,
                              "epan" = 1)
  res <- pdf_kernel_C(x, X, kernel_type_value, bwd)
  return(res)
}
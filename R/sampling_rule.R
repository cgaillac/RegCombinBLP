#' The subsampling rule
#'
#' @param n sample size.
#'
#' @return
#' the subsampling size
#'
sampling_rule_blp_main <- function(n){
    res =   1.5*(0.4*n - 0.2*max(n-5,0) - 0.1*max(n-1000,0) - 0.1*(1-log(4000)/log(n))*max(n-4000,0))
return(res)}

sampling_rule_nonstd <- function(n){
  res =   1*(0.5*n - 0.3*max(n-5,0) - 0.15*max(n-1000,0) - 0.05*(1-log(3000)/log(n))*max(n-3000,0))
return(res)}

sampling_rule_blp_others <- function(n){
  res = 0.5*(0.5*n - 0.3*max(n-5,0) - 0.1*max(n-1000,0) - 0.1*(1-log(4000)/log(n))*max(n-4000,0))
return(res)}


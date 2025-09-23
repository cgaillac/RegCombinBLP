#' Compute Bounds on the BLP
#'
#' Computes bounds on the BLP coefficients, optionally
#' using different types of conditioning variables, weights, bootstrapping, and factor/discrete adjustments.
#'
#' @param sample1 Optional vector of indices for subsampling or bootstrapping (default NULL).
#' @param Xc_x Matrix of conditioning covariates for the right-hand side data.
#' @param Xnc Matrix of  outside regressors (non-commonly observed) (used in regression).
#' @param Xc_y Matrix of conditioning covariates for the left-hand side data.
#' @param Wc_x Optional matrix of additional covariates (weights/factors) for right-hand side data.
#' @param Wc_y Optional matrix of additional covariates (weights/factors) for left-hand side data.
#' @param Y Vector or matrix of outcome/dependent variable values.
#' @param w_x Optional vector of weights for right-hand side data (default NULL).
#' @param w_y Optional vector of weights for left-hand side data (default NULL).
#' @param dimXc Integer: number of conditioning covariates.
#' @param dimW Integer: number of weight/factor covariates.
#' @param dimXnc Integer: number of non-conditioning covariates.
#' @param values Optional values used in bounds computation (default NULL).
#' @param full Logical: whether to use full computation in BLP bounds (default FALSE).
#' @param FW Logical: whether to apply Frisch-Waugh adjustments (default FALSE).
#' @param discrete Logical: whether to treat conditioning variables as discrete (default FALSE).
#' @param unchanged Logical: whether to keep all values when discretizing (default TRUE).
#' @param rule Numeric: sampling rule for bootstrap/subsampling (default 1).
#' @param factor Logical: whether to treat W covariates as factors (default FALSE).
#' @param K_sel Integer: number of clusters for discrete/factor variables (default 0, auto-selected).
#' @param dataset Integer: dataset type flag (default 0).
#' @param O_only Logical, whether to compute bounds only for outside regressors (default FALSE).
#'
#' @return A list or object containing BLP bounds. Returns DGM BLP bounds,#
#'  including hull points and influence function data (if applicable).
#'
#' @details
#' The function proceeds in the following steps:
#' 1. **Resampling**: If `sample1` is provided, performs subsampling or bootstrapping
#'    for predictor (`Xnc`, `Xc_x`, `Wc_x`) and outcome (`Y`, `Xc_y`, `Wc_y`) matrices.
#' 2. **Weight normalization**: Weights (`w_x`, `w_y`) are normalized to sum to 1 for the selected sample.
#' 3. **Bounds computation**: Depending on `DP`:
#'    a. `DP = TRUE`: Calls `bounds_Pacini_Xc()` to compute Pacini-style bounds.
#'    b. `DP = FALSE`: Calls `bounds_BLP()` to compute standard BLP bounds with optional full, discrete, and factor adjustments.
#' 4. Returns the computed bounds and related outputs.
#'
#' @examples
#' \dontrun{
#' # Basic usage with random data
#' results <- compute_BLP(
#'   sample1 = NULL,
#'   Xc_x = matrix(rnorm(100), 50, 2),
#'   Xnc = matrix(rnorm(50), 50, 1),
#'   Xc_y = matrix(rnorm(100), 50, 2),
#'   Wc_x = NULL,
#'   Wc_y = NULL,
#'   Y = rnorm(50),
#'   dimXc = 2,
#'   dimW = 0,
#'   dimXnc = 1
#' )
#' }
#'
#' @export

compute_BLP <- function(sample1 = NULL,Xc_x,Xnc,Xc_y,Wc_x,Wc_y,Y,w_x =NULL,w_y=NULL,dimXc,
                        dimW,dimXnc,
                        values=NULL,
                        full=FALSE,FW=FALSE,
                        discrete=FALSE,
                        unchanged = TRUE, rule=1,
                        factor=FALSE,K_sel=0,
                        dataset=0,
                        O_only = FALSE){



  exact=2
  weights_x=w_x
  weights_y=w_y

  if(!is.null(sample1)){


    n_x = dim(Xnc)[1]
    n_y = dim(Y)[1]
    n_xy = min(n_x,n_y)
    T_xy  = (n_y/(n_x+n_y))*n_x

    #
    if(bootstrap==FALSE){
      if(!FW  & !DP & !discrete){
        bs = rule*floor(sampling_rule_blp_main(T_xy))
      }else{
        bs = rule*floor(sampling_rule_blp_others(T_xy))
      }

      if(!is.null( weights_x)){
        bb = sample(1:n_x,bs, replace=FALSE)
      }else{
        bb = sample(1:n_x,bs, replace=FALSE)
      }

      weights_x= rep(1/dim(Xnc)[1],dim(Xnc)[1])

      if(!is.null(Xc_x)){
        Xc_xb = matrix(Xc_x[bb,],bs,dimXc)
      }else{
        Xc_xb =NULL
      }

      if(!is.null(Wc_x)){
        Wc_xb = matrix(Wc_x[bb,],bs,dimW)
      }else{
        Wc_xb =NULL
      }


      Xncb = as.matrix(Xnc[bb,],bs,dimXnc)
      weights_x =  matrix(weights_x[bb],bs,1)
      weights_x = weights_x/sum(weights_x)
      n_x = dim(Xncb)[1]

      if(!is.null( weights_y)){
        bby = sample(1:n_y,bs, replace=FALSE)
      }else{
        bby = sample(1:n_y,bs, replace=FALSE)
      }

      weights_y= rep(1/length(Y),length(Y))


      if(!is.null(Xc_y)){
        Xc_yb = matrix(Xc_y[bby,],bs,dimXc)
      }else{
        Xc_yb =NULL
      }


      if(!is.null(Wc_y)){
        Wc_yb = matrix(Wc_y[bby,],bs,dimW)
      }else{
        Wc_yb =NULL
      }


      Yb = matrix(Y[bby],bs,1)
      weights_y =  matrix(weights_y[bby],bs,1)
      weights_y = weights_y/sum(weights_y)
      n_y = dim(Yb)[1]
    }else{

      ######## boot
      bs = floor(n_x/rule)

      if(!is.null( weights_x)){
        bb = sample(1:n_x,bs, replace=TRUE, prob=  weights_x)
      }else{
        bb = sample(1:n_x,bs, replace=TRUE)
      }

      weights_x= rep(1/dim(Xnc)[1],dim(Xnc)[1])

      if(!is.null(Xc_x)){
        Xc_xb = matrix(Xc_x[bb,],bs,dimXc)
      }else{
        Xc_xb =NULL
      }

      if(!is.null(Wc_x)){
        Wc_xb = matrix(Wc_x[bb,],bs,dimW)
      }else{
        Wc_xb =NULL
      }

      Xncb = matrix(Xnc[bb,],bs,dimXnc)
      weights_x =  matrix(weights_x[bb],bs,1)
      weights_x = weights_x/sum(weights_x)
      n_x = dim(Xncb)[1]

      bs =  floor(n_y/rule)
      if(!is.null( weights_x)){
        bby = sample(1:n_y,bs, replace=TRUE, prob=  weights_y)
      }else{
        bby = sample(1:n_y,bs, replace=TRUE)
      }

      weights_y= rep(1/length(Y),length(Y))


      if(!is.null(Xc_y)){
        Xc_yb = matrix(Xc_y[bby,],bs,dimXc)
      }else{
        Xc_yb =NULL
      }

      if(!is.null(Wc_y)){
        Wc_yb = matrix(Wc_y[bby,],bs,dimW)
      }else{
        Wc_yb =NULL
      }

      Yb = matrix(Y[bby],bs,1)
      weights_y =  matrix(weights_y[bby],bs,1)
      weights_y = weights_y/sum(weights_y)
      n_y = dim(Yb)[1]
    }

  }else{

    ## point estimate
    Wc_xb =Wc_x
    Wc_yb =Wc_y
    Xc_xb =Xc_x
    Xncb = Xnc
    Xc_yb = Xc_y
    Yb = Y

  }

  dataX =NULL
  dataP =NULL

    out <- bounds_BLP(Xncb,Yb,Xc_xb,Xc_yb,Wc_xb,Wc_yb,
                      w_x=weights_x,w_y=weights_y,
                      dimXc, dimW,values,full,FW,discrete, unchanged,factor,K_sel,
                      dataset,  O_only )



  return(out)




}

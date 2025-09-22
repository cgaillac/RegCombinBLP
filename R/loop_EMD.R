#' This function computes the upper and lower bounds of a coefficient using
#' the 1-dimensional Earth Mover's Distance (EMD) for a specific subset of
#' data determined by the value of `Xc_x0` and `Xc_y0`. It is used as a helper
#' function in computing bounds for the BLP method.
#'
#' @param k Integer. Index of the value in `values` for which to compute the bounds.
#' @param dimXc Integer. Dimension of the common covariates (Xc).
#' @param values Numeric matrix. Unique values of covariates used to define subsets.
#' @param eta_d Numeric vector. Residualized independent variable for the subset.
#' @param Y Numeric vector. Outcome variable.
#' @param w_x Numeric vector. Sample weights for `eta_d`.
#' @param w_y Numeric vector. Sample weights for `Y`.
#' @param Xc_x0 Numeric matrix. Covariate subset for the independent variable.
#' @param Xc_y0 Numeric matrix. Covariate subset for the outcome variable.
#' @param Xnc Numeric matrix. Non-common covariates.
#' @param nu_d Numeric vector. Residuals for the independent variable.
#' @param nu_y Numeric vector. Residuals for the outcome variable.
#' @param dataset Integer (default 0). If 1, uses dataset-specific normalization for probability weights.
#'
#' @return A list of length 2:
#' \describe{
#'   \item{[[1]]}{1x2 numeric matrix. Lower and upper EMD-based bounds for the subset.}
#'   \item{[[2]]}{Numeric scalar. Weight/probability associated with the subset.}
#' }
#'
#' @details
#' The function selects the subset of observations corresponding to the value
#' indexed by `k` in `values`. It then computes the 1D EMD between the residualized
#' independent variable (`nu_d`) and the residualized outcome variable (`nu_y`),
#' separately for the positive and negative directions. This yields the upper and
#' lower bounds for the coefficient corresponding to that subset.
#'
#' It uses the `emd_1d_sorted` function for efficient 1D EMD computation.
#'
#' @examples
#' \dontrun{
#' # Example with simple numeric vectors
#' eta_d <- c(1.0, 2.0, 3.0, 4.0)
#' nu_d <- eta_d - mean(eta_d)
#' Y <- c(2.0, 3.0, 1.0, 5.0)
#' nu_y <- Y - mean(Y)
#' Xc_x0 <- matrix(c(1,1,2,2), ncol=1)
#' Xc_y0 <- matrix(c(1,1,2,2), ncol=1)
#' w_x <- rep(0.25,4)
#' w_y <- rep(0.25,4)
#' values <- matrix(c(1,2), ncol=1)
#' loop_EMD(1, dimXc=1, values, eta_d, Y, w_x, w_y, Xc_x0, Xc_y0, NULL, nu_d, nu_y)
#'}
#' @export
#'
loop_EMD <- function(k,dimXc,values,eta_d,Y,
                     w_x,w_y,Xc_x0,Xc_y0,
                     Xnc,nu_d,nu_y,dataset=0){



  boundsxc = matrix(NA,1,2)

  # sum(table(Xc_x))
  if(dimXc==1){
    val =values[k,]
    sel_x = (Xc_x0==val)
    sel_y = (Xc_y0==val)
  }else{
    val = t(as.matrix(values[k,]))
    sel_x = matrix(1,dim(Xc_x0)[1],1)
    sel_y = matrix(1,dim(Xc_y0)[1],1)
    for(ddd in 1:dimXc){
      sel_x =  sel_x & (Xc_x0[,ddd]==val[ddd])
      sel_y =  sel_y & (Xc_y0[,ddd]==val[ddd])
    }
    sel_x = matrix( sel_x,dim(Xc_x0)[1],1)
    sel_y = matrix( sel_y,dim(Xc_y0)[1],1)
  }


  nu_d0 = nu_d[sel_x]
  nu_y0 = nu_y[sel_y]

  prob <-  NA

  n_x = sum(sel_x)
  n_y = sum(sel_y)

  if(n_x >1 & n_y>1){
    ## sort both


    Umat = cbind( nu_d0,w_x[sel_x]/sum(w_x[sel_x]))
    Utilde_p = Umat[order(Umat[,1],decreasing=FALSE),]

    y_c_mat = cbind( nu_y0,w_y[sel_y]/sum(w_y[sel_y]))
    Ytilde =  y_c_mat[order(y_c_mat[,1],decreasing=FALSE),]


    value = emd_1d_sorted( Utilde_p[,2],  Ytilde[,2], Utilde_p[,1], Ytilde[,1], metric = 'sqeuclidean', p = 2.0)
    boundsxc[1,2] =   0.5*(sum(Utilde_p[,2]*Utilde_p[,1]^2) + sum(Ytilde[,2]*Ytilde[,1]^2) -   value$cost)

    Umat = cbind(-nu_d0,w_x[sel_x]/sum(w_x[sel_x]))
    Utilde_m = Umat[order(Umat[,1],decreasing=FALSE),]


    value = emd_1d_sorted(Utilde_m[,2],  Ytilde[,2], Utilde_m[,1], Ytilde[,1], metric = 'sqeuclidean', p = 2.0)
    boundsxc[1,1] =   0.5*(sum(Utilde_m[,2]*Utilde_m[,1]^2) + sum(Ytilde[,2]*Ytilde[,1]^2) -   value$cost)
    ##### version where integrating on P(w)=1
    if(dataset==1){
      prob <-  sum(w_y[sel_y]) #
    }else{
      prob <-  (sum(w_y[sel_y]) +  sum(w_x[sel_x]))/2 #
    }


  }


  output = vector("list")
  output[[1]] <- boundsxc
  output[[2]] <- prob

  return(output)

}

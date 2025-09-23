#' Bounds on the Best Linear Predictor (BLP) using Data Combination
#'
#' Computes bounds for regression coefficients provided in D'Haultfoeuille, Gaillac, Maurel (2024) <doi:10.48550/arXiv.2412.04816>.
#' It optionally allows incorporating conditioning variables, weights, using asymptotic normality results,
#' and information about the support of the regressors.
#'
#' @param Ldata A data frame containing the dependent variable (`out_var`) and covariates.
#' @param Rdata A data frame containing the covariates used for regression (`nc_var`, `c_var`, `w_var`).
#' @param out_var Name of the dependent variable in `Ldata`.
#' @param nc_var Character vector of outside regressors (non-commonly observed) (used in regression).
#' @param c_var Character vector of covariates commonly observed (default NULL).
#' @param w_var Character vector of auxiliary covariates, available in both datasets but not entering the regression  (default NULL).
#' @param w_x Optional weights for `Rdata` observations (default NULL, equal weights used).
#' @param w_y Optional weights for `Ldata` observations (default NULL, equal weights used).
#' @param nbCores Number of CPU cores to use for parallel computations (default 1).
#' @param full Logical, full computation option (default FALSE).
#' @param FW Logical, whether to compute Frisch-Waugh type bounds (default FALSE).
#' @param discrete Logical, whether to treat conditioning variables as discrete (default FALSE).
#' @param pt_est Logical, whether to return point estimates only (default FALSE).
#' @param ASN Logical, whether to compute asymptotic variance using influence function (default TRUE).
#' @param unchanged Logical, whether to keep all values when discretizing (default FALSE).
#' @param factor Logical, whether to treat W variables as factors (default FALSE).
#' @param K_sel Number of clusters if discrete clustering is used to partition the support of W (default 0, automatically selected).
#' @param dataset Integer flag indicating if we use propensity score weighting (default 0).
#' @param O_only Logical, whether to compute bounds only for outside regressors (default FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{hull_point}{Matrix of BLP hull points for each covariate.}
#'   \item{mat_varb_out_unc}{Matrix of upper and lower bounds without constraints.}
#'   \item{mat_varb_out_unc_asn}{Matrix of asymptotic variance estimates.}
#'   \item{output_influence}{List of influence function calculations for asymptotic variance.}
#' }
#'
#' @details
#' The function proceeds in several steps:
#' 1. Prepares matrices of covariates (`Xc_x`, `Xc_y`), outside regressors (`Xnc`), auxiliary variables (`Wc_x`, `Wc_y`).
#' 2. Calls `compute_BLP()` to compute hull points (BLP bounds) based on the input datasets.
#' 3. If `ASN` is TRUE, computes asymptotic variance using influence functions, optionally conditioning on `W`.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' Ldata = data.frame(Y = rnorm(100), C1 = rnorm(100))
#'Rdata = data.frame(NC1 = rnorm(100), C1 = rnorm(100))
#'
#'results <- regCombin_BLP(Ldata,Rdata,
#'                         out_var = "Y",
#'                         nc_var = c("NC1"),
#'                         c_var = c("C1"),
#'                         w_var=NULL,
#'                         w_x = NULL ,w_y= NULL,nbCores=1,
#'                         full=FALSE,FW=FALSE,discrete=TRUE,
#'                         pt_est=FALSE,
#'                         ASN=TRUE,unchanged =FALSE,factor=FALSE,K_sel= 0,
#'                         dataset=0, O_only =FALSE)
#'
#' }
#'
#' @export
#'
#'
#'
regCombin_BLP <- function(Ldata,Rdata,
                          out_var,nc_var,c_var = NULL,w_var =NULL,
                          w_x = NULL ,w_y= NULL,
                          nbCores=1,
                          full=FALSE,FW=FALSE,discrete=FALSE,
                          pt_est=FALSE,
                          ASN=TRUE,unchanged = FALSE,
                          factor=FALSE,K_sel= 0,
                          dataset=0,
                          O_only = FALSE){

  ###############
  exact=2
  expo = 0.2
  values=NULL
  n_x = dim(Rdata)[1]
  n_y = dim(Ldata)[1]

  if(n_x==n_y){
    same_size=TRUE
  }

  ###############
  if(is.null(w_var)){
    dimW =0
  }else{
    dimW = length(w_var)
  }
  dimXc = length(c_var)
  dimXnc = length(nc_var)
  dimW_all=dimW+dimXc

  ###############
  alpha=0.05
  q95 <-  function(x){quantile(x,1-alpha,na.rm=T)}
  q05 <-  function(x){quantile(x,alpha,na.rm=T)}

  ###############
  if(dimXc!=0 & dimW!=0){

    ### dataset 1
    Wc_x = as.matrix(Rdata[,w_var],dim(Rdata)[1],dimW)
    Xc_x = as.matrix(Rdata[,c_var],dim(Rdata)[1],dimXc)
    Xnc = as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)

    ### dataset 2
    Wc_y = as.matrix(Ldata[,w_var],dim(Ldata)[1],dimW)
    Xc_y = as.matrix(Ldata[,c_var],dim(Ldata)[1],dimXc)
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)

  }else if(dimXc!=0 & dimW==0){

    ### dataset 1
    Wc_x = NULL
    Xc_x = as.matrix(Rdata[,c_var],dim(Rdata)[1],dimXc)
    Xnc = as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)

    ### dataset 2
    Wc_y = NULL
    Xc_y = as.matrix(Ldata[,c_var],dim(Ldata)[1],dimXc)
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)

  }else if(dimXc==0 & dimW!=0){

    ### dataset 1
    Wc_x = as.matrix(Rdata[,w_var],dim(Rdata)[1],dimW)
    Xc_x = NULL
    Xnc = as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)

    ### dataset 2
    Wc_y = as.matrix(Ldata[,w_var],dim(Ldata)[1],dimW)
    Xc_y = NULL
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)


  }else{

    ### dataset 1
    Wc_x = NULL
    Xc_x = NULL
    Xnc = as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)
    ### dataset 2
    Wc_y = NULL
    Xc_y = NULL
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)

  }

  output  <- vector("list")
  hull_point=NULL
  hull_sharp=NULL
  sample1 =NULL
  dataX =NULL
  dataP =NULL

  ############### Compute the point estimate
  hull_point <- compute_BLP(sample1 =sample1,Xc_x,Xnc,Xc_y,Wc_x,Wc_y,Y,
                            w_x,w_y,
                            dimXc,dimW,dimXnc,
                            values,full,FW,
                            discrete,unchanged,rule=1,factor,K_sel,
                            dataset,  O_only )

  Bsamp=500

  mat_varb_out_unc = matrix(0,dim(hull_point)[1],2)
  mat_varb_out_unc_asn = matrix(0,dim(hull_point)[1],2)

  if(Bsamp >0 & pt_est==FALSE){


    if(ASN==TRUE){


      index_infl = 1
      output_influence <- vector("list")

      #### based on ASN ####
      ## compute standard error
      sig =  matrix(0,dim(hull_point)[1],2)
      dimW_all = dimXc+dimW

      lambda = n_y/(n_x+n_y)

      if( dimW_all==0 ){

        if(is.null(w_x)){
          w_x = matrix(1,n_x,1)/n_x
        }
        if(is.null(w_y)){
          w_y = matrix(1,n_y,1)/n_y
        }



        for(p in 1:(dim(hull_point)[1])){

          for(j in 1:2){
            if(j == 1){
              ind = -1
            }else{
              ind = 1
            }



            ######
            reg = cbind(1,Xnc[,-c(p)])
            gamma_hat <- solve(t(reg) %*% reg) %*% t(reg) %*% (ind*Xnc[,p])
            eta_hat <- (ind*Xnc[,p]) - reg %*% gamma_hat
            # eta_hat <- eta_hat

            sU = sort(eta_hat,decreasing = FALSE)
            sY = sort(Y,decreasing = FALSE)

            eval_Y <- function(j){return( ceil(j*length(Y)/length(eta_hat)))}
            eval_eta <- function(j){return( ceil(j*length(eta_hat)/length(Y)))}

            # eval_eta(1:length(Y))
            ###
            psi_4 =  c(0,cumsum(sU[eval_eta(2:length(sY))]*diff(sY)))
            d_j3  =  c(0,cumsum(sY[eval_Y(2:length(sU))]*diff(sU)))

            # c(0,cumsum(sY[-c(1)]*diff(sU)))
            ### Sep. 4. ##########################################

            ranks = cbind(seq(1,n_x),eta_hat)
            ranks = cbind(seq(1,n_x),ranks[order(eta_hat,decreasing = FALSE),])
            ### function which to the order of  eta_hat associate the initial one (x,y)
            ff <- approxfun(ranks[,2],ranks[,1])

            A = colSums(reg*(sort(Y)[pmin(ff(seq(1,n_x)),n_x)]))/n_x
            B = ginv(t(reg)%*%reg/n_x)
            AB = A%*%B

            if((dim(AB)[1])>1){
              infl2 = AB%*%t(reg*(eta_hat%*%matrix(1,1,dim(reg)[2]-1)))
              psi_2 = - t(infl2)
            }else{
              infl2 = as.numeric(AB)*eta_hat
              psi_2 = - infl2
            }

            psi_3 <-   d_j3[ff(seq(1,n_x))]

            psi_1 <- -  (eta_hat^2*w_x-sum(eta_hat^2*w_x))*ind*hull_point[p,j]
            psi_x <- psi_1 + psi_2 + psi_3

            ves0 <-  (lambda*(sum(psi_x^2*w_x) - sum(psi_x*w_x)^2)+
                        (1-lambda)*(sum(psi_4^2*w_y) - sum(psi_4*w_y)^2))

            ves <- ves0/c(sum(eta_hat^2*w_x))^2

            sig[p,j] = sqrt(ves)/sqrt(lambda*n_x)


            output_infl <- vector("list")
            output_infl[["lambda"]] <- lambda
            output_infl[["psi_x"]] <-psi_x
            output_infl[["psi_y"]] <- psi_4
            output_infl[["w_x"]] <- w_x
            output_infl[["w_y"]] <- w_y
            output_infl[["psi_n"]] <- NULL
            output_infl[["denom"]] <- c(mean(eta_hat^2))^2
            output_infl[["sig"]] <-  sig[p,j]

            output_influence[[index_infl]] <-output_infl
            index_infl = index_infl+1

          }


        }

      }else{  ##### conditional on W.



        direct= FALSE
        if(n_x==n_y & is.null(w_x) & is.null(w_y) ){
          direct= TRUE
        }

        if(is.null(w_x)){
          w_x = matrix(1,n_x,1)/n_x
        }
        if(is.null(w_y)){
          w_y = matrix(1,n_y,1)/n_y
        }


        Xnc_c= Xnc
        Y_c = Y

        n_x = dim(Xnc )[1]
        n_y= dim(Y)[1]

        W_all_y =cbind(Xc_y,Wc_y)
        W_all_x =cbind(Xc_x,Wc_x)

        c_var0 = paste0("X",1:dimW_all)

        colnames(W_all_y) <-  c_var0
        colnames(W_all_x) <-  c_var0


        W_all_x1=cbind(rep(1,n_x),W_all_x)
        W_all_y1=cbind(rep(1,n_y),W_all_y)
        W_all_all = rbind(W_all_x1, W_all_y1)

        #####  split each coordinantes support for g(w) ####### => get discrete for each.
        ##### discretize if needed #########################################
        # j=1
        # j=1
        Xcurr = rbind(W_all_y,W_all_x)
        colnames(Xcurr) <- c_var0

        if(is.null(w_x)){
          w_x = matrix(1,n_x,1)/n_x
        }
        if(is.null(w_y)){
          w_y = matrix(1,n_y,1)/n_y
        }

        w_xy = c(w_x *  n_x, w_y *  n_y) / ( n_x +  n_y)


        ###################################### FW ##################################################
        if(FW==TRUE){

          if(O_only == TRUE){
            bounds_dim =  dim(Xnc)[2]
          }else{
            bounds_dim =  dim(Xnc)[2] + dimXc
          }


          bounds= matrix(0,bounds_dim,2)

          #### !! here Xc_x and not the discrete one
          if(dimXc!=0){
            Xnc_c= cbind(matrix(1,n_x,1), Xnc,Xc_x )  ### to build the A matrix.
          }else{
            Xnc_c= cbind(matrix(1,n_x,1), Xnc )  ### to build the A matrix.

          }

          #### when using the discretized version as W
          if(factor==TRUE){

            values =  W_all_y[!duplicated( W_all_y[,c_var0]),c_var0]
            # Xc_x0 = matrix(0, dim(Rdatad)[1], 1)
            Xc_y0 = matrix(0, dim(W_all_y)[1], 1)
            ind = 1
            # j=2
            for (j in dim(values)[1]) {
              if (dimW_all == 1) {
                val = values[j,]
                sel_y = (W_all_y[, c_var0] == val)
              } else {
                val = t(as.matrix(values[j,]))
                # sel_x = matrix(1, dim(W_all_y)[1], 1)
                sel_y = matrix(1, dim(W_all_y)[1], 1)
                for (ddd in 1:dimW_all) {
                  # sel_x = sel_x & (Rdatad[,  c_var0[ddd]] == val[ddd])
                  sel_y = sel_y & (W_all_y[,  c_var0[ddd]] == val[ddd])
                }
                # sum(sel_x)
                # sel_x = matrix(sel_x, dim(Rdatad)[1], 1)
                sel_y = matrix(sel_y, dim(W_all_y)[1], 1)
              }
              # Xc_x0[sel_x] = ind
              Xc_y0[sel_y] = ind
              ind = ind + 1
            }



            dataY = as.data.frame(cbind(Y,  Xc_y0, w_y))
            colnames(dataY) <- c("Y", "W", "w_y")
            dataY$W <- as.factor(dataY$W)


            dataY <- dataY %>% group_by(W) %>% mutate( nu_y = Y - sum(w_y*Y)/sum(w_y),
                                                       projy = sum(w_y*Y)/sum(w_y))
            nu_y = dataY$nu_y

            # dataX0 = as.data.frame(cbind(Xc_x0,w_x))
            # colnames(dataX0) <- c("W","w_x")

            if(dataset==1){
              # dataP = data.frame( D=c(rep(1,n_x), rep(0,n_y)), w= c(w_x,w_y)/2)
              dataP <- as.data.frame(cbind(c(rep(0,n_x), rep(1,n_y)),  c(w_x,w_y), rbind(W_all_x,W_all_y)))
              colnames(dataP) <- c("D","w","W")
              # data
              ##poids ?
              # matX = as.matrix(dataP[,-c(1)],n_x+n_y,1+length(c_var0))
              # colnames(matX) <- c("w",c_var0)
              # fit_w <- glm(    dataP[,"D"] ~ W_all_all  , family="binomial",weights = dataP[,"w"] )
              fit_w <- glm( "D ~ W", data=dataP, family="binomial",weights = as.numeric(dataP$w ))
              p_w1 <- predict(fit_w,data= rbind(W_all_x,W_all_y), type="response")
              p_w1 <-   p_w1[1:n_x]/(1-  p_w1[1:n_x])
              p_w1 <-   p_w1*w_x/sum(  p_w1*w_x)
            }else{
              p_w1 = w_x
              # dataP = data.frame( D=c(rep(1,n_x), rep(0,n_y)), w= c(w_x,w_y)/2)
              # dataP <- as.data.frame(cbind(c(rep(1,n_x), rep(0,n_y)),  c(w_x,w_y)/2 , rbind(Xc_x0,Xc_y0)))
              # colnames(dataP) <- c("D","w","W")
              # # colnames(dataP) <- c("D","w",c_var0)
            }



          }else{


            dataY = as.data.frame(cbind(Y, W_all_y,w_y))
            colnames(dataY) <- c("Y",c_var0,"w_y")

            # fit_y = lm(Y ~ W_all_y  , weights=w_y)
            fit_y = lm(paste0("Y ~ ", paste0(c_var0,collapse = "+"))  ,
                       data =    dataY , weights=as.numeric(w_y ))

            nu_y  = fit_y$residuals
            # summary(fit_y)
            delta_y=fit_y$coefficients

            if(dataset==1){
              # dataP = data.frame( D=c(rep(1,n_x), rep(0,n_y)), w= c(w_x,w_y)/2)
              dataP <- as.data.frame(cbind(c(rep(0,n_x), rep(1,n_y)),  c(w_x,w_y) , W_all_all))
              colnames(dataP) <- c("D","w","1",c_var0)
              # data
              ##poids ?
              # matX = as.matrix(dataP[,-c(1)],n_x+n_y,1+length(c_var0))
              # colnames(matX) <- c("w",c_var0)
              # fit_w <- glm(    dataP[,"D"] ~ W_all_all  , family="binomial",weights = dataP[,"w"] )
              fit_w <- glm( paste0("D ~ ", paste0(c_var0,collapse = "+")) , data=dataP,
                            family="binomial",weights = dataP$w )
              p_w1 <- predict(fit_w,data=W_all_all, type="response")
              p_w1 <-   p_w1[1:n_x]/(1-  p_w1[1:n_x])
              p_w1 <-   p_w1*w_x/sum(  p_w1*w_x)
            }else{
              p_w1 = w_x
              # dataP <- as.data.frame(cbind(c(rep(1,n_x), rep(0,n_y)),  c(w_x,w_y)/2 , W_all_all))
              # colnames(dataP) <- c("D","w","1",c_var0)
            }

          }

          d=2


          for (d in 1:(dim(hull_point)[1])){  ### XXXX add a for loop


            #### regression with weights P(W), last case of 2.3.2
            # fit_d0 = lm(Xnc_c[,(d+1)] ~ Xnc_c[,-c(d+1)] -1  , weights= p_w1)
            # fit_d0 = lm(Xnc_c[,(d+1)] ~ Xnc_c[,-c(d+1)] -1  , weights= w_x)
            fit_d0 = lm(Xnc_c[, (d + 1)] ~ Xnc_c[, -c(d + 1)] - 1  , weights =  as.numeric(p_w1))

            eta_d  = fit_d0$residuals



            if(factor==TRUE){




              values =  W_all_x[!duplicated( W_all_x[,c_var0]),c_var0]
              # Xc_x0 = matrix(0, dim(Rdatad)[1], 1)
              Xc_x0 = matrix(0, dim(W_all_x)[1], 1)
              ind = 1
              # j=2
              for (j in dim(values)[1]) {
                if (dimW_all == 1) {
                  val = values[j,]
                  # sel_x = (W_all_y[, c_var0] == val)
                  sel_y = (W_all_x[, c_var0] == val)
                } else {
                  val = t(as.matrix(values[j,]))
                  sel_x = matrix(1, dim(W_all_x)[1], 1)
                  # sel_y = matrix(1, dim(W_all_y)[1], 1)
                  for (ddd in 1:dimW_all) {
                    sel_x = sel_x & (W_all_x[,  c_var0[ddd]] == val[ddd])
                    # sel_y = sel_y & (W_all_y[,  c_var0[ddd]] == val[ddd])
                  }
                  # sum(sel_x)
                  sel_x = matrix(sel_x, dim(W_all_x)[1], 1)
                  # sel_y = matrix(sel_y, dim(W_all_y)[1], 1)
                }
                Xc_x0[sel_x] = ind
                # Xc_y0[sel_y] = ind
                ind = ind + 1
              }


              dataX = as.data.frame(cbind(eta_d , Xc_x0, w_x,     p_w1))
              colnames(dataX) <- c("eta_d", "W", "w_x", "p_w")
              dataX$W <- as.factor(dataX$W)

              if(dataset==1){
                dataX <- dataX %>% group_by(W) %>% mutate( nu_d = eta_d - sum(p_w*eta_d)/sum(p_w))
                nu_d = dataX$nu_d
                dataXs <- dataX[,c("eta_d","W","w_x","nu_d","p_w")] %>% group_by(W) %>% summarise( projd =  sum(p_w*eta_d)/sum(p_w),
                                                                                                   projd2 =  sum(p_w*eta_d^2)/sum(p_w))
              }else{
                dataX <- dataX %>% group_by(W) %>% mutate( nu_d = eta_d - sum(w_x*eta_d)/sum(w_x))
                nu_d = dataX$nu_d
                dataXs <- dataX[,c("eta_d","W","w_x","nu_d")] %>% group_by(W) %>% summarise( projd =  sum(w_x*eta_d)/sum(w_x),
                                                                                             projd2 =  sum(w_x*eta_d^2)/sum(w_x))
              }

              dataYs <- dataY[,c("Y","W","w_y","nu_y")] %>% group_by(W) %>% summarise( projy =  sum(w_y*Y)/sum(w_y))
              dataXss = inner_join(dataXs,  dataYs, by="W")

              dataWsy = inner_join(dataY[,c("W","w_y")],  dataXss, by="W")
              dataWsx = inner_join(dataX[,c("W","w_x","p_w")],  dataXss, by="W")


              # w_xy = c(w_x*dim(dataX)[1],w_y*dim(dataY)[1])/(dim(dataX)[1] + dim(dataY)[1])

              if(dataset==1){
                # off = sum(dataWsy$projy* dataWsy$projd *t(dataWsy$w_y),na.rm=T)
                # denom = sum(dataWsy$projd2 *t(dataWsy$w_y),na.rm=T)

                off = sum(dataWsy$projy* dataWsy$projd *t(dataWsy$w_y),na.rm=T)
                denom = sum(dataWsy$projd2 *t(dataWsy$w_y),na.rm=T)
                # denom = sum(dataWsy$projd2 *t(dataWsy$w_y),na.rm=T)
              }else{
                # off = sum(dataWsy$projy* dataWsy$projd *t(dataWsy$w_y),na.rm=T)
                off = sum(dataWsx$projy* dataWsx$projd * w_xy[1:dim(dataX)[1]],na.rm=T) +
                  sum(dataWsy$projy* dataWsy$projd * w_xy[(dim(dataX)[1]+1):(dim(dataX)[1] + dim(dataY)[1])],na.rm=T)

                denom = sum(dataWsx$projd2 * w_xy[1:dim(dataX)[1]],na.rm=T)+
                  sum(dataWsy$projd2*  w_xy[(dim(dataX)[1]+1):(dim(dataX)[1] + dim(dataY)[1])],na.rm=T)
              }



            }else{

              dataX = as.data.frame(cbind(eta_d, W_all_x,p_w1))
              colnames(dataX) <- c("eta_d",c_var0,"w_x")
              fit_d = lm(paste0("eta_d ~ ", paste0(c_var0,collapse = "+"))  , data =    dataX,
                         weights=dataX$w_x)

              nu_d  = fit_d$residuals
              # summary(fit_y)
              delta_d=fit_d$coefficients

              if(dataset==1){
                V = t(W_all_y1)%*%(W_all_y1*(w_y%*%matrix(1,1,1+dim(W_all_y)[2])))
              }else{
                V = t(W_all_all)%*%(W_all_all* w_xy%*%matrix(1,1,dim(W_all_all)[2]))
              }

              off = c((t(delta_d)%*%V)%*%delta_y)

              if(dataset==1){
                denom = sum(eta_d^2*p_w1 ,na.rm=T)
              }else{
                denom = sum(eta_d^2*w_x ,na.rm=T)
              }

            }

            ####
            j=1
            for(j in 1:2){

              if(j == 1){
                ind = -1

              }else{
                ind = 1

              }

              eta_hat = eta_d *ind
              nu_hat = nu_d *ind



              ###### computation of psi 1
              # if(factor==TRUE){
              #   psi1_1 =  ind*dataX_1$nu_d*dataX_1$projy
              #   A = colSums(Xnc_c[,-c(d+1)]*( dataX_1$projy))/n_x
              # }else{
              predy <- predict(  fit_y,newdata= as.data.frame(W_all_x))
              predy = fit_y$coefficients%*%t(W_all_x1)
              psi1_1 =   t(nu_hat* predy)
              # A=fit_y$coefficients%*%t(mat1

              # preddx <- predict(  fit_d,newdata= as.data.frame(W_all_x))

              # A = colSums(Xnc_c[,-c(d+1)]*predy)/n_x
              # A = colSums(Xnc_c[,-c(d+1)]*((matrix(predy,dim(Xnc_c)[1],1)*w_x)%*%matrix(1,1,dim(Xnc_c)[2]-1)))#/n_x
              # }


              mat1= t(W_all_x1)%*%(Xnc_c[,-c(d+1)]*(p_w1%*%matrix(1,1,dim(Xnc_c)[2]-1)))
              B = ginv(t(Xnc_c[,-c(d+1)])%*%(Xnc_c[,-c(d+1)]*(p_w1%*%matrix(1,1,dim(Xnc_c)[2]-1))))
              # B = ginv(t(Xnc_c[,-c(d+1)])%*%(Xnc_c[,-c(d+1)])/n_x)
              AB = mat1%*%B



              if((dim(AB)[1])>1){
                infl1_2 = AB%*%t(Xnc_c[,-c(d+1)]*(eta_hat%*%matrix(1,1,dim(Xnc_c)[2]-1)))
                infl1_2 = fit_y$coefficients%*%infl1_2
                psi1_2 = - t(infl1_2)
              }else{
                infl1_2 = as.numeric(AB)*eta_hat
                psi1_2 = -infl1_2
              }



              psi_1 = psi1_1 + psi1_2



              if(dataset==1){


                psi_n_x = t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*    p_w1 -
                  sum(t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*   p_w1)
                psi_n_y = t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y -
                  sum(t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y)

              }else{



                psi_n_x = t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*   w_x -
                  sum(t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*   w_x)
                psi_n_y = t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y -
                  sum(t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y)

              }

              predd <- predict(  fit_d,newdata= as.data.frame(W_all_y))

              psi_2 =  ind* predd *nu_y

              projh <- matrix(NA,1,dimW_all+1)
              projh_5 <- matrix(NA,1,dimW_all+1)
              # }

              psi_4 <- matrix(NA,n_x,1)
              psi_6 <- matrix(NA,n_y,1)

              psi_3 <- matrix(NA,n_x,1)
              psi_5 <- matrix(NA,n_y,1)


              Bw = ginv(t(W_all_y1)%*%as.matrix(W_all_y1)/n_y)

              fd=ecdf( nu_hat)
              fy=ecdf( nu_y)


              k=1
              p3=quantile(nu_y,fd(nu_hat))
              if(length(c_var0)>1){
                projh[k,] <- colMeans(p3*W_all_x1)

              }else{
                projh[k,] <- mean(p3*W_all_x1)

              }

              p5 = quantile(nu_hat,fy(nu_y))
              if(length(c_var0)>1){
                projh_5[k,] <- colMeans((p5%*%matrix(1,1,dim(W_all_y1)[2]))*W_all_y1)
              }else{
                projh_5[k,] <- mean(p5*W_all_y1)

              }

              sU = sort( nu_hat,decreasing = FALSE)
              sY = sort(nu_y,decreasing = FALSE)
              ###
              eval_Y <- function(j){return( j/length(nu_y))}
              eval_eta <- function(j){return( j/length(nu_d))}
              # length(diff(sY))
              d_i6k = c(0,cumsum(quantile(nu_hat,eval_Y(2:n_y))*diff(sY)))
              d_j4k = c(0,cumsum(quantile(nu_y,eval_eta(2:n_x))*diff(sU)))

              ranks_d = cbind(seq(1, n_x),nu_hat)
              ranks_d = cbind(seq(1, n_x),ranks_d[order(nu_hat,decreasing = FALSE),])
              ffk <- approxfun(ranks_d[,2],ranks_d[,1])

              ranks_y = cbind(seq(1, n_y),nu_y)
              ranks_y = cbind(seq(1, n_y),ranks_y[order(nu_y,decreasing = FALSE),])
              ffk_y <- approxfun(ranks_y[,2],ranks_y[,1])

              psi_6<-  d_i6k[ffk_y(seq(1,n_y))]

              psi_4 <-  d_j4k[ffk(seq(1,n_x))]


              ABw = projh[k,]%*%Bw
              psi_3  =  ABw%*%t(W_all_x1*nu_hat)

              ma= ABw%*%AB
              if(length(ma)>1){
                infl3_2 =ma%*%t(Xnc_c[,-c(d+1)]*(eta_hat%*%matrix(1,1,dim(Xnc_c)[2]-1)))

                psi_3 = psi_3  - infl3_2
              }else{
                infl3_2 =as.numeric(ma)*(eta_hat)

                psi_3 = psi_3  - infl3_2
              }
              psi_3 = - t(psi_3)


              ABw5= projh_5[k,]%*%Bw
              psi_5 =  - t(ABw5%*%t(W_all_y1*nu_y))




              if(dataset==1){
                psi_10 <-  - (eta_hat^2-denom)*ind*hull_point[d,j]
                psi_x<-  psi_1 + psi_3  + psi_4  + psi_10 +  psi_n_x
                psi_y  <-  psi_2 + psi_5  + psi_6  +  psi_n_y
                ves <-lambda * sum(psi_x^2*p_w1) + (1 - lambda) *sum(psi_y ^ 2*w_y)
                ves <- ves -   (lambda * sum(psi_x*p_w1)^ 2  +  (1 - lambda) * sum(psi_y*w_y)^2)
              }else{
                psi_10 <-  - (eta_hat^2-denom)*ind*hull_point[d,j]
                psi_x <-  psi_1 + psi_3  + psi_4  + psi_10 +  psi_n_x
                psi_y  <-  psi_2 + psi_5  + psi_6 +  psi_n_y
                ves <-lambda * sum(psi_x^2*w_x) + (1 - lambda) * sum(psi_y^2*w_y) # +
                # lambda*(1-lambda)*sum(psi_n^2* w_xy)
                ves <- ves -   (lambda * sum(psi_x*w_x)^ 2  +
                                  (1 - lambda) * sum(psi_y*w_y)^2) #+ lambda*(1-lambda)*sum(psi_n* w_xy)^2)
              }



              ves0 <-    ves / c(denom) ^ 2


              sig[d,j] = sqrt(ves0)/sqrt(lambda*n_x)



              output_infl <- vector("list")
              output_infl[["lambda"]] <- lambda
              output_infl[["psi_x"]] <- psi_x
              output_infl[["psi_y"]] <- psi_y
              output_infl[["w_x"]] <- w_x
              output_infl[["w_y"]] <- w_y
              # output_infl[["psi_n"]] <- psi_n
              output_infl[["denom"]] <- denom
              output_infl[["sig"]] <-  sig[d,j]

              output_influence[[index_infl]] <-output_infl
              index_infl = index_infl+1

            }

          }







        # }else if(discrete == TRUE && exact<2) {



        }else if(discrete == TRUE && exact==2) {
          ## no Xc projection for this one, already done.

          Xnc_c = Xnc
          Y_c = Y
          Xc_y_c = cbind(rep(1, dim(Xc_y)[1], 1), Xc_y, Wc_y)
          Xc_x_c = cbind(rep(1, dim(Xc_x)[1], 1), Xc_x, Wc_x)

          n_x = dim(Xnc)[1]
          n_y = dim(Y)[1]


          W_all_y = cbind(Xc_y, Wc_y)
          W_all_x = cbind(Xc_x, Wc_x)

          c_var0 = paste0("X", 1:dimW_all)

          colnames(W_all_y) <-  c_var0
          colnames(W_all_x) <-  c_var0


          W_all_x1 = cbind(rep(1, n_x), W_all_x)
          W_all_y1 = cbind(rep(1, n_y), W_all_y)
          W_all_all = rbind(W_all_x1, W_all_y1)
          # w_xy = c(w_x * dim(dataX)[1], w_y * dim(dataY)[1]) / (dim(dataX)[1] + dim(dataY)[1])


          #####  split each coordinantes support for g(w) ####### => get discrete for each.
          ##### discretize if needed #########################################
          # j=1
          Xcurr = rbind(W_all_y, W_all_x)
          dim(Xcurr)
          colnames(Xcurr) <- c_var0

          if (K_sel > 0) {
            K = K_sel
          } else{
            K <- max(2, floor(min(n_y,n_x)^expo))
            # K <- max(2, floor(min(n_y,n_x)^0.5))
          }




          # if (dim(values)[1] > 2) {


          if(O_only == TRUE){
            bounds_dim =  dim(Xnc)[2]
          }else{
            bounds_dim =  dim(Xnc)[2] + dimXc
          }

          bounds = matrix(0,   bounds_dim, 2)

          #### !! here Xc_x and not the discrete one
          if (dimXc != 0) {
            Xnc_c = cbind(matrix(1, n_x, 1), Xnc , Xc_x)  ### to build the A matrix.
          } else{
            Xnc_c = cbind(matrix(1, n_x, 1), Xnc)  ### to build the A matrix.
          }





          #### when using the discretized version as W
          if (factor == TRUE) {
            dataY = as.data.frame(cbind(Y, Xc_y0, w_y))
            colnames(dataY) <- c("Y", "W", "w_y")
            dataY$W <- as.factor(dataY$W)

            dataY <-dataY %>% group_by(W) %>% mutate(
              nu_y = Y - sum(w_y * Y) / sum(w_y),
              projy = sum(w_y * Y) / sum(w_y)
            )
            nu_y = dataY$nu_y


            if (dataset == 1) {
              # dataP = data.frame( D=c(rep(1,n_x), rep(0,n_y)), w= c(w_x,w_y)/2)
              dataP <-as.data.frame(cbind(c(rep( 0, n_x ), rep(1, n_y)),  c(w_x,w_y) , rbind(Xc_x0, Xc_y0)))
              colnames(dataP) <- c("D", "w", "W")
              # data
              ##poids ?
              # matX = as.matrix(dataP[,-c(1)],n_x+n_y,1+length(c_var0))
              # colnames(matX) <- c("w",c_var0)
              # fit_w <- glm(    dataP[,"D"] ~ W_all_all  , family="binomial",weights = dataP[,"w"] )
              fit_w <-
                glm(
                  "D ~ W",
                  data = dataP,
                  family = "binomial",
                  weights = as.numeric(dataP$w)
                )
              p_w1 <- predict(fit_w, data = W_all_all, type = "response")
              p_w1 <-   p_w1[1:n_x] / (1 -  p_w1[1:n_x])
              p_w1 <-   p_w1 * w_x / sum(p_w1 * w_x)
            } else{
              p_w1 = w_x

            }


          } else{

            dataY = as.data.frame(cbind(Y, W_all_y, w_y))
            colnames(dataY) <- c("Y", c_var0, "w_y")

            # fit_y = lm(Y ~ W_all_y  , weights=w_y)
            fit_y = lm(paste0("Y ~ ", paste0(c_var0, collapse = "+"))  ,
                       weights = as.numeric(w_y),
                       data =    dataY)

            nu_y  = fit_y$residuals
            # summary(fit_y)
            delta_y = fit_y$coefficients



            if (dataset == 1) {
              # dataP = data.frame( D=c(rep(1,n_x), rep(0,n_y)), w= c(w_x,w_y)/2)
              dataP <-
                as.data.frame(cbind(c(rep(0, n_x ), rep( 1, n_y)),   c(w_x,w_y) , W_all_all))
              colnames(dataP) <- c("D", "w", "1", c_var0)
              # data
              ##poids ?
              # matX = as.matrix(dataP[,-c(1)],n_x+n_y,1+length(c_var0))
              # colnames(matX) <- c("w",c_var0)
              # fit_w <- glm(    dataP[,"D"] ~ W_all_all  , family="binomial",weights = dataP[,"w"] )
              fit_w <-
                glm(
                  paste0("D ~ ", paste0(c_var0, collapse = "+")) ,
                  data = dataP,
                  family = "binomial",
                  weights = dataP$w
                )
              p_w1 <- predict(fit_w, data = W_all_all, type = "response")
              p_w1 <-   p_w1[1:n_x] / (1 -  p_w1[1:n_x])
              p_w1 <-   p_w1 * w_x / sum(p_w1 * w_x)
            } else{
              p_w1 = w_x

            }

          }


          ######### adapted clustering


          abs_nu_y = abs(nu_y)

          ##### Step 1
          if (factor == TRUE) {
            data_abs = as.data.frame(cbind(abs_nu_y, Xc_y0, w_y))
            colnames(data_abs) <- c("anu_y", "W", "w_y")
            data_abs$W <- as.factor(data_abs$W)

            data_abs <-data_abs %>% group_by(W) %>% mutate(
              nu_abs_nu_y = abs_nu_y - sum(w_y * abs_nu_y) / sum(w_y),
              proj_abs_nu_y = sum(w_y * abs_nu_y) / sum(w_y)
            )
            # nu_y = dataY$nu_y




          } else{

            data_abs = as.data.frame(cbind(abs_nu_y, W_all_y, w_y))
            colnames(data_abs) <- c("anu_y", c_var0, "w_y")

            # fit_y = lm(Y ~ W_all_y  , weights=w_y)
            fit_abs = lm(paste0("anu_y ~ ", paste0(c_var0, collapse = "+"))  ,
                         weights = as.numeric(w_y),
                         data =    data_abs)


          }

          #############################################################################

          d = 1
          for (d in 1:bounds_dim) {
            ### XXXX add a for loop


            W_all_y_disc = as.matrix(W_all_y[,1])
            W_all_x_disc =  as.matrix(W_all_x[,1])

            #### regression with weights P(W), last case of 2.3.2
            fit_d0 = lm(Xnc_c[, (d + 1)] ~ Xnc_c[, -c(d + 1)] - 1  , weights =  as.numeric(p_w1))

            eta_d  = fit_d0$residuals


            if (factor == TRUE) {
              dataX = as.data.frame(cbind(eta_d , Xc_x0, w_x,     p_w1))
              colnames(dataX) <- c("eta_d", "W", "w_x", "p_w")
              dataX$W <- as.factor(dataX$W)

              if (dataset == 1) {
                dataX <- dataX %>% group_by(W) %>% mutate(nu_d = eta_d - sum(p_w * eta_d)/sum(p_w))
                nu_d = dataX$nu_d
                dataXs <- dataX[, c("eta_d", "W", "w_x", "nu_d", "p_w")] %>% group_by(W) %>% summarise(
                  projd =  sum(p_w * eta_d) / sum(p_w),
                  projd2 =  sum(p_w *eta_d ^ 2) / sum(p_w)
                )


              } else{
                dataX <-
                  dataX %>% group_by(W) %>% mutate(nu_d = eta_d - sum(w_x * eta_d) / sum(w_x))
                nu_d = dataX$nu_d
                dataXs <-
                  dataX[, c("eta_d", "W", "w_x", "nu_d")] %>% group_by(W) %>% summarise(
                    projd =  sum(w_x * eta_d) / sum(w_x),
                    projd2 =  sum(w_x *
                                    eta_d ^ 2) / sum(w_x)
                  )
              }

              dataYs <- dataY[, c("Y", "W", "w_y", "nu_y")] %>% group_by(W) %>% summarise(projy =  sum(w_y *
                                                                                                         Y) / sum(w_y))
              dataXss = inner_join(dataXs,  dataYs, by = "W")

              dataWsy = inner_join(dataY[, c("W", "w_y")],  dataXss, by =
                                     "W")
              dataWsx = inner_join(dataX[, c("W", "w_x", "p_w")],  dataXss, by =
                                     "W")


              if (dataset == 1) {

                off = sum(dataWsy$projy * dataWsy$projd * t(dataWsy$w_y),
                          na.rm = T)
                denom = sum(dataWsy$projd2 * t(dataWsy$w_y), na.rm = T)

              } else{

                off = sum(dataWsx$projy * dataWsx$projd * w_xy[1:dim(dataX)[1]], na.rm =T) +
                  sum(dataWsy$projy * dataWsy$projd * w_xy[(dim(dataX)[1] + 1):(dim(dataX)[1] + dim(dataY)[1])], na.rm = T)

                denom = sum(dataWsx$projd2 * w_xy[1:dim(dataX)[1]], na.rm= T) +
                  sum(dataWsy$projd2 *  w_xy[(dim(dataX)[1] + 1):(dim(dataX)[1] + dim(dataY)[1])], na.rm = T)
              }



            } else{

              dataX = as.data.frame(cbind(eta_d, W_all_x, p_w1))
              colnames(dataX) <- c("eta_d", c_var0, "w_x")

              fit_d = lm(
                paste0("eta_d ~ ", paste0(c_var0, collapse = "+"))  ,
                data =    dataX,
                weights = dataX$w_x
              )

              nu_d  = fit_d$residuals
              # summary(fit_y)
              delta_d = fit_d$coefficients

              if (dataset == 1) {
                # V = t(cbind(rep(1,n_y),W_all_y))%*%(cbind(rep(1,n_y),W_all_y)*(w_y%*%matrix(1,1,1+dim(W_all_y)[2])))
                V = t(W_all_y1) %*% (W_all_y1 * (w_y %*% matrix(1, 1, 1 +dim(W_all_y)[2])))
              } else{
                V = t(W_all_all) %*% (W_all_all *    w_xy %*% matrix(1, 1, dim(W_all_all)[2]) )
              }
              off = c((t(delta_d) %*% V) %*% delta_y)


              if (dataset == 1) {
                denom = sum(eta_d ^ 2 * p_w1 , na.rm = T)
              } else{
                denom = sum(eta_d ^ 2 * w_x , na.rm = T)
              }
              # }

            }


            ###################################################################3
            ######### adapted clustering
            # if(exact==2){

            abs_nu_d = abs( nu_d)

            ##### Step 1
            if (factor == TRUE) {
              data_abs_d = as.data.frame(cbind(abs_nu_d, W_all_x, w_x))
              colnames(data_abs_d) <- c("anu_d", "W", "w_y")
              data_abs_d$W <- as.factor(data_abs_d$W)

              data_abs_d <-data_abs_d %>% group_by(W) %>% mutate(
                nu_abs_nu_d_x = abs_nu_d - sum(w_x * abs_nu_d) / sum(w_d),
                proj_abs_nu_d_x = sum(w_d * abs_nu_d) / sum(w_d)
              )

            } else{

              data_abs_d = as.data.frame(cbind(abs_nu_d, W_all_x, w_x))
              colnames(data_abs_d) <- c("anu_d", c_var0, "w_x")

              # fit_y = lm(Y ~ W_all_y  , weights=w_y)
              fit_abs_d = lm(paste0("anu_d ~ ", paste0(c_var0, collapse = "+"))  ,
                             weights = as.numeric(w_x),
                             data =    data_abs_d)

              # nu_abs_nu_d_x  = predict(fit_abs_d$residuals

              # summary(fit_y)
              # delta_y = fit_y$coefficients

            }


            nu_abs_nu_y = predict(fit_abs,newdata=as.data.frame(Xcurr))
            nu_abs_nu_d = predict(fit_abs_d,newdata=as.data.frame(Xcurr))

            Xcurr_res = cbind(nu_abs_nu_y, nu_abs_nu_d)
            # dim(Xcurr_res)
            # c_var1 = paste0("X", 1:2)
            # colnames(Xcurr_res) <-  c_var0  # c("nu_y_res","nu_d_res")


            set.seed(2)
            redo =1
            nb_redo = 1
            ##### Step 3
            while (redo==1){
              K = max(ceil(K/1.5^(nb_redo-1)),2)
              (cl <- kmeans(Xcurr_res, K))

              # for (j in 1:dim(W_all_y)[2]) {
              W_all_y_disc[,1] = as.numeric(cl$cluster[1:n_y])
              W_all_x_disc[,1] = as.numeric(cl$cluster[(1 + n_y):(n_y + n_x)])
              # }
              redo =0
              for( kk in 1:K){
                if(sum( W_all_y_disc==kk) <2 | sum( W_all_x_disc==kk)<2){
                  redo =1
                  nb_redo =   nb_redo + 1
                }
              }

            }

            ##########################################################


            Ldatad <- as.data.frame(cbind(Y, W_all_y_disc))
            colnames(Ldatad) <- c("Y", "X1")
            Rdatad <- as.data.frame(cbind(Xnc, W_all_x_disc))
            nc_var0 = paste0("Xnc", 1:dim(Xnc)[2])
            colnames(Rdatad) <- c(nc_var0,  "X1")
            # Xc_x0 = matrix(0, dim(Rdatad)[1], 1)


            ###########################

            values = as.matrix(Rdatad[!duplicated(Rdatad[,"X1"]),c_var0[1]])
            refs0 = NULL
            for (j in 1:dim(values)[1]) {
              if (sum(values[j,] == 0) == (dim(values)[2] - 1)) {
                refs0 = c(refs0, j)
              }
            }


            ##########################################################################################

            prob = matrix(0, 1, dim(values)[1])
            boundsxc = matrix(NA, dim(values)[1], 2)



            # Xc_x0 = W_all_x_disc
            # Xc_y0 = W_all_y_disc
            #
            #####################
            # start.time <- Sys.time()
            res0 <- lapply(
              1:length(values),
              loop_EMD,
              dimXc=1,
              values,
              eta_d,
              Y,
              w_x,
              w_y,
              W_all_x_disc,
              W_all_y_disc,
              Xnc,
              nu_d,
              nu_y,
              dataset
            )
            # end.time <- Sys.time()
            # end.time  - start.time
            #

            for (k in 1:length(values)) {
              prob[k] <- as.numeric(res0[[k]][[2]])
              boundsxc[k, ] <- c(res0[[k]][[1]])
            }
            prob <- unlist(prob)
            prob <- prob / sum(prob, na.rm = T)

            bounds[d, ] = c(-off, off) + colSums(boundsxc * (t(prob) %*%
                                                               matrix(1, 1, 2)), na.rm = T)
            bounds[d, 1] = -bounds[d, 1]
            bounds[d, ] = bounds[d, ] / denom #sum(w_x*eta_d^2)




            #### psi 8, psi 9
            pt_details <- as.data.frame(cbind(boundsxc,t(prob),values))
            colnames(pt_details) <- c("LB","UB","prob","W")
            pt_details$W <- as.factor(pt_details$W)


            dataX_1 = cbind( dataX,  W_all_x_disc)
            colnames( dataX_1)[dim( dataX_1)[2]] <- "W"
            dataX_1$W <- as.factor(dataX_1$W)
            dataX_1 <- left_join(  dataX_1, pt_details, by = "W")

            dataY_1 = cbind( dataY,  W_all_y_disc)
            colnames( dataY_1)[dim( dataY_1)[2]] <- "W"
            dataY_1$W <- as.factor(dataY_1$W)
            dataY_1 <- left_join(  dataY_1, pt_details, by = "W")


            j=1
            for(j in 1:2){

              if(j == 1){
                ind = -1


                psi_8 <- matrix(0,dim(dataX_1)[1],1)
                psi_9 <- matrix(0,dim(dataY_1)[1],1)
                for(k in 1:length(values)){
                  psi_8 <- psi_8 +  pt_details$LB[k]*(1*(dataX_1$W==values[k]) -pt_details$prob[k])
                  psi_9 <- psi_9 +  pt_details$LB[k]*(1*(dataY_1$W==values[k]) -pt_details$prob[k])
                }


              }else{
                ind = 1


                psi_8 <- matrix(0,dim(dataX_1)[1],1)
                psi_9 <- matrix(0,dim(dataY_1)[1],1)
                for(k in 1:length(values)){
                  psi_8 <- psi_8 +  pt_details$UB[k]*(1*(dataX_1$W==values[k]) -pt_details$prob[k])
                  psi_9 <- psi_9 +  pt_details$UB[k]*(1*(dataY_1$W==values[k]) -pt_details$prob[k])
                }

              }
              #


              eta_hat = eta_d *ind
              nu_hat = nu_d *ind

              predy <- predict(  fit_y,newdata= as.data.frame(W_all_x))
              mat1= t(W_all_x1)%*%(Xnc_c[,-c(d+1)]*(p_w1%*%matrix(1,1,dim(Xnc_c)[2]-1)))
              # A=fit_y$coefficients%*%t(mat1
              psi1_1 =   nu_hat* predy

              if(dataset==1){

                psi_n_x = t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*    p_w1 -
                  sum(t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*   p_w1)
                psi_n_y = t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y -
                  sum(t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y)

              }else{


                psi_n_x = t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*   w_x -
                  sum(t(ind *delta_d %*%t(W_all_x1))* t(delta_y %*%t(W_all_x1))*   w_x)
                psi_n_y = t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y -
                  sum(t(ind *delta_d %*%t(W_all_y1))* t(delta_y %*%t(W_all_y1))*   w_y)

              }

              B = ginv(t(Xnc_c[,-c(d+1)])%*%(Xnc_c[,-c(d+1)]*(p_w1%*%matrix(1,1,dim(Xnc_c)[2]-1))))
              AB = mat1%*%B
              #
              if((dim(AB)[1])>1){
                infl1_2 = AB%*%t(Xnc_c[,-c(d+1)]*(eta_hat%*%matrix(1,1,dim(Xnc_c)[2]-1)))
                infl1_2 = delta_y%*%infl1_2
                psi1_2 = - t(infl1_2)
              }else{
                infl1_2 = as.numeric(AB)*eta_hat
                psi1_2 = -infl1_2
              }

              psi_1 = psi1_1 + psi1_2
              predd <- predict(  fit_d,newdata= as.data.frame(W_all_y))
              psi_2 =  ind* predd *nu_y

              projh <- matrix(NA,length(values),dimW_all+1)
              projh_5 <- matrix(NA,length(values),dimW_all+1)
              # }

              psi_6 <- matrix(0,n_y,1)
              psi_4 <- matrix(0,n_x,1)

              psi_3 <- matrix(0,n_x,1)
              psi_5 <- matrix(0,n_y,1)

              Bw = ginv(t(W_all_y1)%*%(as.matrix(W_all_y1)*(w_y%*%matrix(1,1,dim(W_all_y1)[2]))))

              k=6
              for(k in 1:length(values)){

                val =values[k]
                sel_x = (W_all_x_disc==val)
                sel_y = (W_all_y_disc==val)

                nu_dk = nu_hat[sel_x]
                nu_yk = nu_y[sel_y]

                n_dk=sum(sel_x)
                n_yk=sum(sel_y)

                if( n_dk >1 & n_yk>1){
                  fdk=ecdf( nu_dk)
                  fyk=ecdf( nu_yk)
                  #
                  #
                  if(factor==TRUE){
                    projh[k,1] <- mean(sort(nu_yk)[eval_Y (ffk(seq(1, n_dk)))])
                  }else{
                    if(dim(W_all_x1)[2]>1){
                      projh[k,] <- colMeans(quantile(nu_yk,fdk(nu_dk))*W_all_x1[sel_x,])

                    }else{
                      projh[k,] <- mean(quantile(nu_yk,fdk(nu_dk))*W_all_x1[sel_x,])

                    }
                    #
                    if(dim(W_all_y1)[2]>1){
                      projh_5[k,] <- colMeans(quantile(nu_dk,fyk(nu_yk))*W_all_y1[sel_y,])
                    }else{
                      projh_5[k,] <- mean(quantile(nu_dk,fyk(nu_yk))*W_all_y1[sel_y,])

                    }
                  }
                  #
                  sU = sort( nu_dk,decreasing = FALSE)
                  sY = sort(nu_yk,decreasing = FALSE)
                  ###
                  eval_Y <- function(j){return( j/length(nu_yk))}
                  eval_eta <- function(j){return( j/length(nu_dk))}
                  # length(diff(sY))

                  d_j6k = c(0,cumsum(quantile(nu_dk,eval_Y(2:n_yk))*diff(sY)))
                  d_i4k = c(0,cumsum(quantile(nu_yk,eval_eta(2:n_dk))*diff(sU)))

                  ranks_d = cbind(seq(1, n_dk),nu_dk)
                  ranks_d = cbind(seq(1, n_dk),ranks_d[order(nu_dk,decreasing = FALSE),])
                  ffk <- approxfun(ranks_d[,2],ranks_d[,1])

                  ranks_y = cbind(seq(1, n_yk),nu_yk)
                  ranks_y = cbind(seq(1, n_yk),ranks_y[order(nu_yk,decreasing = FALSE),])
                  ffk_y <- approxfun(ranks_y[,2],ranks_y[,1])

                  psi_6k <-  d_j6k[ffk_y(seq(1,n_yk))]
                  psi_6[sel_y] <-  psi_6k

                  psi_4k <- d_i4k[ffk(seq(1,n_dk))]
                  psi_4[sel_x] <-  psi_4k


                  ABw = projh[k,]%*%Bw
                  psi_3 [sel_x] =  ABw%*%t(W_all_x1[sel_x,]*nu_dk)

                  ma = ABw%*%AB
                  if(length(ma)>1){
                    infl3_2 = ma%*%t(Xnc_c[sel_x,-c(d+1)]*(eta_hat[sel_x]%*%matrix(1,1,dim(Xnc_c)[2]-1)))
                    # infl3_2 = fit_y$coefficients%*%infl1_2
                    psi_3[sel_x] = psi_3[sel_x]  - t(infl3_2)
                  }else{
                    infl3_2 =as.numeric(ma)*(eta_hat[sel_x])
                    # infl3_2 = fit_y$coefficients%*%infl1_2
                    psi_3[sel_x] = psi_3[sel_x]  - t(infl3_2)
                  }
                  psi_3[sel_x] = - psi_3[sel_x]

                  ABw5= projh_5[k,]%*%Bw
                  psi_5 [sel_y] =  - (ABw5%*%t(W_all_y1[sel_y,]*nu_yk))
                }
              }


              if(dataset==1){
                # psi_7 <-  - (eta_hat^2-denom)*ind*hull_point[d,j]
                psi_10 <-  - (eta_hat^2-denom)*ind*hull_point[d,j]
                psi_x <- psi_1  +  psi_3  + psi_4 + psi_10 + psi_8 +  psi_n_x
                psi_y  <- psi_2 + psi_5 + psi_6  + psi_9 +  psi_n_y


                ves <-lambda * sum(psi_x ^ 2*p_w1) + (1 - lambda) * sum(psi_y^ 2*w_y)
                ves <- ves -   (lambda * sum(psi_x*p_w1) ^ 2  +  (1 - lambda) * sum(psi_y*w_y)^2)
              }else{

                psi_10 <-  - (eta_hat^2-denom)*ind*hull_point[d,j]
                psi_x <- psi_1  +  psi_3  + psi_4 + psi_10 + psi_8 +  psi_n_x
                psi_y  <- psi_2 + psi_5 + psi_6  + psi_9 +  psi_n_y
                ves <-lambda * sum(psi_x ^ 2*w_x)  + (1 - lambda) * sum(psi_y ^ 2*w_y) # +

                ves <- ves -   (lambda * sum(psi_x*w_x)  ^ 2  +  (1 - lambda) *sum(psi_y*w_y) ^2)
              }
              #
              #


              ves0 <-    ves / c(denom)^2
              sig[d,j] = sqrt(ves0)/sqrt(lambda*n_x)

              output_infl <- vector("list")
              output_infl[["lambda"]] <- lambda
              output_infl[["psi_x"]] <- psi_x
              output_infl[["psi_y"]] <- psi_y
              output_infl[["w_x"]] <- w_x
              output_infl[["w_y"]] <- w_y

              output_infl[["denom"]] <- denom
              output_infl[["sig"]]   <-  sig[d,j]

              output_influence[[index_infl]] <-output_infl
              index_infl = index_infl+1


            }

          }





        }



      }


      ind_p = 1
      for(p in 1:(dim(hull_point)[1])){

        if( dimW_all ==0){
          mat_varb_out_unc_asn[p,2]  =   hull_point[p,2] + qnorm(1-alpha)*sig[p,2]
          mat_varb_out_unc_asn[p,1]  =   hull_point[p,1] - qnorm(1-alpha)*sig[p,1]

        }else{
          #######
          ### estimate the correlation.
          sig_U = sig[p,2]*sqrt(lambda*n_x)
          sig_L = sig[p,1]*sqrt(lambda*n_x)


          output_infl_L =  output_influence[[ind_p]]
          ind_p =  ind_p +1
          output_infl_U =  output_influence[[ind_p]]
          ind_p =  ind_p +1


          psi_xU <- output_infl_U[["psi_x"]]/output_infl_U[["denom"]]
          psi_yU <- output_infl_U[["psi_y"]]/output_infl_U[["denom"]]
          psi_xL <- output_infl_L[["psi_x"]]/output_infl_L[["denom"]]
          psi_yL <- output_infl_L[["psi_y"]]/output_infl_L[["denom"]]

          ves <-lambda * sum(psi_xU * psi_xL*w_x)  + (1 - lambda) * sum(psi_yU * psi_yL*w_y) # +
          # lambda*(1-lambda)*sum(psi_n^2* w_xy)
          ves <- ves -   (lambda * sum(psi_xU *w_x) *sum( psi_xL*w_x)  +  (1 - lambda) *sum(psi_yU *w_y) *sum( psi_yL*w_y)) #+
          # lambda*(1-lambda)*sum(psi_n* w_xy)^2)
          rho = ves/(sig_U *sig_L)

          # critical value
          c_hat = interpolate_critical_value(0.05,   abs(rho))

          theta_star = (hull_point[p,1]*sig_U +  hull_point[p,2]*sig_L)/(sig_U + sig_L)
          sig_star = sig_U *sig_L*sqrt(2+2*rho)/(sig_U + sig_L)

          mat_varb_out_unc_asn[p,2]  =   max( hull_point[p,2] + c_hat*sig[p,2],  theta_star + qnorm(1-alpha)*sig_star/sqrt(lambda*n_x))
          mat_varb_out_unc_asn[p,1]  =   min( hull_point[p,1] - c_hat*sig[p,1],  theta_star - qnorm(1-alpha)*sig_star/sqrt(lambda*n_x))
        }
      }






    }
  }

  output <- vector("list")
  output[["pt"]] <- hull_point
  if(ASN==TRUE){
    output[["ci"]] <-mat_varb_out_unc_asn

    output[["influence"]] <- output_influence
  }else{
    output[["ci"]] <-  mat_varb_out_unc


  }


  return(output)
}

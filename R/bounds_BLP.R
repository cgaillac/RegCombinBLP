#' This function compute the DGM bounds on the BLP for all the different coefficients.
#'
#' @param Xnc Matrix of non-conditioning covariates (predictor variables).
#' @param Y Vector or matrix of outcome/dependent variable values.
#' @param Xc_x Matrix of conditioning covariates for the right-hand side data.
#' @param Xc_y Matrix of conditioning covariates for the left-hand side data.
#' @param Wc_x Optional matrix of additional covariates (weights/factors) for right-hand side data.
#' @param Wc_y Optional matrix of additional covariates (weights/factors) for left-hand side data.
#' @param w_x Optional vector of weights for right-hand side data (default NULL).
#' @param w_y Optional vector of weights for left-hand side data (default NULL).
#' @param dimXc Integer: number of conditioning covariates (default 0).
#' @param dimW Integer: number of weight/factor covariates (default 0).
#' @param values Optional list or matrix of values used for discrete/factor covariates (default NULL).
#' @param full Logical: whether to compute full BLP bounds (default FALSE).
#' @param FW Logical: whether to apply Frisch-Waugh-type adjustment (default FALSE).
#' @param discrete Logical: whether to treat conditioning variables as discrete (default FALSE).
#' @param dK DEPRECATED, Numeric: discretization or smoothing parameter (default 2.5).
#' @param unchanged Logical: whether to retain all values when discretizing (default FALSE).
#' @param factor Logical: whether to treat W covariates as factors (default FALSE).
#' @param K_sel Integer: number of clusters for discretization (default 0, auto-selected).
#' @param dataset Integer flag indicating if we use propensity score weighting (default 0).
#' @param O_only Logical, whether to compute bounds only for outside regressors (default FALSE).
#'
#' @return A matrix of dimension `(number of covariates) x 2`, where each row represents
#' the lower and upper bounds for a corresponding regression coefficient.
#'
#' @details
#' The function performs the following operations:
#' 1. Normalizes weights (`w_x`, `w_y`) and computes combined weights.
#' 2. Determines whether a direct computation is possible (same sample size, no weights).
#' 3. If no conditioning or factor variables exist, performs direct regression-based bounds.
#' 4. If conditioning/factor covariates are present:
#'    a. Applies Frisch-Waugh adjustment if `FW = TRUE`.
#'    b. Discretizes or clusters `W` if `discrete = TRUE`.
#'    c. Handles special cases for `exact = 1` or `exact = 2`.
#'    d. Computes weighted residuals, projections, and Earth Mover's Distance (EMD)  bounds.
#' 5. Returns the computed lower and upper bounds for each covariate.
#'
#' @examples
#' \dontrun{
#' # Simple usage with random data
#' Xnc <- matrix(rnorm(100), 50, 2)
#' Y <- rnorm(50)
#' results <- bounds_BLP(Xnc, Y, dimXc = 0, dimW = 0)
#' print(results)
#' }
#'
#' @export
bounds_BLP <-
  function(Xnc,
           Y,
           Xc_x = NULL,
           Xc_y = NULL,
           Wc_x = NULL,
           Wc_y = NULL,
           w_x = NULL,
           w_y = NULL,
           dimXc = 0,
           dimW = 0,
           values = NULL,
           full = FALSE,
           FW = FALSE,
           discrete = FALSE,
           dK = 2.5,
           unchanged = FALSE,
           factor = FALSE,
           K_sel = 0,
           dataset = 0,
           O_only = FALSE) {

    exact=2
    expo = 0.2
    n_x = dim(Xnc)[1]
    n_y = dim(Y)[1]
    dimXnc = dim(Xnc)[2]
    p = dim(Xnc)[2]

    dimW_all = dimXc + dimW

    bounds = matrix(NA, p + dimXc, 2)

    lim = 3

    direct = FALSE
    if (n_x == n_y & is.null(w_x) & is.null(w_y)) {
      direct = TRUE
    }

    if (is.null(w_x)) {
      w_x = matrix(1, n_x, 1) / n_x
    }
    if (is.null(w_y)) {
      w_y = matrix(1, n_y, 1) / n_y
    }

    w_xy = c(w_x *  n_x, w_y *  n_y) / ( n_x +  n_y)

    dataX = NULL
    dataP = NULL

    if (n_x >  lim & n_y > lim) {
      if (dimW_all == 0) {

        ###### when no common regressors, direct
        Xnc_c = cbind(rep(1, n_x, 1), Xnc)
        Y_c = Y

        for (k in 1:p) {
          ### regression of T_1 on the others, to get eta_d
          delta_d = solve((t(Xnc_c[, -c(k + 1)] * c(w_x)) %*% Xnc_c[, -c(k +1)]) ,
                          t(Xnc_c[, -c(k + 1)] * c(w_x)) %*% Xnc_c[, (k + 1)])
          eta_d = Xnc_c[, c(k + 1)] -  Xnc_c[, -c(k + 1)] %*% delta_d

          ## sort both
          Umat = cbind(eta_d, w_x)
          Utilde = Umat[order(Umat[, 1], decreasing = FALSE), ]

          y_c_mat = cbind(Y_c, w_y)
          Ytilde =  y_c_mat[order(y_c_mat[, 1], decreasing = FALSE), ]

          if (direct) {
            bounds[k, 2] = sum(Utilde[, 1] * Ytilde[, 1]) / n_x

          } else{
            value = emd_1d_sorted(Utilde[, 2],
                                  Ytilde[, 2],
                                  Utilde[, 1],
                                  Ytilde[, 1],
                                  metric = 'sqeuclidean',
                                  p = 2.0)
            bounds[k, 2] =   0.5 * (sum(Utilde[, 2] * Utilde[, 1] ^ 2) + sum(Ytilde[, 2] *
                                                                               Ytilde[, 1] ^ 2) -   value$cost)
          }

          Umat2 = cbind(-eta_d, w_x)
          Utilde = Umat2[order(Umat2[, 1], decreasing = FALSE), ]

          if (direct) {
            bounds[k, 1] = -sum(Utilde[, 1] * Ytilde[, 1]) / n_x

          } else{
            value = emd_1d_sorted(Utilde[, 2],
                                  Ytilde[, 2],
                                  Utilde[, 1],
                                  Ytilde[, 1],
                                  metric = 'sqeuclidean',
                                  p = 2.0)
            bounds[k, 1] = -0.5 * (sum(Utilde[, 2] * Utilde[, 1] ^ 2) + sum(Ytilde[, 2] *
                                                                              Ytilde[, 1] ^ 2) -   value$cost)
          }

          bounds[k, ] <-  bounds[k, ] / sum(c(w_x) * eta_d ^ 2)
        }


      } else{
        ###################### with common regressors

        if (FW == TRUE) {
          ##### crude bound, take into account W if any, but otherwise bound (4) with g(W) = Cst.

          if(O_only == TRUE){
            bounds_dim =  dim(Xnc)[2]
          }else{
            bounds_dim =  dim(Xnc)[2] + dimXc
          }

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

          #####  split each coordinantes support for g(w) ####### => get discrete for each.
          ##### discretize if needed #########################################

          Xcurr = rbind(W_all_y, W_all_x)
          colnames(Xcurr) <- c_var0

          bounds = matrix(0,   bounds_dim , 2)

          #### !! here Xc_x and not the discrete one
          if (dimXc != 0) {
            Xnc_c = cbind(matrix(1, n_x, 1),Xnc,Xc_x)  ### to build the A matrix.
          } else{
            Xnc_c = cbind(matrix(1, n_x, 1), Xnc)  ### to build the A matrix.

          }

          #### when using the discretized version as W
          if (factor == TRUE) {

            values =  W_all_y[!duplicated( W_all_y[,c_var0]),c_var0]
            Xc_y0 = matrix(0, dim(W_all_y)[1], 1)
            ind = 1
            for (j in dim(values)[1]) {
              if (dimW_all == 1) {
                val = values[j,]
                sel_y = (W_all_y[, c_var0] == val)
              } else {
                val = t(as.matrix(values[j,]))
                sel_y = matrix(1, dim(W_all_y)[1], 1)
                for (ddd in 1:dimW_all) {
                  sel_y = sel_y & (W_all_y[,  c_var0[ddd]] == val[ddd])
                }
                sel_y = matrix(sel_y, dim(W_all_y)[1], 1)
              }
              Xc_y0[sel_y] = ind
              ind = ind + 1
            }



            dataY = as.data.frame(cbind(Y,  Xc_y0, w_y))
            colnames(dataY) <- c("Y", "W", "w_y")
            dataY$W <- as.factor(dataY$W)

            dataY <- dataY %>% group_by(W) %>% mutate(nu_y = Y - sum(w_y * Y) / sum(w_y),
                                                      projy = sum(w_y * Y) /
                                                        sum(w_y))
            nu_y = dataY$nu_y

            if (dataset == 1) {

              dataP <-
                as.data.frame(cbind(
                  c(rep(0, n_x), rep(1, n_y)),
                  c(w_x, w_y) / 2 ,
                  rbind(W_all_x, W_all_y)
                ))
              colnames(dataP) <- c("D", "w", "W")

              fit_w <-
                glm(
                  "D ~ W",
                  data = dataP,
                  family = "binomial",
                  weights = as.numeric(dataP$w)
                )
              p_w1 <-
                predict(fit_w,
                        data = rbind(W_all_x, W_all_y),
                        type = "response")
              p_w1 <-   p_w1[1:n_x] / (1 -  p_w1[1:n_x])
              p_w1 <-   p_w1 * w_x / sum(p_w1 * w_x)
            } else{
              p_w1 = w_x

            }



          } else{
            #
            dataY = as.data.frame(cbind(Y, W_all_y, w_y))
            colnames(dataY) <- c("Y", c_var0, "w_y")


            fit_y = lm(paste0("Y ~ ", paste0(c_var0, collapse = "+"))  ,
                       data =    dataY,
                       weights = as.numeric(w_y))

            nu_y  = fit_y$residuals

            delta_y = fit_y$coefficients

            if (dataset == 1) {

              dataP <- as.data.frame(cbind(c(rep(0, n_x), rep(1, n_y)),  c(w_x, w_y) / 2 , W_all_all))
              colnames(dataP) <- c("D", "w", "1", c_var0)
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



          d = 1
          for (d in 1:bounds_dim){

            fit_d0 = lm(Xnc_c[, (d + 1)] ~ Xnc_c[, -c(d + 1)] - 1  , weights =  as.numeric(p_w1))
            eta_d  = fit_d0$residuals


            if (factor == TRUE) {


              values =  W_all_x[!duplicated( W_all_x[,c_var0]),c_var0]

              Xc_x0 = matrix(0, dim(W_all_x)[1], 1)
              ind = 1

              for (j in dim(values)[1]) {
                if (dimW_all == 1) {
                  val = values[j,]

                  sel_y = (W_all_x[, c_var0] == val)
                } else {
                  val = t(as.matrix(values[j,]))
                  sel_x = matrix(1, dim(W_all_x)[1], 1)

                  for (ddd in 1:dimW_all) {
                    sel_x = sel_x & (W_all_x[,  c_var0[ddd]] == val[ddd])

                  }

                  sel_x = matrix(sel_x, dim(W_all_x)[1], 1)

                }
                Xc_x0[sel_x] = ind

                ind = ind + 1
              }


              dataX = as.data.frame(cbind(eta_d , Xc_x0, w_x,     p_w1))
              colnames(dataX) <- c("eta_d", "W", "w_x", "p_w")
              dataX$W <- as.factor(dataX$W)

              if (dataset == 1) {
                dataX <- dataX %>% group_by(W) %>% mutate(nu_d = eta_d - sum(p_w * eta_d) / sum(p_w))
                nu_d = dataX$nu_d
                dataXs <-
                  dataX[, c("eta_d", "W", "w_x", "nu_d", "p_w")] %>% group_by(W) %>% summarise(
                    projd =  sum(p_w * eta_d) / sum(p_w),
                    projd2 =  sum(p_w * eta_d ^ 2) / sum(p_w)
                  )
              } else{
                dataX <-
                  dataX %>% group_by(W) %>% mutate(nu_d = eta_d - sum(w_x * eta_d) / sum(w_x))
                nu_d = dataX$nu_d
                dataXs <-
                  dataX[, c("eta_d", "W", "w_x", "nu_d")] %>% group_by(W) %>% summarise(
                    projd =  sum(w_x * eta_d) / sum(w_x),
                    projd2 =  sum(w_x * eta_d ^ 2) / sum(w_x)
                  )
              }

              dataYs <-
                dataY[, c("Y", "W", "w_y", "nu_y")] %>% group_by(W) %>% summarise(projy =  sum(w_y *
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

                off = sum(dataWsx$projy * dataWsx$projd * w_xy[1:dim(dataX)[1]], na.rm =
                            T) +
                  sum(dataWsy$projy * dataWsy$projd * w_xy[(dim(dataX)[1] +
                                                              1):(dim(dataX)[1] + dim(dataY)[1])], na.rm = T)

                denom = sum(dataWsx$projd2 * w_xy[1:dim(dataX)[1]], na.rm =
                              T) +
                  sum(dataWsy$projd2 *  w_xy[(dim(dataX)[1] + 1):(dim(dataX)[1] + dim(dataY)[1])], na.rm =
                        T)
              }



            } else{
              dataX = as.data.frame(cbind(eta_d, W_all_x, p_w1))
              colnames(dataX) <- c("eta_d", c_var0, "w_x")
              fit_d = lm(paste0("eta_d ~ ", paste0(c_var0, collapse = "+"))  ,
                         data = dataX,
                         weights = dataX$w_x)

              nu_d  = fit_d$residuals

              delta_d = fit_d$coefficients

              if (dataset == 1) {
                V = t(W_all_y1) %*% (W_all_y1 * (w_y %*% matrix(1, 1, 1 + dim(W_all_y)[2])))
              } else{
                V = t(W_all_all) %*% (W_all_all * w_xy %*% matrix(1, 1, dim(W_all_all)[2]))
              }

              off = c((t(delta_d) %*% V) %*% delta_y)

              if (dataset == 1) {
                denom = sum(eta_d ^ 2 * p_w1 , na.rm = T)
              } else{
                denom = sum(eta_d ^ 2 * w_x , na.rm = T)
              }
            }



            Umat = cbind(nu_d, w_x)
            Utilde_p = Umat[order(Umat[, 1], decreasing = FALSE), ]

            y_c_mat = cbind(nu_y, w_y)
            Ytilde =  y_c_mat[order(y_c_mat[, 1], decreasing = FALSE), ]


            value = emd_1d_sorted(Utilde_p[, 2],
                                  Ytilde[, 2],
                                  Utilde_p[, 1],
                                  Ytilde[, 1],
                                  metric = 'sqeuclidean',
                                  p = 2.0)
            bounds[d, 2] <-
              0.5 * (sum(Utilde_p[, 2] * Utilde_p[, 1] ^ 2) + sum(Ytilde[, 2] * Ytilde[, 1] ^
                                                                    2) -   value$cost)

            Umat2 = cbind(-nu_d, w_x)
            Utilde_m = Umat2[order(Umat2[, 1], decreasing = FALSE), ]

            value = emd_1d_sorted(Utilde_m[, 2],
                                  Ytilde[, 2],
                                  Utilde_m[, 1],
                                  Ytilde[, 1],
                                  metric = 'sqeuclidean',
                                  p = 2.0)
            bounds[d, 1] =  0.5 * (sum(Utilde_m[, 2] * Utilde_m[, 1] ^ 2) + sum(Ytilde[, 2] *
                                                                                  Ytilde[, 1] ^ 2) -   value$cost)

            bounds[d, ] = c(-off, off) + bounds[d, ]
            bounds[d, 1] = -bounds[d, 1]
            bounds[d, ] = bounds[d, ] / denom


          }


        } else if (discrete == TRUE && exact<2) {
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

          }




          W_all_y_disc = W_all_y
          W_all_x_disc = W_all_x


          ##################### previous clustering based on Kmeans, directly on the covariates.
          if(exact==1){

            set.seed(2)
            (cl <- kmeans(Xcurr, K))

            for (j in 1:dim(W_all_y)[2]) {
              W_all_y_disc[, j] = as.numeric(cl$cluster[1:n_y])
              W_all_x_disc[, j] = as.numeric(cl$cluster[(1 + n_y):(n_y + n_x)])
            }

          }


          ##########################################################

          c_var0 = paste0("X", 1:dimW_all)

          Ldatad <- as.data.frame(cbind(Y, W_all_y_disc))
          colnames(Ldatad) <- c("Y", c_var0)
          Rdatad <- as.data.frame(cbind(Xnc, W_all_x_disc))
          nc_var0 = paste0("Xnc", 1:dim(Xnc)[2])
          colnames(Rdatad) <- c(nc_var0, c_var0)
          Xc_x0 = matrix(0, dim(Rdatad)[1], 1)



          ########################### get only the ones with sufficiently enougth data #######
          # unchanged = TRUE

          if (dimW_all != 1) {

            values = Rdatad[!duplicated(Rdatad[,c_var0]),c_var0]
            refs0 = NULL
            for (j in 1:dim(values)[1]) {
              if (sum(values[j,] == 0) == (dim(values)[2] - 1)) {
                refs0 = c(refs0, j)
              }
            }
          } else {

            values = Rdatad[!duplicated(Rdatad[,c_var0]),c_var0]
            values = as.matrix( values, length(values),1)

            refs0 = (1:length(values))[values > 0]
          }


          if (unchanged == FALSE) {
            if (var(w_x) != 0) {
              nb_s = 10 ^ 6
              sample_R = sample(1:dim(Rdatad)[1],
                                nb_s,
                                replace = TRUE,
                                prob = w_x)
              Rdatad_sel <- Rdatad[sample_R, ]
            } else{
              Rdatad_sel = Rdatad
            }

            if (var(w_y) != 0) {
              nb_s = 10 ^ 6
              sample_L = sample(1:dim(Ldatad)[1],
                                nb_s,
                                replace = TRUE,
                                prob = w_y)
              Ldatad_sel <- Ldatad[sample_L, ]
            } else{
              Ldatad_sel = Ldatad
            }

            if (dimW_all == 1) {
              Xc_pool <- c(Rdatad_sel[, c_var0], Ldatad_sel[, c_var0])
            } else {
              Xc_pool <- rbind(Rdatad_sel[, c_var0], Ldatad_sel[, c_var0])
            }

          } else{
            if (dimW_all == 1) {
              Xc_pool <- c(Rdatad[, c_var0], Ldatad[, c_var0])
            } else {
              Xc_pool <- rbind(Rdatad[, c_var0], Ldatad[, c_var0])
            }

            Rdatad_sel = Rdatad
            Ldatad_sel = Ldatad
          }

          if (dimW_all != 1) {
            nv = dim(values)[1]
          }else{
            nv = length(values)
          }

          select_values = NULL
          unselect_values = NULL
          if (unchanged == TRUE) {
            select_values = 1:nv
            unselect_values = NULL
          } else {
            values_tab = matrix(0, nv, 3)
            for (k in 1:nv) {
              values_tab[k, 1] <-
                tabulate_values(k, values, Rdatad_sel [, c_var0], dimW_all)
              values_tab[k, 2] <-
                tabulate_values(k, values, Ldatad_sel [, c_var0],  dimW_all)
              values_tab[k, 3] <-
                tabulate_values(k, values, Xc_pool, dimW_all)
            }



            nb_min = 5 #########################
            for (j in 1: nv) {
              if (values_tab[j, 1]  >= nb_min & values_tab[j, 2] >= nb_min) {
                select_values = c(select_values, j)
              } else {
                unselect_values = c(unselect_values, j)
              }
            }
          }
          values_sel = vector("list")
          values_sel[["selected"]] <- values[select_values, ]
          values_sel[["old"]] <- values

          Xc_x0 = matrix(0, dim(Rdatad)[1], 1)
          Xc_y0 = matrix(0, dim(Ldatad)[1], 1)
          ind = 1
          # j=2
          for (j in select_values[-c(1)]) {
            if (dimW_all == 1) {
              val = values[j,]
              sel_x = (Rdatad[, c_var0] == val)
              sel_y = (Ldatad[, c_var0] == val)
            } else {
              val = t(as.matrix(values[j,]))
              sel_x = matrix(1, dim(Rdatad)[1], 1)
              sel_y = matrix(1, dim(Ldatad)[1], 1)
              for (ddd in 1:dimW_all) {
                sel_x = sel_x & (Rdatad[,  c_var0[ddd]] == val[ddd])
                sel_y = sel_y & (Ldatad[,  c_var0[ddd]] == val[ddd])
              }
              # sum(sel_x)
              sel_x = matrix(sel_x, dim(Rdatad)[1], 1)
              sel_y = matrix(sel_y, dim(Ldatad)[1], 1)
            }
            Xc_x0[sel_x] = ind
            Xc_y0[sel_y] = ind
            ind = ind + 1
          }

          colnames(Xc_y0) <- "Xc"
          colnames(Xc_x0) <- "Xc"
          Rdata0 <- Rdatad
          Ldata0 <- Ldatad
          Rdatad <-
            as.data.frame(cbind(Rdatad[, 1:dim(Xnc)[2]], Xc_x0))
          colnames(Rdatad) <- c(nc_var0, "Xc")
          Ldatad <- as.data.frame(cbind(Ldatad[, 1], Xc_y0))
          colnames(Ldatad) <- c("Y", "Xc")
          dimXc_old = dimXc
          dimWc_old = dimW
          c_var_old = c_var0
          values_old = values
          dimXc = 1
          c_var_g = "Xc"
          groups = 1
          values = matrix(create_values(dimXc, c_var_g, Rdatad))
          refs0 = (1:length(values))[values > 0]
          ##########################################################################################




          if (dim(values)[1] > 2) {


            if(O_only == TRUE){
              bounds_dim =  dim(Xnc)[2]
            }else{
              bounds_dim =  dim(Xnc)[2] + dimXc_old
            }

            bounds = matrix(0,   bounds_dim, 2)
            prob = matrix(0, 1, dim(values)[1])
            boundsxc = matrix(NA, dim(values)[1], 2)

            #### !! here Xc_x and not the discrete one
            if (dimXc_old != 0) {
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

                dataP <-as.data.frame(cbind(c(rep( 0, n_x ), rep(1, n_y)),  c(w_x, w_y) / 2 , rbind(Xc_x0, Xc_y0)))
                colnames(dataP) <- c("D", "w", "W")

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


              fit_y = lm(paste0("Y ~ ", paste0(c_var0, collapse = "+"))  ,
                         weights = as.numeric(w_y),
                         data =    dataY)

              nu_y  = fit_y$residuals

              delta_y = fit_y$coefficients




              if (dataset == 1) {
                dataP <-
                  as.data.frame(cbind(c(rep(0, n_x ), rep( 1, n_y)),  c(w_x, w_y) / 2 , W_all_all))
                colnames(dataP) <- c("D", "w", "1", c_var0)

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


            for (d in 1:bounds_dim) {

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

                dataYs <-
                  dataY[, c("Y", "W", "w_y", "nu_y")] %>% group_by(W) %>% summarise(projy =  sum(w_y *
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
                  off = sum(dataWsx$projy * dataWsx$projd * w_xy[1:dim(dataX)[1]], na.rm =
                              T) +
                    sum(dataWsy$projy * dataWsy$projd * w_xy[(dim(dataX)[1] +
                                                                1):(dim(dataX)[1] + dim(dataY)[1])], na.rm = T)

                  denom = sum(dataWsx$projd2 * w_xy[1:dim(dataX)[1]], na.rm =
                                T) +
                    sum(dataWsy$projd2 *  w_xy[(dim(dataX)[1] + 1):(dim(dataX)[1] + dim(dataY)[1])], na.rm =
                          T)
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
                  V = t(W_all_y1) %*% (W_all_y1 * (w_y %*% matrix(1, 1, 1 +dim(W_all_y)[2])))
                } else{
                  V = t(W_all_all) %*% (W_all_all * w_xy%*% matrix(1, 1, dim(W_all_all)[2]))
                }

                off = c((t(delta_d) %*% V) %*% delta_y)

                if (dataset == 1) {
                  denom = sum(eta_d ^ 2 * p_w1 , na.rm = T)
                } else{
                  denom = sum(eta_d ^ 2 * w_x , na.rm = T)
                }
                # }

              }


              res0 <- lapply(
                1:length(values),
                loop_EMD,
                dimXc,
                values,
                eta_d,
                Y,
                w_x,
                w_y,
                Xc_x0,
                Xc_y0,
                Xnc,
                nu_d,
                nu_y,
                dataset
              )


              for (k in 1:length(values)) {
                prob[k] <- as.numeric(res0[[k]][[2]])
                boundsxc[k, ] <- c(res0[[k]][[1]])
              }
              prob <- unlist(prob)
              prob <- prob / sum(prob, na.rm = T)

              bounds[d, ] = c(-off, off) + colSums(boundsxc * (t(prob) %*%
                                                                 matrix(1, 1, 2)), na.rm = T)
              bounds[d, 1] = -bounds[d, 1]
              bounds[d, ] = bounds[d, ] / denom

            }


          }else{
            bounds <- bounds_BLP(Xnc, Y, Xc_x = NULL, Xc_y = NULL, w_x, w_y)
          }
          #


        } else if(discrete == TRUE && exact==2) {
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


          #####  split each coordinantes support for g(w) ####### => get discrete for each.
          ##### discretize if needed #########################################
          # j=1
          Xcurr = rbind(W_all_y, W_all_x)
          colnames(Xcurr) <- c_var0

          if (K_sel > 0) {
            K = K_sel
          } else{
            K <- max(2, floor(min(n_y,n_x)^expo))
          }

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

              dataP <-as.data.frame(cbind(c(rep( 0, n_x ), rep(1, n_y)),  c(w_x, w_y) / 2 , rbind(Xc_x0, Xc_y0)))
              colnames(dataP) <- c("D", "w", "W")

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

            fit_y = lm(paste0("Y ~ ", paste0(c_var0, collapse = "+"))  ,
                       weights = as.numeric(w_y),
                       data =    dataY)

            nu_y  = fit_y$residuals

            delta_y = fit_y$coefficients



            if (dataset == 1) {
              dataP <-
                as.data.frame(cbind(c(rep(0, n_x ), rep( 1, n_y)),  c(w_x, w_y) / 2 , W_all_all))
              colnames(dataP) <- c("D", "w", "1", c_var0)

              fit_w <-
                glm(
                  paste0("D ~ ", paste0(c_var0, collapse = "+")) ,
                  data = dataP,
                  family = "binomial",
                  weights = NULL
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

              delta_d = fit_d$coefficients

              if (dataset == 1) {
                V = t(W_all_y1) %*% (W_all_y1 * (w_y %*% matrix(1, 1, 1 +dim(W_all_y)[2])))
              } else{
                V = t(W_all_all) %*% (W_all_all *w_xy %*% matrix(1, 1, dim(W_all_all)[2]) )
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

              fit_abs_d = lm(paste0("anu_d ~ ", paste0(c_var0, collapse = "+"))  ,
                             weights = as.numeric(w_x),
                             data =    data_abs_d)

            }


            nu_abs_nu_y = predict(fit_abs,newdata=as.data.frame(Xcurr))
            nu_abs_nu_d = predict(fit_abs_d,newdata=as.data.frame(Xcurr))

            Xcurr_res = cbind(nu_abs_nu_y, nu_abs_nu_d)

            ##### Step 3
            set.seed(2)
            redo =1
            nb_redo = 1
            ##### Step 3
            while (redo==1){

              K = max(ceil(K/1.5^(nb_redo-1)),2)
              (cl <- kmeans(Xcurr_res, K))

              W_all_y_disc[,1] = as.numeric(cl$cluster[1:n_y])
              W_all_x_disc[,1] = as.numeric(cl$cluster[(1 + n_y):(n_y + n_x)])

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

            #####################

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

            for (k in 1:length(values)) {
              prob[k] <- as.numeric(res0[[k]][[2]])
              boundsxc[k, ] <- c(res0[[k]][[1]])
            }
            prob <- unlist(prob)
            prob <- prob / sum(prob, na.rm = T)

            bounds[d, ] = c(-off, off) + colSums(boundsxc * (t(prob) %*%
                                                               matrix(1, 1, 2)), na.rm = T)
            bounds[d, 1] = -bounds[d, 1]
            bounds[d, ] = bounds[d, ] / denom

          }


        } else{
              #### continuous case, see other code.

        }
      }
    }
    return(bounds)
  }

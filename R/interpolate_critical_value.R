#' Interpolate Critical Values for Given Alpha and Correlation
#'
#' Returns the interpolated critical value (z-value) based on a pre-defined table of critical
#' values for different significance levels (alpha) and correlation coefficients (rho).
#'
#' @param alpha Numeric: significance level. Supported values are 0.1, 0.05, 0.01.
#' @param rho Numeric: correlation coefficient (between 0 and 1) at which to interpolate the critical value.
#'
#' @return Numeric. The interpolated critical value corresponding to the specified `alpha` and `rho`.
#'
#' @details
#' The function defines a small grid of critical z-values for combinations of `alpha` and `rho`.
#' It then uses linear interpolation to return the approximate critical value for any input `rho`.
#' Values outside the predefined `rho` range are clamped to the nearest boundary (`rule = 2`).
#'
#' @examples
#' \dontrun{
#' # Interpolate critical value for alpha = 0.05 and rho = 0.92
#' interpolate_critical_value(0.05, 0.92)
#'
#' # Interpolate for alpha = 0.01 and rho = 0.99
#' interpolate_critical_value(0.01, 0.99)
#'}
#' @export
interpolate_critical_value <- function(alpha, rho) {
  # Define the grid

  alphas <- c(0.1, 0.05, 0.01)
  rhos <- c(0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 1.0)

  # Corresponding z-values from your table
  z_table <- matrix(c(
    1.28, 1.29, 1.31, 1.36, 1.44, 1.54, 1.64,   # alpha = 0.1
    1.64, 1.65, 1.65, 1.70, 1.76, 1.81, 1.96,   # alpha = 0.05
    2.33, 2.33, 2.33, 2.34, 2.40, 2.43, 2.58    # alpha = 0.01
  ), nrow = 3, byrow = TRUE)
  rownames(z_table)<-alphas

  return( approx(rhos, y = z_table[rownames(z_table)==alpha,], xout= rho,rule = 2)$y)
}





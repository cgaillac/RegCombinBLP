#' Compute 1D Earth Mover's Distance (EMD) for Sorted Distributions
#'
#' This function computes the 1-dimensional Earth Mover's Distance (EMD) between two
#' sorted distributions with associated weights. The distributions `u` and `v`
#' must be sorted in ascending order. It supports different distance metrics.
#'
#' @param u_weights Numeric vector of weights for the `u` distribution. Must be same length as `u`.
#' @param v_weights Numeric vector of weights for the `v` distribution. Must be same length as `v`.
#' @param u Numeric vector of sorted points for the first distribution.
#' @param v Numeric vector of sorted points for the second distribution.
#' @param metric Character. Distance metric to use. Options:
#'   1. `'sqeuclidean'` (default): squared Euclidean distance.
#'   2. `'euclidean'` or `'cityblock'`: absolute distance.
#'   3. `'minkowski'`: Minkowski distance with order `p`.
#' @param p Numeric. Power parameter for the Minkowski metric (default is 2.0).
#'
#' @return A list with components:
#' \describe{
#'   \item{G}{Numeric vector of transported weights along the optimal flow.}
#'   \item{indices}{Matrix of indices indicating source (`u`) and target (`v`) positions for each transported weight.}
#'   \item{cost}{Numeric scalar. The total transportation cost (EMD) between `u` and `v`.}
#' }
#'
#' @details
#' This function assumes that `u` and `v` are sorted in ascending order. The algorithm
#' implements a simple greedy approach for 1D distributions:
#' it iteratively moves weight from the current source to the current target until
#' all mass is transported. Different distance metrics can be specified to compute
#' the cost, including squared Euclidean, Manhattan, and Minkowski distances.
#'
#' @note
#' Adapted from the Python Optimal Transport (POT) library:
#' Flamary R., Vincent-Cuaz C., Courty N., Gramfort A., et al. POT Python Optimal Transport (version 0.9.5).
#' URL: https://github.com/PythonOT/POT
#'
#'
#' @examples
#' \dontrun{
#' # Simple example with uniform weights
#' u <- c(1, 2, 3)
#' v <- c(2, 3, 4)
#' u_weights <- c(0.3, 0.4, 0.3)
#' v_weights <- c(0.5, 0.2, 0.3)
#' emd_1d_sorted(u_weights, v_weights, u, v)
#'
#' # Using the Minkowski distance
#' emd_1d_sorted(u_weights, v_weights, u, v, metric = 'minkowski', p = 3)
#'}
#' @export
#'
emd_1d_sorted <- function(u_weights, v_weights, u, v, metric = 'sqeuclidean', p = 2.0) {
  # Check if the input parameters are consistent
  if (length(u_weights) != length(u) || length(v_weights) != length(v)) {
    stop("Lengths of weights and corresponding points must be the same.")
  }

  # Initialize variables
  cost <- 0.0
  n <- length(u_weights)
  m <- length(v_weights)
  i <- 0
  w_i <- u_weights[1]
  j <- 0
  w_j <- v_weights[1]
  m_ij <- 0.0

  G <- numeric(n + m - 1)  # Transportation matrix
  indices <- matrix(0, nrow = n + m - 1, ncol = 2)  # Indices matrix

  cur_idx <- 0

  # Main loop for solving the transportation problem
  repeat {
    # Calculate the distance based on the given metric
    if (metric == 'sqeuclidean') {
      m_ij <- (u[i + 1] - v[j + 1])^2
    } else if (metric == 'cityblock' || metric == 'euclidean') {
      m_ij <- abs(u[i + 1] - v[j + 1])
    } else if (metric == 'minkowski') {
      m_ij <- abs(u[i + 1] - v[j + 1])^p
    } else {
      stop("Unsupported metric")
    }

    # Decide how to move between source and target distributions
    if (w_i < w_j || j == m - 1) {
      cost <- cost + m_ij * w_i
      G[cur_idx + 1] <- w_i
      indices[cur_idx + 1, 1] <- i
      indices[cur_idx + 1, 2] <- j
      i <- i + 1
      if (i == n) {
        break
      }
      w_j <- w_j - w_i
      w_i <- u_weights[i + 1]
    } else {
      cost <- cost + m_ij * w_j
      G[cur_idx + 1] <- w_j
      indices[cur_idx + 1, 1] <- i
      indices[cur_idx + 1, 2] <- j
      j <- j + 1
      if (j == m) {
        break
      }
      w_i <- w_i - w_j
      w_j <- v_weights[j + 1]
    }

    cur_idx <- cur_idx + 1
  }

  cur_idx <- cur_idx + 1
  return(list(G = G[1:cur_idx], indices = indices[1:cur_idx, ], cost = cost))
}

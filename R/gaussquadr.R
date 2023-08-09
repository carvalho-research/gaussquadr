#' gaussquadr: Gaussian Quadrature
#'
#' The \code{gaussquadr} package provides routines for multivariate Gaussian
#' quadrature.
#'
#' @section \code{gaussquadr} class:
#' The main object is \code{GaussQuad}, with linear transformation methods
#' \code{shift}, \code{rotate}, and \code{scale} and the main quadrature method
#' \code{apply}.
#'
#' @docType package
#' @name gaussquadr
#' @useDynLib gaussquadr, .registration = TRUE, .fixes = "G_"
NULL


gauss_kinds <- list(
  # classic versions:
  "legendre" = 1L, # w(x) = 1 on (-1, 1)
  "chebyshev1" = 2L, # w(x) = 1 / sqrt(1-x^2) on (-1, 1)
  "chebyshev2" = 3L, # w(x) = sqrt(1-x^2) on (-1, 1)
  "hermite" = 4L, # w(x) = exp(-x^2) on (-inf, inf)
  "jacobi" = 5L, # w(x) = (1-x)^alpha * (1+x)^beta on (-1, 1)
  "laguerre" = 6L, # w(x) = x^alpha * exp(-x) on (0, inf)
  # probabilist versions:
  "uniform" = 1L, # w(x) = 1 on (0, 1)
  "wigner" = 3L, # w(x) = sqrt(1-x^2) / (pi / 2) on (-1, 1)
  "gaussian" = 4L, # w(x) = exp(-x^2 / 2) / sqrt(2 * pi) on (-inf, inf)
  "beta" = 5L, # w(x) = x^(alpha-1) * (1-x)^beta / B(alpha, beta) on (0, 1)
  "gamma" = 6L # w(x) = x^(alpha-1) * exp(-beta*x) /
               #          (beta^alpha / Gamma(alpha)) on (0, inf)
)

gauss_quad <- function (kind, n, alpha = 0., beta = 0., endpts = NULL) {
  n <- as.integer(n)
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  endpts <- as.double(endpts)
  kpts <- min(length(endpts), 2L)

  buffer <- numeric(n)
  res <- list(nodes = numeric(n), weights = numeric(n))
  res$ierr <- .Call(G_gauss_quad, kind, alpha, beta, kpts, endpts, buffer,
                    res$nodes, res$weights, PACKAGE = "gaussquadr")
  res
}

# TODO: filter by threshold
# TODO: iterator
multi_grid <- function (x, n = length(x)) {
  if (n == 1) matrix(x, ncol = 1) else {
    l <- length(x); g <- multi_grid(x, n - 1)
    cbind(kronecker(rep(1, l), g), kronecker(x, rep(1, l ^ (n - 1))))
  }
}

apply_integral <- function (f) {
  gf <- if (is.vector(self$nodes)) sapply(self$nodes, f) else
    apply(self$nodes, 1, f)
  if (is.vector(gf)) sum(gf * self$weights) else
    rowSums(sweep(gf, 2, self$weights, `*`))
}

#' GaussQuad: multivariate Gauss quadrature
#'
#' R6 class with methods for linear transformations and integration.
#' @export
GaussQuad <- R6::R6Class("GaussQuad", public = list(
  #' @field type Type of Gauss quadrature
  type = NA,

  #' @field nodes Grid of nodes
  nodes = NULL,

  #' @field weights Quadrature weights
  weights = NULL,

  #' @field alpha Parameter for Jacobi, Laguerre, beta and gamma quadratures
  alpha = 0.,

  #' @field beta Parameter for Jacobi and beta quadratures
  beta = 0.,

  #' @description
  #' Create a new \code{GaussQuad} object.
  #' @param type Type of Gauss quadrature.
  #' @param n Number of nodes in one dimension.
  #' @param d Number of dimensions.
  #' @param alpha Parameter for Jacobi, Laguerre, beta and gamma quadratures.
  #' @param beta Parameter for Jacobi and beta quadratures.
  #' @param endpts Endpoints to be fixed in quadrature.
  #' @param prune_coef Coefficient controling pruning: only weights above
  #'        \code{exp(d * ((1 - pc) * log(min(w)) + pc * log(max(w))))} are kept.
  #'        Defaults to \code{1 - 1 / d}. Use \code{prune_coef = 0} to keep all
  #'        weights.
  #' @param prune_eps Prune weights below \code{eps(sum(weights))}? Defaults to TRUE.
  #' @return A new \code{GaussQuad} object.
  initialize = function (type, n, d = 1, alpha = 0., beta = 1., endpts = NULL,
                         prune_coef = 1 - 1 / d, prune_eps = TRUE) {
    kind <- gauss_kinds[[type]]
    if (is.null(kind)) stop("invalid type: `", type, "`")
    if (n <= 0) stop("invalid number of nodes: ", n)
    if (d <= 0) stop("invalid number of dimensions: ", d)
    if ((kind == 5 || kind == 6) && alpha < -1)
      stop("invalid alpha exponent: ", alpha)
    if ((type == "jacobi" || type == "laguerre") && alpha < -1)
      stop("invalid alpha exponent: ", alpha)
    if ((type == "gamma" || type == "beta") && alpha < 0)
      stop("invalid alpha exponent: ", alpha)
    if (type == "jacobi" && beta < -1) stop("invalid beta exponent: ", beta)
    if ((type == "gamma" || type == "beta") && beta < 0)
      stop("invalid beta exponent: ", beta)
    if (prune_coef < 0 || prune_coef >= 1)
      stop("invalid prune parameter: ", prune_coef)

    gq <- if (type == "gamma" || type == "beta")
      gauss_quad(kind, n, alpha - 1, beta - 1, endpts)
    else
      gauss_quad(kind, n, alpha, beta, endpts)
    switch(type, # correct for probabilist versions
      "uniform" = {
        gq$nodes <- .5 * (1 + gq$nodes)
        gq$weights <- gq$weights / 2
      },
      "wigner" = { gq$weights <- gq$weights / (pi / 2) },
      "gaussian" = {
        gq$nodes <- gq$nodes * sqrt(2)
        gq$weights <- gq$weights / sqrt(pi)
      },
      "gamma" = {
        gq$nodes <- gq$nodes / beta
        lz <- lgamma(alpha)
        gq$weights <- exp(log(gq$weights) - lz)
      },
      "beta" = {
        gq$nodes <- rev(.5 * (1 - gq$nodes))
        lz <- lbeta(alpha, beta) + (alpha + beta - 1) * log(2)
        gq$weights <- rev(exp(log(gq$weights) - lz))
      }
    )

    if (prune_coef > 0) {
      th <- exp(d * ((1 - prune_coef) * log(gq$weights[1]) +
                     prune_coef * log(gq$weights[(n + 1) / 2])))
    }
    ng <- multi_grid(gq$nodes, d)
    w <- apply(multi_grid(gq$weights, d), 1, prod)
    if (prune_coef > 0) { ng <- ng[w >= th, ]; w <- w[w >= th] }
    if (prune_eps) {
      eps <- .Machine$double.eps * sum(w)
      ng <- ng[w >= eps, ]; w <- w[w >= eps]
    }
    self$type <- type; self$nodes <- ng; self$weights <- w
    self$alpha <- alpha; self$beta <- beta
  },

  #' @description
  #' Prints object information.
  #' @param short Suppress type? Defaults to FALSE.
  print = function (short = FALSE) {
    nd <- dim(self$nodes)
    if (length(nd) == 0) {
      n <- length(self$nodes); d <- 1
    } else {
      n <- nd[1]; d <- nd[2]
    }
    cat(class(self)[1])
    if (!short) cat(paste0(" of type `", self$type, "`"))
    if (self$type == "jacobi" || self$type == "laguerre" ||
        self$type == "gamma" || self$type == "beta")
      cat(paste0(", alpha = ", self$alpha))
    if (self$type == "jacobi" || self$type == "gamma" || self$type == "beta")
      cat(paste0(", beta = ", self$beta))
    cat(paste0(" [ ", n, " nodes"))
    if (d > 1) cat(paste0(" x ", d, " dims"))
    cat(" ]\n")
    invisible(self)
  },

  #' @description
  #' Scale the node grid of a \code{GaussQuad} object.
  #' @param scale Scale to apply to each dimension of the grid.
  scale = function (scale) {
    self$nodes <- if (is.vector(self$nodes)) self$nodes * scale else
      sweep(self$nodes, 2, scale, `*`)
    invisible(self)
  },

  #' @description
  #' Rotates the node grid of a \code{GaussQuad} object.
  #' @param rotation Rotation to apply to each node in the grid.
  rotate = function (rotation) {
    self$nodes <- if (is.vector(self$nodes)) self$nodes * c(rotation) else
      tcrossprod(self$nodes, rotation)
    invisible(self)
  },

  #' @description
  #' Shift the node grid of a \code{GaussQuad} object.
  #' @param shift Shift to apply to each dimension of the grid.
  shift = function (shift) {
    self$nodes <- if (is.vector(self$nodes)) self$nodes + shift else
      sweep(self$nodes, 2, shift, `+`)
    invisible(self)
  },

  #' @description
  #' Computes the quadrature of the integrand \code{f}.
  #' @param f Integrand function.
  #' @return The computed quadrature.
  integrate = apply_integral,
  E = apply_integral) # alias for probabilist types
)


#' GaussianIntegrator: multivariate Gaussian quadrature
#'
#' R6 class for multivariate Gaussian expectation.
#' @export
GaussianIntegrator <- R6::R6Class("GaussianIntegrator",
  inherit = GaussQuad, public = list(
    #' @field mean Mean of Gaussian distribution
    mean = NULL,

    #' @field variance Variance of Gaussian distribution
    variance = NULL,

    #' @description
    #' Create a new \code{GaussianIntegrator} object.
    #' @param n Number of nodes in one dimension or \code{GaussQuad} object of
    #'          type "gaussian" to be copied.
    #' @param mean Mean of Gaussian distribution.
    #' @param variance Variance of Gaussian distribution. Can be an
    #'                 eigendecomposition, of class \code{eigen}.
    #' @param ... Remaining parameters passed to \code{GaussQuad} initializer.
    #' @return A new \code{GaussianIntegrator} object.
    initialize = function (n, mean, variance, ...) {
      if (class(n)[1] == "GaussQuad") {
        if (n$type != "gaussian") stop("invalid quadrature type: ", n$type)
        self$type <- "gaussian";
        self$nodes <- n$nodes; self$weights <- n$weights
      } else {
        super$initialize("gaussian", n, length(mean), ...)
      }
      ev <- if (any(class(variance) == "eigen")) ev else
        eigen(variance, symmetric = TRUE)
      self$mean <- mean; self$variance <- ev
      self$scale(sqrt(ev$values))$rotate(ev$vectors)$shift(mean)
    },

    #' @description
    #' Prints object information.
    print = function () {
      super$print(short = TRUE)
      invisible(self)
    })
)

GaussQuad <- GaussQuad$new
GaussianIntegrator <- GaussianIntegrator$new

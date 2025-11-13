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

#' Compute robust eigendecomposition of symmetric PD matrix \code{A}.
#' @param Symmetric positive-definite matrix to be decomposed
#' @param rank_tol Tolerance for rank determination: only eigenvalues > plus(\code{rank_tol} * max(eigenvalues)) are kept. Defaults to sqrt of machine precision.
#' @return \code{eigen} object with extra field \code{valid} indicating valid eigenvalues
#' @export
symm_eigen <- function (A, rank_tol = sqrt(.Machine$double.eps)) {
  ea <- eigen(A, symmetric = TRUE)
  ea$valid <- ea$values > max(rank_tol * ea$values[1], 0)
  if (!all(ea$valid)) {
    ea$vectors <- ea$vectors[, ea$valid, drop = FALSE]
    ea$values <- ea$values[ea$valid]
  }
  ea
}


#' Efficiently and safely computes the aggregate log-sum-exp of the vector \code{x}.
#' @param x Numeric vector.
#' @return \code{log(sum(exp(x)))}.
#' @export
lse <- function (x, na_rm = FALSE) {
  .Call(G_lse, as.double(x), na_rm, PACKAGE = "gaussquadr")
}


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
  "beta" = 5L, # w(x) = x^(alpha-1) * (1-x)^(beta-1) / B(alpha, beta) on (0, 1)
  "gamma" = 6L # w(x) = x^(alpha-1) * exp(-beta*x) /
  #                       (beta^alpha / Gamma(alpha)) on (0, inf)
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
  ord <- order(res$weights, decreasing = TRUE)
  res$nodes <- res$nodes[ord]; res$weights <- res$weights[ord]
  res
}

# TODO: filter by threshold
# TODO: iterator
multi_grid <- function (x, n = length(x)) {
  if (n == 1) {
    matrix(x, ncol = 1)
  } else {
    l <- length(x); g <- multi_grid(x, n - 1)
    cbind(kronecker(rep(1, l), g), kronecker(x, rep(1, l ^ (n - 1))))
  }
}

apply_integral <- function (f) {
  if (missing(f)) {
    gf <- if (is.vector(self$nodes)) self$nodes else t(self$nodes)
  } else {
    if (is.numeric(f)) {
      if (length(f) == 1) return(f * sum(self$weights))
      gf <- f
    } else { # function?
      gf <- if (is.vector(self$nodes)) sapply(self$nodes, f) else
        apply(self$nodes, 1, f)
    }
  }
  if (is.vector(gf)) sum(gf * self$weights, na.rm = TRUE) else
    rowSums(sweep(gf, 2, self$weights, `*`), na.rm = TRUE)
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

  #' @field nodes_changed By location-scale transform using \code{axpy}?
  nodes_changed = FALSE,

  #' @field weights_changed By importance weighting using \code{reweight}?
  weights_changed = FALSE,

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
  initialize = function (type, n, d = 1, alpha = NA, beta = NA, endpts = NULL,
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
      "wigner" = {
        gq$weights <- gq$weights / (pi / 2)
      },
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
        gq$nodes <- .5 * (1 - gq$nodes)
        lz <- lbeta(alpha, beta) + (alpha + beta - 1) * log(2)
        gq$weights <- exp(log(gq$weights) - lz)
      }
    )

    if (prune_coef > 0) {
      th <- exp(d * ((1 - prune_coef) * log(gq$weights[1]) +
                     prune_coef * log(gq$weights[(n + 1) / 2])))
    }
    ng <- matrix(multi_grid(gq$nodes, d), ncol = d)
    w <- apply(multi_grid(gq$weights, d), 1, prod)
    if (prune_coef > 0) {
      ng <- ng[w >= th, ]; w <- w[w >= th]
    }
    if (prune_eps) {
      eps <- .Machine$double.eps * sum(w)
      ng <- ng[w >= eps, , drop = FALSE]; w <- w[w >= eps]
    }
    self$type <- type; self$nodes <- ng; self$weights <- w
    self$alpha <- alpha; self$beta <- beta
  },

  #' @description
  #' Prints object information.
  print = function () {
    nd <- dim(self$nodes)
    if (length(nd) == 0) {
      n <- length(self$nodes); d <- 1
    } else {
      n <- nd[1]; d <- nd[2]
    }
    cat(paste0(class(self)[1], " of type `", self$type, "`"))
    if (self$type == "jacobi" || self$type == "laguerre" ||
        self$type == "gamma" || self$type == "beta")
      cat(paste0(", alpha = ", self$alpha))
    if (self$type == "jacobi" || self$type == "gamma" || self$type == "beta")
      cat(paste0(", beta = ", self$beta))
    cat(paste0(" [ ", n, " nodes"))
    if (d > 1) cat(paste0(" x ", d, " dims"))
    cat(" ]")
    if (self$nodes_changed) cat("[G]")
    if (self$weights_changed) cat("[W]")
    cat("\n")
    invisible(self)
  },

  #' @description
  #' Applies a location-scale transformation \code{A * x + y} to each node
  #' \code{x} in the grid.
  #' @param location Shift representing the location change.
  #' @param scale Square matrix, vector, or \code{symm_eigen} decomposition representing the scale.
  #' @param squared Is the scale squared? For probabilist kernels, it means
  #' that the scale is a variance parameter. In this case, the "sqrt" of the
  #' scale will be applied first.
  location_scale = function (location = NULL, scale = NULL, squared = FALSE) {
    tg <- self$nodes
    if (!is.null(scale)) {
      if (is.vector(scale)) {
        if (squared) scale <- sqrt(scale)
        sw <- abs(prod(scale))
      } else {
        if (!any(class(scale) == "eigen")) # eigen-decompose?
          scale <- symm_eigen(scale)
        if (squared) scale$values <- sqrt(scale$values)
        valid <- scale$valid
        sw <- abs(prod(scale$values)) # Jacobian: abs det
        scale <- sweep(scale$vectors, 2, scale$values, `*`) # restore
      }
      tg <- if (is.vector(tg)) {
        tg * c(scale)
      } else {
        if (is.vector(scale)) sweep(tg, 2, scale, `*`) else
          tcrossprod(if (all(valid)) tg else tg[, valid], scale)
      }
      self$weights <- self$weights * sw # reweight from Jacobian
    }
    if (!is.null(location)) {
      tg <- if (is.vector(tg)) tg + location else
        sweep(tg, 2, location, `+`)
    }
    self$nodes <- tg; self$nodes_changed <- TRUE
    invisible(self)
  },

  #' @description
  #' Reweights by multiplying importance weights and optionally normalizes.
  #' @param iw Importance weights; if a function, apply to nodes first
  #' @param log_weights Are the importance weights in log scale?
  #' @param normalize Should the new weights be normalized to sum to one?
  reweight = function (iw, ..., log_weights = FALSE, normalize = FALSE) {
    if (is.function(iw)) iw <- iw(self$nodes, ...)
    if (log_weights) {
      lw <- iw + log(self$weights)
      if (normalize) lw <- lw - lse(lw)
      w <- exp(lw)
    } else {
      w <- iw * self$weights
      if (normalize) w <- w / sum(w)
    }
    self$weights <- w; self$weights_changed <- TRUE
    invisible(self)
  },

  #' @description
  #' Computes the quadrature of the integrand \code{f}.
  #' @param f Integrand function.
  #' @return The computed quadrature.
  integrate = apply_integral,
  E = apply_integral) # alias for probabilist types
)


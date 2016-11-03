##' Approximate a function by a spline, refined to a given tolerance
##'
##' @title Approximate function with a spline
##' @param target A function to approximate
##'
##' @param a The lower bound of the interpolation
##'
##' @param b The upper bound of the interpolation
##'
##' @param ... Additional arguments passed through to \code{target}
##'
##' @param method The method of interpolation to used; passed through
##'   to \code{\link{splinefun}}
##'
##' @param n_base The number of equally spaced points to start with
##'
##' @param max_depth The number of times the base points will be
##'   subdivided, at most.
##'
##' @param atol The absolute tolerance; used to determine how the
##'   approximation is doing for values with small absolute value.
##'
##' @param rtol The relative tolerance' used to determine how the
##'   approximation is doing for values with large absolute value.
##'
##' @param verbose Be verbose when fitting the spline?  This might be
##'   useful on expensive functions.
##'
##' @param tol The default value for \code{atol} and \code{rtol}; use
##'   this to set both at once.
##'
##' @param target_vectorised Flag indicating if \code{target} can
##'   accept a vector of \code{x} values; the default, \code{TRUE},
##'   assumes that it can.  If your function requires each \code{x}
##'   value separately set this to \code{FALSE}.
##'
##' @param inverse Indicates if the \emph{inverse} approximate
##'   function is required.  This is useful if you can use
##'   \code{target} to map x to y but you want the inverse mapping
##'   back.  \code{a} and \code{b} will be the domain of \code{x}
##'   still (so the range of \code{y}, and the domain of the returned
##'   function).
##'
##' @export
##' @author Rich FitzJohn
##' @importFrom stats splinefun
approximate <- function(target, a, b,
                        ...,
                        method = "fmm",
                        n_base = 17L, max_depth = 16L,
                        atol = tol, rtol = tol,
                        verbose = FALSE,
                        tol = sqrt(.Machine$double.eps),
                        target_vectorised = TRUE, inverse = FALSE) {
  if (!is.finite(a) || !is.finite(b)) {
    stop("Bounds must be finite and non-missing")
  }
  if (a >= b) {
    stop("Impossible bounds")
  }
  if (!target_vectorised) {
    target_scalar <- target
    target <- function(x, ...) {
      vapply(x, target_scalar, numeric(1), ...)
    }
  }

  dx <- (b - a) / (n_base - 1)
  dx_min <- dx / (2^max_depth)

  make_spline <- function(x, y) {
    if (inverse) {
      spline <- splinefun(yy, xx, method = method)
    } else {
      spline <- splinefun(xx, yy, method = method)
    }
  }

  ## Kick off:
  xx <- seq(a, b, length.out = n_base)
  yy <- target(xx, ...)
  zz <- rep(c(FALSE, TRUE), c(1, n_base - 1L))

  spline <- make_spline(xx, yy)

  while (any(zz)) {
    if (verbose) {
      message(sprintf("refining %d intervals", sum(zz)))
    }
    dx <- dx / 2
    if (dx < dx_min) {
      stop("function is as refined as currently possible")
    }

    x_mid <- xx[zz] - dx
    y_mid <- target(x_mid, ...)

    if (inverse) {
      p_mid <- spline(y_mid)
      z_mid <- !check_err(x_mid, p_mid, atol, rtol)
    } else {
      p_mid <- spline(x_mid)
      z_mid <- !check_err(y_mid, p_mid, atol, rtol)
    }
    zz[zz] <- z_mid

    ## Might be worth working out what order this should be in here;
    ## that might help speed up the spline calculations.
    xx <- c(xx, x_mid)
    yy <- c(yy, y_mid)
    zz <- c(zz, z_mid)

    spline <- make_spline(xx, yy)
  }

  spline
}

check_err <- function(y_true, y_pred, atol, rtol) {
  abs(y_true - y_pred) < atol | abs(1 - y_pred / y_true) < rtol
}

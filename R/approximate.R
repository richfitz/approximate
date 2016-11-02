approximate <- function(target, a, b,
                        ...,
                        method = "fmm",
                        n_base = 17L, max_depth = 16L,
                        atol = tol, rtol = tol,
                        verbose = TRUE,
                        tol = sqrt(.Machine$double.eps),
                        target_vectorised = TRUE) {
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

  ## Kick off:
  xx <- seq(a, b, length.out = n_base)
  yy <- target(xx)
  zz <- rep(c(FALSE, TRUE), c(1, n_base - 1L))

  spline <- splinefun(xx, yy, method = method)

  while (any(zz)) {
    if (verbose) {
      message(sprintf("refining %d intervals", sum(zz)))
    }
    dx <- dx / 2
    if (dx < dx_min) {
      stop("function is as refined as currently possible")
    }

    x_mid <- xx[zz] - dx
    y_mid <- target(x_mid)
    p_mid <- spline(x_mid)

    ## plot(xx, yy)
    ## points(x_mid, y_mid, col="red")
    ## points(x_mid, p_mid, pch=19, cex = .5)
    ## segments(x_mid, y_mid, x_mid, p_mid, col="grey")

    z_mid <- !check_err(y_mid, p_mid, atol, rtol)
    zz[zz] <- z_mid

    ## Might be worth working out what order this should be in here;
    ## that might help speed up the spline calculations.
    xx <- c(xx, x_mid)
    yy <- c(yy, y_mid)
    zz <- c(zz, z_mid)

    spline <- splinefun(xx, yy, method = method)
  }

  spline
}

check_err <- function(y_true, y_pred, atol, rtol) {
  abs(y_true - y_pred) < atol | abs(1 - y_pred / y_true) < rtol
}

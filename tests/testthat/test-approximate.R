context("approximate")

test_that("approximate", {
  target <- function(x) sin(2 * x)
  r <- c(0, 2 * pi)

  s <- approximate(target, r[1], r[2], tol = 1e-6)

  expect_is(s, "function")

  e <- environment(s)
  x <- e$z$x
  y <- e$z$y
  expect_equal(length(x), 241)

  expect_identical(y, target(x))

  x_mid <- (x[-1] + x[-length(x)]) / 2
  y_mid <- target(x_mid)
  z_mid <- s(x_mid)

  expect_equal(z_mid, y_mid, tolerance = 1e-7)
  err <- pmax(abs(z_mid - y_mid), abs(1 - z_mid / y_mid))
  expect_true(all(err < 1e-6))
})

test_that("non-vectorised", {
  target <- function(x) sin(2 * x)
  target2 <- function(x) {
    if (length(x) != 1L) {
      stop("I require scalar x")
    }
    target(x)
  }

  r <- c(0, 2 * pi)

  s1 <- approximate(target, r[1], r[2], tol = 1e-6)
  expect_error(approximate(target2, r[1], r[2], tol = 1e-6),
               "I require scalar x")
  s2 <- approximate(target, r[1], r[2], tol = 1e-6, target_vectorised = FALSE)

  xx <- seq(r[[1]], r[[2]], length.out = 1000)
  expect_equal(s1(xx), s2(xx))
})

test_that("error cases", {
  target <- runif
  expect_error(approximate(target, 1, 0), "Impossible bounds")
  expect_error(approximate(target, 0, NA), "Bounds must be finite")
  expect_error(approximate(target, NA, 1), "Bounds must be finite")
  expect_error(approximate(target, 0, Inf), "Bounds must be finite")
  expect_error(approximate(target, Inf, 1), "Bounds must be finite")
  expect_error(approximate(target, 0, 1, max_depth = 3),
               "function is as refined")
})

test_that("inverse", {
  f <- function(x) x * x
  g <- approximate(f, 0, 5, inverse = TRUE, verbose = TRUE,
                   atol = 1e-4, rtol = 1e-4)
  xx <- seq(0, f(5), length.out = 101)
  yy <- g(xx)
  expect_lt(max(abs(yy - sqrt(xx))), 1e-5)
})

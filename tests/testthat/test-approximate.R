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

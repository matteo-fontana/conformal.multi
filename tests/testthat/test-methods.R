make_data = function(n = 40, q = 2, seed = 1) {
  set.seed(seed)
  x = matrix(stats::rnorm(n * 2), ncol = 2)
  y = cbind(x[, 1] + stats::rnorm(n), -x[, 2] + 2 * stats::rnorm(n))
  if (q == 3) y = cbind(y, x[, 1] - x[, 2] + stats::rnorm(n))
  list(x = x, y = y, x0 = matrix(c(0, 0, 1, -1), ncol = 2, byrow = TRUE))
}

manual_split = function(d, training, alpha, score, s_type) {
  fun = lm_multi()
  calibration = setdiff(seq_len(nrow(d$x)), training)
  l = length(calibration)
  out = fun$train.fun(d$x[training, ], d$y[training, ])
  res = d$y - fun$predict.fun(out, d$x)
  s = computing_s_regression(res[training, ], type = s_type, alpha = alpha, tau = 1)
  resc = t(t(res[calibration, ]) / s)
  if (score == "l2") rho = sqrt(rowSums(resc^2))
  if (score == "max") rho = apply(abs(resc), 1, max)
  if (score == "mahalanobis") {
    sigma = stats::cov(t(t(res[training, ]) / s))
    rho = sqrt(stats::mahalanobis(resc, rep(0, ncol(resc)), sigma))
  }
  k = sort(rho)[ceiling((l + 1) * (1 - alpha))]
  hw = k * s
  if (score == "mahalanobis") hw = hw * sqrt(diag(sigma))
  pred = fun$predict.fun(out, d$x0)
  list(lo = pred - matrix(hw, nrow(pred), length(hw), byrow = TRUE),
       up = pred + matrix(hw, nrow(pred), length(hw), byrow = TRUE), k = k)
}

test_that("split conformal boxes match the definition of each score", {
  d = make_data()
  fun = lm_multi()
  training = 1:20
  for (score in c("l2", "max", "mahalanobis")) {
    for (s_type in c("identity", "st-dev")) {
      out = conformal.multidim.split(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun,
                                     alpha = 0.1, split = training,
                                     score = score, s_type = s_type)
      ref = manual_split(d, training, 0.1, score, s_type)
      expect_equal(unname(out$lo), unname(ref$lo))
      expect_equal(unname(out$up), unname(ref$up))
      expect_equal(out$k_s, ref$k)
    }
  }
})

test_that("split conformal accepts mad functions", {
  d = make_data()
  fun = lm_multi()
  out = conformal.multidim.split(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun,
                                 seed = 1, mad.train.fun = fun$train.fun,
                                 mad.predict.fun = fun$predict.fun)
  expect_equal(dim(out$lo), c(2, 2))
  expect_true(all(out$lo <= out$up))
})

test_that("interval.build returns the convex hull of the retained set", {
  yyy = c(0, 0.5, 3, 1, 1.5, 4) # lower bounds, then upper bounds
  expect_equal(interval.build(yyy, B = 3, tr = 0.5), c(0, 4))
  expect_equal(interval.build(yyy, B = 3, tr = 1.5), c(0.5, 1))
  expect_true(all(is.na(interval.build(yyy, B = 3, tr = 2.5))))
})

test_that("multi split, jackknife+ and full conformal return regions", {
  d = make_data()
  fun = lm_multi()
  for (aggregation in c("componentwise", "depth")) {
    out = conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun,
                                    B = 5, seed = 1, aggregation = aggregation)
    expect_equal(dim(out$lo), c(2, 2))
    expect_true(all(out$lo <= out$up))
  }
  jack = conformal.multidim.jackplus(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun)
  expect_equal(dim(jack$lo), c(2, 2))
  expect_true(all(jack$lo <= jack$up))
  full = conformal.multidim.full(d$x, d$y, d$x0[1, , drop = FALSE], fun$train.fun,
                                 fun$predict.fun, num.grid.pts.dim = 10)
  expect_length(full$valid_points, 1)
  expect_true(nrow(full$valid_points[[1]]) > 0)
})

test_that("bikeMi data are available", {
  expect_equal(dim(conformalInference.multi::bikeMi), c(41, 6))
})

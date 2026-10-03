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
                                    B = 20, seed = 1, aggregation = aggregation)
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

test_that("depth aggregation returns non-degenerate boxes and rejects too few vertices", {
  d = make_data(n = 60)
  fun = lm_multi()
  out = conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun,
                                  B = 20, seed = 1, aggregation = "depth")
  expect_true(all(out$up > out$lo))
  expect_error(conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun,
                                         fun$predict.fun, B = 5, seed = 1,
                                         aggregation = "depth"))
})

test_that("multi split: scalar training_size is recycled, seeds are reproducible", {
  d = make_data(n = 60)
  fun = lm_multi()
  a = conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun,
                                B = 4, seed = 1, training_size = 0.7)
  b = conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun,
                                B = 4, seed = 1, training_size = rep(0.7, 4))
  expect_equal(a$lo, b$lo)
  set.seed(3)
  c1 = conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun, B = 4)
  set.seed(3)
  c2 = conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun, B = 4)
  expect_equal(c1$up, c2$up)
})

test_that("mad functions with negative predictions give valid boxes", {
  d = make_data()
  fun = lm_multi()
  out = conformal.multidim.split(d$x, d$y, rbind(d$x0, c(10, 10)), fun$train.fun,
                                 fun$predict.fun, seed = 1,
                                 mad.train.fun = fun$train.fun,
                                 mad.predict.fun = function(out, newx)
                                   -abs(fun$predict.fun(out, newx)) - 0.1)
  expect_true(all(out$lo <= out$up))
})

test_that("jackknife+ works with elastic net and data frames", {
  skip_if_not_installed("glmnet")
  d = make_data()
  fun = elastic.funs()
  out = conformal.multidim.jackplus(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun)
  expect_equal(dim(out$up), c(2, 2))
  lfun = lm_multi()
  dfx = as.data.frame(d$x)
  pf = function(out, newx) lfun$predict.fun(out, as.matrix(newx))
  tf = function(x, y) lfun$train.fun(as.matrix(x), y)
  s1 = conformal.multidim.split(dfx, d$y, as.data.frame(d$x0), tf, pf, split = 1:20)
  s2 = conformal.multidim.split(d$x, d$y, d$x0, lfun$train.fun, lfun$predict.fun,
                                split = 1:20)
  expect_equal(unname(s1$up), unname(s2$up))
})

test_that("full conformal with q = 1 can be plotted, s_type is the last argument", {
  d = make_data()
  fun = lm_multi()
  full = conformal.multidim.full(d$x, d$y[, 1, drop = FALSE], d$x0[1, , drop = FALSE],
                                 fun$train.fun, fun$predict.fun, num.grid.pts.dim = 10)
  expect_s3_class(plot_multidim(full)[[1]], "ggplot")
  expect_equal(tail(names(formals(conformal.multidim.full)), 1), "s_type")
})

test_that("randomized split returns the whole space when the quantile index exceeds l", {
  d = make_data(n = 20)
  fun = lm_multi()
  found = FALSE
  for (st in 1:200) {
    tau = {set.seed(st); stats::runif(1)}
    if (ceiling(10 + tau - 11 * 0.08) > 10) { found = TRUE; break }
  }
  skip_if_not(found)
  expect_warning(out <- conformal.multidim.split(d$x, d$y, d$x0, fun$train.fun,
                                                 fun$predict.fun, alpha = 0.08,
                                                 split = 1:10, randomized = TRUE,
                                                 seed_tau = st))
  expect_true(all(is.infinite(out$up)))
})

test_that("bikeMi interaction column equals we * rain", {
  expect_equal(conformalInference.multi::bikeMi$we_rain,
               conformalInference.multi::bikeMi$we * conformalInference.multi::bikeMi$rain)
})

test_that("jackknife+ with method = 'max' equals the jackknife+ interval for q = 1", {
  d = make_data(n = 30)
  fun = lm_multi()
  y1 = d$y[, 1, drop = FALSE]
  out = conformal.multidim.jackplus(d$x, y1, d$x0, fun$train.fun, fun$predict.fun,
                                    alpha = 0.1, method = "max")
  n = nrow(d$x)
  loo = lapply(1:n, function(i) fun$train.fun(d$x[-i, ], y1[-i, , drop = FALSE]))
  R = sapply(1:n, function(i) abs(y1[i, 1] - fun$predict.fun(loo[[i]], d$x[i, , drop = FALSE])))
  mu0 = sapply(1:n, function(i) fun$predict.fun(loo[[i]], d$x0))
  lo = apply(mu0, 1, function(m) sort(m - R)[floor(0.1 * (n + 1))])
  up = apply(mu0, 1, function(m) sort(m + R)[ceiling(0.9 * (n + 1))])
  expect_equal(c(out$lo), lo)
  expect_equal(c(out$up), up)
  out2 = conformal.multidim.jackplus(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun,
                                     method = "max", scale = c(1, 2))
  expect_true(all(out2$lo < out2$up))
  expect_error(conformal.multidim.jackplus(d$x, d$y, d$x0, fun$train.fun,
                                           fun$predict.fun, method = "max", scale = 1))
})

test_that("the prediction functions do not change the kind of random number generator", {
  d = make_data()
  fun = lm_multi()
  kind = RNGkind()
  conformal.multidim.jackplus(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun)
  expect_identical(RNGkind(), kind)
  conformal.multidim.msplit(d$x, d$y, d$x0, fun$train.fun, fun$predict.fun, B = 4, seed = 1)
  expect_identical(RNGkind(), kind)
  conformal.multidim.full(d$x, d$y, d$x0[1, , drop = FALSE], fun$train.fun,
                          fun$predict.fun, num.grid.pts.dim = 5)
  expect_identical(RNGkind(), kind)
})

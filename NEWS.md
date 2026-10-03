# conformalInference.multi 1.1.3

## Restored

* `conformal.multidim.jackplus()` and the `bikeMi` data set, which were
  available in version 1.1.1 and were not included in version 1.1.2. The
  jackknife+ function uses the argument conventions of version 1.1.2 and
  retains the `ceiling(2n(1 - alpha/2))` deepest leave-one-out vectors, as in
  `conformalInference.fd` (version 1.1.1 retained one more vector). No
  finite-sample coverage guarantee is available for this extension; its
  coverage can be below `1 - alpha`, in particular when `q` is large relative
  to `n`.

## New features

* `conformal.multidim.msplit()` gains the argument `aggregation`.
  `"componentwise"` (default) is the aggregation of version 1.1.2;
  `"depth"` is the depth-based aggregation of version 1.1.1, which requires
  `floor(2 * tau * B) >= 2`. The default of `tau` now depends on
  `aggregation` (`1 - (B + 1)/(2B)` and `0.1` respectively).
  `mad.train.fun` and `mad.predict.fun` are passed to the split conformal
  function again. A scalar `training_size` is used for every split.
* `conformal.multidim.full()` gains the argument `s_type` (`"st-dev"`,
  default, or `"identity"`), as in version 1.1.1; it is the last argument, so
  that positional calls written for version 1.1.2 still work.
* `plot_multidim()` accepts the output of every prediction function: it
  dispatches to `plot_multidim_full()` for full conformal output and displays
  interval midpoints when no point prediction is available.
  `plot_multidim_full()` also handles a univariate response.

## Changes in behaviour

* `conformal.multidim.msplit()`: with the default `seed = FALSE` the seeds of
  the `B` splits are drawn from the current random number stream (so that
  `set.seed()` makes the result reproducible), instead of resetting the
  stream with `set.seed(0)`. With a numeric `seed`, the `b`-th split uses
  `seed + b` as before; `seed = 0` reproduces the default of version 1.1.2.
* The prediction functions no longer change the option
  `future.rng.onMisuse`, and use parallel-safe random numbers
  (`future.seed = TRUE`).
* `conformal.multidim.split()` and `conformal.multidim.jackplus()` pass `x`
  and `x0` to `train.fun` and `predict.fun` as they are given (matrices or
  data frames), and accept predictions returned as data frames.

## Bug fixes

* `conformal.multidim.split()`: with `score = "l2"` the half-widths of the
  prediction box were computed from the squared Euclidean norm instead of the
  norm. With `score = "mahalanobis"` the box was not the bounding box of the
  Mahalanobis ball; the covariance matrix is now estimated on the training
  set and the box is `pred -/+ k_s * s * sqrt(diag(Sigma))`.
* `conformal.multidim.split()`: user-supplied `split` failed because the size
  of the calibration set was undefined; `seed_tau` was not checked;
  `mad.train.fun`/`mad.predict.fun` could not be used (the absolute value of
  the predicted scale is now used). In the randomized version, a draw of tau
  for which the region is the whole space now returns `-Inf`/`Inf` with a
  warning instead of an error.
* `conformal.multidim.full()` failed with the default value of `score` and
  with a data frame `x0`; a `train.fun` or `mad.train.fun` using its `out`
  argument failed or received a model fitted on different data (the models
  are now fitted from scratch for each trial value).
* `conformal.multidim.jackplus()` failed with `elastic.funs()`,
  `lasso.funs()` and `ridge.funs()` when glmnet was not loaded: the
  leave-one-out models are now fitted and evaluated on the same worker.
* `conformal.multidim.msplit()` with `aggregation = "componentwise"`: when the
  set of values covered by more than `tau * B` intervals is not an interval,
  its convex hull is returned (previously only the last interval was kept);
  an empty set is returned as `NA`. A warning is given when `lambda` exceeds
  `tau * B - floor(tau * B)`, for which the coverage guarantee does not hold.
* `lm_multi()`: the coefficients of aliased columns (rank-deficient training
  set, e.g. an indicator that is constant on the training half of a split)
  are set to zero, as in `predict.lm()`, instead of producing `NA`
  predictions.
* `bikeMi`: the column `we_rain` differed from `we * rain` on 26 February,
  28 February and 4 March 2016.
* `plot_multidim()`: `same.scale = TRUE` reversed the y axis.
* `plot_multidim_full()` failed when the response had column names.
* The prediction functions restore the user's `future` plan on exit.

# conformalInference.multi 1.1.2

* Argument names harmonised (`training_size`, `s_type`, `seed_tau`).

# conformalInference.multi 1.1.1

* Unified plots in `plot_multidim`.
* Added example data with BikeMi.

# conformalInference.multi 1.1.0

* Added `conformal.multidim.jackplus()`.
* Modified the input of `conformal.multidim.split()`,
  `conformal.multidim.msplit()` and `conformal.multidim.full()` to make it
  consistent.

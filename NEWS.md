# conformalInference.multi 1.1.3

## Restored

* `conformal.multidim.jackplus()` and the `bikeMi` data set, which were
  available in version 1.1.1 and were not included in version 1.1.2. The
  jackknife+ function now uses the argument conventions of version 1.1.2 and
  retains the `ceiling(2n(1 - alpha/2))` deepest leave-one-out vectors, as in
  `conformalInference.fd`.

## New features

* `conformal.multidim.msplit()` gains the argument `aggregation`.
  `"componentwise"` (default) is the aggregation of version 1.1.2;
  `"depth"` is the depth-based aggregation of version 1.1.1. The default of
  `tau` now depends on `aggregation` (`1 - (B + 1)/(2B)` and `0.1`
  respectively). `mad.train.fun` and `mad.predict.fun` are passed to the split
  conformal function again.
* `conformal.multidim.full()` gains the argument `s_type`
  (`"st-dev"`, default, or `"identity"`), as in version 1.1.1.
* `plot_multidim()` accepts the output of every prediction function: it
  dispatches to `plot_multidim_full()` for full conformal output and displays
  interval midpoints when no point prediction is available.

## Bug fixes

* `conformal.multidim.split()`: with `score = "l2"` the half-widths of the
  prediction box were computed from the squared Euclidean norm instead of the
  norm. With `score = "mahalanobis"` the box was not the bounding box of the
  Mahalanobis ball; the covariance matrix is now estimated on the training
  set and the box is `pred -/+ k_s * s * sqrt(diag(Sigma))`.
* `conformal.multidim.split()`: user-supplied `split` failed because the size
  of the calibration set was undefined; `seed_tau` was not passed to the
  checks; `mad.train.fun`/`mad.predict.fun` could not be used.
* `conformal.multidim.full()` failed with the default value of `score`.
* `conformal.multidim.msplit()` with `aggregation = "componentwise"`: when the
  set of values covered by more than `tau * B` intervals is not an interval,
  its convex hull is returned (previously only the last interval was kept);
  an empty set is returned as `NA`.
* `lm_multi()`: the coefficients of aliased columns (rank-deficient training
  set, e.g. an indicator that is constant on the training half of a split)
  are set to zero, as in `predict.lm()`, instead of producing `NA`
  predictions.
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

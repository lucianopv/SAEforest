# Closed-form smearing indicators ---------------------------------------------
#
# The smearing set of a domain is the Minkowski sum {preds[k] + res[j]} over all
# population units k in the domain and all sample residuals j. Materialising it
# costs O(N_i * n) time and memory; Mean, Hcr and Pgap can instead be read off
# the sorted residual vector in O(N_i log n) time and O(N_i + n) memory.
#
# Every other indicator (Gini, the quantiles, Qsr) needs the sorted smearing set
# itself, and custom indicators are arbitrary functions of it, so those still
# materialise -- via the exact same construction as before, preserving parity.

# Indicators computable without materialising the Minkowski sum: Mean/Hcr/Pgap
# in closed form, the quantiles by order-statistic selection. Gini and Qsr need
# the sorted set itself and are deliberately absent.
smear_quantile_probs <- c(Quant10 = 0.10, Quant25 = 0.25, Median = 0.50,
                          Quant75 = 0.75, Quant90 = 0.90)
smear_fast_indicators <- c("Mean", "Hcr", "Pgap", names(smear_quantile_probs))

# Selection costs ~100 counting passes per order statistic, so it only overtakes
# materialise-and-sort once the smearing set is large. Measured break-even is
# around 1e5 cells (1.05x there, 12.8x at 2e6, 29.3x at 2e7).
smear_select_min_cells <- 1e5

smear_indicators <- function(preds, res, thresh, custom = NULL,
                             select.indicator = NULL,
                             transformation = c("none", "log")) {

  transformation <- match.arg(transformation)
  if (transformation == "log") {
    return(smear_indicators_log(preds, res, thresh, custom, select.indicator))
  }

  wants_quant <- any(select.indicator %in% names(smear_quantile_probs))
  use_fast <- is.null(custom) &&
    !is.null(select.indicator) &&
    all(select.indicator %in% smear_fast_indicators) &&
    (!wants_quant || length(preds) * length(res) >= smear_select_min_cells)

  if (!use_fast) {
    # Fallback: verbatim the original construction, so the flattened element
    # order (and hence floating-point summation order) is unchanged.
    smear_i <- matrix(rep(res, length(preds)), nrow = length(preds),
                      ncol = length(res), byrow = TRUE)
    smear_i <- smear_i + preds
    return(calc_indicat(c(smear_i), threshold = thresh, custom = custom,
                        select.indicator = select.indicator))
  }

  N <- length(preds)
  rs <- sort(res)
  n <- length(rs)

  # Hcr = mean(Y < thresh). For each k, count residuals strictly below
  # (thresh - preds[k]); findInterval(left.open = TRUE) returns #{rs < x}.
  # NOTE: this tests res[j] < thresh - preds[k] rather than
  # preds[k] + res[j] < thresh. The two can disagree by one ULP for a cell
  # sitting within rounding distance of the threshold; for continuous data
  # that is a measure-zero event, but it is a deviation, not an identity.
  n_below <- findInterval(thresh - preds, rs, left.open = TRUE)

  # Pgap = mean over cells with y >= 0 of (y < thresh) * (thresh - y)/thresh.
  # calc_indicat() maps negative cells to NA and averages with na.rm = TRUE, so
  # they leave the denominator as well as the numerator. For unit k the cells
  # that count are those with rs[j] in [-preds[k], thresh - preds[k]).
  cs <- c(0, cumsum(rs))                                  # prefix sums of rs
  n_neg <- findInterval(-preds, rs, left.open = TRUE)      # #{rs < -preds[k]}
  hi <- pmax(n_below, n_neg)                              # guard empty bands
  cnt_band <- hi - n_neg
  sum_band <- cs[hi + 1L] - cs[n_neg + 1L]
  n_valid <- sum(n - n_neg)                               # cells with y >= 0
  Pgap <- if (n_valid == 0) NaN else {
    sum(cnt_band * (thresh - preds) - sum_band) / thresh / n_valid
  }

  indicators <- cbind(
    Mean = mean(preds) + mean(rs),
    Hcr = sum(as.numeric(n_below)) / (N * n),
    Pgap = Pgap
  )

  # Quantiles by order-statistic selection, reproducing calc_indicat()'s
  # interpolation (type 7) on the same order statistics it would have indexed
  # out of the fully sorted smearing set.
  if (wants_quant) {
    ps <- sort(preds)
    m <- N * n
    qsel <- intersect(select.indicator, names(smear_quantile_probs))
    qvals <- vapply(smear_quantile_probs[qsel], function(prob) {
      h <- (m - 1) * prob + 1
      lo <- floor(h); hi <- ceiling(h)
      a <- minkowski_order_stat(ps, rs, lo)
      b <- if (hi == lo) a else minkowski_order_stat(ps, rs, hi)
      a + (h - lo) * (b - a)
    }, numeric(1))
    indicators <- cbind(indicators, t(qvals))
  }

  indicators[, select.indicator]
}

# r-th order statistic of the Minkowski sum {A[k] + B[j]} ----------------------
#
# A and B must be sorted ascending. The implied matrix M[k, j] = A[k] + B[j] is
# sorted along both axes, so #{cells <= v} is computable in O(|A| log |B|)
# without materialising it. Bisect on value to bracket the answer, then snap to
# the largest cell not exceeding the upper bound -- that cell IS the r-th
# smallest once the bracket has narrowed to adjacent doubles.
minkowski_order_stat <- function(A, B, r) {
  nA <- length(A); nB <- length(B)

  # #{cells <= v}. as.numeric() guards the integer sum: findInterval() returns
  # integers and nA * nB can exceed .Machine$integer.max for a large domain.
  count_le <- function(v) sum(as.numeric(findInterval(v - A, B)))
  # largest cell value not exceeding v
  snap <- function(v) {
    idx <- findInterval(v - A, B)
    keep <- idx >= 1L
    max(A[keep] + B[idx[keep]])
  }

  lo <- A[1L] + B[1L]
  hi <- A[nA] + B[nB]
  if (r <= 1L || count_le(lo) >= r) return(lo)
  if (r >= nA * nB) return(hi)

  # invariant: count_le(lo) < r <= count_le(hi)
  repeat {
    mid <- lo + (hi - lo) / 2
    if (mid <= lo || mid >= hi) break        # adjacent doubles, no progress
    if (count_le(mid) >= r) hi <- mid else lo <- mid
  }
  snap(hi)
}

# Log-path (multiplicative) smearing -------------------------------------------
#
# Cells are z[k, j] = exp(preds[k]) * exp(res[j]) = a[k] * b[j], with a, b > 0.
# point_nonLin() maps non-finite cells to NA before reducing, so calc_indicat()'s
# na.rm denominators count only the FINITE cells. The closed forms below assume
# every cell is finite; the guard below proves that cheaply (a, b > 0 and finite
# means max(a) * max(b) bounds every product) and otherwise falls back.
#
# Note z >= 0 always, so calc_indicat()'s negative-to-NA rule for Pgap is inert
# here -- unlike on the additive path.
smear_indicators_log <- function(preds, res, thresh, custom = NULL,
                                 select.indicator = NULL) {
  a <- exp(preds)
  b <- exp(res)

  fast_ok <- is.null(custom) &&
    !is.null(select.indicator) &&
    all(select.indicator %in% c("Mean", "Hcr", "Pgap")) &&
    all(is.finite(a)) && all(is.finite(b)) && all(a > 0) &&
    is.finite(max(a) * max(b))

  if (!fast_ok) {
    val <- c(outer(a, b, "*"))
    val[!is.finite(val)] <- NA
    return(calc_indicat(val, threshold = thresh, custom = custom,
                        select.indicator = select.indicator))
  }

  N <- length(a)
  bs <- sort(b)
  n <- length(bs)
  cs <- c(0, cumsum(bs))                       # prefix sums of sorted exp(res)

  # z < thresh  <=>  b[j] < thresh / a[k]   (a[k] > 0, so this preserves order).
  # As on the additive path this is a rearrangement, not a bit-identical test:
  # b[j] < thresh/a[k] and a[k]*b[j] < thresh can disagree by one ULP for a cell
  # within rounding distance of the threshold.
  n_below <- findInterval(thresh / a, bs, left.open = TRUE)
  sum_below <- cs[n_below + 1L]

  indicators <- cbind(
    Mean = mean(a) * mean(b),
    Hcr = sum(as.numeric(n_below)) / (N * n),
    Pgap = sum(n_below * thresh - a * sum_below) / thresh / (N * n)
  )
  indicators[, select.indicator]
}

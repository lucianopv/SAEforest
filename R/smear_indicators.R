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

# Indicators with a closed form over the Minkowski sum.
smear_fast_indicators <- c("Mean", "Hcr", "Pgap")

smear_indicators <- function(preds, res, thresh, custom = NULL,
                             select.indicator = NULL) {

  use_fast <- is.null(custom) &&
    !is.null(select.indicator) &&
    all(select.indicator %in% smear_fast_indicators)

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
    Hcr = sum(n_below) / (N * n),
    Pgap = Pgap
  )
  indicators[, select.indicator]
}

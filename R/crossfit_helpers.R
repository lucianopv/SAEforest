# Helpers for re_method = "crossfit" (MERFranger). Internal.

# The grouping variable named by a `random` formula fragment such as "(1|area)"
# or "(1|district)". Only single, non-nested groupings are supported: cross-fitting
# by area has no obvious meaning for crossed effects, so those are refused rather
# than silently folded on the first term.
grouping_factor_ <- function(data, random) {
  inside <- sub(".*\\|\\s*", "", sub("\\).*$", "", random))
  nm <- trimws(inside)
  if (grepl("[/:]", nm)) {
    stop('re_method = "crossfit" supports a single grouping factor; got "', nm, '".')
  }
  if (!nm %in% names(data)) {
    stop('re_method = "crossfit" cannot find the grouping variable "', nm,
         '" in `data`.')
  }
  factor(data[[nm]])
}

# Assign whole AREAS to folds -- never units, or an area would appear in both the
# training and the held-out side and the forest could memorise it after all,
# which is the very thing this is meant to prevent.
#
# Deterministic given the level order, so a fit is reproducible without consuming
# RNG (MERFranger's own seed governs the forests).
assign_area_folds_ <- function(grp, folds) {
  lv <- levels(grp)
  k <- max(2L, min(as.integer(folds), length(lv)))
  stats::setNames(rep_len(seq_len(k), length(lv)), lv)
}

# ranger's num.threads if the caller passed one, else 1. predict.ranger defaults
# to NULL, which ranger reads as "every core" -- unusable inside a forked worker.
dots_threads_ <- function(...) {
  dots <- list(...)
  if (!is.null(dots$num.threads)) dots$num.threads else 1
}

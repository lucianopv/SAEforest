# FROZEN reference, generated ONCE at 3250638 (2026-10-22), BEFORE the
# degenerate-case guards in MERFranger.R and adjust_ErrorSD_weights.R (a clamped
# EM iteration, never-out-of-bag rows, NA OOB predictions). Those guards are
# meant to change only the degenerate fits they target; every non-degenerate
# var.adjust fit must stay byte-identical, which test-MERFranger.R checks against
# this file. Do not regenerate unless a change is MEANT to move these fits.
devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))
saveRDS(varadjust_guard_fits(),
        testthat::test_path("fixtures", "varadjust_fits_ref.rds"))

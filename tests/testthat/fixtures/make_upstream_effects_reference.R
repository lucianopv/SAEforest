# FROZEN reference, generated ONCE at e261df1 (2026-10-22), BEFORE var.adjust
# predictions left out the effect model's intercept (Task 13a). That change must
# touch only the var.adjust path: the upstream effect model (Target ~ -1 + (1 |
# area)) has no fixed intercept, so every var.adjust = FALSE output must stay
# byte-identical, which test-ranef-intercept.R checks against this file. Do not
# regenerate unless a change is MEANT to move the upstream path.
devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))
saveRDS(upstream_effect_outputs(),
        testthat::test_path("fixtures", "upstream_effects_ref.rds"))

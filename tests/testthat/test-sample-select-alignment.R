# Regression suite for the sample_select() domain/size misalignment.
# See GridSAE docs/superpowers/specs/2026-08-24-sample-select-fix-and-hpc-efficiency-design.md

drawn_counts <- function(x, dName, doms) {
  as.integer(table(factor(as.character(x[[dName]]), levels = doms)))
}

test_that("the unsorted fixture really does diverge (guards the test itself)", {
  d <- tiny_saef_data_unsorted()
  first_appearance <- as.character(unique(d$pop[[d$dName]]))
  level_order <- names(split(seq_len(nrow(d$pop)), d$pop[[d$dName]]))
  expect_false(identical(first_appearance, level_order))
  # and the sample sizes must not all be equal, or a permutation is invisible
  expect_gt(length(unique(as.integer(table(d$smp[[d$dName]])))), 1L)
})

test_that("sample_select draws each domain's OWN sample size when pop is unsorted", {
  d <- tiny_saef_data_unsorted()
  doms <- sort(unique(as.character(d$smp[[d$dName]])))

  set.seed(11)
  got <- SAEforest:::sample_select(d$pop, smp = d$smp, dName = d$dName)

  expect_equal(drawn_counts(got, d$dName, doms),
               drawn_counts(d$smp, d$dName, doms))
})

test_that("a population domain with no sampled units draws nothing", {
  d <- tiny_saef_data_unsorted()
  doms <- sort(unique(as.character(d$smp[[d$dName]])))
  drop_dom <- doms[[1]]
  d$smp <- d$smp[as.character(d$smp[[d$dName]]) != drop_dom, , drop = FALSE]

  set.seed(12)
  got <- SAEforest:::sample_select(d$pop, smp = d$smp, dName = d$dName)

  expect_equal(sum(as.character(got[[d$dName]]) == drop_dom), 0L)
  expect_equal(nrow(got), nrow(d$smp))
})

# The pre-fix implementation, vendored verbatim so this parity test keeps its
# meaning after R/support_functions.R is edited.
sample_select_prefix_reference <- function(pop, smp, dName) {
  smpSizes <- table(smp[dName])
  smpSizes <- data.frame(
    smpidD = as.character(names(smpSizes)), n_smp = as.numeric(smpSizes),
    stringsAsFactors = FALSE
  )
  smpSizes <- dplyr::left_join(
    data.frame(idD = as.character(unique(pop[[dName]])), stringsAsFactors = FALSE),
    smpSizes, by = c("idD" = "smpidD")
  )
  smpSizes$n_smp[is.na(smpSizes$n_smp)] <- 0
  splitPop <- split(pop, pop[[dName]])
  stratSamp <- function(dfList, ns) {
    do.call(rbind, mapply(dfList, ns, FUN = function(df, n) {
      popInd <- seq_len(nrow(df))
      sel <- base::sample(popInd, n, replace = FALSE)
      df[sel, ]
    }, SIMPLIFY = F))
  }
  stratSamp(splitPop, smpSizes$n_smp)
}

test_that("on SORTED pop the fix is byte-identical to the pre-fix implementation", {
  # A fixture that is genuinely aligned (pop ordered by domain, so
  # first-appearance order == split()'s factor-level order) AND has unequal
  # per-domain sample sizes. Both matter: on aligned data the alignment fix is
  # a no-op, so any difference here is the index rewrite changing numbers --
  # and unequal sizes are what make a mis-ordering observable at all.
  # (tiny_saef_data() is NOT aligned, and its equal sizes of 8 would mask a
  # mis-ordering entirely, which is why it is not used here.)
  d <- tiny_saef_data_unsorted()
  # Sort on the factor itself, not as.character(): order() on a factor uses its
  # integer codes, i.e. split()'s own grouping order, independent of locale.
  # order(as.character(...)) instead sorts by the session's LC_COLLATE, which
  # need not agree with the factor's baked-in level order (observed on this
  # machine's en_US.UTF-8: "Bludenz" sorts before "Braunau am Inn" as text,
  # but the district factor's levels place "Braunau am Inn" first).
  d$pop <- d$pop[order(d$pop[[d$dName]]), , drop = FALSE]
  rownames(d$pop) <- NULL
  stopifnot(identical(as.character(unique(d$pop[[d$dName]])),
                      names(split(seq_len(nrow(d$pop)), d$pop[[d$dName]]))))

  set.seed(4242)
  ref <- sample_select_prefix_reference(d$pop, smp = d$smp, dName = d$dName)
  set.seed(4242)
  got <- SAEforest:::sample_select(d$pop, smp = d$smp, dName = d$dName)

  rownames(ref) <- NULL
  rownames(got) <- NULL
  expect_equal(got, ref)
})

test_that("a precomputed `groups` gives the same result as computing it inline", {
  d <- tiny_saef_data_unsorted()
  grp <- split(seq_len(nrow(d$pop)), d$pop[[d$dName]])

  set.seed(77)
  a <- SAEforest:::sample_select(d$pop, smp = d$smp, dName = d$dName)
  set.seed(77)
  b <- SAEforest:::sample_select(d$pop, smp = d$smp, dName = d$dName, groups = grp)

  rownames(a) <- NULL
  rownames(b) <- NULL
  expect_equal(b, a)
})

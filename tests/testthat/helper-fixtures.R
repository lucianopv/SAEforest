# Tiny deterministic subset of eusilcA for fast, reproducible tests.
tiny_saef_data <- function(n_dom = 4, per_dom = 8, seed = 2026) {
  testthat::skip_if_not_installed("emdi")
  set.seed(seed)
  data("eusilcA_smp", package = "emdi", envir = environment())
  data("eusilcA_pop", package = "emdi", envir = environment())
  doms <- sort(unique(as.character(eusilcA_smp$district)))[seq_len(n_dom)]

  smp <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(eusilcA_smp$district == d)
    eusilcA_smp[sample(rows, min(per_dom, length(rows))), ]
  }))
  pop <- eusilcA_pop[eusilcA_pop$district %in% doms, ]
  # cap population per domain to keep smearing cheap
  pop <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(pop$district == d)
    pop[rows[seq_len(min(40, length(rows)))], ]
  }))

  smp <- droplevels(smp); pop <- droplevels(pop)
  list(
    Y = smp$eqIncome,
    X = smp[, c("gender", "eqsize", "cash", "self_empl", "unempl_ben")],
    smp = smp[, c("district", "gender", "eqsize", "cash",
                  "self_empl", "unempl_ben", "eqIncome")],
    pop = pop[, c("district", "gender", "eqsize", "cash",
                  "self_empl", "unempl_ben")],
    dName = "district"
  )
}

# Same as tiny_saef_data(), plus a "subunit" column (2 subunits per district)
# nested within district, added to BOTH pop and smp -- for testing aggregate_to
# in either direction: finer than dName (dName="district", aggregate_to="subunit")
# or coarser than dName (dName="subunit", aggregate_to="district").
tiny_saef_data_with_subunit <- function(n_dom = 4, per_dom = 8, seed = 2026) {
  base <- tiny_saef_data(n_dom = n_dom, per_dom = per_dom, seed = seed)
  base$pop$subunit <- paste0(base$pop$district, "_", rep(1:2, length.out = nrow(base$pop)))
  base$smp$subunit <- paste0(base$smp$district, "_", rep(1:2, length.out = nrow(base$smp)))
  base
}

# Tiny deterministic aggData=TRUE fixture (aggregated-covariate / ELM-calibration path),
# using SAEforest's own bundled data (no emdi dependency). Restricted to in-sample-only
# domains so the out-of-sample augmentation branch of point_meanAGG() is never exercised.
tiny_agg_data <- function(seed = 2026, per_dom = 15) {
  data("eusilcA_smp", "eusilcA_popAgg", "popNsize", package = "SAEforest", envir = environment())

  doms <- c("Graz (Stadt)", "Linz (Stadt)", "Salzburg-Umgebung")
  covars <- c("eqsize", "cash", "self_empl", "unempl_ben")

  set.seed(seed)
  smp <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(as.character(eusilcA_smp$district) == d)
    eusilcA_smp[sample(rows, min(per_dom, length(rows))), ]
  }))
  smp <- droplevels(smp)

  # popAgg and popN must have identical row counts AND identical dName ordering
  # (error_checkFunctions.R checks all.equal(popnsize[[dName]], pop_data[[dName]])).
  popAgg <- eusilcA_popAgg[match(doms, as.character(eusilcA_popAgg$district)), ]
  popAgg <- droplevels(popAgg)
  rownames(popAgg) <- NULL

  popN <- popNsize[match(doms, as.character(popNsize$district)), ]
  rownames(popN) <- NULL

  list(
    Y = smp$eqIncome,
    X = smp[, covars],
    smp = smp[, c("district", covars, "eqIncome")],
    popAgg = popAgg[, c("district", covars)],
    popN = popN,
    dName = "district"
  )
}

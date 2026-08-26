# Wrapper function for point and MSE estimates for domain-level nonlinear indicators

SAEforest_nonLin <- function(Y,
                             X,
                             dName,
                             smp_data,
                             pop_data,
                             smearing = TRUE,
                             MSE = "none",
                             importance = "none",
                             initialRandomEffects = 0,
                             ErrorTolerance = 0.0001,
                             MaxIterations = 25,
                             B = 100,
                             var.adjust = FALSE,
                             B_adj = 100,
                             adj_tol = 0,
                             adj_tol_mse = NULL,
                             mse_tol = 0,
                             B_MC = 100,
                             transformation = c("none", "log"),
                             threshold = NULL,
                             custom_indicator = NULL,
                             aggregate_to = NULL,
                             na.rm = TRUE,
                             select.indicator = NULL,
                             out_call,
                             cores = 1,
                             worker.threads = 1L,
                             n_smear_residuals = NULL,
                             dedup_by = NULL,
                             adjust_mse = c("refit", "plugin"),
                             ...) {

  transformation <- match.arg(transformation, c("none", "log"))
  adjust_mse <- match.arg(adjust_mse, c("refit", "plugin"))

  if (na.rm == TRUE) {
    comp_smp <- complete.cases(smp_data)
    smp_data <- smp_data[comp_smp, ]
    Y <- Y[comp_smp]
    X <- X[comp_smp, ]
    pop_data <- pop_data[complete.cases(pop_data[,c(colnames(X),dName)]), ]
  }

  # make domain variable to character and sort data-sets
  smp_data[[dName]] <- factor(smp_data[[dName]], levels = unique(smp_data[[dName]]))
  pop_data[[dName]] <- factor(pop_data[[dName]], levels = unique(pop_data[[dName]]))

  # order data according to factors to ease MSE estimation
  order_in <- order(smp_data[[dName]])
  smp_data <- smp_data[order_in, ]
  X <- X[order_in, ]
  Y <- Y[order_in]
  pop_data <- pop_data[order(pop_data[[dName]]), ]

  # Validate the dedup key ONCE, up front, on the final pop_data (the stable
  # sort above preserves level grouping). A violated key must fail loudly here,
  # not silently produce wrong predictions downstream.
  if (!is.null(dedup_by)) check_dedup_by(pop_data, X, dedup_by)


  # Point and MSE estimates for domain-level indicators and unit-level data (smearing) ------
  if (smearing == TRUE) {

    # point estimation
    nonLin_preds <- point_nonLin(
      Y = Y,
      X = X,
      dName = dName,
      threshold = threshold,
      smp_data = smp_data,
      pop_data = pop_data,
      initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance,
      MaxIterations = MaxIterations,
      importance = importance,
      custom_indicator = custom_indicator,
      aggregate_to = aggregate_to,
      var.adjust = var.adjust,
      B_adj = B_adj,
      adj_tol = adj_tol,
      transformation = transformation,
      select.indicator = select.indicator,
      dedup_by = dedup_by,
      cores = cores,
      n_smear_residuals = n_smear_residuals,
      ...
    )

    # Plug-in K for the MSE replicates: the point fit's per-iteration K values,
    # estimated exactly as before, are handed to every replicate fit instead of
    # being re-estimated with B_adj inner forests per iteration. The point fit
    # itself is never plugged.
    K_fixed_mse <- NULL
    if (adjust_mse == "plugin" && MSE != "none") {
      K_fixed_mse <- unlist(nonLin_preds[[2]]$K)
      stopifnot(is.numeric(K_fixed_mse), length(K_fixed_mse) >= 1)
    }

    if(is.null(aggregate_to)){
      data_specs <- sae_specs(dName = dName, cns = pop_data, smp = smp_data)
    } else{
      data_specs <- sae_specs(dName = aggregate_to, cns = pop_data, smp = smp_data)
    }

    if (MSE == "none") {
      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = data_specs$dName),
        MSE_Estimates = NULL,
        AdjustedSD = NULL
      )

      class(result) <- c("SAEforest_nonLin", "SAEforest")
      return(result)
    }

    # MSE estimation
    if (MSE != "none") {
      if (isTRUE(var.adjust)) {
        # ErrorSD is already Mendez-Lohr bias-corrected inside the MERF loop
        # (var.adjust path); a second correction here would double-count the bias.
        adj_SD <- nonLin_preds[[2]]$ErrorSD
      } else {
        message("Error SD Bootstrap started:")
        adj_SD <- adjust_ErrorSD(Y = Y, X = X, smp_data = smp_data, mod = nonLin_preds[[2]], B = B_adj, ...)
      }
      message(paste("MSE Bootstrap with", B, "rounds started:"))
    }

    if (MSE == "wild") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = TRUE,
        MC = FALSE,
        B_point = B_MC,
        aggregate_to = aggregate_to,
        var.adjust = var.adjust,
        B_adj = B_adj,
        adj_tol = adj_tol,
        adj_tol_mse = adj_tol_mse,
        mse_tol = mse_tol,
        transformation = transformation,
        select.indicator = select.indicator,
        dedup_by = dedup_by,
        K_fixed = K_fixed_mse,
        cores = cores,
        worker.threads = worker.threads,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = data_specs$dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = data_specs$dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("SAEforest_nonLin", "SAEforest")
      return(result)
    }

    if (MSE == "nonparametric") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = FALSE,
        MC = FALSE,
        B_point = B_MC,
        aggregate_to = aggregate_to,
        var.adjust = var.adjust,
        B_adj = B_adj,
        adj_tol = adj_tol,
        adj_tol_mse = adj_tol_mse,
        mse_tol = mse_tol,
        transformation = transformation,
        select.indicator = select.indicator,
        dedup_by = dedup_by,
        K_fixed = K_fixed_mse,
        cores = cores,
        worker.threads = worker.threads,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = data_specs$dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = data_specs$dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("SAEforest_nonLin", "SAEforest")
      return(result)
    }
  }


  # Point and MSE estimates for domain-level indicators and unit-level data (MC version) ----

  if (smearing == FALSE) {
    nonLin_preds <- point_MC_nonLin(
      Y = Y,
      X = X,
      dName = dName,
      threshold = threshold,
      smp_data = smp_data,
      pop_data = pop_data,
      initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance,
      MaxIterations = MaxIterations,
      importance = importance,
      custom_indicator = custom_indicator,
      B_point = B_MC,
      aggregate_to = aggregate_to,
      var.adjust = var.adjust,
      B_adj = B_adj,
      adj_tol = adj_tol,
      select.indicator = select.indicator,
      dedup_by = dedup_by,
      ...
    )

    # Plug-in K for the MSE replicates: the point fit's per-iteration K values,
    # estimated exactly as before, are handed to every replicate fit instead of
    # being re-estimated with B_adj inner forests per iteration. The point fit
    # itself is never plugged.
    K_fixed_mse <- NULL
    if (adjust_mse == "plugin" && MSE != "none") {
      K_fixed_mse <- unlist(nonLin_preds[[2]]$K)
      stopifnot(is.numeric(K_fixed_mse), length(K_fixed_mse) >= 1)
    }

    if(is.null(aggregate_to)){
      data_specs <- sae_specs(dName = dName, cns = pop_data, smp = smp_data)
    } else{
      data_specs <- sae_specs(dName = aggregate_to, cns = pop_data, smp = smp_data)
    }

    if (MSE == "none") {
      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = data_specs$dName),
        MSE_Estimates = NULL,
        AdjustedSD = NULL
      )

      class(result) <- c("MC_MERF_nonLin", "SAEforest")
      return(result)
    }

    # MSE estimation
    if (MSE != "none") {
      if (isTRUE(var.adjust)) {
        # ErrorSD is already Mendez-Lohr bias-corrected inside the MERF loop
        # (var.adjust path); a second correction here would double-count the bias.
        adj_SD <- nonLin_preds[[2]]$ErrorSD
      } else {
        message(paste("Error SD Bootstrap started:"))
        adj_SD <- adjust_ErrorSD(Y = Y, X = X, smp_data = smp_data, mod = nonLin_preds[[2]], B = B_adj, ...)
      }
      message(paste("MSE Bootstrap with", B, "rounds started:"))
    }

    if (MSE == "wild") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = TRUE,
        MC = TRUE,
        B_point = B_MC,
        aggregate_to = aggregate_to,
        var.adjust = var.adjust,
        B_adj = B_adj,
        adj_tol = adj_tol,
        adj_tol_mse = adj_tol_mse,
        mse_tol = mse_tol,
        select.indicator = select.indicator,
        dedup_by = dedup_by,
        K_fixed = K_fixed_mse,
        cores = cores,
        worker.threads = worker.threads,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = data_specs$dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = data_specs$dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("MC_MERF_nonLin", "SAEforest")
      return(result)
    }

    if (MSE == "nonparametric") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = FALSE,
        MC = TRUE,
        B_point = B_MC,
        aggregate_to = aggregate_to,
        var.adjust = var.adjust,
        B_adj = B_adj,
        adj_tol = adj_tol,
        adj_tol_mse = adj_tol_mse,
        mse_tol = mse_tol,
        select.indicator = select.indicator,
        dedup_by = dedup_by,
        K_fixed = K_fixed_mse,
        cores = cores,
        worker.threads = worker.threads,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = data_specs$dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = data_specs$dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("MC_MERF_nonLin", "SAEforest")
      return(result)
    }
  }
}

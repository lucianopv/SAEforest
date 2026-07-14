# Wrapper Functions for Point Estimates
#
# Point estimates for mean and unit-level data---------------------------------------------

point_mean <- function(Y,
                       X,
                       dName,
                       smp_data,
                       pop_data,
                       initialRandomEffects,
                       ErrorTolerance,
                       MaxIterations,
                       importance = "none",
                       aggregate_to,
                       ...) {

  random <- paste0(paste0("(1|", dName), ")")

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    ...
  )

  unit_preds <- predict(unit_model$Forest, pop_data)$predictions +
    predict(unit_model$EffectModel, pop_data, allow.new.levels = TRUE)

  if(is.null(aggregate_to)){
    unit_preds_ID <- cbind(pop_data[dName], unit_preds)
    f0 <- as.formula(paste0("unit_preds ", " ~ ", dName))
    mean_preds <- aggregate(f0,
                            data = unit_preds_ID,
                            FUN = mean
    )
    colnames(mean_preds) <- c(dName, "Mean")
    }else{
      unit_preds_ID <- cbind(pop_data[aggregate_to], unit_preds)
      f0 <- as.formula(paste0("unit_preds ", " ~ ", aggregate_to))
      mean_preds <- aggregate(f0,
                              data = unit_preds_ID,
                              FUN = mean
      )
      colnames(mean_preds) <- c(aggregate_to, "Mean")
  }


  out_ob <- vector(mode = "list", length = 2)

  out_ob[[1]] <- mean_preds
  out_ob[[2]] <- unit_model

  return(out_ob)
}


# Point estimates for nonlinear indicators and unit-level data (smearing)------------------

point_nonLin <- function(Y,
                         X,
                         dName,
                         threshold,
                         smp_data,
                         pop_data,
                         initialRandomEffects,
                         ErrorTolerance,
                         MaxIterations,
                         importance = "none",
                         custom_indicator,
                         aggregate_to,
                         var.adjust = FALSE,
                         B_adj = 100,
                         adj_tol = 0,
                         transformation = c("none", "log"),
                         select.indicator = NULL,
                         cores = 1,
                         n_smear_residuals = NULL,
                         ...) {

  transformation <- match.arg(transformation, c("none", "log"))

  random <- paste0(paste0("(1|", dName), ")")
  if(is.null(aggregate_to)){
    domains <- names(table(pop_data[[dName]]))
    popSize <- as.numeric(table(pop_data[[dName]]))
  } else{
    domains <- names(table(pop_data[[aggregate_to]]))
    popSize <- as.numeric(table(pop_data[[aggregate_to]]))
  }
  thresh <- get_thresh(Y, threshold = threshold)

  if (transformation == "log" && any(Y <= 0, na.rm = TRUE)) {
    stop("transformation='log' requires strictly positive Y.")
  }
  if (transformation == "log") Y <- log(Y)

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    var.adjust = var.adjust,
    B_adj = B_adj,
    adj_tol = adj_tol,
    ...
  )

  unit_preds <- predict(unit_model$Forest, pop_data)$predictions +
    predict(unit_model$EffectModel, pop_data, allow.new.levels = TRUE)


  #  smearing step
  if(is.null(aggregate_to)){
    dName = dName
  } else{
    dName = aggregate_to
  }

  residuals_to_use <- unit_model$OOBresiduals
  if (!is.null(n_smear_residuals)) {
    if (!is.numeric(n_smear_residuals) || length(n_smear_residuals) != 1L ||
        is.na(n_smear_residuals) || n_smear_residuals < 1) {
      stop("n_smear_residuals must be a single positive number or NULL.", call. = FALSE)
    }
    if (n_smear_residuals < length(residuals_to_use)) {
      residuals_to_use <- sample(residuals_to_use, size = n_smear_residuals)
    }
  }

  exp_residuals <- if (transformation == "log") exp(residuals_to_use) else NULL

  group_idx <- split(seq_along(unit_preds), pop_data[[dName]])

  compute_domain <- function(i) {
    if (cores > 1L) pin_blas_threads()
    domain_preds <- unit_preds[group_idx[[domains[i]]]]
    if (transformation == "log") {
      val_i <- c(outer(exp(domain_preds), exp_residuals, "*"))
      val_i[!is.finite(val_i)] <- NA
    } else {
      smear_i <- matrix(rep(residuals_to_use, popSize[i]), nrow = popSize[i],
                         ncol = length(residuals_to_use), byrow = TRUE)
      smear_i <- smear_i + domain_preds
      val_i <- c(smear_i)
    }
    calc_indicat(val_i, threshold = thresh, custom = custom_indicator,
                 select.indicator = select.indicator)
  }

  # Point estimates involve no randomness (unlike the MSE bootstrap), so no RNG
  # scoping is needed here -- unlike with_parallel_rng, cores just picks pbapply's
  # backend. pbapply::pblapply(cl = <integer>) dispatches to parallel::mclapply
  # (fork) on non-Windows, matching the existing MSE-bootstrap cores mechanism's
  # "ignored on Windows" behavior, and to cl = NULL (serial, byte-identical to the
  # original for loop) when cores == 1. Each domain does no ranger/lme4 call of its
  # own (the Forest/EffectModel predict() calls already happened once, above, for
  # the whole population), so there is no thread-oversubscription risk to guard
  # against here -- force_serial_threads is not needed for this loop.
  cl <- if (cores > 1L) cores else NULL
  smear_list <- pbapply::pblapply(seq_along(domains), compute_domain, cl = cl)

  indicators <- as.data.frame(do.call(rbind, smear_list))
  indicators_out <- cbind(domains, indicators)
  names(indicators_out)[1] <- dName

  out_ob <- vector(mode = "list", length = 2)

  out_ob[[1]] <- indicators_out
  out_ob[[2]] <- unit_model

  return(out_ob)
}


# Point estimates for mean and aggregated covariate data-----------------------------------

point_meanAGG <- function(Y,
                          X,
                          dName,
                          smp_data,
                          Xpop_agg,
                          initialRandomEffects,
                          ErrorTolerance,
                          MaxIterations,
                          OOsample_obs,
                          ADDsamp_obs,
                          w_min,
                          wSet = NULL,
                          importance = "none",
                          verbose = TRUE,
                          samp_seed = 0712,
                          ...) {

  random <- paste0(paste0("(1|", dName), ")")
  groupNames <- as.character(unique(smp_data[[dName]]))
  groupNamesCens <- as.character(unique(Xpop_agg[[dName]]))
  OOsamp <- groupNamesCens[!groupNamesCens %in% groupNames]


  # similarity for out-of-sample
  similarXcens <- Xpop_agg[, colnames(Xpop_agg) != dName]
  rownames(similarXcens) <- groupNamesCens
  simXcensMatrix <- as.matrix(dist(similarXcens))
  diag(simXcensMatrix) <- NA

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    ...
  )

  unit_preds <- predict(unit_model$Forest, smp_data)$predictions
  u_ij <- predict(unit_model$EffectModel, smp_data, allow.new.levels = TRUE)

  smp_data$forest_preds <- unit_preds + u_ij
  smp_data$u_ij <- u_ij

  joint_smp_data <- smp_data

  # order vars by simularity
  if (is.null(wSet)) {
    wSet <- names(sort(ranger::importance(unit_model$Forest), decreasing = TRUE))
  }

  wSet <- wSet[wSet %in% names(Xpop_agg)]

  if (length(OOsamp) != 0) {
    sim_groups <- apply(simXcensMatrix[!(groupNamesCens %in% OOsamp), OOsamp], 2, FUN = which.min)
    sim_groups <- groupNamesCens[!(groupNamesCens %in% OOsamp)][sim_groups]

    # sampling and incorporating OOsample data
    smp_oos <- vector(mode = "list", length = length(OOsamp))

    for (i in seq(length(OOsamp))) {
      samp_from <- smp_data[as.character(smp_data[[dName]]) == sim_groups[i], ]
      set.seed(samp_seed)
      return_oos <- dplyr::sample_n(samp_from, OOsample_obs, replace = TRUE)
      return_oos[dName] <- OOsamp[i]
      smp_oos[[i]] <- return_oos
    }

    OOs_smp_data <- do.call(rbind.data.frame, smp_oos)

    # out-of-sample observations
    unit_preds_add <- predict(unit_model$Forest, OOs_smp_data)$predictions
    u_ij <- 0

    OOs_smp_data$forest_preds <- unit_preds_add
    OOs_smp_data$u_ij <- u_ij

    joint_smp_data <- rbind(smp_data, OOs_smp_data)
  }

  # find weights and adjust for failure
  smp_weightsIncluded <- vector(mode = "list", length = length(groupNamesCens))
  smp_weightsNames <- vector(mode = "list", length = length(groupNamesCens))

  for (i in groupNamesCens) {
    pos <- which(i == groupNamesCens)

    X_input_elm <- as.matrix(joint_smp_data[wSet])[joint_smp_data[dName] == i, ]
    mu_input_elm <- as.matrix(Xpop_agg[wSet])[Xpop_agg[dName] == i, ]

    X_input_elm <- X_input_elm[, colSums(X_input_elm != 0) > 0]
    mu_input_elm <- mu_input_elm[colnames(X_input_elm)]

    ELMweight <- elm_wrapper(X_input_elm, mu_input_elm)
    sum_w <- round(sum(ELMweight$prob), digits = 7)

    w_smp_data <- joint_smp_data[joint_smp_data[dName] == i, ]

    if (sum_w == 1) {
      w_smp_data$weights <- ELMweight$prob
      rownames(w_smp_data) <- NULL
      smp_weightsIncluded[[pos]] <- w_smp_data
      smp_weightsNames[[pos]] <- colnames(X_input_elm)
    } else {
      mod_smp_data <- joint_smp_data[joint_smp_data[dName] == i, ]
      rownames(mod_smp_data) <- NULL

      if (ADDsamp_obs != 0) {
        similarXcens <- Xpop_agg[, colnames(Xpop_agg) != dName]
        rownames(similarXcens) <- groupNamesCens
        simXcensMatrix <- as.matrix(dist(similarXcens))
        diag(simXcensMatrix) <- NA

        sim_group <- which.min(simXcensMatrix[, pos])
        samp_add <- joint_smp_data[joint_smp_data[dName] == groupNamesCens[sim_group], ]
        set.seed(samp_seed)
        return_add <- dplyr::sample_n(samp_add, ADDsamp_obs, replace = TRUE)
        return_add[dName] <- i

        mod_smp_data <- rbind(
          joint_smp_data[joint_smp_data[dName] == i, ],
          return_add
        )

        rownames(mod_smp_data) <- NULL
      }

      # recalculation
      X_input_elm <- as.matrix(mod_smp_data[wSet])
      mu_input_elm <- as.matrix(Xpop_agg[wSet])[Xpop_agg[dName] == i, ]

      X_input_elm <- X_input_elm[, colSums(X_input_elm != 0) > 0]
      mu_input_elm <- mu_input_elm[colnames(X_input_elm)]

      ELMweight <- elm_wrapper(X_input_elm, mu_input_elm)
      sum_w <- round(sum(ELMweight$prob), digits = 7)

      if (sum_w == 1) {
        mod_smp_data$weights <- ELMweight$prob
        rownames(mod_smp_data) <- NULL
        smp_weightsIncluded[[pos]] <- mod_smp_data
        smp_weightsNames[[pos]] <- colnames(X_input_elm)
      } else {
        sum_w <- 0

        while ((sum_w != 1) & (dim(X_input_elm)[2] > w_min)) {
          X_input_elm <- X_input_elm[, -dim(X_input_elm)[2]]
          mu_input_elm <- mu_input_elm[-dim(X_input_elm)[2]]

          ELMweight <- elm_wrapper(X_input_elm, mu_input_elm)
          sum_w <- round(sum(ELMweight$prob), digits = 7)
        }
      }

      if (sum_w == 1) {
        mod_smp_data$weights <- ELMweight$prob
        rownames(mod_smp_data) <- NULL
        smp_weightsIncluded[[pos]] <- mod_smp_data
        smp_weightsNames[[pos]] <- colnames(X_input_elm)
      } else {
        if (verbose == TRUE) {
          message(paste("Calculation of weights failed for area:", i))
        }
        rownames(w_smp_data) <- NULL
        w_smp_data$weights <- 1 / dim(w_smp_data)[1]
        smp_weightsIncluded[[pos]] <- w_smp_data
        smp_weightsNames[[pos]] <- NA
      }
    }
  }

  final_smp_data <- do.call(dplyr::bind_rows, smp_weightsIncluded)

  # estimate final means
  final_smp_data$W_mean <- with(final_smp_data, forest_preds * weights)
  f0 <- as.formula(paste("W_mean", paste(dName), sep = " ~ "))
  Mean_preds <- aggregate(f0, data = final_smp_data, FUN = sum)
  colnames(Mean_preds)[2] <- c("Mean")

  # Prepare return object
  return(list(
    Indicators = Mean_preds,
    MERFmodel = unit_model,
    ModifiedSet = final_smp_data,
    ADDsamp_obs = ADDsamp_obs,
    OOsample_obs = OOsample_obs,
    wSet = wSet, w_min = w_min,
    wAreaInfo = smp_weightsNames
  ))
}


# Point estimates for nonlinear indicators and unit-level data (MC version) ---------------

point_MC_nonLin <- function(Y,
                            X,
                            dName,
                            threshold,
                            smp_data,
                            pop_data,
                            initialRandomEffects,
                            B_point,
                            ErrorTolerance,
                            MaxIterations,
                            importance = "none",
                            custom_indicator,
                            aggregate_to,
                            var.adjust = FALSE,
                            B_adj = 100,
                            adj_tol = 0,
                            select.indicator = NULL,
                            ...) {

  domains <- names(table(pop_data[[dName]]))
  random <- paste0(paste0("(1|", dName), ")")

  popSize_N <- data.frame(table(pop_data[[dName]]))
  popSize_n <- data.frame(table(smp_data[[dName]]))
  colnames(popSize_N) <- c(dName, "N_i")
  colnames(popSize_n) <- c(dName, "n_i")
  popSize <- dplyr::left_join(popSize_N, popSize_n, by = dName)
  popSize[, 3][is.na(popSize[, 3])] <- 0

  thresh <- get_thresh(Y, threshold = threshold)

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    var.adjust = var.adjust,
    B_adj = B_adj,
    adj_tol = adj_tol,
    ...
  )

  unit_preds <- predict(unit_model$Forest, pop_data)$predictions +
    predict(unit_model$EffectModel, pop_data, allow.new.levels = TRUE)

  # preparing data for MC step
  ran_obj <- ran_comp(Y = Y, smp_data = smp_data, mod = unit_model, ADJsd = unit_model$ErrorSD, dName = dName)
  ran_effs <- ran_obj$ran_effs
  forest_res <- ran_obj$forest_res
  smp_data <- ran_obj$smp_data

  pred_val <- matrix(unit_preds,
    ncol = B_point,
    nrow = length(unit_preds), byrow = FALSE
  )

  y_hat <- pred_val + sample(forest_res, size = length(pred_val), replace = TRUE)

  gamm_i <- (unit_model$RanEffSD^2) / (unit_model$RanEffSD^2 + (unit_model$ErrorSD^2 / popSize$n_i))

  v_i <- apply(y_hat, 2, function(x) {
    rep(sample(ran_effs, size = length(popSize$N_i), replace = TRUE), popSize$N_i)
  })
  gamm_1 <- rep((1 - gamm_i), popSize$N_i)
  y_star <- y_hat + gamm_1 * v_i

  if(is.null(aggregate_to)){
    indi_agg <- rep(1:length(popSize$N_i), popSize$N_i)
    comb <- function(x) {
      matrix(unlist(x), nrow = length(popSize$N_i), byrow = TRUE)
    }
  } else{
    N_i_agg <- as.numeric(table(pop_data[[aggregate_to]]))
    indi_agg <- rep(1:length(N_i_agg), N_i_agg)
    comb <- function(x) {
      matrix(unlist(x), nrow = length(N_i_agg), byrow = TRUE)
    }
  }

  my_agg <- function(x) {
    tapply(x, indi_agg, calc_indicat, threshold = thresh, custom = custom_indicator,
           select.indicator = select.indicator, simplify = FALSE)
  }
  tau_star <- apply(y_star, my_agg, MARGIN = 2, simplify = FALSE)

  col_names <- colnames(tau_star[[1]]$`1`)
  if (is.null(col_names)) col_names <- select.indicator

  tau_star <- sapply(tau_star, comb, simplify = FALSE)

  indicators <- Reduce("+", tau_star) / length(tau_star)
  colnames(indicators) <- col_names

  # preparing final output
  result <- list(
    if(is.null(aggregate_to)){
    Indicators = cbind(popSize[dName], indicators)
    } else {
      Indicators = cbind(unique(pop_data[aggregate_to]), indicators)
    }
    ,
    model = unit_model
  )

  return(result)
}

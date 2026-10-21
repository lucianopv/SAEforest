#' Main function for unit-level MERF
#'
#' This function enables the use of Mixed Effects Random Forests (MERFs) by effectively
#' combining a random forest from \pkg{ranger} with a model capturing random effects from
#' \pkg{lme4}. The MERF algorithm is an algorithmic procedure reminiscent of an EM-algorithm
#' (see Details). The function is the base-function for the wrapping function (\code{\link{SAEforest_model}}
#' and should not be directly used by the ordinary user. Recommended exceptions are applications exceeding
#' the scope of existing wrapper functions or further research. The function \code{MERFranger}
#' allows to model complex patterns of structural relations (see Examples). The function returns
#' an object of class \code{MERFranger}, which can be used to produce unit-level predictions. In contrast to
#' the wrapping functions, this function does not directly provide SAE estimates on domain-specific indicators.
#'
#' @param Y Continuous input value of target variable.
#' @param X Matrix of predictive covariates.
#' @param random Specification of random effects terms following the syntax of \link[lme4]{lmer}.
#' Random effect terms are specified by vertical bars \code{(|)} separating expressions for design matrices
#' from grouping factors. For further details see \link[lme4]{lmer} and the example below.
#' @param data data.frame of sample data including the specified elements of \code{Y} and
#' \code{X}.
#' @param initialRandomEffects Numeric value or vector of initial estimate of random effects.
#' Defaults to 0.
#' @param ErrorTolerance Numeric value to monitor the MERF algorithm's convergence. Defaults to 1e-04.
#' @param MaxIterations Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm. Defaults to 25.
#' @param importance Variable importance mode processed by the
#' random forest from the \pkg{ranger}. Must be 'none', 'impurity', 'impurity_corrected',
#' 'permutation'. For further details see \link[ranger]{ranger}.
#' @param na.rm Logical. Whether missing values should be removed. Defaults to \code{TRUE}.
#' @param var.adjust Logical. If \code{TRUE}, applies a per-iteration variance bias correction
#' (Krennmair et al., 2026) to the random-effects residuals before refitting, and the returned
#' object additionally carries \code{K}, \code{sigma_e}, \code{sigma_u}, \code{RanEffSD_unadj},
#' \code{ErrorSD_unadj}, \code{sd_Mendez_Lohr_naive_unad} and \code{clamped}. Defaults to \code{FALSE}, which
#' reproduces the original, unadjusted MERF algorithm exactly.
#' @param adjust_method Character. How the bias correction \code{K} is computed when
#' \code{var.adjust = TRUE}. \code{"refit"} (default) follows Krennmair & Schmid (2022),
#' refitting a random forest for every bootstrap draw, and is byte-identical to previous
#' behaviour. \code{"weights"} follows Mendez & Lohr (2011, Procedure 2): the forest's
#' out-of-bag weight matrix is held fixed and no new trees are grown, which is the variant
#' those authors validate and which they report the refitting version can overcorrect.
#' The forest is then grown with \code{keep.inbag = TRUE} automatically. Dense in memory,
#' so \code{O(n^2)}; intended for moderate \code{n}.
#' @param re_method Character. How the random-effects step sees the forest.
#' \code{"insample"} (default) uses the forest's own out-of-bag predictions and is
#' byte-identical to previous behaviour. \code{"crossfit"} refits the forest on
#' other areas and predicts each area from a model that never saw it, which stops a
#' flexible learner from absorbing the area effect through area-constant covariates.
#' Costs \code{re_folds} extra forest fits per iteration.
#' @param re_folds Number of area folds used when \code{re_method = "crossfit"}.
#' Defaults to 5.
#' @param B_adj Numeric value specifying the number of bootstrap replicates used within the
#' variance bias correction when \code{var.adjust = TRUE}. Defaults to 100. Ignored when
#' \code{var.adjust = FALSE}.
#' @param adj_tol Numeric \eqn{\ge 0}. Relative Monte-Carlo standard-error tolerance
#' for the adaptive early-stopping of the residual-variance bias correction when
#' \code{var.adjust = TRUE}. The inner bootstrap fits forests in batches and stops
#' once \code{sd(g)/sqrt(m)/mean(g) < adj_tol} (after at least \code{min(20, B_adj)}
#' forests, capped at \code{B_adj}). \code{adj_tol = 0} (or \code{NULL}) disables early-stopping and
#' always fits the full \code{B_adj} replicates. Defaults to \code{0.05}. Ignored when
#' \code{var.adjust = FALSE}.
#' @param K_fixed Optional finite numeric vector of per-iteration bias terms.
#' When supplied (requires \code{var.adjust = TRUE}), iteration \code{i} uses
#' \code{K_fixed[[min(i, length(K_fixed))]]} instead of estimating K with the
#' inner bootstrap. Used by the MSE bootstrap's plug-in mode
#' (\code{adjust_mse = "plugin"}), which passes the point fit's per-iteration
#' K values. Defaults to \code{NULL} (estimate K, current behaviour).
#' @param seed Integer value used to seed R's random number generator once at the start of the
#' function, enabling reproducibility for direct, standalone calls of \code{MERFranger}. If
#' \code{NULL} (the default), no seed is set. Note that \code{\link{SAEforest_model}} does not
#' pass its own \code{seed} argument down to this function; each MERF fit called within a
#' bootstrap therefore continues to draw its own randomness from the caller's seeded RNG stream.
#' @param ... Additional parameters are directly passed to the random forest \link[ranger]{ranger}.
#' Most important parameters are for instance \code{mtry} (number of variables to possibly split at
#' in each node), or \code{num.trees} (number of trees). For further details on possible parameters
#' see \link[ranger]{ranger} and the example below.
#'
#' @return An object of class MERFranger includes the following elements:
#'
#' \item{\code{Forest}}{A random forest of class \link[ranger]{ranger} modelling fixed effects
#' of the model.}
#' \item{\code{EffectModel}}{A model of random effects of class \code{\link[lme4]{merMod}} capturing
#' structural components of MERFs and modeling random components.}
#' \item{\code{RandomEffects}}{List element containing the values of random intercepts from \code{EffectModel}.}
#' \item{\code{RanEffSD}}{Numeric value of the standard deviation of random intercepts.}
#' \item{\code{ErrorSD}}{Numeric value of standard deviation of unit-level errors.}
#' \item{\code{VarianceCovariance}}{VarCorr matrix from \code{EffectModel}.}
#' \item{\code{LogLik}}{Vector with numerical entries showing the loglikelihood of the MERF algorithm.}
#' \item{\code{IterationsUsed}}{Numeric number of iterations used until convergence of the MERF algorithm.}
#' \item{\code{OOBresiduals}}{Vector of OOB-residuals.}
#' \item{\code{Random}}{Character specifying the random intercept in the random effects model.}
#' \item{\code{ErrorTolerance}}{Numerical value to monitor the MERF algorithm's convergence.}
#' \item{\code{initialRandomEffects}}{Numeric value or vector of initial specification of random effects.}
#' \item{\code{MaxIterations}}{Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm.}
#'
#' If \code{var.adjust = TRUE}, the following additional elements are returned:
#' \item{\code{K}}{List of the per-iteration Mendez and Lohr (2011) bootstrap bias terms that are
#' subtracted from the naive residual variance.}
#' \item{\code{sigma_e}}{List of the per-iteration residual standard deviations of the random effects
#' model fitted on the bias-corrected residuals.}
#' \item{\code{sigma_u}}{List of the per-iteration random intercept standard deviations.}
#' \item{\code{RanEffSD_unadj}}{Numeric unadjusted standard deviation of random intercepts, taken from a
#' reference model fitted on the uncorrected residuals.}
#' \item{\code{ErrorSD_unadj}}{Numeric unadjusted residual standard deviation, taken from a reference model
#' fitted on the uncorrected residuals.}
#' \item{\code{sd_Mendez_Lohr_naive_unad}}{Numeric naive residual standard deviation
#' \code{sd(Target - forest_preds)} at the final iteration, before bias correction.}
#' \item{\code{clamped}}{Integer number of EM iterations in which \code{K} exceeded the naive
#' residual variance, so the adjusted SD was clamped to 0 (each one also raises a warning).
#' 0 when none.}
#'
#' @details
#' There exists a generic function for \code{predict} for objects obtained by \code{MERFranger}.
#'
#' The MERF algorithm iteratively optimizes two separate steps: a) the random forest
#' function, assuming the random effects term to be correct and b) estimates the random
#' effects part, assuming the OOB-predictions from the forest to be correct. Overall convergence
#' of the algorithm is monitored by the log-likelihood of a joint model of both components. For
#' further details see Krennmair & Schmid (2022) or Hajjem et al. (2014).
#'
#' Note that the \code{MERFranger} object is a composition of elements from a random forest of class
#' \code{ranger} and a random effects model of class \code{\link[lme4]{merMod}}. Thus, all generic functions are
#' applicable to corresponding objects. For further details on generic functions see \code{\link[ranger]{ranger}}
#' and \code{\link[lme4]{lmer}} as well as the examples below.
#'
#' @references
#' Hajjem, A., Bellavance, F., & Larocque, D. (2014). Mixed-Effects Random Forest for Clustered
#' Data. Journal of Statistical Computation and Simulation, 84 (6), 1313–1328.
#'
#' Krennmair, P., & Schmid, T. (2022). Flexible Domain Prediction Using Mixed Effects
#' Random Forests. Journal of Royal Statistical Society: Series C (Applied Statistics) (forthcoming).
#'
#' @seealso \code{\link{SAEforest}}, \code{\link[ranger]{ranger}}, \code{\link[lme4]{lmer}},
#' \code{\link{SAEforest_model}}
#'
#' @export
#' @import ranger
#' @import lme4
#'
#' @examples
#' \donttest{
#' # Load Data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]
#'
#' # Example 1:
#' # Calculating general model used in wrapper functions
#'
#' model1 <- MERFranger(Y = income, X = X_covar, random = "(1|district)",
#'                      data = eusilcA_smp, num.trees=50)
#'
#' # get individual predictions:
#'
#' ind_pred <- predict(model1, eusilcA_pop)
#'}
MERFranger <- function(Y, X, random, data,
                       importance = "none",
                       initialRandomEffects = 0,
                       ErrorTolerance = 0.0001,
                       MaxIterations = 25,
                       na.rm = TRUE,
                       var.adjust = FALSE,
                       B_adj = 100,
                       adj_tol = 0.05,
                       adjust_method = c("refit", "weights"),
                       re_method = c("insample", "crossfit"),
                       re_folds = 5L,
                       K_fixed = NULL,
                       seed = NULL,
                       ...) {
  adjust_method <- match.arg(adjust_method)
  re_method <- match.arg(re_method)

  if (na.rm == TRUE) {
    comp_smp <- complete.cases(data)
    data <- data[comp_smp, ]
    Y <- Y[comp_smp]
    X <- X[comp_smp, ]
  }

  input_checks_MERF(
    Y = Y, X = X, data = data, initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
    importance = importance, na.rm = na.rm
  )

  if (isTRUE(var.adjust) &&
      !(is.numeric(B_adj) && length(B_adj) == 1 && !is.na(B_adj) && B_adj >= 1)) {
    stop("When var.adjust = TRUE, B_adj must be a single numeric value >= 1.")
  }

  # adjust_method is only meaningful when the correction runs at all. It is
  # IGNORED rather than rejected when var.adjust = FALSE, matching how B_adj and
  # adj_tol behave -- and necessarily so: point_mean() and the aggregated-data
  # path call MERFranger without var.adjust while forwarding their dots, so a
  # SAEforest_model(var.adjust = TRUE, adjust_method = "weights") call reaches
  # them carrying adjust_method but not var.adjust. Erroring here made the whole
  # MSE path unusable with this route.
  if (!isTRUE(var.adjust)) adjust_method <- "refit"

  # Same reasoning for the plugin path: K_fixed means K is not estimated at all
  # in this fit, so how it WOULD have been estimated is moot. Erroring here would
  # make adjust_mse = "plugin" unusable with this route, which is precisely the
  # configuration the MSE bootstrap needs to stay affordable.
  if (!is.null(K_fixed)) adjust_method <- "refit"

  if (identical(adjust_method, "weights") &&
      "keep.inbag" %in% names(list(...)) && !isTRUE(list(...)$keep.inbag)) {
    stop('adjust_method = "weights" needs keep.inbag = TRUE; it is set ',
         "automatically, so pass either nothing or keep.inbag = TRUE.")
  }

  if (identical(re_method, "crossfit") &&
      !(is.numeric(re_folds) && length(re_folds) == 1 && re_folds >= 2)) {
    stop('re_method = "crossfit" needs re_folds to be a single number >= 2.')
  }

  if (!is.null(K_fixed)) {
    if (!isTRUE(var.adjust)) {
      stop("K_fixed requires var.adjust = TRUE.")
    }
    if (!(is.numeric(K_fixed) && length(K_fixed) >= 1 && all(is.finite(K_fixed)))) {
      stop("K_fixed must be a finite numeric vector of length >= 1.")
    }
  }

  if (!is.null(seed)) set.seed(seed)

  Target <- Y
  ContinueCondition <- TRUE
  iterations <- 0
  AdjustedTarget <- Target - initialRandomEffects
  oldLogLik <- 0

  K_list <- list(); sigma_e <- list(); sigma_u <- list()
  r_ij_new <- NULL; naive_unadj <- NA_real_
  clamped <- 0L   # EM iterations whose adjusted SD was clamped to 0

  while (ContinueCondition) {
    iterations <- iterations + 1
    rf <- if (identical(adjust_method, "weights") &&
              !("keep.inbag" %in% names(list(...)))) {
      # The weights are read off the bootstrap multiplicities, which ranger only
      # records when asked. Added here rather than pushed onto the user.
      ranger::ranger(x = X, y = AdjustedTarget, importance = importance,
                     keep.inbag = TRUE, ...)
    } else {
      ranger::ranger(x = X, y = AdjustedTarget, importance = importance, ...)
    }
    forest_preds <- rf$predictions

    # Cross-fitting by AREA. The random-effect step is the only consumer of
    # forest_preds, and out-of-bag is not enough for it: an observation held out
    # of a tree is still predicted from its neighbours in the same area, so a
    # forest that can identify the area from area-constant covariates reproduces
    # the area effect anyway and the variance step sees nothing left to explain
    # (absorption; simulations/PREDICTION_crossfit.md). Refitting on other areas
    # removes that route. The returned forest is still the full-data one, so only
    # the variance components and the BLUPs change.
    if (identical(re_method, "crossfit")) {
      grp <- grouping_factor_(data, random)
      folds <- assign_area_folds_(grp, re_folds)
      cf <- rep(NA_real_, length(AdjustedTarget))
      for (k in sort(unique(folds))) {
        hold <- folds[as.character(grp)] == k
        if (all(hold) || !any(hold)) next
        rf_k <- ranger::ranger(x = X[!hold, , drop = FALSE],
                               y = AdjustedTarget[!hold],
                               importance = "none", ...)
        cf[hold] <- stats::predict(rf_k, X[hold, , drop = FALSE],
                                   num.threads = dots_threads_(...))$predictions
      }
      # Any observation the folding could not cover keeps its OOB prediction.
      forest_preds <- ifelse(is.finite(cf), cf, forest_preds)
    }

    if (!var.adjust) {
      ## ---- upstream path (KEEP byte-identical for parity) ----
      f0 <- as.formula(paste0("Target ~ -1+", random))
      lmefit <- lme4::lmer(f0, data = data, REML = FALSE, offset = forest_preds)
      newLogLik <- as.numeric(stats::logLik(lmefit))
      ContinueCondition <- (abs(abs(newLogLik - oldLogLik[iterations]) / oldLogLik[iterations]) > ErrorTolerance &
                              iterations < MaxIterations)
      oldLogLik <- c(oldLogLik, newLogLik)
      AllEffects <- predict(lmefit)
      AdjustedTarget <- Target - (AllEffects - forest_preds)
    } else {
      ## ---- variance-adjusted path (Krennmair et al., 2026) ----
      # na.rm throughout: forest_preds are OOB predictions, and ranger returns NA
      # for an observation in-bag in every tree. Without na.rm one such value makes
      # naive_unadj NA, and the comparison below then throws rather than fitting.
      r_ij <- Target - forest_preds
      r_ij <- r_ij - mean(r_ij, na.rm = TRUE)
      naive_unadj <- sd(Target - forest_preds, na.rm = TRUE)

      if (is.null(K_fixed)) {
        res_K <- if (identical(adjust_method, "weights")) {
          # Mendez & Lohr Procedure 2: reuse this forest's weights, grow nothing.
          adjust_ErrorSD_weights_(
            Y = AdjustedTarget, X = X, rf = rf, B = B_adj,
            num.threads = if ("num.threads" %in% names(list(...))) list(...)$num.threads else 1
          )
        } else {
          adjust_ErrorSD_(
            Y = AdjustedTarget, X = X, smp_data = data, rf = rf,
            B = B_adj, adj_tol = adj_tol, ...
          )
        }
        K <- res_K$K
        # K_list deliberately records the RAW estimate (pre !is.finite defence
        # below); sanitize_K_fixed() relies on this ordering when plugging
        # these values into MSE replicates.
        K_list[[iterations]] <- res_K$K
      } else {
        # Plug-in K (adjust_mse = "plugin"): reuse the point fit's
        # per-iteration K instead of re-estimating it with B_adj inner
        # forests. Iterations beyond length(K_fixed) reuse the last value.
        # Consumes no RNG, so plugin replicate streams differ from refit
        # streams by design; everything else on this path is unchanged
        # (own naive_unadj, clamp warning, rescale, REML lmer).
        K <- K_fixed[[min(iterations, length(K_fixed))]]
        K_list[[iterations]] <- K
      }
      # Defence in depth. adjust_ErrorSD_ now drops missing OOB cells, so K should
      # always be finite; if some future path still yields NA/NaN, fall back to no
      # correction for this iteration instead of dying in the comparison below --
      # `if (NA > x)` is exactly the "missing value where TRUE/FALSE needed" that
      # killed a 90-minute GridSAE fit 19 minutes in.
      if (!is.finite(K)) {
        warning("Variance bias correction K is not finite at iteration ", iterations,
                "; skipping the correction for this iteration.")
        K <- 0
      }
      if (K > naive_unadj^2) {
        clamped <- clamped + 1L
        warning("Variance bias correction K exceeds naive residual variance; ",
                "clamping adjusted SD to 0 at iteration ", iterations, ".")
      }
      sd.adjust <- sqrt(max(0, naive_unadj^2 - K))
      r_ij_new <- (r_ij / sd(r_ij, na.rm = TRUE)) * sd.adjust

      lmefit <- lme4::lmer(as.formula(paste0("r_ij_new ~", random)),
                           data = data, REML = TRUE)
      newLogLik <- as.numeric(stats::logLik(lmefit))
      # isTRUE: a clamped iteration hands lmer an all-zero response, for which
      # lme4 returns logLik = +Inf (sigma 0, ranef 0). The clamped iteration's
      # own check is Inf > tol and continues; the next one divides by that Inf,
      # gets NaN, and `while (NA)` would throw. NaN now ends the loop on the
      # fit just made. Elsewhere the comparison is never NA, so unchanged.
      ContinueCondition <- (isTRUE(abs(abs(newLogLik - oldLogLik[iterations]) / oldLogLik[iterations]) > ErrorTolerance) &
                              iterations < MaxIterations)
      oldLogLik <- c(oldLogLik, newLogLik)
      AllEffects <- predict(lmefit) - lme4::fixef(lmefit)
      AdjustedTarget <- Target - AllEffects
      sigma_e[[iterations]] <- stats::sigma(lmefit)
      sigma_u[[iterations]] <- as.data.frame(lme4::VarCorr(lmefit))$sdcor[1]
    }
  }

  data$forest_preds <- NULL
  # OOBresiduals is the pool the smearing point estimator draws from
  # (pointestim_wrapper.R), so it must represent the UNIT-LEVEL error, which has
  # mean zero. On the upstream path that holds automatically: predict(lmefit)
  # carries the forest_preds offset, so Target - predict(lmefit) is already
  # centred. On the variance-adjusted path it does NOT. r_ij_new is centred, but
  # AllEffects is the per-OBSERVATION random effect, and lme4's intercepts sum to
  # zero over DOMAINS, not over observations. With unbalanced domains the
  # observation-weighted BLUP mean is non-zero whenever domain size correlates
  # with the domain effect, and that mean lands in the pool as a pure level
  # error in every indicator computed from it.
  #
  # Measured on GridSAE (localidade anchor, 709 domains, 1-555 households each,
  # 2026-10-02): pool mean 0.00000 absorbed vs -0.10856 cross-fitted on the
  # wealth index, which accounted for 96.8% of a -0.1122 shift in the district
  # means; +0.01342 on log consumption. The bug needs sigma_v > 0 to bite, so it
  # was invisible wherever the forest absorbed the random effect.
  residuals <- if (!var.adjust) {
    Target - predict(lmefit)
  } else {
    res_ij <- r_ij_new - AllEffects
    res_ij - mean(res_ij, na.rm = TRUE)
  }

  result <- list(
    Forest = rf,
    EffectModel = lmefit,
    RandomEffects = lme4::ranef(lmefit),
    RanEffSD = as.data.frame(lme4::VarCorr(lmefit))$sdcor[1],
    ErrorSD = stats::sigma(lmefit),
    VarianceCovariance = lme4::VarCorr(lmefit),
    LogLik = oldLogLik,
    IterationsUsed = iterations,
    OOBresiduals = residuals,
    Random = random,
    ErrorTolerance = ErrorTolerance,
    initialRandomEffects = initialRandomEffects,
    MaxIterations = MaxIterations
  )

  if (var.adjust) {
    # Unadjusted reporting fit: fit ONCE after convergence on the final forest's
    # raw marginal residuals. This was previously refit every EM iteration with only
    # the last retained; lme4::lmer uses no RNG, so hoisting it is byte-identical and
    # saves IterationsUsed-1 lmer fits.
    Input_lmefit_ <- Target - forest_preds
    lmefit_ <- lme4::lmer(as.formula(paste0("Input_lmefit_ ~", random)),
                          data = data, REML = TRUE)
    result$K <- K_list
    result$sigma_e <- sigma_e
    result$sigma_u <- sigma_u
    result$RanEffSD_unadj <- as.data.frame(lme4::VarCorr(lmefit_))$sdcor[1]
    result$ErrorSD_unadj <- stats::sigma(lmefit_)
    result$sd_Mendez_Lohr_naive_unad <- naive_unadj
    result$clamped <- clamped
  }

  class(result) <- "MERFranger"
  return(result)
}

# Generic predict function for MERFranger ------------------------------------------------
#' @export
predict.MERFranger <- function(object, ...) {
  retval <- predict(object$Forest, ...)$predictions +
    predict(object$EffectModel, allow.new.levels = TRUE, ...)

  return(retval)
}

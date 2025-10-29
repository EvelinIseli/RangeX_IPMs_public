# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DATA ANALYSIS SCRIPT: FUNCTIONS FOR MODELLING --------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Author              : Evelin Iseli
### Data used           : none
### Date last modified  : 29.11.2024
### Purpose             : define functions used for modelling all the vital rates for all species plus model checks, function for setting up IPM kernels and do basic IPM analyses

### code adapted from: 
### - Lyu, S. and J. M. Alexander (2022). "Competition contributes to both warm and cool range edges." Nature Communications 13(1): 2502.
### - Nomoto, H. and Alexander, J.M. (2021). "Drivers of local extinction risk in alpine plants under warming climate." Ecology Letters 24: 1157-1166.
### - Ellner et al. (2016). Data-Driven Modelling of Structured Populations : A Practical Guide to the Integral Projection Model. Cham : Springer.



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#rm(list = ls())

### packages etc. --------------------------------------------------------------

# basic packages
library(dplyr); library(tidylog); library(janitor); library(tidyverse) # data manipulation
library(stringr) # working with regex

# task-specific packages (include short description of what it is used for)
library(doBy) # used in making size classes based on quantiles


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MODEL DIAGNOSTICS ------------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### PLOTS: linear models (growth & fecundity) ----------------------------------

## for supplying raw data 
#plot.diag.lm <- function(x, y, data, fc) {
#  # x: column name of size t column (as character)
#  # y: vital rate to test (as character, e.g. "biomass_comb_2023_log")
#  # data: data including size and vital rate (growth/ fecundity) variables
#  # fc: name of species
#  
#  e <- order(data[[x]])
#  data <- data[e, ]
#  
#  # lm model 
#  m.lm = lm(as.formula(paste(y, "~", x)) , data = data)
#  
#  par(mfrow=c(2,2))
#  
#  # panel a: residuals versus fitted
#  # Ellner 2016, p.42: residuals should have constant variance and no trend in mean
#  zhat <- fitted(m.lm)
#  resid <- residuals(m.lm)
#  plot(zhat, resid, xlab = "Fitted values", ylab = "Residuals")
#  mtext(fc, outer = TRUE, cex = 1.5, line = -2)
#  gam.resid <- gam(resid ~ s(zhat), method = "REML") # add smooth line 
#  rhat <- predict(gam.resid, type = "response")
#  points(zhat, rhat, type = "l")
#  
#  # panel b: normal qq-plot
#  # Ellner 2016, p.43: residuals should be Gaussian (< 5% out of the band)
#  sresid <- rstandard(m.lm)
#  qqPlot(sresid, main = "", xlab = "Normal quantiles", ylab = "Stndrdized residual quantiles", 
#         col.lines = "black", lwd = 1)
#  
#  # panel c: scale-location plot (absolute residuals vs. fitted)
#  # Ellner 2016, p.44: standardized residuals should have constant variance and no trend in mean
#  plot(zhat, sqrt(abs(sresid)), xlab = "Fitted values", ylab = "sqrt(|Std Residuals|)")
#  gam.sresid <- gam(sqrt(abs(sresid)) ~ s(zhat), method = "REML")
#  rhat <- predict(gam.sresid, type = "response")
#  points(zhat, rhat, type = "l")
#  
#  # panel d: compare linear to non-linear fit (gam)
#  # Ellner 2016, p. 44: check out potential non-linearity found in panel a
#  gam.grow <- gam(as.formula(paste(y, "~s(", x, ")")), data = data, method = "REML")
#  AIC(gam.grow, m.lm)
#  gam.grow.fitted <- predict(gam.grow, type = "response")
#  matplot(data[[x]], cbind(fitted(m.lm), gam.grow.fitted), type = "l", 
#          lty = c(1, 2), lwd = 2, xlab = "Size t", ylab = "Fitted size t+1")
#  
#  par(mfrow=c(1,1))
#}



### PLOTS: generalized linear models (survival & flowering) --------------------

plot.diag.glm <- function(x, y, data, fc, mod.family) { # no expectation of Gaussian residuals
  # data: data including size and vital rate (survival/ flowering) variables
  # x: column name of size t column (as character)
  # y: vital rate to test (as character, e.g. "survival_2023")
  # fc: name of species
  # mod.family: binomial or poisson
  
  e <- order(data[[x]])
  data <- data[e, ]
  
  # fit models to the reduced data set
  m.glm = glm(as.formula(paste0(y, "~", x)), data = data, family = mod.family)
  m.gam = gam(as.formula(paste0(y, "~s(", x, ")")), data = data, family = mod.family, method = "REML")
  
  par(mfrow=c(1,2))
  
  # Ellner 2016, p. 45: compare model predictions of the glm (solid line) and a GAM (red points) fit and compare them 
  # with survival rates within size classes (open points)
  sizerank <- seq_along(data[[x]])/nrow(data)
  data$sizeclass <- round(9 * sizerank) + 1
  data[[y]] <- as.numeric(as.character(data[[y]])) # make sure the binary variable is treated as numeric (0 and 1), not factors (1 and 2)
  glm.ps <- summaryBy(as.formula(paste0(x, "+", y, "~ sizeclass")), data = data, na.rm = TRUE)
  plot(glm.ps[[paste0(x, ".mean")]], glm.ps[[paste0(y, ".mean")]], xlab = "Size z", ylab = "Event probability", 
       pch = 1, xlim = range(data[[x]]), ylim = range(c(0, glm.ps[[paste0(y, ".mean")]]))) # c(0,1)
  mtext(fc, outer = TRUE, cex = 1.5, line = -2)
  points(data[[x]], fitted(m.glm), type = "l", lty = 1, lwd = 2)
  svals <- seq(min(data[[x]]), max(data[[x]]), length = 30)
  ghat <- predict(m.gam, newdata = setNames(data.frame(svals), x), type = "response")
  points(svals, ghat, type = "p", col = "red", pch = 16)
  
  # Ellner 2016, p. 45: compare the linear model with a nonlinear model using gam
  # y axis: estimated effect of of the predictor on the response
  # edf: complexity of the spine, values close to 1 indicate linear relationship, higher values non-linearity
  edf_value <- summary(m.gam)$s.table[, "edf"]
  plot(m.gam, seWithMean = TRUE, xlab = "Size z", ylab = "Spline(z)")
  mtext(paste0("edf = ", round(edf_value, 2)),
        side = 3, adj = 1)
  
}

# for supplying a model
plot.diag.glm.mod <- function(mod.test, x, y, data, fc, site, mod.family) { # no expectation of Gaussian residuals
  # data: data including size and vital rate (survival/ flowering) variables
  # x: column name of size t column (as character)
  # y: vital rate to test (as character, e.g. "survival_2023")
  # fc: name of species
  # mod.test: model to test
  # site: site of model
  # mod.family: binomial or poisson
  
  cat("Processing species:", fc, "at site:", site, "\n")
  
  # delete NAs in relevant columns
  if (grepl("seed", y)) {
    data <- data %>%
      dplyr::filter(!is.na(.data[[y]]) & flower_status != 0 & !is.na(size_t0_log))
  } else {
    data <- data %>%
      dplyr::filter(!is.na(.data[[y]]) & !is.na(size_t0_log))
  }
  
  # fitted values (population-level for mixed models)  # <<< NEW
  is_mer <- inherits(mod.test, "merMod")
  fitted_values <- if (is_mer) {
    predict(mod.test, type = "response", re.form = NA)
  } else {
    predict(mod.test, type = "response")
  }
  
  # order by fitted to avoid zig-zag
  ord <- order(fitted_values)
  fitted_values <- fitted_values[ord]
  data_ordered  <- data[ord, ]
  
  # simple data checks for GAM fit
  response_distribution <- table(data_ordered[[y]])
  skip_gam <- FALSE
  if (any(response_distribution < 10) && y == "survival") {
    skip_gam <- TRUE
    cat("Skipping GAM fit due to class imbalance in", y,
        "for species:", fc, "at site:", site, "\n")  # <<< NEW (fc/site)
  }
  if (nrow(data_ordered) < 5 && y == "survival") {
    skip_gam <- TRUE
    cat("Skipping GAM fit due to insufficient data for species:",
        fc, "at site:", site, "\n")                  # <<< NEW (fc/site)
  }
  
  # optional: nonparametric smoother for comparison (skip if less observations than the default spline basis size n = 10)
  m.gam <- NULL
  if (!skip_gam) {
    m.gam <- try(
      mgcv::gam(as.formula(paste0(y, " ~ s(", x, ")")),
                data = data_ordered, family = mod.family, method = "REML"),
      silent = TRUE
    )
    if (inherits(m.gam, "try-error")) m.gam <- NULL
  }
  
  par(mfrow = c(1, 2))
  
  # Ellner 2016 p.45: binned empirical rates vs model & GAM
  sizerank <- seq_along(data_ordered[[x]]) / nrow(data_ordered)
  data_ordered$sizeclass <- round(9 * sizerank) + 1
  data_ordered[[y]] <- as.numeric(as.character(data_ordered[[y]]))
  glm.ps <- doBy::summaryBy(
    stats::as.formula(paste0(x, "+", y, "~ sizeclass")),
    data = data_ordered, na.rm = TRUE
  )
  plot(glm.ps[[paste0(x, ".mean")]], glm.ps[[paste0(y, ".mean")]],
       xlab = "Size z", ylab = "Event probability",
       pch = 1, xlim = range(data_ordered[[x]]),
       ylim = range(c(0, glm.ps[[paste0(y, ".mean")]])))
  mtext(fc, outer = TRUE, cex = 1.5, line = -2)
  mtext(site, outer = TRUE, cex = 1.5, line = -4)
  
  points(data_ordered[[x]], fitted_values,
         type = "p", lty = 0.5, lwd = 1, cex = 0.5, col = "blue")
  
  if (!skip_gam && !is.null(m.gam)) {
    svals <- seq(min(data_ordered[[x]]), max(data_ordered[[x]]), length = 30)
    ghat  <- predict(m.gam, newdata = setNames(data.frame(svals), x), type = "response")
    points(svals, ghat, type = "p", col = "red", pch = 16)
  }
  
  # Residuals vs fitted (use response-scale residuals consistently)  # <<< NEW
  resid_values <- data_ordered[[y]] - fitted_values
  plot(fitted_values, resid_values,
       xlab = "Fitted Values", ylab = "Response residuals")
  
  if (!skip_gam && !is.null(m.gam)) {
    uf <- unique(fitted_values)
    if (length(uf) > 2) {
      gam_fit <- mgcv::gam(resid_values ~ s(fitted_values, k = min(10, length(uf) - 1)))
      pred_values <- predict(gam_fit, newdata = list(fitted_values = fitted_values))
      lines(fitted_values[order(fitted_values)], pred_values[order(fitted_values)],
            col = "red", lwd = 1)
    } else {
      cat("Not enough unique values to fit a GAM.\n")
    }
  }


  
}

### LINEARITY: linear models (growth & fecundity) ------------------------------

#linearity.lm <- function(x, y, data, fc) {
#  # data: data including size and vital rate (size at t1/ potentially fecundity) variables
#  # x: column name of size t column (as character)
#  # y: vital rate to test (as character, e.g. "survival_2023")
#  # fc: name of species
#
#  
#  # linear or quadratic
#  m.linear = lm(as.formula(paste(y, "~", x)) , data = data); AIC.linear = AICc(m.linear)
#  m.quadratic =  lm(as.formula(paste(y, "~", x, "+ I(", x,"^2)")), data = data); AIC.quadratic = AICc(m.quadratic)
#  m.gam = gam(as.formula(paste(y, "~ s(", x, ")")), data = data, method = "REML"); AIC.gam = AICc(m.gam)
#  
#  data.frame(species = fc, n = nrow(m.linear$model), AIC.linear = AIC.linear, AIC.quadratic = AIC.quadratic, AIC.gam = AIC.gam, AIC.delta = AIC.linear - AIC.quadratic)
#}



### LINEARITY: generalized linear models (survival & flowering) ----------------

#linearity.glm <- function(x, y, data, fc, mod.family) {
#  # data: data including size and vital rate (survival/ flowering) variables
#  # x: column name of size t column (as character)
#  # y: vital rate to test (as character, e.g. "survival_2023")
#  # fc: name of species
#  
#  # linear or quadratic
#  m.linear = glm(as.formula(paste(y, "~", x)) , data=data, family = mod.family); AIC.linear = AICc(m.linear)
#  m.quadratic =  glm(as.formula(paste(y, "~", x, "+ I(", x,"^2)")), data=data, family = mod.family); AIC.quadratic = AICc(m.quadratic)
#  m.gam = gam(as.formula(paste(y, "~ s(", x, ")")), data = data, family = mod.family, method = "REML"); AIC.gam = AICc(m.gam)
#  
#  data.frame(species = fc, n = nrow(m.linear$model), AIC.linear = AIC.linear, AIC.quadratic = AIC.quadratic, AIC.gam = AIC.gam, AIC.delta = AIC.linear - AIC.quadratic)
#}


### NORMALITY: linear models (growth & fecundity) ------------------------------

# --> not applicable to glms

#normality.lm <- function(x, y, data, fc) {
#  # data: data including size and vital rate (size at t1/ potentially fecundity) variables
#  # x: column name of size t column (as character)
#  # y: vital rate to test (as character, e.g. "biomass_comb_2023_log")
#  # fc: name of species
#  
#  m.linear = lm(as.formula(paste(y, "~", x)) , data = data)
#  sresid = rstandard(m.linear)
#  test = shapiro.test(sresid) # p < 0.05 if residueals are not normal
#  data.frame(species = fc, test.statistic = test$statistic, p.value = test$p.value)
#}


### CONSTANT VARIANCE: linear models (growth & fecundity) ----------------------

# --> not applicable to glms

#constant.var.lm <- function(x, y, data, fc) {
#  # data: data including size and vital rate (size at t1/ potentially fecundity) variables
#  # x: column name of size t column (as character)
#  # y: vital rate to test (as character, e.g. "biomass_comb_2023_log")
#  # fc: name of species
#  
#  m.linear = lm(as.formula(paste(y, "~", x)) , data = data)
#  zhat <- fitted(m.linear)
#  sresid = rstandard(m.linear)
#  test = cor.test(zhat, sqrt(abs(sresid)), method = "k") # positive correlation indicates that variance increases with increasing fitted values (i.e. heteroscedasticity)
#  data.frame(species = fc, test.statistic = test$statistic, p.value = test$p.value) # p < 0.05 indicates heteroscedasticity
#}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MODEL SELECTION --------------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# model selection function to inlcude random effects (both block and plot or only plot)

modelcmp.glm.mixed <- function(data, response, candidate, mod.family) {
  # data: data containing all the model variables
  # response: response variable as a character, e.g. "survival_2023"
  # candidate: candidate models as a string vector, e.g. "biomass_comb_2021_log * site * warm_treat"
  # mod.family: gaussian for linear models
  # weights conditionally included for survival and flowering! --> no
  
  # 1) data preparation
  
  # prepare the data, delete NAs in necessary columns
  if (grepl("seed", response)) {
    data <- data %>% dplyr::filter(!is.na(.data[[response]]) & flower_status != 0)
  } else {
    data <- data %>% dplyr::filter(!is.na(.data[[response]]))
  }
  N_size <- nrow(data)
  
  # ensure that random effects are factors if present
  if ("block_ID"  %in% names(data))  data$block_ID  <- as.factor(data$block_ID)
  if ("plot_ID_original" %in% names(data)) data$plot_ID_original <- as.factor(data$plot_ID_original)
  
  # set the family of the models (depending on the vital rate) --> transform strings into real family objects
  fam <- if (is.character(mod.family)) { # if the family is passed as a string, then...
    switch(mod.family, # ...transform it
           "gaussian" = gaussian(),
           "binomial" = binomial(),
           "poisson"  = poisson(),
           stop("Unsupported mod.family: ", mod.family)
    )
  } else mod.family # if it's already passed as a family object, just keep it
  
  # 2) prepare the model formula of the full model
  
  # create the fixed effect formula
  fixed_form <- stats::as.formula(paste(response, "~", candidate))
  
  # check which random effect terms are even there and possible to use (exist & have >1 level)?
  can_block  <- "block_ID"  %in% names(data)  && dplyr::n_distinct(data$block_ID)  > 1 # i.e. can block be used as random effect at all?
  can_plot   <- "plot_ID_original" %in% names(data) && dplyr::n_distinct(data$plot_ID_original) > 1
  
  # function: construct the model structure (lmer or glmer)
  fit_once <- function(fixed_formula, use_block, use_plot) {
    # build the random-effects part as text, depending on whether use_block and use_plot are are true
    re_rhs <- ""
    if (use_block) re_rhs <- paste0(re_rhs, " + (1|block_ID)") # use block or not?
    if (use_plot)  re_rhs <- paste0(re_rhs,  " + (1|plot_ID_original)")
    
    # glue fixed + random parts into a single formula object
    form <- stats::as.formula(
      paste0(deparse(fixed_formula), re_rhs)
    )
    
    # choose the right function
    #    - Gaussian -> lmer (LMM)
    #    - Binomial/Poisson -> glmer (GLMM)
    if (identical(fam$family, "gaussian")) { # if Gaussian, use lmer from lm4
      #lme4::lmer(form, data = data, REML = FALSE, na.action = na.fail)
      lme4::lmer(form, data = data, REML = FALSE, na.action = na.fail,
                 control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) # use bobyqa ratehr than default optimizer in lmer (same as in glmer)
    } else { # if not, use glmer from lme4
      lme4::glmer(form, data = data, family = fam, na.action = na.fail,
                  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)), nAGQ = 0)
    }
  }
  
  # try both random effects (if feasible) -> if singular drop plot -> if still singular drop to fixed effects
  fit_status  <- "ok"
  fit_message <- ""
  re_used     <- "none"
  fit <- NULL
  
  try({
    # derive what to USE from what we CAN fit (never plot-only)
    use_block0 <- can_block
    use_plot0  <- can_block && can_plot   # only allow plot if block is also used
    
    fit <- fit_once(fixed_form, use_block0, use_plot0)
    
    re_used <- if (use_block0 && use_plot0) {
      "block+plot"
    } else if (use_block0) {
      "block"
    } else {
      "none"
    }
  }, silent = TRUE)
  
  # 3) determine the random structure of the model
  
  # does the mixed effect model fit work at all? if not --> use fixed effects only
  if (is.null(fit)) {
    fe_form <- fixed_form # use the above defined fixed_form as model formula if the fit with random structure is NULL
    fit <- tryCatch(
      stats::glm(fe_form, data = data, family = fam, na.action = na.fail),
      error = function(e) { fit_status <<- "failed"; fit_message <<- paste0("fe_fit_error: ", e$message); NULL } # note that using random structures failed completely
    )
    re_used <- "none"
  }
  
  # does the mixed effect model work, but is singular? --> use less random effects or even none at all 
  if (!is.null(fit) && inherits(fit, "merMod") && lme4::isSingular(fit, tol = 1e-4)) { # if the model is not singular, this block doesn't run at all, if the fit is not NULL, but singular --> run
    if (re_used == "block+plot") { # block + plot were used as random effects, but the model is singular
      # drop plot, keep block
      fit2 <- tryCatch(fit_once(fixed_form, use_block = TRUE, use_plot = FALSE), error = function(e) NULL) # fit again, but this time use_plot = FALSE
      if (!is.null(fit2)) { # if the second fit is not NULL, use second fit and save a few infos
        fit <- fit2
        re_used <- "block"
        fit_status  <- "singular_stepdown"
        fit_message <- paste(c(fit_message, "dropped plot_ID_original"), collapse = " | ")
      }
      # if block-only is still singular, fall back to fixed effects only
      if (inherits(fit, "merMod") && lme4::isSingular(fit, tol = 1e-4)) { # if the new fit (i.e. fit2) is still singular, also run this bit
        fit3 <- tryCatch(stats::glm(fixed_form, data = data, family = fam, na.action = na.fail), error = function(e) NULL) # fit again, but this time both use_plot and use_block = FALSE
        if (!is.null(fit3)) { # if the thrid fit is not NULL, use third fit and save a few infos
          fit <- fit3
          re_used <- "none"
          fit_message <- paste(c(fit_message, "dropped block_ID → FE"), collapse = " | ")
        }
      }
    } else if (re_used == "block") { # if there's only block to start with (because levels are missing or similar) and the model is singular...
      # already block-only and singular → FE
      fit2 <- tryCatch(stats::glm(fixed_form, data = data, family = fam, na.action = na.fail), error = function(e) NULL) # ...fit with only fixed effects
      if (!is.null(fit2)) { # if this second fixed-effects-only fit is not NULL, use it and save a few infos
        fit <- fit2
        re_used <- "none"
        fit_status  <- "singular_stepdown"
        fit_message <- paste(c(fit_message, "block_ID singular → FE"), collapse = " | ")
      }
    }
  }
  
  # if the model still NULL after all these possibilities (because data empty etc.), add placeholders to the output list
  if (is.null(fit)) {
    dredge.df <- data.frame(model = NA_character_, AICc = NA_real_, delta = NA_real_, n = N_size)
    return(list(dredge_output=dredge.df, best_model=NULL, chosen_delta_aicc=NA_real_,
                chosen_aicc=NA_real_, N_no=N_size, re_used="none", fit_status=fit_status, fit_message=fit_message))
  }
  
  # determine fixed effects structure with model selection
  
  # now that the random structure is determined, use dredge to do the fixed effect model selection
  dredge.output <- tryCatch(MuMIn::dredge(fit, rank = "AICc", trace = FALSE), error = function(e) NULL) # if derdge() throws an error, NULL is returned
  if (is.null(dredge.output) || nrow(as.data.frame(dredge.output)) == 0) { # if the dredge.output is NULL, just add placeholders to the output list
    dredge.df <- data.frame(model = NA_character_, AICc = NA_real_, delta = NA_real_, n = N_size)
    return(list(dredge_output=dredge.df, best_model=NULL, chosen_delta_aicc=NA_real_,
                chosen_aicc=NA_real_, N_no=N_size, re_used=re_used, fit_status="dredge_failed", fit_message=fit_message))
  }
  
  dredge.df <- as.data.frame(dredge.output) %>% dplyr::mutate(n = N_size)
  
  # select best model
  top.model    <- MuMIn::get.models(dredge.output, 1)[[1]] # best model
  aicc2.models <- MuMIn::get.models(dredge.output, subset = delta < 2) # other models within 2 AICc
  pars.model <- aicc2.models[[which.min(sapply(aicc2.models, function(m) length(stats::coef(m))))]] # select most parsimonious model within 2 AICc
  
  # get some metrics
  chosen_aicc <- MuMIn::AICc(pars.model)
  best_aicc   <- MuMIn::AICc(top.model)
  delta_aicc  <- chosen_aicc - best_aicc
  
  # save output in a list
  list(
    dredge_output      = dredge.df,
    best_model         = pars.model,
    chosen_delta_aicc  = delta_aicc,
    chosen_aicc        = chosen_aicc,
    N_no               = N_size,
    re_used            = re_used,
    fit_status         = fit_status,
    fit_message        = fit_message
  )
}


# model selection function for binomial models with a matrix response (mixed)  -
modelcmp.glm.mixed_matrix <- function(data, response_matrix, candidate, mod.family) {
  # data: data containing all the model variables
  # response_matrix: a matrix with two columns (successes, failures)
  # candidate: candidate models as a string vector, e.g. "treat_comp * treat_warm"
  # mod.family: binomial for logistic regression (string or family object)
  
  # 1) data preparation
  
  # add success and failure from the provided response_matrix & delete NAs
  data <- data %>%
    dplyr::mutate(
      success = response_matrix[, 1],
      failure = response_matrix[, 2]
    ) %>%
    dplyr::filter(!is.na(success) & !is.na(failure))
  
  # number of rows after filtering
  N_size <- nrow(data)
  
  # ensure that random effects are factors if present
  if ("block_ID" %in% names(data)) data$block_ID <- as.factor(data$block_ID)
  if ("plot_ID_original" %in% names(data)) data$plot_ID_original <- as.factor(data$plot_ID_original)
  
  # set the family of the models (string -> family object if needed)
  fam <- if (is.character(mod.family)) {
    switch(mod.family,
           "gaussian" = gaussian(),
           "binomial" = binomial(),
           "poisson"  = poisson(),
           stop("Unsupported mod.family: ", mod.family))
  } else mod.family
  
  # 2) prepare the model formula of the full model 
  
  # construct fixed-effects formula with cbind(success, failure) on the LHS
  fixed_form <- reformulate(
    termlabels = candidate,
    response = as.call(list(as.name("cbind"), as.name("success"), as.name("failure")))
  ) |> stats::as.formula()
  
  # check which random effect terms are even there and possible to use (exist & >1 level)?
  can_block <- "block_ID" %in% names(data) && dplyr::n_distinct(data$block_ID) > 1
  can_plot  <- "plot_ID_original" %in% names(data) && dplyr::n_distinct(data$plot_ID_original) > 1
  
  # function: construct the model structure (lmer or glmer)
  fit_once <- function(fixed_formula, use_block, use_plot) {
    # build the random-effects part as text, depending on whether use_block and use_plot are true
    re_rhs <- ""
    if (use_block) re_rhs <- paste0(re_rhs, " + (1|block_ID)")
    if (use_plot)  re_rhs <- paste0(re_rhs, " + (1|plot_ID_original)")
    
    # glue fixed + random parts into a single formula object
    form <- stats::as.formula(paste0(deparse(fixed_formula), re_rhs))
    
    # choose the right function
    #    - Gaussian -> lmer (LMM)
    #    - Binomial/Poisson -> glmer (GLMM)
    if (identical(fam$family, "gaussian")) {
      lme4::lmer(form, data = data, REML = FALSE, na.action = na.fail)
    } else {
      lme4::glmer(form, data = data, family = fam, na.action = na.fail,
                  control = lme4::glmerControl(optimizer = "bobyqa"), nAGQ = 0)
    }
  }
  
  # 3) determine the random structure of the model 
  
  # try both random effects (if feasible) -> if singular drop plot -> if still singular drop to fixed effects
  fit_status  <- "ok"
  fit_message <- ""
  re_used     <- "none"
  fit <- NULL
  
  # initial fit per policy (never plot-only)
  try({
    use_block0 <- can_block
    use_plot0  <- can_block && can_plot
    fit <- fit_once(fixed_form, use_block0, use_plot0)
    re_used <- if (use_block0 && use_plot0) "block+plot" else if (use_block0) "block" else "none"
  }, silent = TRUE)
  
  # if initial mixed model couldn't be fit at all -> fall back to fixed effects glm
  if (is.null(fit)) {
    fit <- tryCatch(
      stats::glm(formula = fixed_form, data = data, family = fam, na.action = na.fail),
      error = function(e) { fit_status <<- "failed"; fit_message <<- paste0("fe_fit_error: ", e$message); NULL }
    )
    re_used <- "none"
  }
  
  # if a mixed model fit exists but is singular -> step-down
  if (!is.null(fit) && inherits(fit, "merMod") && lme4::isSingular(fit, tol = 1e-4)) {
    if (re_used == "block+plot") {
      # drop plot, keep block
      fit2 <- tryCatch(fit_once(fixed_form, use_block = TRUE, use_plot = FALSE), error = function(e) NULL)
      if (!is.null(fit2)) {
        fit <- fit2
        re_used <- "block_ID"
        fit_status  <- "singular_stepdown"
        fit_message <- paste(c(fit_message, "dropped plot_ID_original"), collapse = " | ")
      }
      # if block-only is still singular, fall back to fixed effects only
      if (inherits(fit, "merMod") && lme4::isSingular(fit, tol = 1e-4)) {
        fit3 <- tryCatch(stats::glm(formula = fixed_form, data = data, family = fam, na.action = na.fail), error = function(e) NULL)
        if (!is.null(fit3)) {
          fit <- fit3
          re_used <- "none"
          fit_message <- paste(c(fit_message, "dropped block_ID → FE"), collapse = " | ")
        }
      }
    } else if (re_used == "block") {
      # already block-only and singular → FE
      fit2 <- tryCatch(stats::glm(formula = fixed_form, data = data, family = fam, na.action = na.fail), error = function(e) NULL)
      if (!is.null(fit2)) {
        fit <- fit2
        re_used <- "none"
        fit_status  <- "singular_stepdown"
        fit_message <- paste(c(fit_message, "block_ID singular → FE"), collapse = " | ")
      }
    }
  }
  
  # if still NULL (because data empty or errors), bail with placeholders
  if (is.null(fit)) {
    dredge.df <- data.frame(model = NA_character_, AICc = NA_real_, delta = NA_real_, n = N_size)
    return(list(dredge_output = dredge.df, best_model = NULL, chosen_delta_aicc = NA_real_,
                chosen_aicc = NA_real_, N_no = N_size, re_used = "none",
                fit_status = fit_status, fit_message = fit_message))
  }
  
  # determine fixed effects structure with model selection 
  
  # dredge over fixed effects only (random structure is fixed now)
  dredge.output <- tryCatch(MuMIn::dredge(fit, rank = "AICc", trace = FALSE), error = function(e) NULL)
  if (is.null(dredge.output) || nrow(as.data.frame(dredge.output)) == 0) {
    dredge.df <- data.frame(model = NA_character_, AICc = NA_real_, delta = NA_real_, n = N_size)
    return(list(dredge_output = dredge.df, best_model = NULL, chosen_delta_aicc = NA_real_,
                chosen_aicc = NA_real_, N_no = N_size, re_used = re_used,
                fit_status = "dredge_failed", fit_message = fit_message))
  }
  
  dredge.df   <- as.data.frame(dredge.output) %>% dplyr::mutate(n = N_size)
  top.model   <- MuMIn::get.models(dredge.output, 1)[[1]]
  aicc2.models<- MuMIn::get.models(dredge.output, subset = delta < 2)
  pars.model  <- aicc2.models[[which.min(sapply(aicc2.models, function(m) length(stats::coef(m))))]]
  
  chosen_aicc <- MuMIn::AICc(pars.model)
  best_aicc   <- MuMIn::AICc(top.model)
  delta_aicc  <- chosen_aicc - best_aicc
  
  # return results 
  list(
    dredge_output      = dredge.df,
    best_model         = pars.model,
    chosen_delta_aicc  = delta_aicc,
    chosen_aicc        = chosen_aicc,
    N_no               = N_size,
    re_used            = re_used,
    fit_status         = fit_status,
    fit_message        = fit_message
  )
}



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BOOTSTRAPPING ----------------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# function to bootstrap model coefficients (for all models, sizedependent and independent)
bootstrap.npara <- function(model, data, nboot = 5000) {
  # get model family (poisson, binomial, gaussian)
  family_type <- family(model)$family
  
  # make matrix to store bootstrapped results and a vector to store sigma
  boot_results <- matrix(NA, nrow = nboot, ncol = length(coef(model)))
  sigma_results <- numeric(nboot) 
  
  for (b in 1:nboot) { # for each bootstrapping iteration...
    set.seed(449 + b)
    resampled_data <- data[sample(nrow(data), replace = TRUE), ] # sample some data
    
    # fit model to randomly sampled data
    boot_model <- glm(formula(model),
                      data = resampled_data,
                      family = family_type)
    
    # store coefficients and sigma of this iteration
    boot_results[b, ] <- coef(boot_model)
    sigma_results[b] <- sigma(boot_model) 
  }
  
  # combine bootstrapped coefficients and sigma into a matrix with an additional column
  boot_results_sigma <- cbind(boot_results, sigma = sigma_results)
  
  return(boot_results_sigma)
}


# function for parametric bootstrapping (simulating new data based on the fitted model)
# works for glmerMod, lmerMod, glm --> branches out internally (bootMer for mixed models, manual simulate - refit for glm)
bootstrap.para <- function(model, data, nboot = 5000) {
  
  # 1) Mixed models -> use bootMer (simulate + refit)
  if (inherits(model, c("glmerMod", "lmerMod"))) { # if it's a mixed model, proceed here
    # define what to extract from each refit:
    statfun <- function(fit) {
      # fixed effects:
      fe <- lme4::fixef(fit)
      # sigma only for lmer (Gaussian LMM); NA for glmer:
      sig <- if (inherits(fit, "lmerMod")) sigma(fit) else NA_real_
      c(fe, sigma = sig)
    }
    
    # run parametric bootstrap with bootMer
    # (use.u = TRUE --> simulate using fitted random effects; type = "parametric" --> simulate from fitted model)
    set.seed(449)  # make reproducible
    b <- lme4::bootMer(model, FUN = statfun, nsim = nboot,
                       use.u = TRUE, type = "parametric", parallel = "no")
    res <- b$t # extract numeric matrix of size nboot × k, where k is the length of the vector statfun returns (i.e. nboot x extracted coefficients)
    colnames(res) <- c(names(lme4::fixef(model)), "sigma") # name columns
    return(res)
  }
  
  # 2) Simple models (glm) -> simulate + refit manually with glm 
  if (inherits(model, "glm")) {
    family_obj <- stats::family(model)# get family (binomial/poisson/gaussian)
    lhs        <- as.character(stats::formula(model)[[2]])  # get response variables for model formula
    
    # prepare data copy for simulate/refit 
    fun_dat <- data
    if (lhs %in% names(fun_dat)) {
      # ensure numeric for gaussian/binomial 0/1 responses
      fun_dat[[lhs]] <- as.numeric(as.character(fun_dat[[lhs]]))
    } else {
      # if LHS is not a simple column name (e.g., cbind(...) ), stop — current function is for size-dependent GLMs
      stop("bootstrap.para(): non-scalar response (e.g., cbind) detected in GLM. ",
           "Use a matrix-aware bootstrap for size-independent models.")
    }
    
    fun_model <- stats::glm(stats::formula(model), data = fun_dat, family = family_obj) # fit the model
    
    # store the output
    p <- length(stats::coef(fun_model))
    boot_results  <- matrix(NA_real_, nrow = nboot, ncol = p) # matrix to store results
    sigma_results <- rep(NA_real_, nboot)
    
    for (b in seq_len(nboot)) { # for each bootstrapping iteration... (1:nboot)
      set.seed(449 + b)
      
      # simulate new response data based on the model's coefficients
      sim <- stats::simulate(fun_model, nsim = 1)
      fun_dat[[lhs]] <- sim[[1]] # add simulated response
      
      # refit model to simulated data
      new_model <- stats::glm(stats::formula(fun_model), data = fun_dat, family = family_obj)
      
      # store the coefficients
      boot_results[b, ] <- as.vector(stats::coef(new_model))
      sigma_results[b] <- if (identical(family_obj$family, "gaussian")) sigma(new_model) else NA_real_ # sigma only meaningful for gaussian glm, therefore extract only then
    }
    
    # store the results of either mixed or simple models, with and without sigma
    out <- cbind(boot_results, sigma = sigma_results)
    colnames(out) <- c(names(stats::coef(fun_model)), "sigma")
    return(out)
  }
  
  stop("bootstrap.para(): unsupported model class: ", paste(class(model), collapse = "/"))
}



# function for parametric bootstrapping with a matrix response variable (only binomial models)
# works for glmerMod andglm --> branches out internally (bootMer for mixed models, manual simulate - refit for glm)
bootstrap.para_matrix <- function(model, data, nboot = 5000) {
  
  # check whether family really is binomial (only works for binomial models)
  fam <- tryCatch(stats::family(model), error = function(e) NULL)
  if (is.null(fam) || !identical(fam$family, "binomial")) {
    stop("bootstrap.para_matrix(): model must be binomial.")
  }
  
  # 1) Mixed models -> use bootMer (simulate + refit)
  if (inherits(model, "glmerMod")) {
    statfun <- function(fit) {
      # fixed effects only for binomial GLMMs
      lme4::fixef(fit)
    }
    set.seed(449)
    b <- lme4::bootMer(
      model, FUN = statfun, nsim = nboot,
      use.u = TRUE, type = "parametric", parallel = "no"
    )
    res <- b$t
    colnames(res) <- names(lme4::fixef(model))
    return(res)
  }
  
  # 2) Simple models (glm) -> simulate new cbind(success, failure) + refit manually with glm 
  if (inherits(model, "glm")) {
    # pull the two response names from cbind(success, failure) ~ ...
    # all.vars(formula) starts with the response cbind variables for a binomial matrix response
    vars <- all.vars(stats::formula(model))
    if (length(vars) < 2)
      stop("bootstrap.para_matrix(): couldn't recover matrix response names from formula.")
    response_success <- vars[1]
    response_failure <- vars[2]
    
    # make a clean copy of the data and ensure numeric columns
    fun_dat <- data
    fun_dat[[response_success]] <- as.numeric(as.character(fun_dat[[response_success]]))
    fun_dat[[response_failure]] <- as.numeric(as.character(fun_dat[[response_failure]]))
    
    fun_model <- stats::glm(stats::formula(model), data = fun_dat, family = fam)
    
    p <- length(stats::coef(fun_model))
    boot_results <- matrix(NA_real_, nrow = nboot, ncol = p)
    
    for (b in seq_len(nboot)) {
      set.seed(449 + b)
      sim <- stats::simulate(fun_model, nsim = 1)
      # simulate.glm returns a matrix for binomial cbind(), accessible as sim[[1]]
      fun_dat[[response_success]] <- sim[[1]][, 1]
      fun_dat[[response_failure]] <- sim[[1]][, 2]
      
      new_model <- stats::glm(stats::formula(fun_model), data = fun_dat, family = fam)
      boot_results[b, ] <- as.vector(stats::coef(new_model))
    }
    
    colnames(boot_results) <- names(stats::coef(fun_model))
    return(boot_results)
  }
  
  stop("bootstrap.para_matrix(): unsupported model class: ", paste(class(model), collapse = "/"))
}



## function to bootstrap mean and sd for recruit size (non-paramtric)
#bootstrap.npara.mean.sd <- function(x, nboot = 5000) {
#  results <- data.frame(mean = numeric(nboot), sd = numeric(nboot))
#  
#  # randomly sample, then get mean and sd for each random sample
#  for (b in 1:nboot) {
#    set.seed(449 + b)
#    
#    boot_sample <- sample(x, length(x), replace = TRUE)
#    results$mean[b] <- mean(boot_sample, na.rm = TRUE)
#    results$sd[b] <- sd(boot_sample, na.rm = TRUE)
#    results$N_size[b] <- length(boot_sample)
#  }
#  
#  return(results)
#}

# function to bootstrap mean and sd for recruit size (paramtric)
bootstrap.para.mean.sd <- function(x, nboot = 5000) {
  results <- data.frame(mean = numeric(nboot), sd = numeric(nboot), N_size = numeric(nboot))
  
  # calculate mean and standard deviation of the original data
  original_mean <- mean(x, na.rm = TRUE)
  original_sd <- sd(x, na.rm = TRUE)
  
  # generate parametric bootstrapping samples
  for (b in 1:nboot) {
    set.seed(449 + b)
    
    # simulate new data based on the original mean and sd
    boot_sample <- rnorm(length(x), mean = original_mean, sd = original_sd)
    results$mean[b] <- mean(boot_sample, na.rm = TRUE)
    results$sd[b] <- sd(boot_sample, na.rm = TRUE)
    results$N_size[b] <- length(boot_sample) 
  }
  
  return(results)
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CALCULATE COEFFICIENTS -------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# function to calculate the intercept based on the structure of the best model saved in params
calculate.intercepts.parallel <- function(params) {
  cl <- makeCluster(4)  
  clusterEvalQ(cl, library(stringr))  # load stringr on all workers!
  registerDoParallel(cl)
  
  # initialize a list to store results
  results <- foreach(i = 1:nrow(params), .combine = rbind) %dopar% {
    best_model <- params$best_model_predictors[i]
    intercept <- params$Intercept[i]
    
    # initialize intercepts for combinations
    intercepts <- data.frame(
      intercept.bare.ambi = intercept,
      intercept.bare.warm = intercept,
      intercept.vege.ambi = intercept,
      intercept.vege.warm = intercept,
      stringsAsFactors = FALSE
    )
    
    # logic based on the presence of treat_comp and treat_warm
    # case 1: only treat_comp included
    if (str_detect(best_model, "treat_comp") && !str_detect(best_model, "treat_warm")) {
      if (!is.na(params$treat_compvege[i])) {
        intercepts$intercept.vege.ambi <- intercept + params$treat_compvege[i]
        intercepts$intercept.vege.warm <- intercept + params$treat_compvege[i]  # no treat_warm adjustment, because treat_warm not included
      }
    }
    
    # case 2: only treat_warm included
    if (str_detect(best_model, "treat_warm") && !str_detect(best_model, "treat_comp")) {
      if (!is.na(params$treat_warmwarm[i])) {
        intercepts$intercept.bare.warm <- intercept + params$treat_warmwarm[i]
        intercepts$intercept.vege.warm <- intercept + params$treat_warmwarm[i]  # no treat_comp adjustment, because treat_comp is not included
      }
    }
    
    # case 3: both treat_comp and treat_warm are included, but no interaction
    if (str_detect(best_model, "treat_comp") && str_detect(best_model, "treat_warm") &&
        !str_detect(best_model, "treat_comp:treat_warm")) {
      if (!is.na(params$treat_compvege[i])) {
        intercepts$intercept.vege.ambi <- intercept + params$treat_compvege[i] # ambi, so only add vege effect
      }
      if (!is.na(params$treat_warmwarm[i])) {
        intercepts$intercept.bare.warm <- intercept + params$treat_warmwarm[i] # bare, so only add warm effect
        intercepts$intercept.vege.warm <- intercept + params$treat_compvege[i] + params$treat_warmwarm[i] # add both
      }
    }
    
    # case 4: both treat_comp and treat_warm are included, with their interaction
    if (str_detect(best_model, "treat_comp:treat_warm")) {
      if (!is.na(params$treat_compvege[i])) {
        intercepts$intercept.vege.ambi <- intercept + params$treat_compvege[i] # ambi, so only add vege effect
      }
      if (!is.na(params$treat_warmwarm[i])) {
        intercepts$intercept.bare.warm <- intercept + params$treat_warmwarm[i] # bare, so only add warm effect
      }
      if (!is.na(params$`treat_compvege:treat_warmwarm`[i])) {
        intercepts$intercept.vege.warm <- intercept + 
          params$treat_compvege[i] + 
          params$treat_warmwarm[i] + 
          params$`treat_compvege:treat_warmwarm`[i] # add both plus interaction
      }
    }
    
    return(intercepts)
  }
  
  stopCluster(cl)
  results_df <- cbind(params, results)
  return(results_df)
}


# same function, but for slopes (same logic)
calculate.slopes.parallel <- function(params) {
  cl <- makeCluster(4)  # adjust the number of cores as needed
  clusterEvalQ(cl, library(stringr))  # load stringr on all workers!
  registerDoParallel(cl)
  
  # initialize a list to store results
  results <- foreach(i = 1:nrow(params), .combine = rbind) %dopar% {
    best_model <- params$best_model_predictors[i]
    slope <- params$size_t0_log[i]
    
    # initialize slopes for combinations
    slopes <- data.frame(
      slope.bare.ambi = slope,
      slope.bare.warm = slope,
      slope.vege.ambi = slope,
      slope.vege.warm = slope,
      stringsAsFactors = FALSE
    )
    
    # logic based on the presence of treat_comp and treat_warm
    if (best_model != "") {
      # case 1: only treat_comp is included
      if (str_detect(best_model, "size_t0_log:treat_comp") && 
          !str_detect(best_model, "size_t0_log:treat_warm")) {
        if (!is.na(params$`size_t0_log:treat_compvege`[i])) {
          slopes$slope.vege.ambi <- slope + params$`size_t0_log:treat_compvege`[i]
          slopes$slope.vege.warm <- slope + params$`size_t0_log:treat_compvege`[i] # no treat_warm adjustment, because treat_warm not included
        }
      }
      
      # case 2: only treat_warm is included
      if (str_detect(best_model, "size_t0_log:treat_warm") && 
          !str_detect(best_model, "size_t0_log:treat_comp")) {
        if (!is.na(params$`size_t0_log:treat_warmwarm`[i])) {
          slopes$slope.bare.warm <- slope + params$`size_t0_log:treat_warmwarm`[i]
          slopes$slope.vege.warm <- slope + params$`size_t0_log:treat_warmwarm`[i] # no treat_comp adjustment, because treat_comp is not included
        }
      }
      
      # case 3: both treat_comp and treat_warm, but no interaction
      if (str_detect(best_model, "size_t0_log:treat_comp") && 
          str_detect(best_model, "size_t0_log:treat_warm") && 
          !str_detect(best_model, "size_t0_log:treat_comp:treat_warm")) {
        if (!is.na(params$`size_t0_log:treat_compvege`[i])) {
          slopes$slope.vege.ambi <- slope + params$`size_t0_log:treat_compvege`[i] # ambi, so only add vege effect
        }
        if (!is.na(params$`size_t0_log:treat_warmwarm`[i])) {
          slopes$slope.bare.warm <- slope + params$`size_t0_log:treat_warmwarm`[i] # bare, so only add warm effect
          slopes$slope.vege.warm <- slope + params$`size_t0_log:treat_compvege`[i] + params$`size_t0_log:treat_warmwarm`[i] # add both
        }
      }
      
      # case 4: both treat_comp and treat_warm with interaction
      if (str_detect(best_model, "size_t0_log:treat_comp:treat_warm")) {
        if (!is.na(params$`size_t0_log:treat_compvege`[i])) {
          slopes$slope.vege.ambi <- slope + params$`size_t0_log:treat_compvege`[i] # ambi, so only add vege effect
        }
        if (!is.na(params$`size_t0_log:treat_warmwarm`[i])) {
          slopes$slope.bare.warm <- slope + params$`size_t0_log:treat_warmwarm`[i] # bare, so only add warm effect
        }
        if (!is.na(params$`size_t0_log:treat_compvege:treat_warmwarm`[i])) {
          slopes$slope.vege.warm <- slope + 
            params$`size_t0_log:treat_compvege`[i] + 
            params$`size_t0_log:treat_warmwarm`[i] + 
            params$`size_t0_log:treat_compvege:treat_warmwarm`[i] # add both plus interaction
        }
      }
    }
    
    return(slopes)
  }
  
  stopCluster(cl)
  results_df <- cbind(params, results)
  return(results_df)
}



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SET UP VITAL RATE FUNCTIONS --------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. probability of surviving (logistic)
survival.t0 <- function(z, params) {
  if (is.na(params$surv.slope)) { # check whether the mdoel is intercept only
    # if yes, define the linear predictor using only the intercept
    u <- params$surv.int
  } else {
    # if not, define the linear predictor using both intercept and slope
    u <- params$surv.int + params$surv.slope * z
  }

  # inverse logit transformation
  survival.p <- exp(u)/(1 + exp(u))
  survival.p[u>700] <- 1 # exp(710) gives Inf values
  # return the probability of survival
  return(survival.p)
}

# 2.1 growth function (gaussian)
growth.t0t <- function(z1, z, params) {
  # mean size at t based on size at t0
  if (is.na(params$growth.slope)) {
    mu <- params$growth.int  # use only  intercept
  } else {
    mu <- params$growth.int + params$growth.slope * z  # use both intercept and slope
  }

  # sd around mean
  sig <- params$growth.sigma
  
  ## debugging
  #print(paste("z:", z, "mu:", mu, "sig:", sig))
  
  # probability density of new biomass at time t for current biomass at t0
  pd <- dnorm(z1, mean = mu, sd = sig)
  # return probability density
  return(pd)
}


# 2.2 growth function (ceiling)
growth.t0t_ceiling <- function(z1, z, params) {
  # mean size at t based on size at t0
  if (is.na(params$growth.slope)) {
    mu <- params$growth.int  # use only  intercept
  } else {
    mu <- params$growth.int + params$growth.slope * z  # use both intercept and slope
  }
  
  # ceiling
  mu[mu > params$U] = params$U
  print(sum(mu > params$U)) # only mention if ceiling was used
  
  # sd around mean
  sig <- params$growth.sigma
  
  ## debugging
  #print(paste("z:", z, "mu:", mu, "sig:", sig))
  
  # probability density of new biomass at t for current biomass at t0
  pd <- dnorm(z1, mean = mu, sd = sig)
  # return probability density
  return(pd)
}

# 2.3 growth function (gaussian)
# to predict size
growth.t0 <- function(z, params) {
  # mean size at t based on size at t0
  if (is.na(params$growth.slope)) {
    size_t1_pred <- params$growth.int  # use only  intercept
  } else {
    size_t1_pred <- params$growth.int + params$growth.slope * z  # use both intercept and slope
  }
  return(size_t1_pred)
}


# 3. flowering probability (logistic)
flowering.t0 <- function(z, params) {
  if (is.na(params$flow.slope)) { # check whether the model is intercept only
    # if yes, define the linear predictor using only the intercept
    u <- params$flow.int
  } else {
    # if not, define the linear predictor using both intercept and slope
    u <- params$flow.int + params$flow.slope * z
  }
  
  # probability of flowering, inverse logit
  flowering.p <- exp(u)/(1+exp(u))
  flowering.p[u>700] <- 1 # exp(710) gives Inf values
  # return probability of flowering
  return(flowering.p) 
}

# 4.1 seed production (poisson)
seeds.t0 <- function(z, params) {
  # mean seed production at t based on size at t0
  if (is.na(params$seed.slope)) {
    u <- params$seed.int  # use only  intercept
  } else {
    u <- params$seed.int + params$seed.slope * z  # use both intercept and slope
  }
  
  # convert linear predictor into mean seed production
  mu <- exp(u) 
  
  # generate mean seed count based on Poisson distribution (only integers and no negatives)
  seeds <- pmax(0, rpois(length(mu), lambda = mu))  
  # return seed production
  return(seeds)

}

# 4.2 seed production (ceiling)
seeds.t0_ceiling <- function(z, params) {
  # mean seed production at t based on size at t0
  if (is.na(params$seed.slope)) {
    u <- params$seed.int  # use only  intercept
  } else {
    u <- params$seed.int + params$seed.slope * z  # use both intercept and slope
  }
  
  # convert linear predictor into mean seed production
  mu <- exp(u) 
  # generate mean seed count based on Poisson distribution (only integers and no negatives)
  seeds <- pmax(0, rpois(length(mu), lambda = mu))  
  
  # add a ceiling
  seeds[seeds > params$max_fecundity] <- params$max_fecundity
  # return seed production
  return(seeds)
}

# 4.3 mean seed production (poisson)
seeds.t0_mean <- function(z, params) {
  # mean seed production at t based on size at t0
  if (is.na(params$seed.slope)) {
    u <- params$seed.int  # use only  intercept
  } else {
    u <- params$seed.int + params$seed.slope * z  # use both intercept and slope
  }
  
  # convert linear predictor into mean seed production
  mu <- exp(u) 
  
  return(mu)
}

# 5. probability of seed germination (logistic)
germination <- function(params) {
  # define linear predictor (intercept)
  u <- params$germ.int
  
  # inverse logit transformation
  germination.p <- exp(u) / (1 + exp(u))
  
  # return probability
  return(germination.p)
}

# 6. probability of seedling establishment (logistic)
establishment <- function(params) {
  # define linear predictor (intercept)
  u <- params$est.int
  
  # inverse logit transformation
  establishment.p <- exp(u) / (1 + exp(u))
  
  # return probability
  return(establishment.p)
}


# 7. seedling size (constant)
seedling.t <- function(z1, params) {
  # probability density of recruits
  rpd <- dnorm(z1, mean = params$r_mean, sd = params$r_sd)
  # return probability density
  return(rpd)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SET UP IPM KERNELS -----------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# create parameter list our of data frame for each species * site * treatment combination
mk.params.list <- function(data_vrcols, data_sstcombi) {
  
  # initialize empty list to store the parameters 
  params_list <- list()
  
  # loop to create a list of params for each species-site-treatment combination
  for (i in 1:nrow(data_sstcombi)) {
    
    # set current combi
    species <- data_sstcombi$species[i]
    site <- data_sstcombi$site[i]
    treat_combi <- data_sstcombi$treat_combi[i]
    
    # filter for the parameters of the current species * site * treatment combination
    dat_current <- data_vrcols %>%
      filter(species == !!species, site == !!site, treat_combi == !!treat_combi)
    
    # create a list of parameters
    params <- list(
      r_mean = dat_current$r_mean[dat_current$vital_rate == "recruit_size"],
      r_sd = dat_current$r_sd[dat_current$vital_rate == "recruit_size"],
      
      surv.int = dat_current$intercept[dat_current$vital_rate == "survival"],
      surv.slope = dat_current$slope[dat_current$vital_rate == "survival"],
      
      growth.int = dat_current$intercept[dat_current$vital_rate == "size_t1_log"],
      growth.slope = dat_current$slope[dat_current$vital_rate == "size_t1_log"],
      growth.sigma = dat_current$sigma[dat_current$vital_rate == "size_t1_log"],
      
      flow.int = dat_current$intercept[dat_current$vital_rate == "flower_status"],
      flow.slope = dat_current$slope[dat_current$vital_rate == "flower_status"],
      
      seed.int = dat_current$intercept[dat_current$vital_rate == "number_seeds"],
      seed.slope = dat_current$slope[dat_current$vital_rate == "number_seeds"],
      
      germ.int = dat_current$intercept[dat_current$vital_rate == "germ_rate"],
      
      est.int = dat_current$intercept[dat_current$vital_rate == "establishment"],
      
      # L and U values (identical for each vital rate)
      L = dat_current$L[dat_current$vital_rate == "survival"],
      U = dat_current$U[dat_current$vital_rate == "survival"]
    )
    
    # store the params list in the params_list with a unique name for each combination
    params_list[[paste(species, site, treat_combi, sep = "_")]] <- params
  }
  
  # return the complete list of parameters
  return(params_list)
}

# growth kernel: survival-growth
G.t0t <- function(z1, z, params) {
  # combine survival and growth
  surv_growth <- survival.t0(z, params) * growth.t0t(z1, z, params)
  # return
  return(surv_growth)
}

# growth kernel: survival-growth, but with ceiling
G.t0t_ceiling <- function(z1, z, params) {
  # combine survival and growth
  surv_growth <- survival.t0(z, params) * growth.t0t_ceiling(z1, z, params)
  # return
  return(surv_growth)
}

# reproduction kernel: reproduction
R.t0t <- function(z1, z, params) {
  # calculate fecundity kernel
  fecundity <- flowering.t0(z, params) * 
    seeds.t0_mean(z, params) * 
    germination(params) *
    establishment(params) *
    seedling.t(z1, params)
  # return
  return(fecundity)
}

# reproduction kernel: reproduction, but with seed production ceiling
R.t0t_ceiling <- function(z1, z, params) {
  # calculate fecundity kernel
  fecundity <- flowering.t0(z, params) * 
    seeds.t0_ceiling(z, params) * 
    germination(params) *
    establishment(params) *
    seedling.t(z1, params)
  # return
  return(fecundity)
}

# full kernel
full.t0t <- function(z1, z, params) {
  # combine growth and reproduction kernel
  G.t0t(z1, z, params) + R.t0t(z1, z, params)
}


# make the IPM kernel
mk.kernel <- function(L, U, par, n, ceiling_growth = FALSE, ceiling_seeds = FALSE) {
  # n: number of meshpoints
  # par: vital rates
  # L, U: lower and upper limit of size
  # fun: full kernel 
  # mesh points 
  h <- (U - L)/n
  meshpts <- L + ((1:n) - 1/2) * h
  
  # no ceiling
  if(ceiling_growth == TRUE & ceiling_seeds == TRUE) {
    # ceiling
    G <- h * (outer(meshpts, meshpts, G.t0t_ceiling, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t_ceiling, params = par))
  } else if (ceiling_growth == TRUE & ceiling_seeds == FALSE){
    G <- h * (outer(meshpts, meshpts, G.t0t_ceiling, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t, params = par))
  } else if (ceiling_growth == FALSE & ceiling_seeds == TRUE) {
    G <- h * (outer(meshpts, meshpts, G.t0t, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t_ceiling, params = par))
  } else {
    G <- h * (outer(meshpts, meshpts, G.t0t, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t, params = par))
  }
  K <- G+R
  return(list(meshpts = meshpts, G=G, R=R, K = K))
}

## make the IPM kernel, but include correction for size eviction
#mk.kernel.sizeev <- function(L, U, par, n, ceiling_growth = FALSE, ceiling_seeds = FALSE) {
#  # n: number of meshpoints
#  # par: vital rates
#  # L, U: lower and upper limit of size
#  # fun: full kernel 
#  # mesh points 
#  h <- (U - L)/n
#  meshpts <- L + ((1:n) - 1/2) * h
#  
#  # no ceiling
#  if(ceiling_growth == TRUE & ceiling_seeds == TRUE) {
#    # ceiling
#    G <- h * (outer(meshpts, meshpts, G.t0t_ceiling, params = par))
#    R <- h * (outer(meshpts, meshpts, R.t0t_ceiling, params = par))
#  } else if (ceiling_growth == TRUE & ceiling_seeds == FALSE){
#    G <- h * (outer(meshpts, meshpts, G.t0t_ceiling, params = par))
#    R <- h * (outer(meshpts, meshpts, R.t0t, params = par))
#  } else if (ceiling_growth == FALSE & ceiling_seeds == TRUE) {
#    G <- h * (outer(meshpts, meshpts, G.t0t, params = par))
#    R <- h * (outer(meshpts, meshpts, R.t0t_ceiling, params = par))
#  } else {
#    G <- h * (outer(meshpts, meshpts, G.t0t, params = par))
#    R <- h * (outer(meshpts, meshpts, R.t0t, params = par))
#  }
#  
#  # prepare for size eviction correction
#  surv <- survival.t0(meshpts, params = par)
#  P <- G
#  
#  # fix eviction of offspring
#  for (i in 1:(n / 2)) {  # loops through first half of G as this is where the offspring is likely to be found
#    G[1, i] <- G[1, i] + (1 - sum(G[, i]))  # correct the first row for evicted individuals (add back any evicted individuals)
#    P[,i] <- G[,i] * surv[i]
#  }
#  # ...and adults
#  for (i in (n / 2 + 1):n) {  # loops through second half of G as this is where the offspring is likely to be found
#    G[n, i] <- G[n, i] + (1 - sum(G[, i]))  # correct the last row for evicted individuals
#    P[,i] <- G[,i] * surv[i]
#  }
#  
#  K <- P + R
#  
#  return(list(meshpts = meshpts, G=P, R=R, K=K))
#}

# plot the IPM kernel

plot.kernel <- function(x,y,k,... ) {
  # x,y: meshpoints
  # k: IPM kernel
  image(x, y, t(k), ...)
}



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SIZE EVICTION ----------------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## set up a function to create a plot showing survival probability vs. colsums of the survivaö/growth kernel
#plot.size.eviction <- function(species_site_treatment, G_colsums, params, size_summary) {
#  
#  # extract variables of current species * site *treatment combination
#  species_curr <- str_sub(name, 1, 6)
#  site_curr <- str_sub(name, 8, 9)
#  treat_comp_curr <- str_sub(name, 11, 14)
#  treat_warm_curr <- str_sub(name, 16, 19)
#  
#  # get the min and max biomass for the current combination
#  maxmin_biomass <- size_summary %>%
#    filter(species == species_curr & site == site_curr & 
#             treat_comp == treat_comp_curr & treat_warm == treat_warm_curr) 
#  
#  # define size range
#  size_range <- seq(maxmin_biomass$min_biomass, maxmin_biomass$max_biomass, length.out = 500)
#  
#  # calculate survival probabilities
#  survival_probs <- survival.t0(size_range, params[[name]])
#  
#  # make plot
#  plot(size_range, G_colsums, type = "p", col = "red", pch = 19, cex = 1,
#       xlab = "Size", ylab = "Size Distribution and Survival Probability",
#       main = paste("Survival Probability vs Size for", name),
#       ylim = c(0, 1.1))
#  
#  #plot(size_range, survival_probs, type = "l", lwd = 2,
#  #     xlab = "Size", ylab = "Survival Probability",
#  #     main = paste("Survival Probability vs Size for", species_curr))
#  
#  # overlay survival probabilities (horizontal line only if survival probability is size-independent)
#  if (length(unique(survival_probs)) == 1) {
#    # if survival_probs has only one unique value (intercept-only model)...
#    abline(h = unique(survival_probs), col = "blue", lwd = 2) # ...add a horizontal line
#  } else {
#    # otherwise plot the survival probabilities as a line
#    lines(size_range, survival_probs, lwd = 2, col = "black")
#  }
#  
#  # optionally save the plot as a PNG
#  # dev.copy(png, filename = paste0("plot_", species_curr, "_", site_curr, "_", treat_comp_curr, "_", treat_warm_curr, ".png"))
#  # dev.off()
#}
#
#
## function calculate the propability of size eviction for recruits and adults (how many indivudals fall outside the specified size limits?)
#prob.size.eviction <- function(params = NA, L = NA, U = NA, n, growth = FALSE, ...) {
#  # params: vital rates parameters
#  # L, U: the lower and upper bounds implemented in the IPM
#  # n: number of bins implemented in the IPM
#  # ellipsis: can be used for rel.tol, which defines the precision of of approximation of the integral
#  
#  # probability of recruit eviction
#  p_eviction_recruit_lower <- pnorm(L, mean = params$r_mean, sd = params$r_sd)
#  p_eviction_recruit_upper <- pnorm(U, mean = params$r_mean, sd = params$r_sd,lower.tail = FALSE)
#  p_eviction_recruit <- data.frame(p_eviction_recruit_lower = p_eviction_recruit_lower, p_eviction_recruit_upper = p_eviction_recruit_upper)
#  
#  # probability of growth eviction
#  if(growth == TRUE) {
#    # mesh points
#    h <- (U - L)/n
#    meshpts <- L + ((1:n) - 1/2) * h
#    
#    # calculate the probability of eviction for each initial size z
#    p_eviction_growth.full = 1 - sapply(meshpts, function(z) integrate(function(u) growth.t0t(u, z, params), L, U, ...)$value); 
#    #p_eviction_growth_upper = sapply(meshpts,function(z) integrate(function(u) growth.t0t(u,z,params), U, U*10, ...)$value); 
#    #p_eviction_growth_lower = sapply(meshpts,function(z) integrate(function(u) growth.t0t(u,z,params), -L*10, L, ...)$value); 
#    p_eviction_growth_upper = sapply(meshpts,function(z) integrate(function(u) growth.t0t(u,z,params), U, Inf, ...)$value); 
#    p_eviction_growth_lower = sapply(meshpts,function(z) integrate(function(u) growth.t0t(u,z,params), -Inf, L, ...)$value); 
#    p_eviction_growth <- data.frame(z = meshpts,
#                                    p_eviction_growth.full = p_eviction_growth.full,
#                                    p_eviction_growth_lower = p_eviction_growth_lower,
#                                    p_eviction_growth_upper = p_eviction_growth_upper)
#  } else {
#    p_eviction_growth <- data.frame(z = NA,
#                                    p_eviction_growth.full = NA,
#                                    p_eviction_growth_lower = NA,
#                                    p_eviction_growth_upper = NA)
#  }
#  
#  # output
#  list(p_eviction_recruit = p_eviction_recruit,
#       p_eviction_growth = p_eviction_growth)
#}
#
#
#
## function to estimate population growth rate after correcting size eviction by iteration
## (eviction measures for "floor-ceiling" solution to eviction, using iteration to compute dominant eigenvalue/vectors --> computes approximation 
## to the effect of expanding the size range from (minsize, maxsize) by applying a demographic floor at minsize and demographic ceiling at maxsize)
#
## code adapted from: Williams, J.L., Miller, T.E.X. & Ellner, S.P. (2012). Avoiding unintentional eviction from integral projection models. Ecology, 93, 2008-2014.
#
#eviction_delta.lambda_iter=function(growthKernel = NA, kernel = NA, survivalFunction = NA,
#                                    minsize = NA, maxsize = NA, params = NA, n.big.matrix = NA) {
#  
#  # growthKernel: growth kernel function g(new.size, old.size, params)
#  # kernel: complete kernel function K(new.size,old.size,params) 
#  # minsize,maxsize: size limits in the model, aka [L,U] or [xmin,xmax] 
#  # params: parameter vector passed to the kernels (params must be defined and used in the call even if it is not used by kernel or growthKernel)
#  # survivalFunction: size-dependent survival function s(x).
#  # OPTIONAL ARGUMENTS
#  # n.big.matrix: linear dimension of approximating matrix for the kernel 
#  
#  # converting provided fnctions into useable format
#  growthKernel=match.fun(growthKernel); 
#  survivalFunction=match.fun(survivalFunction);
#  kernel=match.fun(kernel); 
#  
#  # calculate size range and meshpoints
#  h = (maxsize-minsize)/n.big.matrix;
#  y = minsize + h*c(1:n.big.matrix)-(h/2);
#  
#  # define meshpts for corrected full krnel (needs two more rows)
#  y_corr <- c(minsize, y, maxsize)
#  
#  # construct kernel, v, w for uncorrected model 
#  Kmat = matrix(0, n.big.matrix, n.big.matrix);
#  for(j in 1:n.big.matrix) {
#    Kmat[,j] = h * kernel(y, y[j], params); 
#  }
#  
#  # calculate uncorrected lambda
#  domeigK = lambda.iter(Kmat); 
#  lambda = domeigK$lambda; w = domeigK$w; 
#  v = lambda.iter(t(Kmat))$w; 
#  
#  # calculate the probability of size eviction (size eviction rate)
#  eps = 1-sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), minsize, maxsize)$value); 
#  eps.U = sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), maxsize, Inf)$value); 
#  eps.L = sapply(y,function(z) integrate(function(u) growthKernel(u,z,params), -Inf, minsize)$value); 
#  
#  sx = survivalFunction(y, params);
#  rho = eps*sx; rho.U = eps.U*sx; rho.L = eps.L*sx; 
#  
#  # construct kernel the corrected model (add big and small classes as 2 columns at the right )
#  Kmat2 = cbind(Kmat, kernel(y, maxsize, params), kernel(y, minsize, params));
#  
#  # add the bottom rows: evictees are sent to the small or large class
#  Kmat2 = rbind(Kmat2, c(h*rho.U, rho.U[n.big.matrix], rho.U[1])); 
#  Kmat2 = rbind(Kmat2, c(h*rho.L, rho.L[n.big.matrix], rho.L[1])); 
#  
#  # calculate corrected lambda
#  lambda2 = lambda.iter(Kmat2)$lambda; 
#  
#  vnew = v+(1/lambda)*(v[n.big.matrix]*rho.U+v[1]*rho.L);  
#  dlambdaU = vnew[n.big.matrix]*sum(rho.U*w)/sum(vnew*w); 
#  dlambdaL = vnew[1]*sum(rho.L*w)/sum(vnew*w); 
#  
#  return(list(evict = eps, evict.U = eps.U, evict.L = eps.L, lambda = lambda,
#              lambda2 = lambda2, dlambda = lambda2-lambda,
#              dlambdaU = dlambdaU, dlambdaL = dlambdaL, fullK_uncorr = Kmat, fullK_corr = Kmat2, meshpts_uncorr = y, meshpts_corr = y_corr))
#}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# POPULATION GROWTH RATES ------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# function to project an IPM kernel forward (output: matrix (if n > 1) with the number of individuals in different size classes over time)
project.ipm <- function(n0, k, nstep){
  # k: value of mk.kernel, with meshpoints and kernel
  # nstep: number of step to project
  # n0: number of individuals in different size classes
  
  # project only one step
  if(nstep==1) { 
    nt <- k %*% n0
    re <- nt
  }
  # project n steps
  else if(nstep >=2) {
    ntmat <- matrix(0, nrow=length(n0), ncol = nstep)
    ntmat[,1] <- n0; # each column is nt in one year
    for(t in 2:nstep) { 
      ntmat[,t] <- k %*% ntmat[,t-1] 
      }
    re <- ntmat
  }
  return(re)
}

## function to calculate how many plants would need to be supplemented to maintain a population of a certain size
#calculate_supplementation <- function(pop_size, lambda) {
#  # pop_size: wanted population size
#  # lambda: population growth rate
#
#  # calculate the number of plants to add
#  S_plants <- pop_size * (1 - lambda)
#  
#  
#  return(S_plants)
#}


# function to get lambda via an iterative method (faster for 1000x1000 or bigger matrices)
lambda.iter = function(k, tol = 1e-8) {
  # k: IPM full kernel
  # tol: tolerance, i.e. track convergence while iteratively getting lambda
  
  if(sum(k == 0) == length(k)) { return(list(lambda = NA, w = NA)) }
  else {
    qmax = 10*tol
    lam = 1
    x = rep(1,nrow(k))   
    k = Matrix(k)
    while(qmax > tol) {
      x1 = k%*%x
      qmax = sum(abs(x1-lam*x))  
      lam = sum(x1)
      # lambda is 0 when population is 0
      if(lam == 0) { x = 0; qmax = 0 }
      else x = x1/lam # one individual in total
    } 
    # lambda: converged population growth rate
    # w: stablized size distribution
    return(list(lambda = lam,w = x/sum(x)))
  }
} 


## function to get lambda via the eigenvalue (faster for 250x250 or smaller matrices)
#lambda.k <- function(k, only.lambda = TRUE) {
#  # k: IPM full kerbel
#  # only.lambda: whether or not calcaulte lamdab only or also stable size distribution
#  
#  # lambda is NA when the kernel is null
#  if(sum(k == 0) == length(k)) { return(list(lambda = NA, w=NA)) }
#  else {
#    eigen.k = eigen(k, only.values = only.lambda)
#    lam = abs(eigen.k$values[1]) 
#    if(only.lambda) return(list(lambda=lam, w=NA))
#    else { 
#      # calculate stable size distribution if indicated in call
#      eigen.vec = eigen.k$vectors[,1]
#      w = abs(eigen.vec)/sum(abs(eigen.vec))
#      return(list(lambda=lam, w=w)) 
#    }
#  }
#}


# function to determine the number of bins (min. number bins so that lambda still converges)
number.bins <- function(params, start = 100, end = 5000, step.n = 100, tol=1e-8) {
  # k: IPM full kernel
  # start/ end: range of mesh points
  # step.n: step size
  # tol: tolerance, i.e. track convergence while iteratively getting lambda
  
  # seq of n to try
  nn <- seq(start, end, by = step.n)
  ll <- rep(NA, length(nn))
  k1 <- mk.kernel(L = as.numeric(params$L), U = as.numeric(params$U), par = params, n = nn[1])
  ll[1] <- lambda.iter(k1$K)$lambda 
  
  for(i in 2:length(nn)) {
    ki = mk.kernel(L = as.numeric(params$L), U = as.numeric(params$U), par=params, n=nn[i])
    ll[i] = lambda.iter(ki$K)$lambda # store the lambda of this iteration
    dif.lam = ll[i] - ll[i-1] # compare it to first lambda 
    # check
    if(dif.lam < tol) { # check whether difference is small enough to assume convergence
      print("IPM converged")
      out <- list(n.converged = nn[i], 
                  data.frame(number.bins = nn, lambdas = ll))
      break
    }
  }
  
  # check after loop
  if(nn[i] == end) {
    print("IPM failed to converge")
    out <- list(n.converged = NA, 
                data.frame(number.bins = nn, lambdas = ll))
  }
  return(out)
} 

## central finite difference sensitivities at the *midpoint* params, using RELATIVE step sizes (no log-steps for positive-only params).
## (adapted from Ellner et al. (2016), LTRE analysis of Carlina, chapter 7.6)
#sensitivity.analysis.finitediff <- function(params, eps = 1e-4, ceil_growth = FALSE, ceil_seeds  = FALSE, n = 500) { # , pos_params = c("sigma","r_sd")
#  # eps: relative step factor (try 1e-4; probe 1e-5 and 1e-3 for robustness)
#  # pos_params: parameters that must stay > 0 (perturbed on log scale)
#  
#  # helper function for relative perturbation
#  rel_step <- function(theta, eps) eps * pmax(1, abs(theta))
#  
#  # set some variables
#  n.params <- length(params)  # number of parameters (e.g., survival, growth, etc.)
#  #lambda_original <- NULL
#  L_current <- params$L
#  U_current <- params$U
#  
#  # calculate the original/baseline lambda with the current set of parameters
#  K0 <- mk.kernel(L_current, U_current, params, n = n,
#                  ceiling_growth = ceil_growth,
#                  ceiling_seeds  = ceil_seeds)$K
#  lambda0 <- lambda.iter(K0)$lambda
#  
#  # initialize storage
#  nm   <- names(params)
#  sens <- rep(NA_real_, length(nm))   # sensitivities 
#  lam0 <- rep(lambda0, length(nm))    # baseline/ original lambda (same for all rows)
#  
#  # loop through each parameter to compute its sensitivity
#  for (i in 1:n.params) {
#    
#    pname <- nm[i] # which parameter to perturb?
#    
#    # skip non-perturbable parameters
#    if (pname %in% c("L","U")) next
#    
#    # extract the parameter value
#    theta <- params[[pname]]
#    
#    # calculate the relative perturbation size for this parameter
#    d <- rel_step(theta, eps)
#    
#    # create perturbed parameter sets
#    par_plus  <- params; par_plus[[pname]]  <- theta + d   # θ + Δθ
#    par_minus <- params; par_minus[[pname]] <- theta - d   # θ - Δθ
#    denom <- 2 * d                                         # central difference denominator
#    
#    # build kernels and compute lambda for perturbed parameters
#    Kp <- mk.kernel(L_current, U_current, par_plus,  n = n,
#                    ceiling_growth = ceil_growth,
#                    ceiling_seeds  = ceil_seeds)$K
#    Km <- mk.kernel(L_current, U_current, par_minus, n = n,
#                    ceiling_growth = ceil_growth,
#                    ceiling_seeds  = ceil_seeds)$K
#    
#    lam_p <- lambda.iter(Kp)$lambda
#    lam_m <- lambda.iter(Km)$lambda
#
#    # compute the sensitivity via the central finite-difference formula (dλ/dθ ≈ (λ_up - λ_down) / (2 * Δθ))
#    sens[i] <- (lam_p - lam_m) / denom
#      
#  }
#  
#  # return results 
#  out <- data.frame(
#    parameter       = nm,
#    sensitivity     = sens,
#    original_lambda = lam0,
#    stringsAsFactors = FALSE)
#  
#  return(out)
#}


# function to bootstrap and get a confidence interval around the lambdas (used on the cluster)
lambda.boot <- function(a, params_list, return_kernel = TRUE) {
  # a = index of the params_list
  da <- params_list[[a]]  # get parameters for the specific combination
  
  # generate the kernel with an appropriate number of mesh points (n)
  n_meshpoints <- 500 
  ipm <- mk.kernel(L = as.numeric(da$L), U = as.numeric(da$U), par = da, n = n_meshpoints)
  #print(ipm$K)
  
  # check for NA values in the kernel (shouldn't have any)
  if (sum(is.na(ipm$K)) > 0) { 
    lambda.a <- list(lambda = "null", kernel = NULL)
    print(paste(a, "_null")) 
  } else {
    # calculate lambda using the lambda.k function
    cal.a <- lambda.k(ipm$K)
    lam.a <- cal.a$lambda
    
    # return lambda and optionally the kernel
    if (return_kernel) {
      return(list(lambda = lam.a, K = ipm$K))
    } else {
      return(lam.a)
    }
  }
}

## function to bootstrap and get a confidence interval around the lambdas as well as find the starting population size with starting population  n = 1
#lambda.startpop.boot <- function(a, params_list, start_pop = TRUE) {
#  # a = index of the params_list
#  
#  da <- params_list[[a]]  # get parameters for the specific combination
#  
#  # generate the kernel with an appropriate number of mesh points (n)
#  n_meshpoints <- 500 
#  ipm <- mk.kernel(L = as.numeric(da$L), U = as.numeric(da$U), par = da, n = n_meshpoints)
#  
#  # check for NA values in the kernel (shouldn't have any)
#  if (sum(is.na(ipm$K)) > 0) { 
#    lambda.a <- list(lambda = "null", kernel = NULL)
#    print(paste(a, "_null")) 
#  } else {
#    # calculate lambda using the lambda.k function
#    cal.a <- lambda.k(ipm$K)
#    lam.a <- cal.a$lambda
#    
#    # calculate starting population size only for lambda > 1
#    if (lam.a > 1) {
#      n0 <- find_starting_population(target_pop_size = 10000, kernel = ipm$K, nstep = 10)
#    } else {
#      n0 <- NA  # No need to calculate for lambda <= 1
#    }
#    
#    # return lambda and required starting population (without kernel)
#    if (start_pop == TRUE) {
#      return(list(lambda = lam.a, n0 = n0))  # don't return the kernel --> needs to much memory
#    } else {
#      return(lam.a)
#    }
#  }
#}
  
  
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STARTING POPULATION SIZE -----------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## define function to find starting population with starting population  n = 1 in smallest size class
#find_starting_population_sc1 <- function(target_pop_size, kernel, nstep, tolerance = 1e-3) {
#	# target_pop_size: wanted population size
#	# kernel: IPM kernel K
#	# nstep: number of steps/ years of wanted projection
#	# tolerance: level of difference between pop.size of an iteration and target_pop_size accepted to stop the iterations
#  
#  # initial guess for n0
#  n0 <- c(rep(1, 1), rep(0, 499)) # CHANGE IF MESHPOINTS NUMBER IS CHANGING
#  
#  # iterate to find the right starting population size
#  while (TRUE) {
#    # project IPM
#    projected_kernel <- project.ipm(n0 = n0, k = kernel, nstep = nstep)
#    
#    # calculate total population size at the target time step
#    total_population <- sum(projected_kernel[, nstep])
#    
#    # check if total population is close to the target
#    if (abs(total_population - target_pop_size) < tolerance) { 
#      return(sum(n0))  # return the total starting population size
#    }
#    
#    # adjust n0 based on the scaling factor (adjustment of n0 based on difference between current and wanted pop. size to speed up convergence - there's downscaling if the population overshoots)
#    scaling_factor <- target_pop_size / total_population
#    n0 <- n0 * scaling_factor
#  }
#}

# define function to calculate the required number of seeds for a target pop. size after nstep years
required_seeds_from_params <- function(par, n_mesh, target_pop_size = 5000, nstep = 10, sc1_index = 1) {
  # par: list of vital rate parameters
  # target_pop_size: wanted population size
  # nstep: number of steps/ years of wanted projection
  # sc1_index: which position in state vector corresponds to smallest size class
  # n_mesh: number of meshpoints
  
  K_list <- mk.kernel(L = par$L, U = par$U, par = par, n = n_mesh)
  out_sc1 <- required_sc1_simple(target_pop_size = target_pop_size, kernel = K_list$K,
                                 nstep = nstep, sc1_index = sc1_index)
  
  p_seed_to_sc1 <- plogis(par$germ.int) * plogis(par$est.int) # seed->sc1 probability implied by the same params (logit intercepts)
  
  tibble(
    gT            = out_sc1$gT,
    required_sc1  = out_sc1$required_sc1,
    required_seeds = out_sc1$required_sc1 / p_seed_to_sc1,
    p_seed_to_sc1 = p_seed_to_sc1
  )
}

# define function to calculate the number of individuals in the smallest size class required to reach a target pop. size after nstep years
required_sc1_simple <- function(target_pop_size, kernel, nstep, sc1_index = 1) { # target_pop_size: wanted population size 
  # kernel: IPM kernel K 
  # nstep: number of steps/ years of wanted projection 
  # sc1_index: which position in state vector corresponds to smallest size class 
  
  m <- nrow(kernel)  # infer number of meshpoints 
  n0 <- c(1, rep(0, m-1))  # 1 individual in smallest class 
  
  proj <- project.ipm(n0 = n0, k = kernel, nstep = nstep) 
  gT <- sum(proj[, nstep]) # per-1-sc1 yield at year nstep 
  tibble(gT = gT, required_sc1 = target_pop_size / gT) 
  }



# define function to bump vital rate parameters by a small amount to calculate elasticities (how and how much they get changed depends on their scale etc.)
bump_param <- function(par, par_name, step = 0.01) {
  # par: vital rate parameter to bump
  # par_name: name of the vital rate parameter
  # step: how much to bump
  
  par2 <- par # original parameter value
  
  # define bumper function for changes on logit scale
  .bump_prob <- function(p, step, eps = 1e-6) {
    p2 <- p * (1 + step) # relative bump on the PROBABILITY scale
    p2 <- max(min(p2, 1 - eps), eps)  # clamp inside (eps, 1-eps) to avoid 0/1
    p2
  }
  
  # bump depending on scale
  if (par_name %in% c("surv.int","flow.int","germ.int","est.int")) {
    # logit intercepts: bump probability by +step, back-transform to intercept
    p  <- plogis(par[[par_name]])
    p2 <- .bump_prob(p, step)
    par2[[par_name]] <- qlogis(p2)
  } else if (par_name %in% c("growth.sigma","r_sd")) {
    # positive-only scales
    par2[[par_name]] <- par[[par_name]] * (1 + step)
  } else {
    # default: multiplicative bump (works for positive or negative)
    par2[[par_name]] <- par[[par_name]] * (1 + step)
  }
  par2
}


## define function to calculet elasticities for the required seed numbers
#elasticity_required_seeds <- function(par, par_name, n_mesh,
#                                      target_pop_size = 5000, nstep = 10,
#                                      sc1_index = 1, step = 0.01) {
#  # par: vital rate parameter to bump
#  # par_name: name of the parameter
#  # target_pop_size: wanted population size
#  # nstep: number of steps/ years of wanted projection
#  # sc1_index: which position in state vector corresponds to smallest size class
#  # n_mesh: number of meshpoints
#  # step: how much to bump
#  
#  # calculate "basis" no. required seeds
#  base <- required_seeds_from_params(par, n_mesh, target_pop_size, nstep, sc1_index)
#  
#  # bump one parameter and recalculate the number of required seeds
#  par2 <- bump_param(par, par_name, step)
#  alt  <- required_seeds_from_params(par2, n_mesh, target_pop_size, nstep, sc1_index)
#  
#  # guard against numerics
#  if (!is.finite(base$required_seeds) || !is.finite(alt$required_seeds)) {
#    return(tibble(par_name = par_name, elasticity = NA_real_))
#  }
#  
#  # calculate elasticity
#  E <- (log(alt$required_seeds) - log(base$required_seeds)) / log(1 + step)
#  tibble(par_name = par_name, elasticity = as.numeric(E)) # save in tibble
#}


## define function to find starting population with starting population n = 1 in 10 smallest size classes
#find_starting_population_sc10 <- function(target_pop_size, kernel, nstep, tolerance = 1e-3) {
#  # target_pop_size: wanted population size (could be minimal viable population MVP)
#  # kernel: IPM kernel K
#  # nstep: number of steps/ years of wanted projection
#  # tolerance: level of difference between pop.size of an iteration and target_pop_size accepted to stop the iterations
#  
#  # initial guess for n0
#  n0 <- c(rep(0.1, 10), rep(0, 490)) # CHANGE IF MESHPOINTS NUMBER IS CHANGING
#  
#  # iterate to find the right starting population size
#  while (TRUE) {
#    # project IPM
#    projected_kernel <- project.ipm(n0 = n0, k = kernel, nstep = nstep)
#    
#    # calculate total population size at the target time step
#    total_population <- sum(projected_kernel[, nstep])
#    
#    # check if total population is close to the target
#    if (abs(total_population - target_pop_size) < tolerance) { 
#      return(sum(n0))  # return the total starting population size
#    }
#    
#    # adjust n0 based on the scaling factor (adjustment of n0 based on difference between current and wanted pop. size to speed up convergence - there's downscaling if the population overshoots)
#    scaling_factor <- target_pop_size / total_population
#    n0 <- n0 * scaling_factor
#  }
#}


## function to find required n0 for a specific single size class
#find_starting_population_single_class <- function(target_pop_size, kernel, size_class, nstep, tolerance = 1e-3) {
#  # target_pop_size: wanted population size
#	# kernel: IPM kernel K
#	# nstep: number of steps/ years of wanted projection
#	# tolerance: level of difference between pop.size of an iteration and target_pop_size accepted to stop the iterations
#	# size_class: size class the 1 starting individual should be in
#	
#  # create starting population vector with 1 individual in the specified size class
#  n0 <- rep(0, 3000)
#  n0[size_class] <- 1
#  
#  # iterate to find the right scaling factor for this size class
#  while (TRUE) {
#    # project IPM
#    projected_kernel <- project.ipm(n0 = n0, k = kernel, nstep = nstep)
#    
#    # calculate total population size at the target time step
#    total_population <- sum(projected_kernel[, nstep])
#    
#    # check if total population is close to the target
#    if (abs(total_population - target_pop_size) < tolerance) {
#      return(sum(n0))  # return the total starting population size
#    }
#    
#    # adjust n0 based on the scaling factor
#    scaling_factor <- target_pop_size / total_population
#    n0 <- n0 * scaling_factor
#  }
#}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ANALYSIS FUNCTIONS -----------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# function to calculate bias corrected bootstrap CIs
bcpi <- function(real_lambda, t, alpha) {
  # real_lambda: "real" lambda of this specific site * species * treatment combination
  # t: vector of bootstrap replicates of lambda
  # alpha: significance level
  B <- length(t)
  z0 <- qnorm(mean(t < real_lambda))
  a1 <- pnorm(2 * z0 + qnorm(alpha/2))
  a2 <- pnorm(2 * z0 + qnorm(1 - alpha/2))
  c1 <- quantile(t, a1)
  c2 <- quantile(t, a2)
  return(as.numeric(c(c1, c2)))}


# hierarchical bootstrapping to get competition effect 
hier_bootstrap_lrr <- function(d_pairs, B = 10000, seed = 2025) {
  # d_pairs: data frame containing the paired bootstrap estimates
  # B: number of bootstrap iterations (default = 10000)
  # seed: random seed for reproducibility
  
  # set random seed so results are reproducible
  set.seed(seed)
  
  # prepare a list of bootstrap draws for each species × warming treatment (get all bootstraps into a vector column)
  draws_tbl <- d_pairs %>%
    group_by(species, site_warm) %>%
    summarise(draws = list(lrr_lambda), .groups = "drop")
  
  # list of unique species and how many there are
  species_vec <- draws_tbl %>% pull(species) %>% unique() %>% sort()
  nS <- length(species_vec)
  
  # index each species to the rows in draws_tbl
  idx_by_species <- split(seq_len(nrow(draws_tbl)), draws_tbl$species) # save row numbers taht each species has in draws_tbl
  
  # one iteration of the hierarchical bootstrap:
  # - resample species without replacement (the species composition is not random, so no need to generate inter-species variation!)
  # - for each resampled species, randomly select one lrr draw per warming treatment
  # - compute the mean lrr across species for each warming treatment
  one_iter <- function() {
    
    # resample species (with replacement)
    samp_species <- sample(species_vec, size = nS, replace = FALSE) # this will actually always be a vector of the ten species, but if for some reason species have to be sampleled with replacement, this will be relevant
    
    # for each resampled species, sample one lrr draw per warming treatment
    sampled <- map_dfr(samp_species, function(sp) {
      idxs <- idx_by_species[[sp]]
      map_dfr(idxs, function(i) {
        v <- sample(draws_tbl$draws[[i]], 1, replace = TRUE)
        tibble(site_warm = draws_tbl$site_warm[[i]], value = v)
      })
    })
    
    # compute mean lrr across all sampled species for each warming treatment
    sampled %>%
      group_by(site_warm) %>%
      summarise(mean_lrr = mean(value), .groups = "drop")
  }
  
  # repeat the hierarchical bootstrap B times
  boot_results <- map_dfr(seq_len(B), ~ one_iter() %>% mutate(iter = .x))
  
  return(boot_results)
}


SL.mgcv = function(Y, X, newX, family,
                    obsWeights = NULL, id = NULL,
                    degree = 3,#2,#-1,  # xtan
                    continuous_values = 5,
                    smooth_fn = "s",
                    exclude_spline = NULL,
                    method = "REML",
                    select = TRUE,  # xtan
                    gamma = 1,
                    formula_rhs = NULL,
                    ...) {
  
  #SuperLearner:::.SL.require("mgcv")
  
  reload_gam = FALSE
  
  # Check for GAM package's namespace and unload to ensure proper formula handling.
  # This is not ideal because it adds extra unnecessary computation, and may
  # not be 100% reversible.
  if (F && isNamespaceLoaded("gam") && "package:gam" %in% search()) {
    #gam_ns = getNamespace("gam")
    reload_gam = TRUE
    # Unload gam namespace to avoid any function name conflicts.
    unloadNamespace("gam")
  }
  
  if (is.null(formula_rhs)) {
    
    # Identify covariates that are considered continuous based on # of unique values.
    possible_vars = setdiff(colnames(X), exclude_spline)
    spline_vars = names(which(apply(X[, possible_vars, drop = FALSE], 2,
                                    function(x) {
                                      length(unique(x)) >= continuous_values
                                    })))
    
    standard_vars = setdiff(colnames(X), spline_vars)
    
    formula_spline = NULL
    if (length(spline_vars) > 0) {
      # Formulas don't allow specifying the package of a function :/
      formula_spline = paste(paste0(#"mgcv::",
        smooth_fn,
        "(", spline_vars,
        ", k = ", degree,
        ")"),
        # ", bs='tp' )"),  # xtan
        collapse = " + ")
    }
    
    formula_std = NULL
    # Check if we have any non-continuous vars.
    if (length(standard_vars) > 0) {
      formula_std = paste(standard_vars, collapse = " + ")
    }
    
    formula_rhs = paste(c(formula_spline, formula_std), collapse = " + ")
  }
  
  # Combine with a + between formulas, which is skipped if one side is NULL.
  # Specify outcome name to be _Y so that a covariate can be named "Y".
  # NOTE: mgcv doesn't seem to support `_Y` in its formula interface :/
  formula_str = paste("Y ~ ", formula_rhs)
  
  exports = c("s", "te", "ti", "t2", "gam.control", "gam.fit")
  
  # Copy these functions into this environment.
  # This does seem to resolve potential conflicts with package:gam.
  importIntoEnv(environment(), exports, asNamespace("mgcv"), exports)
  
  formula = as.formula(formula_str)
  
  # NOTE: this can generate an error if there are too few observations per covariate.
  fit = mgcv::gam(formula,
                  #data = cbind(`_Y` = Y, X),
                  data = X,
                  family = family,
                  weights = obsWeights,
                  method = method,
                  select = select,
                  gamma = gamma)
  #})
  
  pred = mgcv::predict.gam(fit, newdata = newX, type = "response")
  fit = list(object = fit)
  out = list(pred = pred, fit = fit)
  class(out$fit) = c("SL.mgcv")
  
  #if (reload_gam) {
  #  loadNamespace("gam")
  #}
  
  return(out)
}


predict.SL.mgcv = function(object, newdata, ...) {
  # SuperLearner:::.SL.require('mgcv')
  pred = mgcv::predict.gam(object = object$object, newdata = newdata,
                           type = "response")
  return(pred)
}

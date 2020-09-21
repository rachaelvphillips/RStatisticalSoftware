

#' @export
make_super_learner <- function(candidates, tmle_task, likelihood, loss_generator = make_eff_loss) {
  loss_function <- loss_generator(tmle_task, likelihood)
  Lrnr_sl$new(candidates,  metalearner = Lrnr_cv_selector$new(loss_function = loss_function))
}
#' @export
make_RR_task <- function(tmle_task, likelihood, for_cv = T) {
  revere_task <- make_revere(tmle_task, likelihood, "gen")
  if(!for_cv){
    revere_task <- revere_task$revere_fold_task("full")
  }
  return(revere_task)
}


#' @export
make_npsem <- function(baseline_covariates, treatment_variable, outcome_variable) {
  npsem <- list(define_node("W", baseline_covariates),
                define_node("A", trtment_variable, c("W")),
                define_node("R", outcome_variable, c("A", "W")),
                define_node("RR", outcome_variable, c("W")))

  return(npsem)
}
#' @export
make_tmle3_task <- function(data, npsem, ...) {
  tmle3_Task(data, npsem, long_format = F, ...)
}
#' @export
make_likelihood <- function(tmle_task, treatment_learner, outcome_learner = NULL, meta_outcome_learner = "default", undersmooth_type = "fourier", ...) {

  factor_list <- list(LF_emp$new("W"),
                      LF_fit$new("A", treatment_learner)  )
  likelihood <- Likelihood$new(factor_list)
  likelihood <- likelihood$train(tmle_task)

  if(undersmooth_type == "fourier"){
    lrnr_derived <- make_learner(Lrnr_undersmooth, basis_generator = fourier_basis(...), offset_learner = outcome_learner,
                                 stratify_by = "A",  ...)

  } else if (undersmooth_type == "hal") {
    lrnr_derived <- make_learner(Lrnr_undersmooth, basis_generator = hal_basis(...),   ...)
  } else {
    message("No undersmoothing will be done.")
    lrnr_derived <- outcome_learner
  }
  lrnr_derived <- make_learner(Lrnr_sl, lrnr_derived, metalearner = Lrnr_cv_selector$new())
  if(undersmooth_type %in% c("fourier", "hal")) {
    lf <- LF_derived$new("R", lrnr_derived, likelihood$clone(), derived_generator())
  } else {
    lf <- LF_fit$new("R", lrnr_derived)
  }
  likelihood$add_factors(lf)
  return(likelihood)
}





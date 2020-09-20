
make_super_learner <- function(candidates, training_task, likelihood, loss_generator = make_eff_loss) {
  loss_function <- loss_generator(training_task, likelihood)
  Lrnr_sl$new(candidates,  metalearner = Lrnr_cv_selector$new(loss_function = loss_function))
}



#' @export
cv_risk <- function(learner, loss_fun, coefs = NULL) {
  # warning(paste("cv_risks are for demonstration purposes only.",
  # "Don't trust these for now."))
  if (!("cv" %in% learner$properties)) {
    stop("learner is not cv-aware")
  }

  task <- learner$training_task
  preds <- learner$predict_fold(task, "validation")
  task <- task$revere_fold_task("validation")
  if (!is.data.table(preds)) {
    preds <- data.table(preds)
    setnames(preds, names(preds), learner$name)
  }
  losses <- preds[, lapply(.SD, loss_fun, task$Y)]
  # multiply each loss (L(O_i)) by the weights (w_i):
  losses_by_id <- losses[, lapply(.SD, function(loss) {
    task$weights *
      loss
  })]
  # for clustered data, this will first evaluate the mean weighted loss
  # within each cluster (subject) before evaluating SD
  losses_by_id <- losses_by_id[, lapply(.SD, function(loss) {
    mean(loss, na.rm = TRUE)
  }), by = task$id]
  losses_by_id[, "task" := NULL]

  # n_obs for clustered data (person-time observations), should be equal to
  # number of independent subjects
  n_obs <- nrow(losses_by_id)
  # evaluate risk SE for each learner incorporating: a) weights and b) using
  # the number of independent subjects
  se <- unlist((1 / sqrt(n_obs)) * losses_by_id[, lapply(
    .SD, sd,
    na.rm = TRUE
  )])

  # get fold specific risks
  validation_means <- function(fold, losses, weight) {
    risks <- lapply(
      origami::validation(losses), weighted.mean,
      origami::validation(weight)
    )
    return(as.data.frame(risks))
  }

  fold_risks <- lapply(
    task$folds,
    validation_means,
    losses,
    task$weights
  )

  fold_risks <- rbindlist(fold_risks)
  fold_mean_risk <- apply(fold_risks, 2, mean)
  fold_min_risk <- apply(fold_risks, 2, min)
  fold_max_risk <- apply(fold_risks, 2, max)
  fold_SD <- apply(fold_risks, 2, sd)

  learner_names <- names(preds)

  risk_dt <- data.table::data.table(
    learner = learner_names,
    coefficients = NA * 0.0,
    mean_risk = fold_mean_risk,
    SE_risk = se,
    fold_SD = fold_SD,
    fold_min_risk = fold_min_risk,
    fold_max_risk = fold_max_risk
  )
  if (!is.null(coefs)) {
    set(risk_dt, , "coefficients", coefs)
  }
  return(risk_dt)
}

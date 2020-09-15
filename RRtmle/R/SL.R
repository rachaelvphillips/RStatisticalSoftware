
make_super_learner <- function(candidates, training_task, likelihood, loss_generator = make_eff_loss) {
  loss_function <- loss_generator(training_task, likelihood)
  make_learner(Pipeline, Lrnr_cv$new(candidates), Lrnr_cv_selector$new(loss_function))
}

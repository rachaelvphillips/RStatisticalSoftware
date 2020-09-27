
make_thresh_npsem_survival <- function(node_list, max_time, data_adaptive = F) {
  baseline_covariates <- node_list[["W"]]
  marker_covariates <- node_list[["A"]]
  outcome_covariate <- node_list[["Nt"]]
  censoring_indicator <- node_list[["At"]]

  nodes <- c("At", "Nt")
  competing_risks_indicator <- function(data, time, args, cols) {
    all_risks <- setdiff(cols, "t")
    parents <- args$parents


    past_jump_id <- data[t<time, last(.SD), .SDcols = all_risks, by = id]

    past_jump_id <- past_jump_id$id[rowSums(past_jump_id[, ..all_risks]) == 0]
    if(length(parents) > 0){
      cur_jump_id <- data[id %in% past_jump_id, last(.SD), .SDcols = parents, by = id]
      cur_jump_id <- cur_jump_id$id[rowSums(cur_jump_id[, ..parents]) == 0]
    } else {
      cur_jump_id <- past_jump_id
    }
    set(data, , "keep", as.numeric(data$id %in% cur_jump_id) )
    return(data[, c("id", "keep")])
  }
  risk_set_map_Nt <- Summary_measure$new(c(nodes, "t"), competing_risks_indicator, args_to_pass = list(parents = c()), group_by_id = F)
  risk_set_map_At <- Summary_measure$new(c(nodes, "t"), competing_risks_indicator, args_to_pass = list(parents = c("Nt")), group_by_id = F)


  times <- 1:max_time

  if(!data_adaptive) {
    npsem <- list(define_node("W", baseline_covariates, c(), time = 0),
                  define_node("A", marker_covariates, "W", time = 0),
                  define_node("Nt", outcome_covariate, c("W", "A"), risk_set_map = risk_set_map_Nt, time = times),
                  define_node("At", censoring_indicator, c("W", "A"), risk_set_map = risk_set_map_At, time = times))
  } else {
    npsem <- list(define_node("W", baseline_covariates, c()),
                  define_node("A_learned", outcome_covariate, c("A")),
                  define_node("A", marker_covariates, "W"),
                  define_node("Nt", outcome_covariate, c("W", "A"), risk_set_map = risk_set_map_Nt, time = times),
                  define_node("At", censoring_indicator, c("W", "A"), risk_set_map = risk_set_map_At, time = times))
  }
  return(npsem)

}

make_thresh_task_survival <- function(data, npsem, ...) {
  tmle3_Task$new(data, npsem, long_format = T, ...)
}

#' @export
make_thresh_likelihood_survival <- function(tmle_task, learner_list,
                                   cutoffs= function(A) {as.vector(quantile(A, seq(0.05, 0.95, length.out = 10)))},
                                   bins = 10, cv = T, marker_learner = NULL) {
  # covariates
  W_factor <- define_lf(LF_emp, "W")

  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }
  data_adaptive <- !is.null(marker_learner)
  if(data_adaptive) {
    learned_marker_node <- "A_learned"
    if(cv){
      short_lik <- Likelihood$new(LF_fit$new("A_learned", Lrnr_sl$new(marker_learner), type = "mean"))
    } else {
      short_lik <- Likelihood$new(LF_fit$new("A_learned", marker_learner, type = "mean"))
    }
    short_lik <- short_lik$train(tmle_task)
    Aval <- short_lik$get_likelihood(tmle_task, "A_learned")
    mnode <- "A"
  } else {
    learned_marker_node <- "A"
    short_lik <- NULL
    Aval <- tmle_task$get_tmle_node("A")
    mnode <- "A"
  }

  if(is.function(cutoffs)) {
    cutoffs <- cutoffs(Aval)
  }


  generator_A <- learner_marker_task_generator(learned_marker_node = learned_marker_node, learned_marker_var = "A", marker_node = "A", node = "A", data_adaptive = data_adaptive)

  generator_Nt <- learner_marker_task_generator(learned_marker_node = learned_marker_node, learned_marker_var = "A", marker_node = "A", node = "Nt", data_adaptive = data_adaptive, is_time_variant = T)
  generator_At <- learner_marker_task_generator(learned_marker_node = learned_marker_node, learned_marker_var = "A", marker_node = "A", node = "At", data_adaptive = data_adaptive, is_time_variant = T)


  A_factor <- define_lf(LF_derived, "A", learner = Lrnr_CDF$new(learner_list[["A"]], bins, cutoffs, cv = cv),short_lik, generator_A,  type = "mean", bound = A_bound)



  # outcome
  Nt_factor <- LF_derived2$new("Nt", Lrnr_thresh$new(learner_list[["Nt"]], mnode, cutoffs =cutoffs, cv = cv ),short_lik, generator_Nt, type = "mean")


  At_factor <- LF_derived2$new("At", learner_list[["At"]],short_lik, generator_At, type = "mean")

  # construct and train likelihood
  factor_list <- c(list(W_factor, A_factor, At_factor, Nt_factor),short_lik$factor_list)





  if (!is.null(tmle_task$npsem[["A"]]$censoring_node)) {
    stop("A is subject to censoring, this isn't supported yet")
  }

  if (!is.null(tmle_task$npsem[["W"]]$censoring_node)) {
    stop("W is subject to censoring, this isn't supported yet")
  }

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}

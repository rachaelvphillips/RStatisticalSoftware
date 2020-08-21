
#A general sampler from likelihood objects
#' @importFrom AR AR.Sim
#' @export

Sampler <- R6Class(
  classname = "Sampler",
  portable = TRUE,
  class = TRUE,
  active = list(
  params = function(){
    private$.params
  },
  time_ordering = function(){
    self$params$time_ordering
  },
  likelihood = function(){
    self$params$likelihood
  }
  ),
  public = list(
    initialize = function(likelihood, time_ordering, use_lf = c()){
      params <- sl3::args_to_list()
      private$.params <- params
    },
    sample = function(tmle_task, start_node, end_node){
      start_index <- which(self$time_ordering == start_node)
      end_index <- which(self$time_ordering == end_node)
      if(end_index < start_index){
        stop("Start and end node do not satisfy time ordering.")
      }
      for(i in start_index:end_index){
        node <- self$time_ordering[[i]]
        tmle_task <- self$next_in_chain(tmle_task, node)
      }
      return(tmle_task)

    },
    compute_conditional_mean = function(tmle_task, start_node, end_node, num_iter = 100){
      # Computes conditional mean of end_node given the (strict) past of start_node via monte carlo simulation
      # The sampling begins with start_node.
      outcomes = matrix(nrow = length(unique(tmle_task$id)), ncol = num_iter)
      for(i in 1:num_iter){
        new_task <- self$sample(tmle_task, start_node, end_node)
        outcomes[,i] <- new_task$get_tmle_node(end_node)[,end_node,with=F][[1]]
      }
      return(rowMeans(outcomes))
    },
    next_in_chain = function(tmle_task, node){
      #Takes a task and node name and then generated a new task
      # with the node values replaced with sampled values.
      #TODO if there are time pooled nodes then this doesnt work.
      samples <- self$sample_from_node(tmle_task, node, 1)
      samples <- as.data.table(samples)
      setnames(samples, node)
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), samples)
      return(cf_task)
    },
    sample_from_node = function(tmle_task, node, n_samples = 1, force_time_value = "NO", use_LF_factor = node %in% self$params$use_lf, fold_number = "full"){
      #Samples from a node. Returns matrix of n by n_samples of values.

      if(use_LF_factor){
        return(self$likelihood$factor_list[[node]]$sample(tmle_task, n_samples, fold_number))
      }
      outcome_type <- tmle_task$npsem[[node]]$variable_type

      times_to_pool <- tmle_task$npsem[[node]]$times_to_pool

      num_id <- length(unique(tmle_task$id))
      #num_rows <- ifelse(is.null(times_to_pool), num_id, times_to_pool*num_id)
      reg_task <- tmle_task$get_regression_task(node, drop_censored = F, force_time_value = force_time_value)

      node_variable <- reg_task$nodes$outcome
      at_risk_id <- reg_task$get_data(, c("id", "t", "at_risk", paste0("last_val_", reg_task$nodes$outcome)))
      # Only sample from those who are at_risk (so that their value may change)

      num_rows <- nrow(at_risk_id)

      if(outcome_type$type == "binomial"){
        # TODO If node is pooled across time and at_risk set changes then this must be done sequentially.
        # However using force_time_value arg, this can be handled.
        cf_outcome <- data.table(t= at_risk_id$t, id = at_risk_id$id, rep(1,num_rows))
        setnames(cf_outcome, c("t", "id", node_variable))

        cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cf_outcome)
        p <- self$likelihood$get_likelihood(cf_task, node)
        if(is.numeric(force_time_value)){
          p <- p[which(cf_task$get_regression_task(node)$data$t == force_time_value),]
        }
        p <- p[match(p$id, at_risk_id$id), node, with = F][[1]]
        values <- do.call(cbind, lapply(1:num_rows, function(i) {
          as.vector(rbinom(1, 1, p[i]))
        }))

      }
      else if (outcome_type$type == "categorical"){
        levels <- outcome_type$levels
        cf_tasks <- lapply(levels, function(level){
          cf_outcome <- data.table(rep(level,num_rows))
          cf_outcome <- data.table(t= at_risk_id$t, id = at_risk_id$t, cf_outcome)
          setnames(cf_outcome, c("t", "id", node))
          cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cf_outcome)
        })
        probs <- as.matrix(lapply(cf_tasks, function(cf_task){
          p <- self$likelihood$get_likelihood(cf_task, node)[, node, with = F][[1]]
        }))
        values <- apply(probs, 1, function(p){
            apply(
              rmultinom(n_samples, 1, p) == 1, 2,
              function(onehots) levels[which(onehots)]
            )
          })
      }
      else if (outcome_type$type == "continuous") {
        #TODO those who aren't at risk have degenerate likelihood
        if(!is.null(times_to_pool)){
          stop("Sampling from continuous variables is not supported for pooled time.")
        }
        #outcome <- tmle_task$get_tmle_node(node)
        ids <- at_risk_id$id

        values <- matrix(nrow = n_samples, ncol = num_rows)
        for (i in seq_along(ids)) {
          id <- ids[i]
          if(at_risk_id$at_risk[i]==0){
            #If this individual is not at risk then sample is just the last value
            samples <- rep(n_samples, at_risk_id$last_val[i] )
            values[, i] <- samples
            next
          } else {
            subject <- tmle_task[which(tmle_task$id == id)]
            f_X <- function(a) {

              #TODO might be more efficient to pass the regression task to likelihood so we dont recompute
              cf_data <- data.table(a)

              setnames(cf_data, names(cf_data), node)
              subject_a <- subject$generate_counterfactual_task(UUIDgenerate(), cf_data)
              likelihood <- self$likelihood$get_likelihood(subject_a, node)[, node, with = F][[1]]
              return(likelihood)
            }
            samples <- AR::AR.Sim(n_samples, f_X, Y.dist = "norm", Y.dist.par = c(mean(outcome_stripped), var(outcome_stripped)),
                                  xlim = c(min(outcome_stripped), max(outcome_stripped))
            )
            values[, i] <- samples
          }

        }

      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      values <- t(values)
      return(values)

    }
  ),
  private = list(
    .params = NULL
  )
)


#A general sampler from likelihood objects

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
  },
  ),
  public = list(
    initialize = function(likelihood, time_ordering){
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
      samples <- data.table(self$sample_from_node(tmle_task, node, 1))
      setnames(samples, node)
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), samples)
      return(cf_task)
    },
    sample_from_node = function(tmle_task, node, n_samples, use_LF_factor = F, fold_number = "full"){
      #Samples from a node. Returns matrix of n by n_samples of values.
      if(use_LF_factor){
        return(self$likelihood$factor_list[[node]]$sample(tmle_task, n_samples, fold_number))
      }
      outcome_type <- tmle_task$npsem[[node]]$variable_type
      times_to_pool <- tmle_task$npsem[[node]]$times_to_pool

      num_id <- length(unique(tmle_task$id))
      num_rows <- ifelse(is.null(times_to_pool), num_id, times_to_pool*num_id)
      if(outcome_type$type == "binomial"){
        cf_outcome <- data.table(rep(1,num_rows))
        setnames(cf_outcome, node)
        cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cf_outcome)
        p <- self$likelihood$get_likelihood(cf_task, node)[, node, with = F][[1]]
        values <- as.matrix(lapply(1:num_rows, function(i) {
          rbinom(n_samples, 1, p[i])
        }))
      }
      else if (outcome_type$type == "categorical"){
        levels <- outcome_type$levels
        cf_tasks <- lapply(levels, function(level){
          cf_outcome <- data.table(rep(level,num_rows))
          setnames(cf_outcome, node)
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
        if(!is.null(times_to_pool)){
          stop("Sampling from continuous variables is not supported for pooled time.")
        }
        outcome <- tmle_task$get_tmle_node(node)
        ids <- outcome$id
        values <- matrix(nrow = n_samples, ncol = num_rows)
        for (id in ids) {
          subject <- tmle_task[which(tmle_task$id == id)]
          f_X <- function(a) {
            #TODO might be more efficient to pass the regression task to likelihood so we dont recompute
            cf_data <- data.table(a)
            setnames(cf_data, names(cf_data), self$name)
            subject_a <- subject$generate_counterfactual_task(UUIDgenerate(), cf_data)
            likelihood <- self$likelihood$get_likelihood(subject_a, node)[, node, with = F][[1]]
            return(likelihood)
          }
          samples <- AR.Sim(n_samples, f_X,
                            xlim = c(min(outcome[, node, with = F][[1]]), max(outcome[, node, with = F][[1]]))
          )
          values[, i] <- samples
        }

      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      values <- t(values)
      return(values)

    }
  ),
  private = list(
    private$.params = NULL
  )
)

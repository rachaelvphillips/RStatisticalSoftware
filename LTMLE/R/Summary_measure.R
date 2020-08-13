#' @importFrom R6 R6Class
#' @import data.table

#' @export
Summary_measure <- R6Class(
  classname = "Summary_measure",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(column_names, summary_function, name = "Summary"){
        summary_function_wrap <- function(data){
          result <- summary_function(data)
          if(!is.data.table(result)){
            result <- data.table(matrix(result, nrow =1))

          }
          return(result)
        }
        params <- sl3::args_to_list()
        params$summary_function <- summary_function_wrap
        private$.params <- params
    },
    summarize = function(data, add_id = T){
      data <- private$.process_data(data, NULL)
      func <- private$.params$summary_function
      # Needed since  pass by promise would break next line
      data <- data[,]


      reduced_data <- data[,func(.SD), by = id,
                           .SDcols = self$params$column_names]

     #  num_sample <- length(unique(reduced_data$id))
     #  num_summary_vars <- nrow(reduced_data) / num_sample
     #  reduced_data$summary_id <- c(1:num_summary_vars, num_sample)
     #  reduced_data <- reshape(reduced_data, idvar = "id", timevar = "summary_id", direction = "wide")

       assertthat::assert_that(is.null(self$params$name) | ncol(reduced_data)-1 == length(self$params$name),
                              msg = "The summary measure names does not match length of summary measure function output.")
      if(!is.null(self$params$name)){
        colnames(reduced_data) <- c("id", self$params$name)
      }
      if(!add_id){
        reduced_data$id = NULL
      }
      return(reduced_data)
    }
  ),
  active = list(
    column_names = function(){
      self$params$column_names
    },
    name = function(){
      self$params$name
    },
    params = function(){
      private$.params
    }
  ),
  private = list(
    .params = NULL,
    .process_data = function(data, row_index){
      assertthat::assert_that("id" %in% colnames(data), msg = "Error: Column 'id' not found in data.")
      if(!is.data.table(data)){
        data = as.data.table(data)
      }

      if(is.null(row_index)){
        return(data)
      }
      return(data[row_index,])
    }
  )
)

make_summary_measure_NULL <- function(column_names = ""){
  name =  NULL
  summary_function <- function(data){
    return(data.table(NULL))

  }

  return(Summary_measure$new(column_names, summary_function, name))
}


make_summary_measure_FULL <- function(column_names){
  column_names <- union("t", column_names)
  name =  NULL

  summary_function <- function(data){
    t <- data$t
    data$t <- NULL
    data <- do.call(cbind, lapply(1:ncol(data), function(i){
      dat <- data.table::transpose(data[,i, with =F])
      colnames(dat) <- paste(colnames(data)[i], t, sep = "_")

      return(dat)}))
    return(data)

  }

  return(Summary_measure$new(column_names, summary_function, name))
}

make_summary_measure_baseline <- function(column_names){
  name = paste( column_names, "baseline", sep = "_")

  summary_function <- function(data){
    return(data[1,])
  }

  return(Summary_measure$new(column_names, summary_function, name))
}


make_summary_measure_last_value <- function(column_names){
  name = paste(column_names, "most_recent", sep = "_")

  summary_function <- function(data){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    return(data[nrow(data),])

  }

  return(Summary_measure$new(column_names, summary_function, name))
}


make_summary_measure_apply <- function(column_names, FUN){
  summary_function <- function(data){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    data <- as.data.table(t(apply(data, 2, FUN)))
    colnames(data) <- as.character(1:ncol(data))

    return(data)

  }
  return(Summary_measure$new(column_names, summary_function, paste(column_names, as.character(substitute(FUN)), sep = "_")))

}


make_summary_measure_running_average <- function(column_names){
  name = paste(column_names, "avg", sep = "_")
  return(make_summary_measure_apply(column_names, mean))
}

make_summary_measure_running_median <- function(column_names){
  name = paste(column_names, "median", sep = "_")
  return(make_summary_measure_apply(column_names,name ))
}
make_summary_measure_relative_difference_from_t0 <- function(column_names){
  summary_function <- function(data){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    return(data[nrow(data),] - data[1,])

  }
  return(Summary_measure$new(column_names, summary_function, paste(column_names, "rel_diff_t0")))

}

make_summary_measure_relative_difference_from_last_t <- function(column_names){
  name = paste(column_names, "rel_diff_last_t", sep = "_")

  summary_function <- function(data){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    data <- data[nrow(data) - data[nrow(data)-1,],]
    data <- data.table(matrix(data,nrow=1))
  }
  return(Summary_measure$new(column_names, summary_function, name))

}

make_summary_measure_slope <- function(column_names){
  name = paste(column_names, "slope_in_t", sep = "_")

  summary_function <- function(data){

    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
   if("t" %in% colnames(data)){
     t = data[,"t",with = F]
   }
    else{
      t = 1:nrow(data)
      data = cbind(t, data)
    }


    slopes = sapply(colnames(data)[-1], function(name){
      return(as.vector(coef(lm(as.formula(paste(name, "~ t")), data.frame(data)))[2]))
    })


  }
  return(Summary_measure$new(column_names, summary_function, name))
}

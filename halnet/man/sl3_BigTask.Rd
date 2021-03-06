% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sl3_BigTask.R
\docType{class}
\name{sl3_BigTask}
\alias{sl3_BigTask}
\title{Define a Machine Learning Task}
\format{
\code{\link{R6Class}} object.
}
\value{
\code{sl3_Task} object
}
\description{
An increasingly less thin wrapper around a \code{data.table} containing the
data. Contains metadata about the particular machine learning problem,
including which variables are to be used as covariates and outcomes.
}
\keyword{data}
\section{Super class}{
\code{\link[sl3:sl3_Task]{sl3::sl3_Task}} -> \code{sl3_BigTask}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{sl3_BigTask$new()}}
\item \href{#method-add_interactions}{\code{sl3_BigTask$add_interactions()}}
\item \href{#method-add_columns}{\code{sl3_BigTask$add_columns()}}
\item \href{#method-next_in_chain}{\code{sl3_BigTask$next_in_chain()}}
\item \href{#method-subset_task}{\code{sl3_BigTask$subset_task()}}
\item \href{#method-get_data}{\code{sl3_BigTask$get_data()}}
\item \href{#method-has_node}{\code{sl3_BigTask$has_node()}}
\item \href{#method-get_node}{\code{sl3_BigTask$get_node()}}
\item \href{#method-offset_transformed}{\code{sl3_BigTask$offset_transformed()}}
\item \href{#method-print}{\code{sl3_BigTask$print()}}
\item \href{#method-revere_fold_task}{\code{sl3_BigTask$revere_fold_task()}}
\item \href{#method-clone}{\code{sl3_BigTask$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$new(
  data,
  covariates,
  outcome = NULL,
  bm = NULL,
  bm_columns = NULL,
  outcome_type = NULL,
  outcome_levels = NULL,
  id = NULL,
  weights = NULL,
  offset = NULL,
  time = NULL,
  nodes = NULL,
  column_names = NULL,
  row_index = NULL,
  folds = NULL,
  flag = TRUE,
  drop_missing_outcome = FALSE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add_interactions"></a>}}
\if{latex}{\out{\hypertarget{method-add_interactions}{}}}
\subsection{Method \code{add_interactions()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$add_interactions(interactions, warn_on_existing = TRUE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add_columns"></a>}}
\if{latex}{\out{\hypertarget{method-add_columns}{}}}
\subsection{Method \code{add_columns()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$add_columns(new_data, column_uuid = uuid::UUIDgenerate())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-next_in_chain"></a>}}
\if{latex}{\out{\hypertarget{method-next_in_chain}{}}}
\subsection{Method \code{next_in_chain()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$next_in_chain(
  bm_columns = NULL,
  covariates = NULL,
  outcome = NULL,
  id = NULL,
  weights = NULL,
  offset = NULL,
  time = NULL,
  folds = NULL,
  column_names = NULL,
  new_nodes = NULL,
  ...
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-subset_task"></a>}}
\if{latex}{\out{\hypertarget{method-subset_task}{}}}
\subsection{Method \code{subset_task()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$subset_task(row_index, drop_folds = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_data"></a>}}
\if{latex}{\out{\hypertarget{method-get_data}{}}}
\subsection{Method \code{get_data()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$get_data(rows = NULL, columns, expand_factors = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-has_node"></a>}}
\if{latex}{\out{\hypertarget{method-has_node}{}}}
\subsection{Method \code{has_node()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$has_node(node_name)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_node"></a>}}
\if{latex}{\out{\hypertarget{method-get_node}{}}}
\subsection{Method \code{get_node()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$get_node(node_name, generator_fun = NULL, expand_factors = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-offset_transformed"></a>}}
\if{latex}{\out{\hypertarget{method-offset_transformed}{}}}
\subsection{Method \code{offset_transformed()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$offset_transformed(link_fun = NULL, for_prediction = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-revere_fold_task"></a>}}
\if{latex}{\out{\hypertarget{method-revere_fold_task}{}}}
\subsection{Method \code{revere_fold_task()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$revere_fold_task(fold_number)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{sl3_BigTask$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

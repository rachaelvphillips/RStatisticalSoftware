% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Likelihood_cache.R
\docType{class}
\name{Likelihood_cache}
\alias{Likelihood_cache}
\title{Cache Likelihood values, update those values
NOTE: Added 'node' argument. This is neccessary if likelihood_factor is a pooled LF
so that likelihood factors do not uniquely characterize nodes, and a single LF can be associated with multiple nodes.
E.g. generate multiple predictions.}
\description{
Cache Likelihood values, update those values
NOTE: Added 'node' argument. This is neccessary if likelihood_factor is a pooled LF
so that likelihood factors do not uniquely characterize nodes, and a single LF can be associated with multiple nodes.
E.g. generate multiple predictions.

Cache Likelihood values, update those values
NOTE: Added 'node' argument. This is neccessary if likelihood_factor is a pooled LF
so that likelihood factors do not uniquely characterize nodes, and a single LF can be associated with multiple nodes.
E.g. generate multiple predictions.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Likelihood_cache$new()}}
\item \href{#method-tasks_at_step}{\code{Likelihood_cache$tasks_at_step()}}
\item \href{#method-get_update_step}{\code{Likelihood_cache$get_update_step()}}
\item \href{#method-key}{\code{Likelihood_cache$key()}}
\item \href{#method-set_values}{\code{Likelihood_cache$set_values()}}
\item \href{#method-get_values}{\code{Likelihood_cache$get_values()}}
\item \href{#method-cache_lf}{\code{Likelihood_cache$cache_lf()}}
\item \href{#method-cache_task}{\code{Likelihood_cache$cache_task()}}
\item \href{#method-set_best}{\code{Likelihood_cache$set_best()}}
\item \href{#method-update_best}{\code{Likelihood_cache$update_best()}}
\item \href{#method-clone}{\code{Likelihood_cache$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-tasks_at_step"></a>}}
\if{latex}{\out{\hypertarget{method-tasks_at_step}{}}}
\subsection{Method \code{tasks_at_step()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$tasks_at_step(current_step)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_update_step"></a>}}
\if{latex}{\out{\hypertarget{method-get_update_step}{}}}
\subsection{Method \code{get_update_step()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$get_update_step(
  likelihood_factor,
  tmle_task,
  fold_number,
  node = ""
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-key"></a>}}
\if{latex}{\out{\hypertarget{method-key}{}}}
\subsection{Method \code{key()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$key(likelihood_factor, tmle_task, fold_number, node = "")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_values"></a>}}
\if{latex}{\out{\hypertarget{method-set_values}{}}}
\subsection{Method \code{set_values()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$set_values(
  likelihood_factor,
  tmle_task,
  update_step = 0,
  fold_number,
  values,
  node = ""
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_values"></a>}}
\if{latex}{\out{\hypertarget{method-get_values}{}}}
\subsection{Method \code{get_values()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$get_values(
  likelihood_factor,
  tmle_task,
  fold_number,
  node = ""
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cache_lf"></a>}}
\if{latex}{\out{\hypertarget{method-cache_lf}{}}}
\subsection{Method \code{cache_lf()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$cache_lf(likelihood_factor)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cache_task"></a>}}
\if{latex}{\out{\hypertarget{method-cache_task}{}}}
\subsection{Method \code{cache_task()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$cache_task(task)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_best"></a>}}
\if{latex}{\out{\hypertarget{method-set_best}{}}}
\subsection{Method \code{set_best()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$set_best()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-update_best"></a>}}
\if{latex}{\out{\hypertarget{method-update_best}{}}}
\subsection{Method \code{update_best()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$update_best()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Likelihood_cache$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

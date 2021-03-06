% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LF_collapsed.R
\name{LF_collapsed}
\alias{LF_collapsed}
\title{Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.}
\format{
\code{\link{R6Class}} object.
}
\value{
\code{LF_base} object
}
\description{
Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.

Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
}
\seealso{
Other Likelihood objects: 
\code{\link{CF_Likelihood_pooled}},
\code{\link{LF_emp}},
\code{\link{LF_fit_pooled}},
\code{\link{LF_static}},
\code{\link{Targeted_Likelihood_pooled}}
}
\concept{Likelihood objects}
\keyword{data}
\section{Super class}{
\code{\link[tmle3:LF_base]{tmle3::LF_base}} -> \code{LF_collapsed}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{LF_collapsed$new()}}
\item \href{#method-get_mean}{\code{LF_collapsed$get_mean()}}
\item \href{#method-get_density}{\code{LF_collapsed$get_density()}}
\item \href{#method-clone}{\code{LF_collapsed$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="tmle3" data-topic="LF_base" data-id="cf_values">}\href{../../tmle3/html/LF_base.html#method-cf_values}{\code{tmle3::LF_base$cf_values()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="tmle3" data-topic="LF_base" data-id="delayed_train">}\href{../../tmle3/html/LF_base.html#method-delayed_train}{\code{tmle3::LF_base$delayed_train()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="tmle3" data-topic="LF_base" data-id="get_likelihood">}\href{../../tmle3/html/LF_base.html#method-get_likelihood}{\code{tmle3::LF_base$get_likelihood()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="tmle3" data-topic="LF_base" data-id="print">}\href{../../tmle3/html/LF_base.html#method-print}{\code{tmle3::LF_base$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="tmle3" data-topic="LF_base" data-id="sample">}\href{../../tmle3/html/LF_base.html#method-sample}{\code{tmle3::LF_base$sample()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="tmle3" data-topic="LF_base" data-id="train">}\href{../../tmle3/html/LF_base.html#method-train}{\code{tmle3::LF_base$train()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LF_collapsed$new(name, factor_list, ..., type = "density")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_mean"></a>}}
\if{latex}{\out{\hypertarget{method-get_mean}{}}}
\subsection{Method \code{get_mean()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LF_collapsed$get_mean(tmle_task, fold_number)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_density"></a>}}
\if{latex}{\out{\hypertarget{method-get_density}{}}}
\subsection{Method \code{get_density()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LF_collapsed$get_density(tmle_task, fold_number)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LF_collapsed$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

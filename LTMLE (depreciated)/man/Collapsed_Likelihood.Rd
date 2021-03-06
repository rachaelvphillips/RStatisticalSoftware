% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Collapsed_Likelihood.R
\name{Collapsed_Likelihood}
\alias{Collapsed_Likelihood}
\title{A likelihood that represents a collapsed likelihood,
meaning that some factors of the underlying likelihood are collapsed into one factor by
multiplication of probabilities.}
\description{
A likelihood that represents a collapsed likelihood,
meaning that some factors of the underlying likelihood are collapsed into one factor by
multiplication of probabilities.

A likelihood that represents a collapsed likelihood,
meaning that some factors of the underlying likelihood are collapsed into one factor by
multiplication of probabilities.
}
\section{Super classes}{
\code{\link[sl3:Lrnr_base]{sl3::Lrnr_base}} -> \code{\link[LTMLE:Likelihood_pooled]{LTMLE::Likelihood_pooled}} -> \code{Collapsed_Likelihood}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Collapsed_Likelihood$new()}}
\item \href{#method-get_hidden_nodes}{\code{Collapsed_Likelihood$get_hidden_nodes()}}
\item \href{#method-clone}{\code{Collapsed_Likelihood$clone()}}
}
}
\if{html}{
\out{<details ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="assert_trained">}\href{../../sl3/html/Lrnr_base.html#method-assert_trained}{\code{sl3::Lrnr_base$assert_trained()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="base_chain">}\href{../../sl3/html/Lrnr_base.html#method-base_chain}{\code{sl3::Lrnr_base$base_chain()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="base_predict">}\href{../../sl3/html/Lrnr_base.html#method-base_predict}{\code{sl3::Lrnr_base$base_predict()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="chain">}\href{../../sl3/html/Lrnr_base.html#method-chain}{\code{sl3::Lrnr_base$chain()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="custom_chain">}\href{../../sl3/html/Lrnr_base.html#method-custom_chain}{\code{sl3::Lrnr_base$custom_chain()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="get_outcome_range">}\href{../../sl3/html/Lrnr_base.html#method-get_outcome_range}{\code{sl3::Lrnr_base$get_outcome_range()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="get_outcome_type">}\href{../../sl3/html/Lrnr_base.html#method-get_outcome_type}{\code{sl3::Lrnr_base$get_outcome_type()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="predict">}\href{../../sl3/html/Lrnr_base.html#method-predict}{\code{sl3::Lrnr_base$predict()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="predict_fold">}\href{../../sl3/html/Lrnr_base.html#method-predict_fold}{\code{sl3::Lrnr_base$predict_fold()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="set_train">}\href{../../sl3/html/Lrnr_base.html#method-set_train}{\code{sl3::Lrnr_base$set_train()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="subset_covariates">}\href{../../sl3/html/Lrnr_base.html#method-subset_covariates}{\code{sl3::Lrnr_base$subset_covariates()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="train">}\href{../../sl3/html/Lrnr_base.html#method-train}{\code{sl3::Lrnr_base$train()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="sl3" data-topic="Lrnr_base" data-id="train_sublearners">}\href{../../sl3/html/Lrnr_base.html#method-train_sublearners}{\code{sl3::Lrnr_base$train_sublearners()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="add_factors">}\href{../../LTMLE/html/Likelihood_pooled.html#method-add_factors}{\code{LTMLE::Likelihood_pooled$add_factors()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="base_train">}\href{../../LTMLE/html/Likelihood_pooled.html#method-base_train}{\code{LTMLE::Likelihood_pooled$base_train()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="get_likelihood">}\href{../../LTMLE/html/Likelihood_pooled.html#method-get_likelihood}{\code{LTMLE::Likelihood_pooled$get_likelihood()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="get_likelihoods">}\href{../../LTMLE/html/Likelihood_pooled.html#method-get_likelihoods}{\code{LTMLE::Likelihood_pooled$get_likelihoods()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="get_possible_counterfactuals">}\href{../../LTMLE/html/Likelihood_pooled.html#method-get_possible_counterfactuals}{\code{LTMLE::Likelihood_pooled$get_possible_counterfactuals()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="print">}\href{../../LTMLE/html/Likelihood_pooled.html#method-print}{\code{LTMLE::Likelihood_pooled$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="sample">}\href{../../LTMLE/html/Likelihood_pooled.html#method-sample}{\code{LTMLE::Likelihood_pooled$sample()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="LTMLE" data-topic="Likelihood_pooled" data-id="validate_task">}\href{../../LTMLE/html/Likelihood_pooled.html#method-validate_task}{\code{LTMLE::Likelihood_pooled$validate_task()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Collapsed_Likelihood$new(hidden_likelihood, collapse_by, cache = NULL, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_hidden_nodes"></a>}}
\if{latex}{\out{\hypertarget{method-get_hidden_nodes}{}}}
\subsection{Method \code{get_hidden_nodes()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Collapsed_Likelihood$get_hidden_nodes(node)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Collapsed_Likelihood$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

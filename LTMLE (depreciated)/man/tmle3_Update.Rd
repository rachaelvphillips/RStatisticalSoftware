% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ltmle3_Update.R
\name{tmle3_Update}
\alias{tmle3_Update}
\title{Defines an update procedure (submodel+loss function)}
\description{
Current Limitations:
loss function and submodel are hard-coded (need to accept arguments for these)
}
\section{Constructor}{

  \code{define_param(maxit, cvtmle, one_dimensional, constrain_step, delta_epsilon, verbose)}

  \describe{
    \item{\code{maxit}}{The maximum number of update iterations
    }
    \item{\code{cvtmle}}{If \code{TRUE}, use CV-likelihood values when
       calculating updates.
    }
    \item{\code{one_dimensional}}{If \code{TRUE}, collapse clever covariates
       into a one-dimensional clever covariate scaled by the mean of their
       EIFs.
    }
    \item{\code{constrain_step}}{If \code{TRUE}, step size is at most
       \code{delta_epsilon} (it can be smaller if a smaller step decreases
       the loss more).
    }
    \item{\code{delta_epsilon}}{The maximum step size allowed if
       \code{constrain_step} is \code{TRUE}.
    }
    \item{\code{convergence_type}}{The convergence criterion to use: (1)
       \code{"scaled_var"} corresponds to sqrt(Var(D)/n)/logn (the default)
       while (2) \code{"sample_size"} corresponds to 1/n.
    }
    \item{\code{fluctuation_type}}{Whether to include the auxiliary covariate
       for the fluctuation model as a covariate or to treat it as a weight.
       Note that the option \code{"weighted"} is incompatible with a
       multi-epsilon submodel (\code{one_dimensional = FALSE}).
    }
    \item{\code{use_best}}{If \code{TRUE}, the final updated likelihood is set to the
       likelihood that minimizes the ED instead of the likelihood at the last update
       step.
    }
    \item{\code{verbose}}{If \code{TRUE}, diagnostic output is generated
       about the updating procedure.
    }
    }
}

\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{tmle3_Update$new()}}
\item \href{#method-collapse_covariates}{\code{tmle3_Update$collapse_covariates()}}
\item \href{#method-update_step}{\code{tmle3_Update$update_step()}}
\item \href{#method-key_to_node_bundle}{\code{tmle3_Update$key_to_node_bundle()}}
\item \href{#method-node_bundle_to_key}{\code{tmle3_Update$node_bundle_to_key()}}
\item \href{#method-generate_submodel_data}{\code{tmle3_Update$generate_submodel_data()}}
\item \href{#method-fit_submodel}{\code{tmle3_Update$fit_submodel()}}
\item \href{#method-apply_submodel}{\code{tmle3_Update$apply_submodel()}}
\item \href{#method-debundle_submodel}{\code{tmle3_Update$debundle_submodel()}}
\item \href{#method-apply_update}{\code{tmle3_Update$apply_update()}}
\item \href{#method-check_convergence}{\code{tmle3_Update$check_convergence()}}
\item \href{#method-update_best}{\code{tmle3_Update$update_best()}}
\item \href{#method-update}{\code{tmle3_Update$update()}}
\item \href{#method-register_param}{\code{tmle3_Update$register_param()}}
\item \href{#method-clone}{\code{tmle3_Update$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$new(
  maxit = 100,
  cvtmle = TRUE,
  one_dimensional = FALSE,
  constrain_step = FALSE,
  delta_epsilon = 1e-04,
  convergence_type = c("scaled_var", "sample_size"),
  fluctuation_type = c("standard", "weighted"),
  optim_delta_epsilon = TRUE,
  use_best = FALSE,
  verbose = FALSE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-collapse_covariates"></a>}}
\if{latex}{\out{\hypertarget{method-collapse_covariates}{}}}
\subsection{Method \code{collapse_covariates()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$collapse_covariates(ED, clever_covariates)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-update_step"></a>}}
\if{latex}{\out{\hypertarget{method-update_step}{}}}
\subsection{Method \code{update_step()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$update_step(
  likelihood,
  tmle_task,
  fold_number = "full",
  update_nodes = NULL
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-key_to_node_bundle"></a>}}
\if{latex}{\out{\hypertarget{method-key_to_node_bundle}{}}}
\subsection{Method \code{key_to_node_bundle()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$key_to_node_bundle(key)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-node_bundle_to_key"></a>}}
\if{latex}{\out{\hypertarget{method-node_bundle_to_key}{}}}
\subsection{Method \code{node_bundle_to_key()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$node_bundle_to_key(node_bundle)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-generate_submodel_data"></a>}}
\if{latex}{\out{\hypertarget{method-generate_submodel_data}{}}}
\subsection{Method \code{generate_submodel_data()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$generate_submodel_data(
  likelihood,
  tmle_task,
  fold_number = "full",
  update_node_key = "Y",
  drop_censored = FALSE,
  for_fitting = F
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit_submodel"></a>}}
\if{latex}{\out{\hypertarget{method-fit_submodel}{}}}
\subsection{Method \code{fit_submodel()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$fit_submodel(submodel_data)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-apply_submodel"></a>}}
\if{latex}{\out{\hypertarget{method-apply_submodel}{}}}
\subsection{Method \code{apply_submodel()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$apply_submodel(submodel_data, epsilon, submodel)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-debundle_submodel"></a>}}
\if{latex}{\out{\hypertarget{method-debundle_submodel}{}}}
\subsection{Method \code{debundle_submodel()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$debundle_submodel(bundle, node)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-apply_update"></a>}}
\if{latex}{\out{\hypertarget{method-apply_update}{}}}
\subsection{Method \code{apply_update()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$apply_update(
  tmle_task,
  likelihood,
  fold_number,
  new_epsilon,
  update_node
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-check_convergence"></a>}}
\if{latex}{\out{\hypertarget{method-check_convergence}{}}}
\subsection{Method \code{check_convergence()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$check_convergence(tmle_task, fold_number = "full")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-update_best"></a>}}
\if{latex}{\out{\hypertarget{method-update_best}{}}}
\subsection{Method \code{update_best()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$update_best(likelihood)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-update"></a>}}
\if{latex}{\out{\hypertarget{method-update}{}}}
\subsection{Method \code{update()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$update(likelihood, tmle_task, update_spec = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-register_param"></a>}}
\if{latex}{\out{\hypertarget{method-register_param}{}}}
\subsection{Method \code{register_param()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$register_param(new_params)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{tmle3_Update$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

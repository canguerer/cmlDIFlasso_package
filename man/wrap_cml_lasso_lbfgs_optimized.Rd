% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{wrap_cml_lasso_lbfgs_optimized}
\alias{wrap_cml_lasso_lbfgs_optimized}
\title{Title
Wrapper function for use in simulation studied, taking all parameters,
passing them to sim_scenario_DIF_data and computing results via
cml_lasso_win_lbfgs_optimized}
\usage{
wrap_cml_lasso_lbfgs_optimized(
  n,
  I,
  nvars,
  seed,
  nlambdas,
  betavalues = rep(0, I),
  deltavalues = matrix(c(-0.5, 0, 0, 0, 0.5), nvars, I),
  thetavalues = rep(0, n),
  covariatevalues = matrix(rnorm(n * nvars), n, nvars),
  no_of_cores_used = "unknown"
)
}
\arguments{
\item{n}{number of persons (only used of if not specified below)}

\item{I}{number of items (only used of if not specified below)}

\item{nvars}{number of covariates (only used if not specified below)}

\item{seed}{seed for data simulation}

\item{betavalues}{Vector with true values of the item parameters.}

\item{deltavalues}{Matrix with the true values of the delta parameters (ncovs x nitems)}

\item{covariatevalues}{Matrix of covariate values of all persons (n x ncovs)}

\item{no_of_cores_used}{Number of cores to be used for parallel computing}
}
\value{
list of results and simulation scenario details.
}
\description{
Title
Wrapper function for use in simulation studied, taking all parameters,
passing them to sim_scenario_DIF_data and computing results via
cml_lasso_win_lbfgs_optimized
}

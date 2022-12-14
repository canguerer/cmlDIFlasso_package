% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{degrees-of-freedom-helpers}
\alias{degrees-of-freedom-helpers}
\alias{df_pen_unequal_zero}
\alias{df_veclength}
\alias{df_lambda}
\title{df calculation for RM-DIF-model with penalization}
\usage{
df_pen_unequal_zero(deltas, DIF_type, tolerance = 0)

df_veclength(deltas, DIF_type, delta_unres)

df_lambda(deltas, DIF_type, delta_unres, tolerance)
}
\arguments{
\item{deltas}{values for DIF-parameters as matrix for certain Lambda
per item (columns) per covariate (rows)}

\item{DIF_type}{for which Lasso-term should the non-nil values be calculated?
options: "items", "variables", "all.interactions"}

\item{tolerance}{Value under which parameters are considered equal
to zero for the calculation of degreees of freedom}

\item{delta_unres}{values for DIF-parameters as matrix form unrestricted cml solution (lambda = 0)}
}
\value{
calculated first (df_pen_unequal_zero) and second part (df_veclength)
of df_delta, or already the complete degrees of freedom (df_lambda) for a
lasso penalized (cml) model according to Yuan and Lin(2006)
 and Zou et al. (2007); remark: 10^(-6) will be zero (default value of argument tolerance)
}
\description{
calculates the first and second parts of df_delta depending on DIF-type given
parameter estimates for delta under a certain lambda and the unrestricted DIF-model.
}
\details{
Title
}

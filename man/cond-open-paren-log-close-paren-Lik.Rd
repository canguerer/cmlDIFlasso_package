% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{cond(log)Lik}
\alias{cond(log)Lik}
\title{conditional(log)Likelihoods}
\arguments{
\item{x}{0/1-vector of pattern of one person}

\item{Y}{0/1-matrix of patterns of all persons (npersons x nitems)}

\item{betas}{vector of item easiness parameters, either for a person,
or a Matrix of calculated itemparameters for all persons already including
DIF-effects (npersons x nitems)}
}
\value{
conditional Likelihood/log conditional Likelihood for the person
or the sample
}
\description{
Calculates the cLik/logcLik for one person or a whole sample given
parameters beta.
}

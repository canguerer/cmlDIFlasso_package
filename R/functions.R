###
### To Do:


##################################
##
#' conditional(log)Likelihoods
#'
#' Calculates the cLik/logcLik for one person or a whole sample given
#' parameters beta.
#'
#' @param x 0/1-vector of pattern of one person
#' @param Y  0/1-matrix of patterns of all persons (npersons x nitems)
#' @param betas vector of item easiness parameters, either for a person,
#' or a Matrix of calculated itemparameters for all persons already including
#' DIF-effects (npersons x nitems)
#'
#' @return conditional Likelihood/log conditional Likelihood for the person
#' or the sample
#' @name cond(log)Lik
NULL

##################################
##
## function:
## name: cL_person_2, cl_person_2
## cL for one person via psychotools

##################################
##
#' conditional(log)Likelihoods
#'
#' Calculates the cLik/logcLik for one person or a whole sample given
#' parameters beta.
#' cL for one person via psychotools
#'
#' @param x 0/1-vector of pattern of one person
#' @param Y  0/1-matrix of patterns of all persons (npersons x nitems)
#' @param betas vector of item easiness parameters, either for a person,
#' or a Matrix of calculated itemparameters for all persons already including
#' DIF-effects (npersons x nitems)
#' @return conditional Likelihood for the person
#' @export

cL_person_2 <- function(x, betas) {

  numerator <- sum((betas * x), na.rm = TRUE)

  betas_symm <- betas[!is.na(x)]

  denominator <- log(elementary_symmetric_functions(-betas_symm, order = 0)[[1]][sum(x, na.rm = TRUE)+1])
#  denominator <- log(psychotools:::elementary_symmetric_functions(-betas_symm, order = 0)[[1]][sum(x, na.rm = TRUE)+1])
  cL <-  exp(numerator-denominator)

  return(cL)
}


##################################
##
#' conditional(log)Likelihoods
#'
#' Calculates the cLik/logcLik for one person or a whole sample given
#' parameters beta.
#' cl for one person via psychotools
#'
#' @param x 0/1-vector of pattern of one person
#' @param Y  0/1-matrix of patterns of all persons (npersons x nitems)
#' @param betas vector of item easiness parameters, either for a person,
#' or a Matrix of calculated itemparameters for all persons already including
#' DIF-effects (npersons x nitems)
#' @return log conditional Likelihood for the person
#' @export

cl_person_2 <- function(x, betas) {

  numerator <- sum((betas * x), na.rm = TRUE)

  betas_symm <- betas[!is.na(x)]

  denominator <- log(elementary_symmetric_functions(-betas_symm, order = 0)[[1]][sum(x, na.rm = TRUE)+1])
#  denominator <- log(psychotools:::elementary_symmetric_functions(-betas_symm, order = 0)[[1]][sum(x, na.rm = TRUE)+1])
  cl <-  numerator-denominator

  return(cl)
}

##################################
##
## function:
## name: cL_sampe_2, cl_sample_2

## adding them all together
## betas: Matrix of all calculated itemparameters already including DIF-effects
## Y: matrix of answers of all persons to all items

##################################
##
#' conditional(log)Likelihoods
#'
#' Calculates the cLik/logcLik for one person or a whole sample given
#' parameters beta.
#' Conditional Likelihood of the whole sample, adding them all together.
#'
#' @param x 0/1-vector of pattern of one person
#' @param Y  0/1-matrix of patterns of all persons (npersons x nitems)
#' @param betas vector of item easiness parameters, either for a person,
#' or a Matrix of calculated itemparameters for all persons already including
#' DIF-effects (npersons x nitems)
#' @return conditional Likelihood/for the sample
#' @export


cL_sample_2 <- function(Y, betas, ...) {
  perperson <- lapply(1:nrow(Y), function(x) cl_person_2(Y[x,],betas[x,] ))
  Liksample <- exp(do.call(sum,perperson))
  return(Liksample)
}


##################################
##
#' conditional(log)Likelihoods
#'
#' Calculates the cLik/logcLik for one person or a whole sample given
#' parameters beta.
#' logarithm of conditional Likelihood of the whole sample, adding them all together.
#'
#' @param x 0/1-vector of pattern of one person
#' @param Y  0/1-matrix of patterns of all persons (npersons x nitems)
#' @param betas vector of item easiness parameters, either for a person,
#' or a Matrix of calculated itemparameters for all persons already including
#' DIF-effects (npersons x nitems)
#' @return logarithm of the conditional Likelihood for the sample
#' @export


cl_sample_2 <- function(Y, betas, ...) {
  perperson <- lapply(1:nrow(Y), function(x) (cl_person_2(Y[x,],betas[x,])))
  logliksample <- do.call(sum, perperson)
  return(logliksample)
}

##################################
##
## function:
## name: penalization

#' penalization term
#'
#' Calculates the certain amount of penalization depeding on DIF-type,
#' lambda and for restriction of intercepts added
#'
#' @param deltas A (ncovars x nitems) matrix of DIF-parameters
#' @param beta_intercept A vector of #items parameters
#' @param DIF_type can be "items", "variables", "all.interactions", "none"
#' @param restr_pen For additional penalization to ensure zero sum of intercepts
#' @param lambda penalization parameter
#' @param lambda_restr penalization parameter for restriction
#'
#' @return A numeric value for penalization of the Loglikelihood
#' @export

penalization <- function(deltas,
                         beta_intercept,
                         DIF_type = c("items", "variables", "all.interactions", "none"),
                         restr_pen = FALSE,
                         lambda = 0,
                         lambda_restr = 0) {

  #  DIF_type <- match.arg(DIF_type)
  #  if (!(DIF_type %in% c("items", "variables", "all.interactions", "none"))) {
  #    stop("Unknown DIF_type")
  #  }

  penalization <- 0
  penalization_restr <- 0

  if (DIF_type == "items") {
    penalization <- sum(sqrt(colSums(deltas^2)))
  }

  if (DIF_type == "variables") {
    penalization <- sum(sqrt(rowSums(deltas^2)))
  }

  if (DIF_type == "all.interactions") {
    penalization <-  sum(abs(deltas))
  }

  if (DIF_type == "none") {
    penalization <-  0
  }

  if (restr_pen == TRUE) {
    penalization_restr <- abs(sum(beta_intercept))
  }

  if (restr_pen == FALSE) {
    penalization_restr <- 0
  }

  penalization_complete <- lambda * penalization + lambda_restr * penalization_restr

  return(penalization_complete)
}



##################################
##
## function:
## name: calc_itemvecs


#' Calculation Item parameters per person
#'
#' Calculates beta vectors per person based on covariate values.
#'
#' @param betas Itemparameter intercepts
#' @param delta DIF-parameter-matrix (nrow = ncovariates, ncol =  nitems)
#' @param covariates Matrix of covariate values of all persons (n x num. of covs)
#'
#' @return Matrix with itemparameters per person.
#' @export
#'

calc_itemvecs <- function( betas, delta, covariates) {
  DIF_part <- covariates %*% delta
  beta_new <- t(betas + t(DIF_part))
  return(beta_new)
}


##################################
##
## function:
## name: expand_grid_JOHN

#' Faster verson of expand.grid for use in pairwise_prods
#'
#' Does the same as expand.grid in base R, just a lot faster.
#' Code is the fastest version taken from Stackoverflow:
#' https://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid
#'
#'
#' @param seq1 first vector
#' @param seq2 second vector
#'
#' @return A matrix with two columns, the rows containing every possible combination
#' of the two vectors seq1 and seq2 ine element each.
#'
#' @export
#'

expand_grid_JOHN <- function(seq1,seq2) {
  cbind(Var1 = rep.int(seq1, length(seq2)),
        Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))))
}




##################################
##
## function: pairwise_prods


#' Pairwise products of two vectors
#'
#' Given vectors 1 and 2, delivers matrix with pairwise products.
#' Necessary for calculation of the score function.
#'
#' @param vec1 first vector
#' @param vec2 second vector
#'
#' @return Matrix with all products (length(vec1) = nrows, length(vec2) = ncols))
#' @export
#'

pairwise_prods <- function(vec1,vec2){
  matrix((expand_grid_JOHN(vec1,vec2)[,1] * expand_grid_JOHN(vec1,vec2)[,2]),
         length(vec1),
         length(vec2))
}

##################################
##
## function: score_func_person



#' Score function for one person
#'
#' Calculates the derivatives of the conditional Loglikelihood for all
#' intercept and delta parameters.
#'
#' @param beta_person Parameter vector of a certain person.
#' @param y_person Pattern of answers of a person.
#' @param covariates_person Covariate values of a person.
#'
#' @return a list containing the values of the score function
#' \item{scoref_betas}{Values of the score function for every intercept.}
#' \item{scoref_delta}{Values of the score function for every delta-parameter (ncovs X nitems)}
#' @export
#'

score_func_person <- function(beta_person,
                              y_person,
                              covariates_person) {
  # calc #correct for person
  r_p <- sum(y_person)

  esf_person <- elementary_symmetric_functions(-beta_person, order = 1)

  scoref_betas <-  y_person - (esf_person[[2]][r_p+1,]  / esf_person[[1]][r_p+1])

  scoref_deltas <- pairwise_prods(scoref_betas, covariates_person)
  ### noch Benennungen für Spalten und zeilen einbauen

  out <- list(scoref_betas = scoref_betas,
              scoref_deltas = scoref_deltas) # vector for betas, matrix for deltas

}

##################################
##
## function: score_func_sample (preferable)



#' Score values of the sample per parameter
#'
#' Calculates the derivatives of the conditional Loglikelihood with respect
#' to each intercept and delta parameter for the full sample.
#'
#' @param beta_intercepts Vecotr of beta intercept values.
#' @param delta A (ncovs x nitems) matrix containing the DIF-parameter values.
#' @param covariates A (npersons x ncovs) matrix containing the covariate values.
#' @param Y (npersons x nitems) matrix of answers of all persons to all items.
#'
#' @return a list containing the values of the score function
#' \item{score_betas}{Values of the score function for every intercept.}
#' \item{score_deltas}{Values of the score function for every delta-parameter (ncovs x nitems)}
#' @export
#'

score_func_sample <- function(beta_intercepts, delta, covariates, Y) {

  betas_per_person <- calc_itemvecs(betas = beta_intercepts,
                                    delta = delta,
                                    covariates = covariates)

  score_func_per_person <- lapply(1:nrow(Y), function(x) score_func_person(betas_per_person[x,],
                                                                           Y[x,],
                                                                           covariates[x,]))
  score_betas_per_person <- lapply(score_func_per_person, `[[`, 1)
  score_deltas_per_person <- lapply(score_func_per_person, `[[`, 2)

  res_score <- list(score_betas = Reduce("+", score_betas_per_person),
                    score_deltas = t(Reduce("+", score_deltas_per_person)) )

  return(res_score)

}

##################################
##
## function:
## name: calc_lambda_max_cml

#' Calculate Lambda max for a sequence
#'
#' Calculates the maximum Lambda-value for the sequence (the minimum Lambda-value, for which all parameters are
#' equal to zero). For L1, this is the maximum-norm of the score function for the full modell evaluated at the
#' cml-estimates of the restricted model (hence all delta-parameters equal to zero). For Group-lasso this is the
#' maximum-norm of a vector containing the Euclidian norms per parameter group of the score vector evaluated
#' at the cml-estimates of the restricted model
#'
#' @param beta_intercepts beta_intercepts from restricted model, eg Rasch model without DIF
#' @param covariates Matrix of covariate values of all persons (n x num. of covs)
#' @param Y matrix of responses

#'
#' @return A list with the maximal Lambda-value for the sequence of Lambdas to be used for
#' optimization with three elements: selection for items, covariates and all interactions.
#' Note that for grouplasso this value is actually lambda_max * sqrt(groupsize), the groupsize
#' being either number of items or covariates, depending on the selection of penalization model.
#' The product is returned for these because group sizes are equal for all groups, so this
#' simplification is justified.
#'
#' @export
#'

calc_lambda_max_cml <- function(beta_intercepts,
                                covariates ,
                                Y){
  n = dim(Y)[1]
  I = dim(Y)[2]
  ncovs = dim(covariates)[2]

  scorevals <- score_func_sample(beta_intercepts = beta_intercepts,
                                 delta = matrix(0, ncovs, I),
                                 covariates = covariates,
                                 Y = Y)

  lambda_max = list(items = max(sqrt(colSums(scorevals$score_deltas^2))),
                    variables = max(sqrt(rowSums(scorevals$score_deltas^2))),
                    all.interactions = max(abs(scorevals$score_deltas)))
  lambda_max
}

##################################
##
## function:
## name: loss_function
## does: calculates logclik with or without penalization for DIF-model
##       given data, parameter values for items, DIF per variable per item
##
## input: data = matrix of responses
##        beta_intercept = values for beta per item
##        delta_matrix = DIF-parameters per item (columns) per covariate (rows)
##        DIFmodel = should DIF be calculated? options: "items", "variables", "all.interactions"
##                   If FALSE a Rasch Model is used
##        lambda = strength of penalization
##
## output: the value of the (penalized) log conditional Likelihood
##
###################################



#' loss-function with penalization
#'
#' calculates logclik with or without penalization for DIF-model
#' given data, parameter values for items, DIF-parameters per variable per item
#' and value of lambda
#'
#' @param data Answers of all persons (0/1, n x I)
#' @param covariates Matrix of covariate values of all persons (n x ncovs)
#' @param beta_intercept vecotr of parameter values for intercepts per item
#' @param delta_matrix Matrix (ncovs x nitems) of DIF-parameters per variable and item
#' @param DIFmodel can be "items", "variables", "all.interactions", "none".
#' And FALSE for Rasch model estimation
#' @param restr_pen For additional penalization to ensure zero sum of intercepts
#' @param lambda penalization parameter
#' @param lambda_restr penalization parameter for restriction
#'
#' @return the value of the (penalized) log conditional Likelihood
#' @export

loss_function <- function(data,
                          covariates,
                          beta_intercept = rep(0, ncol(data)),
                          delta_matrix = matrix(0, ncol(covariates), ncol(data)),
                          DIFmodel = "RM",
                          restr_pen = FALSE,
                          lambda = .5,
                          lambda_restr = 100)   {

  if (DIFmodel == FALSE) {betaperperson <- t(matrix(beta_intercept,
                                                    length(beta_intercept),
                                                    nrow(data)))}


  # if (DIFmodel %in% c("RM","none")) {betaperperson <- t(matrix(beta_intercept,
  #                                                      length(beta_intercept),
  #                                                      nrow(data)))}

  else {betaperperson <- calc_itemvecs(betas = beta_intercept,
                                       delta = delta_matrix,
                                       covariates = covariates)}

  data_log_c_lik <- cl_sample_2(Y = data, betas = betaperperson) -
    penalization(deltas = delta_matrix,
                 beta_intercept = beta_intercept,
                 DIF_type = DIFmodel,
                 restr_pen = restr_pen,
                 lambda = lambda,
                 lambda_restr = lambda_restr)

  return(data_log_c_lik)
}


##################################
##
## function:
## name: RM_can_nlm RM_can_optim
## does: RM_optimization via loss_function and cl_sample_2
##
## input: starting_values: starting values for optimization
##        data: matrix of response data
##        optimization: either "optim" or "nlm"
##
## output: optimized result of nlm()/otpim()-optimization
##
###################################


#' Rasch Model optimization with own loss
#'
#' RM_optimization via loss_function and cl_sample_2, no penalization.
#'
#' @param starting_values starting values for optimization
#' @param data matrix of response data
#' @param restr_pen To conduct penalization for sum=0 restriction set this to TRUE.
#' @param optimization either "optim" or "nlm"
#'
#' @return
#' @export
#'

RM_can_nlm <- function(starting_values, data, restr_pen = FALSE, optimization = "nlm") {

  res <- ("optimization musst be nlm or optim")

  if (optimization == "nlm") {
    res <-   nlm(
      f = function(x) -loss_function(data = data, beta_intercept = x, restr_pen = restr_pen),
      p = starting_values
    )
  }

  if (optimization == "optim") {
    res <- optim(
      starting_values,
      fn = function(x) -loss_function(data = data, beta_intercept = x, restr_pen = restr_pen)
    )
  }

  return(res)

}

##################################
##
## function:
## name: simulate_DIF_data


#' Simulating DIF-data
#'
#' Draws binary data randomly for given item intercepts, DIF-parameters,
#' covariate values and person parameters. (Based on easiness parameters.)
#'
#'
#' @param betas beta_intercepts for all items (vector)
#' @param delta values for DIF-parameters as matrix per item (cols) per covariate (rows)
#' @param covariates Matrix of covariate values of all persons (n x ncovs)
#' @param thetas vector of person parameters
#'
#' @return randomly drawn response matrix for DIF-model and given parameters
#' @export
#'

simulate_DIF_data <- function(betas, delta, covariates, thetas = 0) {

  n <- nrow(covariates)
  i <- length(betas)

  linpred <- matrix(thetas,n,i) + calc_itemvecs(betas = betas,
                                                delta = delta,
                                                covariates = covariates)

  probs <- exp(linpred)/(1+ exp(linpred))
  simdata <- 1* (runif(n*i,0,1) < probs )
  return(simdata)
}

##################################
##
## function:
## name: sim_scenario_DIF_data


#' Simulating DIF-data outputting scenario
#'
#' Draws binary data randomly for given item intercepts, DIF-parameters,
#' covariate values and person parameters. Returns a list with the simulated
#' data and all the parameters used to draw the sample.
#'
#' @param npersons number of persons
#' @param nitems number of items
#' @param ncovs number of covariates
#' @param seed seed to be set for simulation
#' @param betavalues beta_intercepts for all items (vector)
#' @param deltavalues values for DIF-parameters as matrix per item (cols) per covariate (rows)
#' @param thetavalues vector of person parameters
#' @param covariatevalues Matrix of covariate values of all persons (n x ncovs)

#'
#' @return randomly drawn response matrix for DIF-model and given parameters
#' @export
#'

sim_scenario_DIF_data <- function(npersons, nitems, ncovs, seed,
                                  betavalues = rep(0,nitems),
                                  deltavalues = matrix(c(-.5, 0, 0, 0, .5) ,ncovs ,nitems ),
                                  thetavalues = rep(0,npersons),
                                  covariatevalues = matrix(rnorm(npersons*ncovs), npersons, ncovs)) {
  set.seed(seed)
  return(list(data = simulate_DIF_data(betavalues, deltavalues, covariatevalues, thetavalues),
              betavalues = betavalues,
              deltavalues = deltavalues,
              thetavalues = thetavalues,
              covariatevalues = covariatevalues,
              seed = seed))
}

##################################
##
## function:
## name: df_pen_unequal_zero, df_veclength
## does: calculates the first and second parts of df_delta depending on DIF-type
##
## input: deltas:      values for DIF-parameters as matrix for certain Lambda
##                     per item (columns) per covariate (rows)
##        DIF_type:    for which Lasso-term should the non-nil values be calculated?
##                     options: "items", "variables", "all.interactions"
##        delta_unres: values for DIF-parameters as matrix form unrestricted cml solution (lambda = 0)
##                     per item (columns) per covariate (rows)
##
## output: the first and second part of df_delta; remark: 10^(-6) will be zero
##         df_pen_unequal_zero: either the # of delta parameters != 0 or
##                              the # of items or variables with delta parameters != 0
##                              delivering the corresponding part for the calculation of df_pen
##         df_veclength:        either the eucl. length of the deltavector per item or covariate
##                              divided by the corresponding length of the unrestricted cml deltas
##                              or 0 in case of all.interactions.
##
###################################

##################################
##
#' df calculation for RM-DIF-model with penalization
#'
#' calculates the first and second parts of df_delta depending on DIF-type given
#' parameter estimates for delta under a certain lambda and the unrestricted DIF-model.
#'
#' Title
#'
#' @param deltas values for DIF-parameters as matrix for certain Lambda
#' per item (columns) per covariate (rows)
#' @param DIF_type for which Lasso-term should the non-nil values be calculated?
#' options: "items", "variables", "all.interactions"
#' @param delta_unres values for DIF-parameters as matrix form unrestricted cml solution (lambda = 0)
##                     per item (columns) per covariate (rows)
#' @param tolerance Value under which parameters are considered equal
#' to zero for the calculation of degreees of freedom
#'
#' @return calculated first (df_pen_unequal_zero) and second part (df_veclength)
#' of df_delta, or already the complete degrees of freedom (df_lambda) for a
#' lasso penalized (cml) model according to Yuan and Lin(2006)
#'  and Zou et al. (2007); remark: 10^(-6) will be zero (default value of argument tolerance)
#'
#' @name degrees-of-freedom-helpers
NULL

##################################
##
## function:
## name: df_pen_unequal_zero df_veclength

#' @rdname degrees-of-freedom-helpers
#' @export

df_pen_unequal_zero <- function(deltas, DIF_type, tolerance = 0) {

  true_false <- apply(deltas, 1:2, function(x) all.equal(x,0,tolerance = tolerance) ) == "TRUE"

  dfs <- NA

  if (DIF_type == "items") {
    dfs = sum(colSums(true_false) != 0)
  }

  if (DIF_type == "variables") {
    dfs = sum(rowSums(true_false) != 0)
  }

  if (DIF_type == "all.interactions") {
    dfs =  sum(true_false == FALSE)
  }

  if (DIF_type == "none") {
    dfs =  0
  }


  return(dfs)
}

##################################
##
## function:
## name: df_veclength

#' @rdname degrees-of-freedom-helpers
#' @export

df_veclength <- function(deltas, DIF_type, delta_unres) {

  nitems <- ncol(deltas)
  ncovs <- nrow(deltas)

  dfs <- NA

  if (DIF_type == "items") {
    dfs = sum(sqrt(colSums(deltas^2)) / sqrt(colSums(delta_unres^2)))  * (ncovs - 1)
  }

  if (DIF_type == "variables") {
    dfs = sum(sqrt(rowSums(deltas^2)) / sqrt(rowSums(delta_unres^2)))  * (nitems - 1)
  }

  if (DIF_type == "all.interactions") {
    dfs = 0
  }

  if (DIF_type == "none") {
    dfs = 0
  }

  return(dfs)
}



##################################
##
## function:
## name: df_lambda
## does: calculates the degrees of freedom of a lasso-penalized model
##
## input: deltas:      values for DIF-parameters as matrix for certain Lambda
##                     per item (columns) per covariate (rows)
##        DIF_type:    for which Lasso-term should the non-nil values be calculated?
##                     options: "items", "variables", "all.interactions"
##        delta_unres: values for DIF-parameters as matrix form unrestricted cml solution (lambda = 0)
##                     per item (columns) per covariate (rows)
##
## output: sum of df_items, df_df_pen_unequal_zero and df_veclength
##         delivers the degrees of freedom for a lasso penalized (cml) model according to Yuan and Lin(2006)
##         and Zou et al. (2007)
##
###################################

#' @rdname degrees-of-freedom-helpers
#' @export

df_lambda <-  function(deltas, DIF_type, delta_unres, tolerance) {

  dfs = ncol(deltas) - 1 +
    df_pen_unequal_zero(deltas, DIF_type, tolerance) +
    df_veclength(deltas, DIF_type, delta_unres)

  return(dfs)
}

##################################
#
# BIC
#

#' BIC Lasso
#'
#' Calculates the BIC for an L1 penalized model. Needed for picking a
#' hyper parameter lambda.
#'
#' @param Y data
#' @param betas vecotr of intercepts
#' @param deltas (ncov x nitems) matrix of DIF-parameters
#' @param covariates (n x ncovs) matrix of covariate values
#' @param DIF_type can be "items", "variables", "all.interactions", "none"
#' @param delta_unres DIF-parameter estimates from unrestricted estimation
#' @param tolerance Value under which parameters are considered equal
#' to zero for the calculation of degreees of freedom
#'
#' @return Returns the BIC for an estimated DIF model with (or without) penalization term
#' @export
#'

BIC_lasso <- function(Y, betas, deltas, covariates, DIF_type, delta_unres, tolerance) {

  BIC <- -2 * cl_sample_2(Y = Y, betas = calc_itemvecs(betas, deltas, covariates)) +
    df_lambda(deltas, DIF_type, delta_unres, tolerance) * log(nrow(Y) * ncol(Y))

  return(BIC)
}




##################################
##
## function:
## name: RM_DIF_nlm

#' Rasch model with DIF via nlm
#'
#' optimizes loss_function with DIF-part for certain tuning-parameter lambda
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x num. of covs)
#' @param DIF_type should DIF be calculated? options: "items", "variables", "all.interactions"
#' @param lambda strength of penalization
#' @param starting_betas beta_intercepts starting values for all items (vector)
#' @param starting_deltas  starting values for DIF-parameters as matrix (ncovs x nitems)
#' @param restr_pen For additional penalization to ensure zero sum of intercepts
#' @param lambda_restr penalization parameter for restriction
#' @param ... other arguments to be passed
#'
#' @return result of nlm-optimization for DIF-model via cML-loss-function for
#' a certain penalisation-type and lambda value
#' @export
#'

RM_DIF_nlm <- function(data, covariates, DIF_type, lambda,
                       starting_betas,
                       starting_deltas,
                       restr_pen = TRUE,
                       lambda_restr = 100,
                       ...) {

  nitems <- ncol(data)
  ncov <- ncol(covariates)

  starting_values <- c(starting_betas, starting_deltas)

  nlm_results <- nlm(
    f = function(x) -loss_function(data = data,
                                   covariates = covariates,
                                   beta_intercept = x[1:nitems],
                                   delta_matrix = matrix(x[(nitems+1):(nitems+ncov*nitems)],ncov),
                                   DIFmodel = DIF_type,
                                   lambda = lambda,
                                   restr_pen = restr_pen,
                                   lambda_restr = lambda_restr),
    p = starting_values, ...=...
  )

  beta_estimates = nlm_results$estimate[1:nitems]
  delta_estimates = matrix(nlm_results$estimate[(nitems+1):(nitems+ncov*nitems)],ncov,
                           dimnames = list(covariates = paste("covariate",1:ncov, sep = "_"),
                                           items = paste("item",1:nitems, sep = "_")))


  nlm_results$estimatevalues <- list(betas = beta_estimates,
                                     deltas = delta_estimates  #,BIC = BIC_model
  )
  nlm_results$inputs <- list(data = data,
                             covariates = covariates,
                             DIFmodel = DIF_type,
                             lambda = lambda)

  return(nlm_results)
}


##################################
##
## function:
## name: RM_DIF_optim

#' Rasch model with DIF via nlm
#'
#' optimizes loss_function with DIF-part for certain tuning-parameter
#' lambda with optim
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x num. of covs)
#' @param DIF_type should DIF be calculated? options: "items", "variables", "all.interactions"
#' @param lambda strength of penalization
#' @param starting_betas beta_intercepts starting values for all items (vector)
#' @param starting_deltas  starting values for DIF-parameters as matrix (ncovs x nitems)
#' @param restr_pen For additional penalization to ensure zero sum of intercepts
#' @param lambda_restr penalization parameter for restriction
#' @param method optimization method to be past to optim()
#' @param ... other arguments to be passed
#'
#' @return result of optim-optimization for DIF-model via cML-loss-function for
#' a certain penalisation-type and lambda value
#' @export
#'

RM_DIF_optim <- function(data, covariates, DIF_type, lambda,
                         starting_betas,
                         starting_deltas,
                         restr_pen = TRUE,
                         lambda_restr = 100,
                         method = "Nelder-Mead",
                         ...) {

  nitems <- ncol(data)
  ncov <- ncol(covariates)

  starting_values <- c(starting_betas, starting_deltas)

  optim_results <- optim(
    par  = starting_values,
    fn = function(x) -loss_function(data = data,
                                    covariates = covariates,
                                    beta_intercept = x[1:nitems],
                                    delta_matrix = matrix(x[(nitems+1):(nitems+ncov*nitems)],ncov),
                                    DIFmodel = DIF_type,
                                    lambda = lambda,
                                    restr_pen = restr_pen,
                                    lambda_restr = lambda_restr),
    method = method,
    ...=...
  )

  beta_estimates = optim_results$par[1:nitems]
  delta_estimates = matrix(optim_results$par[(nitems+1):(nitems+ncov*nitems)],ncov,
                           dimnames = list(covariates = paste("covariate",1:ncov, sep = "_"),
                                           items = paste("item",1:nitems, sep = "_")))


  optim_results$estimatevalues <- list(betas = beta_estimates,
                                       deltas = delta_estimates  #,BIC = BIC_model
  )
  optim_results$inputs <- list(data = data,
                               covariates = covariates,
                               DIFmodel = DIF_type,
                               lambda = lambda,
                               method = method)

  return(optim_results)
}


##################################
##
## function:
## name: RM_DIF_lbfgs

#' Rasch model with DIF via lbfgs (only L1)
#'
#' optimizes loss_function with DIF-part for certain tuning-parameter
#' lambda with lbfgs
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x num. of covs)
#' @param starting_betas beta_intercepts starting values for all items (vector)
#' @param starting_deltas starting values for DIF-parameters as matrix (ncovs x nitems)
#' @param lambda penalization parameter
#' @param ... other arguments to be passed to lbfgs
#'
#' @return result of optim-optimization for DIF-model via cML-loss-function for
#' a certain penalisation-type (for now only all.interactions) and lambda value
#' @export
#'

RM_DIF_lbfgs <- function(data, covariates,  ## DIF_type, (evtl. später wieder einbauen)
                         starting_betas = rep(0,ncol(data)),
                         starting_deltas = matrix(rep(0,ncol(data)*ncol(covariates)), ncol(covariates), ncol(data)),
                         lambda = 0,
                         ...) {

  nitems <- ncol(data)
  ncov <- ncol(covariates)

  starting_values <- c(starting_betas, starting_deltas)

  objective <- function(x) -cl_sample_2(Y = data,
                                        betas = calc_itemvecs(betas = x[1:nitems],
                                                              delta = matrix(x[(nitems+1):length(x)],ncov),
                                                              covariates = covariates))
  gradient <- function(x) -unlist(score_func_sample(beta_intercepts = x[1:nitems],
                                                    delta = matrix(x[(nitems+1):length(x)],ncov),
                                                    covariates = covariates,
                                                    Y = data))

  lbfgs_results <- lbfgs(
    call_eval = objective,
    call_grad = gradient,
    starting_values,
    orthantwise_c = lambda,
    orthantwise_start = ncol(data),
    ...
  )

  return(list(beta_estimates = lbfgs_results$par[1:nitems] - mean(lbfgs_results$par[1:nitems]),
              delta_estimates = matrix(lbfgs_results$par[(nitems+1):length(lbfgs_results$par)],
                                       ncov,
                                       nitems),
              lambda = lambda,
              lbfgs_res = lbfgs_results))

}

##################################
##
## function:
## name: cml_lasso_win_lbfgs_optimized


#' Title
#' Fitting cmlDIF with L1-penalty via lbfgs (opitmized with lambda_max)
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x ncovs)
#' @param DIF_type can be "items", "variables", "all.interactions", "none"
#' @param nlambdas number of lambda parameters evenly spread between 10^-5 and
#' lambda_max, which is calculated automatically
#' @param tolerance_df Value under which parameters are considered equal
#' to zero for the calculation of degreees of freedom
#'
#' @return list of lbfgs-optimization results for each lambda, the picked model
#' and its parameters, the computing time
#' @export
#'

cml_lasso_win_lbfgs_optimized <- function (data, covariates,
                                           DIF_type = "all.interactions",
                                           nlambdas = 21,
                                           tolerance_df = 0,
                                           ...){
  npersons <- nrow(data)
  nitems <- ncol(data)
  ncov <- ncol(covariates)

  time_unres <- system.time(res_unrestricted <- RM_DIF_lbfgs(data = data,
                                                             covariates = covariates,
                                                             lambda = 0,
                                                             invisible = 1)
  )

  time_res <- system.time(beta_estimates_restricted <- -itempar(raschmodel(data)))

  lambda_max <- calc_lambda_max_cml(beta_intercepts = beta_estimates_restricted,
                                    covariates = covariates,
                                    Y = data)[[DIF_type]]

  lambdas <- seq(10^-5, lambda_max, length.out = nlambdas-1)
  lambdas <- c(lambdas, lambdas[nlambdas-1] +(lambdas[nlambdas-1] - lambdas[nlambdas-2])) # add one lambda to surely include restricted result

  time_lbfgs <- system.time(results <- lapply(1:length(lambdas),function(x) RM_DIF_lbfgs(data = data,
                                                                                         covariates = covariates,
                                                                                         lambda = lambdas[x],
                                                                                         invisible = 1,
                                                                                         ...)))

  # rows: covariates, cols: items, "slices": different lambdas
  delta_est <- array(unlist(lapply(1:length(lambdas), function(x)results[[x]]$delta_estimates)) , c(ncov,nitems,length(lambdas)))
  # cols: items, rows: different lambdas
  beta_est <- t(sapply(1:length(lambdas), function(x)results[[x]]$beta_estimates))

  BIC_list <- NULL
  for(j in 1:length(lambdas)) {
    BIC_list[j] <- results[[j]]$estimatevalues$BIC <-
      BIC_lasso(Y = data,
                betas = results[[j]]$beta_estimates,
                deltas = results[[j]]$delta_estimates,
                covariates = covariates,
                DIF_type = DIF_type,
                delta_unres = res_unrestricted$delta_estimates,
                tolerance = 0)
  }


  # pick model - if several minimum values (usually the first few) exist and , then pick the highest position,
  model_picked <- ifelse(1 %in% which(BIC_list == min(BIC_list)), max(which(BIC_list == min(BIC_list))), which.min(BIC_list))
  res_picked <- results[[model_picked]]

  # collect results into list
  result <- list(results_per_lambda = results,
                 BICs = BIC_list,
                 beta_estimates = beta_est,
                 delta_estimates = delta_est,
                 lambdas = lambdas,
                 times = list(time_unrestricted_model = time_unres,
                              time_penalized_part = time_lbfgs,
                              time_restricted_model = time_res),
                 model_picked = model_picked,
                 lambda_picked = lambdas[model_picked],
                 result_picked = res_picked,
                 result_unrestricted = res_unrestricted,
                 beta_estimates_RM = beta_estimates_restricted)

}




##################################
##
## function:
## name: cml_lasso_win_lbfgs


#' Title
#' cmlDIFlasso calculated with lbfgs
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x ncovs)
#' @param DIF_type can be "items", "variables", "all.interactions", "none"
#' @param lambdas vecotr of penalization parameters to be used
#' @param tolerance_df Value under which parameters are considered equal
#' to zero for the calculation of degreees of freedom
#'
#' @return list of lbfgs-optimization results for each lambda, the picked model
#' and its parameters, the computing time
#' @export
#'

cml_lasso_win_lbfgs <- function (data, covariates,
                                 DIF_type = "all.interactions",
                                 lambdas = 10^(-(-2:7)),
                                 tolerance_df = 0,
                                 ...){
  npersons <- nrow(data)
  nitems <- ncol(data)
  ncov <- ncol(covariates)
  time_unres <- system.time(res_unrestricted <- RM_DIF_lbfgs(data = data,
                                                             covariates = covariates,
                                                             lambda = 0,
                                                             invisible = 1)
  )

  time_lbfgs <- system.time(results <- lapply(1:length(lambdas),function(x) RM_DIF_lbfgs(data = data,
                                                                                         covariates = covariates,
                                                                                         lambda = lambdas[x],
                                                                                         invisible = 1,
                                                                                         ...)))

  # rows: covariates, cols: items, "slices": different lambdas
  delta_est <- array(unlist(lapply(1:length(lambdas), function(x)results[[x]]$delta_estimates)) , c(ncov,nitems,length(lambdas)))
  # cols: items, rows: different lambdas
  beta_est <- t(sapply(1:length(lambdas), function(x)results[[x]]$beta_estimates))

  BIC_list <- NULL
  for(j in 1:length(lambdas)) {
    BIC_list[j] <- results[[j]]$estimatevalues$BIC <-
      BIC_lasso(Y = data,
                betas = results[[j]]$beta_estimates,
                deltas = results[[j]]$delta_estimates,
                covariates = covariates,
                DIF_type = DIF_type,
                delta_unres = res_unrestricted$delta_estimates,
                tolerance = 0)
  }


  # pick model - if several minimum values (usually the first few) exist and , then pick the highest position,
  model_picked <- ifelse(1 %in% which(BIC_list == min(BIC_list)), max(which(BIC_list == min(BIC_list))), which.min(BIC_list))
  res_picked <- results[[model_picked]]

  # collect results into list
  result <- list(results_per_lambda = results,
                 BICs = BIC_list,
                 beta_estimates = beta_est,
                 delta_estimates = delta_est,
                 lambdas = lambdas,
                 times = list(time_unrestricted_model = time_unres,
                              time_penalized_part = time_lbfgs),
                 model_picked = model_picked,
                 lambda_picked = lambdas[model_picked],
                 result_picked = res_picked)

}




####
#
# no parallelisation:
#
####

##################################
##
## function:
## name: cml_lasso_win_lbfgs_optimized


#' Title
#' Fitting cmlDIF with L1-penalty via lbfgs (opitmized with lambda_max)
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x ncovs)
#' @param DIF_type can be "items", "variables", "all.interactions", "none"
#' @param nlambdas number of lambda parameters evenly spread between 10^-5 and
#' lambda_max, which is calculated automatically
#' @param tolerance_df Value under which parameters are considered equal
#' to zero for the calculation of degreees of freedom
#'
#' @return list of lbfgs-optimization results for each lambda, the picked model
#' and its parameters, the computing time
#' @export
#'

cml_lasso_win_lbfgs_optimized <- function (data, covariates,
                                           DIF_type = "all.interactions",
                                           nlambdas = 21,
                                           tolerance_df = 0,
                                           ...){
  npersons <- nrow(data)
  nitems <- ncol(data)
  ncov <- ncol(covariates)

  time_unres <- system.time(res_unrestricted <- RM_DIF_lbfgs(data = data,
                                                             covariates = covariates,
                                                             lambda = 0,
                                                             invisible = 1)
  )

  time_res <- system.time(beta_estimates_restricted <- -itempar(raschmodel(data)))

  lambda_max <- calc_lambda_max_cml(beta_intercepts = beta_estimates_restricted,
                                    covariates = covariates,
                                    Y = data)[[DIF_type]]

  lambdas <- seq(10^-5, lambda_max, length.out = nlambdas-1)
  lambdas <- c(lambdas, lambdas[nlambdas-1] +(lambdas[nlambdas-1] - lambdas[nlambdas-2])) # add one lambda to surely include restricted result

  time_lbfgs <- system.time(results <- lapply(1:length(lambdas),function(x) RM_DIF_lbfgs(data = data,
                                                                                         covariates = covariates,
                                                                                         lambda = lambdas[x],
                                                                                         invisible = 1,
                                                                                         starting_betas = res_unrestricted$beta_estimates,
                                                                                         starting_deltas = res_unrestricted$delta_estimates,
                                                                                         ...)))

  # rows: covariates, cols: items, "slices": different lambdas
  delta_est <- array(unlist(lapply(1:length(lambdas), function(x)results[[x]]$delta_estimates)) , c(ncov,nitems,length(lambdas)))
  # cols: items, rows: different lambdas
  beta_est <- t(sapply(1:length(lambdas), function(x)results[[x]]$beta_estimates))

  BIC_list <- NULL
  for(j in 1:length(lambdas)) {
    BIC_list[j] <- results[[j]]$estimatevalues$BIC <-
      BIC_lasso(Y = data,
                betas = results[[j]]$beta_estimates,
                deltas = results[[j]]$delta_estimates,
                covariates = covariates,
                DIF_type = DIF_type,
                delta_unres = res_unrestricted$delta_estimates,
                tolerance = 0)
  }


  # pick model - if several minimum values (usually the first few) exist and , then pick the highest position,
  model_picked <- ifelse(1 %in% which(BIC_list == min(BIC_list)), max(which(BIC_list == min(BIC_list))), which.min(BIC_list))
  res_picked <- results[[model_picked]]

  # collect results into list
  result <- list(results_per_lambda = results,
                 BICs = BIC_list,
                 beta_estimates = beta_est,
                 delta_estimates = delta_est,
                 lambdas = lambdas,
                 times = list(time_unrestricted_model = time_unres,
                              time_penalized_part = time_lbfgs,
                              time_restricted_model = time_res),
                 model_picked = model_picked,
                 lambda_picked = lambdas[model_picked],
                 result_picked = res_picked,
                 result_unrestricted = res_unrestricted,
                 beta_estimates_RM = beta_estimates_restricted)

}



#####
#
# parallelisation foreach
#
#####


##################################
##
## function:
## name: cml_lasso_win_lbfgs_optimized_foreach


#' Title
#' Fitting cmlDIF with L1-penalty via lbfgs (opitmized with lambda_max)
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x ncovs)
#' @param DIF_type can be "items", "variables", "all.interactions", "none"
#' @param nlambdas number of lambda parameters evenly spread between 10^-5 and
#' lambda_max, which is calculated automatically
#' @param tolerance_df Value under which parameters are considered equal
#' to zero for the calculation of degreees of freedom
#' @param no_of_cores Number of cores to be used for parallel computing
#'
#' @return list of lbfgs-optimization results for each lambda, the picked model
#' and its parameters, the computing time
#' @export
#'

cml_lasso_win_lbfgs_optimized_foreach <- function (data, covariates,
                                                   DIF_type = "all.interactions",
                                                   nlambdas = 21,
                                                   tolerance_df = 0,
                                                   no_of_cores_to_use = ifelse(detectCores()==1,1,detectCores()-1),
                                                   ...){
  npersons <- nrow(data)
  nitems <- ncol(data)
  ncov <- ncol(covariates)

  time_unres <- system.time(res_unrestricted <- RM_DIF_lbfgs(data = data,
                                                             covariates = covariates,
                                                             lambda = 0,
                                                             invisible = 1)
  )

  time_res <- system.time(beta_estimates_restricted <- -itempar(raschmodel(data)))

  lambda_max <- calc_lambda_max_cml(beta_intercepts = beta_estimates_restricted,
                                    covariates = covariates,
                                    Y = data)[[DIF_type]]

  lambdas <- seq(10^-5, lambda_max, length.out = nlambdas-1)
  lambdas <- c(lambdas, lambdas[nlambdas-1] +(lambdas[nlambdas-1] - lambdas[nlambdas-2])) # add one lambda to surely include restricted result

  time_lbfgs <- system.time({

    cl = makeCluster(no_of_cores_to_use)
    registerDoParallel(cl)
    #    clusterEvalQ(cl, {library(cmlDIFlasso2)})

    results <- foreach(i = lambdas,
                       .packages = c("cmlDIFlasso2","psychotools","lbfgs")) %dopar% {
                         RM_DIF_lbfgs(data = data,
                                      covariates = covariates,
                                      lambda = i,
                                      invisible = 1,
                                      starting_betas = res_unrestricted$beta_estimates,
                                      starting_deltas = res_unrestricted$delta_estimates,
                                      ...)}
    stopCluster(cl)
  })

  # rows: covariates, cols: items, "slices": different lambdas
  delta_est <- array(unlist(lapply(1:length(lambdas), function(x)results[[x]]$delta_estimates)) , c(ncov,nitems,length(lambdas)))
  # cols: items, rows: different lambdas
  beta_est <- t(sapply(1:length(lambdas), function(x)results[[x]]$beta_estimates))

  BIC_list <- NULL
  for(j in 1:length(lambdas)) {
    BIC_list[j] <- results[[j]]$estimatevalues$BIC <-
      BIC_lasso(Y = data,
                betas = results[[j]]$beta_estimates,
                deltas = results[[j]]$delta_estimates,
                covariates = covariates,
                DIF_type = DIF_type,
                delta_unres = res_unrestricted$delta_estimates,
                tolerance = 0)
  }


  # pick model - if several minimum values (usually the first few) exist and , then pick the highest position,
  model_picked <- ifelse(1 %in% which(BIC_list == min(BIC_list)), max(which(BIC_list == min(BIC_list))), which.min(BIC_list))
  res_picked <- results[[model_picked]]

  # collect results into list
  result <- list(results_per_lambda = results,
                 BICs = BIC_list,
                 beta_estimates = beta_est,
                 delta_estimates = delta_est,
                 lambdas = lambdas,
                 times = list(time_unrestricted_model = time_unres,
                              time_penalized_part = time_lbfgs,
                              time_restricted_model = time_res),
                 model_picked = model_picked,
                 lambda_picked = lambdas[model_picked],
                 result_picked = res_picked,
                 result_unrestricted = res_unrestricted,
                 beta_estimates_RM = beta_estimates_restricted)

}



##################################
##
## function:
## name: cml_lasso_win_lbfgs_optimized_parLapply


#' Title
#' Fitting cmlDIF with L1-penalty via lbfgs (opitmized with lambda_max)
#'
#' @param data matrix of responses
#' @param covariates Matrix of covariate values of all persons (n x ncovs)
#' @param DIF_type can be "items", "variables", "all.interactions", "none"
#' @param nlambdas number of lambda parameters evenly spread between 10^-5 and
#' lambda_max, which is calculated automatically
#' @param tolerance_df Value under which parameters are considered equal
#' to zero for the calculation of degreees of freedom
#' @param no_of_cores Number of cores to be used for parallel computing
#'
#' @return list of lbfgs-optimization results for each lambda, the picked model
#' and its parameters, the computing time
#' @export
#'

cml_lasso_win_lbfgs_optimized_parLapply <- function (data, covariates,
                                                     DIF_type = "all.interactions",
                                                     nlambdas = 21,
                                                     tolerance_df = 0,
                                                     no_of_cores_to_use = ifelse(detectCores()==1,1,detectCores()-1),
                                                     ...){
  npersons <- nrow(data)
  nitems <- ncol(data)
  ncov <- ncol(covariates)

  time_unres <- system.time(res_unrestricted <- RM_DIF_lbfgs(data = data,
                                                             covariates = covariates,
                                                             lambda = 0,
                                                             invisible = 1)
  )

  time_res <- system.time(beta_estimates_restricted <- -itempar(raschmodel(data)))

  lambda_max <- calc_lambda_max_cml(beta_intercepts = beta_estimates_restricted,
                                    covariates = covariates,
                                    Y = data)[[DIF_type]]

  lambdas <- seq(10^-5, lambda_max, length.out = nlambdas-1)
  lambdas <- c(lambdas, lambdas[nlambdas-1] +(lambdas[nlambdas-1] - lambdas[nlambdas-2])) # add one lambda to surely include restricted result

  time_lbfgs <- system.time({
    cl <- makeCluster(no_of_cores_to_use)
    clusterEvalQ(cl, {library(cmlDIFlasso2)})

    results <- parLapplyLB(cl,1:length(lambdas),function(x) RM_DIF_lbfgs(data = data,
                                                                         covariates = covariates,
                                                                         lambda = lambdas[x],
                                                                         starting_betas = res_unrestricted$beta_estimates,
                                                                         starting_deltas = res_unrestricted$delta_estimates,
                                                                         invisible = 1,
                                                                         ...))
  })
  stopCluster(cl)

  # rows: covariates, cols: items, "slices": different lambdas
  delta_est <- array(unlist(lapply(1:length(lambdas), function(x)results[[x]]$delta_estimates)) , c(ncov,nitems,length(lambdas)))
  # cols: items, rows: different lambdas
  beta_est <- t(sapply(1:length(lambdas), function(x)results[[x]]$beta_estimates))

  BIC_list <- NULL
  for(j in 1:length(lambdas)) {
    BIC_list[j] <- results[[j]]$estimatevalues$BIC <-
      BIC_lasso(Y = data,
                betas = results[[j]]$beta_estimates,
                deltas = results[[j]]$delta_estimates,
                covariates = covariates,
                DIF_type = DIF_type,
                delta_unres = res_unrestricted$delta_estimates,
                tolerance = 0)
  }

  # pick model - if several minimum values (usually the first few) exist and , then pick the highest position,
  model_picked <- ifelse(1 %in% which(BIC_list == min(BIC_list)), max(which(BIC_list == min(BIC_list))), which.min(BIC_list))
  res_picked <- results[[model_picked]]

  # collect results into list
  result <- list(results_per_lambda = results,
                 BICs = BIC_list,
                 beta_estimates = beta_est,
                 delta_estimates = delta_est,
                 lambdas = lambdas,
                 times = list(time_unrestricted_model = time_unres,
                              time_penalized_part = time_lbfgs,
                              time_restricted_model = time_res),
                 model_picked = model_picked,
                 lambda_picked = lambdas[model_picked],
                 result_picked = res_picked,
                 result_unrestricted = res_unrestricted,
                 beta_estimates_RM = beta_estimates_restricted)

}



##################################
##
## function:
## name: session_all


#' Title
#' Session information for simulation and estimation
#'
#' @param no_of_cores_used number of cores used for anlysis.
#' To be completely correct, this should be the number of threads used.
#'
#' @return A list including the output of sessionInfo(),
#' the CPU name, the amount of system RAM,
#' the number of cores used (to be completely correct, this is the number of threads used),
#' the number of threads available and the number of physical cores available
#'
#' @export
#'


session_all <- function(no_of_cores_used = "unknown") {

  list(sessioninfo = sessionInfo(),
       CPU_used = benchmarkme::get_cpu(),
       RAM_available = benchmarkme::get_ram(),
       no_of_cores_used = no_of_cores_used,
       no_of_threads_available = parallel::detectCores(logical = TRUE),
       no_of_physical_cores_available = parallel::detectCores(logical = FALSE))
}

##################################
##
## function:
## name: wrap_cml_lasso_lbfgs_optimized


#' Title
#' Wrapper function for use in simulation studied, taking all parameters,
#' passing them to sim_scenario_DIF_data and computing results via
#' cml_lasso_win_lbfgs_optimized
#'
#' @param n number of persons (only used of if not specified below)
#' @param I number of items (only used of if not specified below)
#' @param nvars number of covariates (only used if not specified below)
#' @param seed seed for data simulation
#' @param betavalues Vector with true values of the item parameters.
#' @param deltavalues Matrix with the true values of the delta parameters (ncovs x nitems)
#' @param covariatevalues Matrix of covariate values of all persons (n x ncovs)
#' @param no_of_cores_used Number of cores to be used for parallel computing
#'
#' @return list of results and simulation scenario details.
#' @export
#'

wrap_cml_lasso_lbfgs_optimized <- function(n, I, nvars, seed, nlambdas,
                                           betavalues = rep(0,I),
                                           deltavalues = matrix(c(-.5, 0, 0, 0, .5), nvars, I ),
                                           thetavalues = rep(0, n),
                                           covariatevalues = matrix(rnorm(n * nvars), n, nvars),
                                           no_of_cores_used = "unknown"){

  scen_1 <- sim_scenario_DIF_data(n, I, nvars, seed,
                                  betavalues = betavalues,
                                  deltavalues = deltavalues,
                                  thetavalues = thetavalues,
                                  covariatevalues = covariatevalues)
  time_lbfgs <- system.time(result <- cml_lasso_win_lbfgs_optimized(data = scen_1$data,
                                                                    covariates = scen_1$covariatevalues,
                                                                    nlambdas = nlambdas))
  list(result = result,
       time_lbfgs = time_lbfgs,
       scenario = scen_1,
       lambdas = result$lambdas,
       session = session_all(no_of_cores_used = no_of_cores_used))
}

##################################
##
## function:
## name: extract_obs_alpha_beta_power


#' Title
#' Helper function to extract information from result list with
#' wrap_cml_lasso_lbfgs_optimized results
#'
#' @param x element for which to calculate
#' @param resultlist list of results from wrap_cml_lasso_lbfgs_optimized
#'
#' @return list of results and simulation scenario details.
#' @export
#'


extract_obs_alpha_beta_power <- function (x, resultlist) {
  tb1 <- table(abs(resultlist[[x]]$result$result_picked$delta_estimates)>.0000001,
               abs(resultlist[[x]]$scenario$deltavalues) != 0)

  tb <- matrix(0,2,2)
  if ("FALSE" %in% rownames(tb1) & "FALSE" %in% colnames(tb1)) {tb[1,1] <- tb1["FALSE","FALSE"]}
  if ("FALSE" %in% rownames(tb1) & "TRUE" %in% colnames(tb1)) {tb[1,2] <- tb1["FALSE","TRUE"]}
  if ("TRUE" %in% rownames(tb1) & "FALSE" %in% colnames(tb1)) {tb[2,1] <- tb1["TRUE","FALSE"]}
  if ("TRUE" %in% rownames(tb1) & "TRUE" %in% colnames(tb1)) {tb[2,2] <- tb1["TRUE","TRUE"]}

  list(obs.alpha = round(prop.table(tb,2)[2,1],3),
       obs.beta = round(prop.table(tb,2)[1,2],3),
       obs.pwr = round(prop.table(tb,2)[2,2],3))
}

##################################
##
## function:
## name: make_df_from_resultlist
#' Title
#' Helper function to extract information from result list with
#' wrap_cml_lasso_lbfgs_optimized results
#'
#' @param resultlist list of results from wrap_cml_lasso_lbfgs_optimized
#'
#' @return Dataframe with information per simulation result, including n, I, nvars, times,
#' observed alpha, power and others.
#'
#' @export
#'
make_df_from_resultlist <- function(resultlist){

  pars <- do.call(rbind,
                  lapply(1:length(resultlist),
                         function(x) cbind(n = nrow(resultlist[[x]]$scenario$data),
                                           I = ncol(resultlist[[x]]$scenario$deltavalues),
                                           nvars = nrow(resultlist[[x]]$scenario$deltavalues))))



  t_and_scen <- lapply(1:nrow(pars),FUN = function(x) unlist(list(time = resultlist[[x]]$time_lbfgs,
                                                                  scenario = pars[x,],
                                                                  extract_obs_alpha_beta_power(x,resultlist))))
  df_t_and_scen <- as.data.frame(do.call(rbind,t_and_scen))

}

##################################
##
## function:
## name:
## does: calculates bla
##
## input: x:     data
##        betas: bla
##
## output: blabla
##
###################################

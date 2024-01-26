## DATASETS
#' @title Data set toy example
#' @description A simulated data set containing 8 variables (5 responses variables and 3 covariates).
#' \itemize{
#'   \item \code{y1} numeric, gaussian response variable
#'   \item \code{y2} numeric, gaussian response variable
#'   \item \code{y3} numeric, gaussian response variable
#'   \item \code{Po1} numeric, response variable with counts simulated from a Poison distribution
#'   \item \code{Be1} factor with two categories, binary response variable
#'   \item \code{X1} continuous covariate
#'   \item \code{X2} continuous covariate
#'   \item \code{X3} categorical covariate with three levels.
#'}
#' @name data_toy
#' @docType data
#' @usage data(data_toy)
#' @format A data frame with 500 rows and 8 variables
NULL

## IMPORTS
#' @importFrom stats family gaussian make.link glm as.formula dpois dbinom dnorm dcauchy model.offset pnorm printCoefmat model.matrix model.weights nlminb setNames
#' @importFrom utils write.table combn
#' @importFrom Formula Formula as.Formula model.part
#' @importFrom statmod gauss.quad
#' @importFrom numDeriv hessian
#' @importFrom mvtnorm dmvnorm
#' @importFrom optimx optimx
#' @importFrom Matrix bdiag


#############################################################################################
#' @title Fitting A Generalized Linear Model For Mixed Outcomes Of Type Gaussian, Poisson and
#'     Binomial Using A Shared Random Effect
#' @description mixoglmm is used to fit a model with shared random effects for Gaussian, Poisson and Binomial random variables.
#' @param formula A formula as in the Formula package?
#' @param data data
#' @param families named list of families. Each element of the list is a family() as in glm. names must correspond to the responses
#' specified in the formula.only gaussian(), binomial() and poisson() are supported so far.
# #' @param gf grouping factor
#' @param cor_struct_gauss correlation structure of the gaussian responses
#' @param constraints.beta list of constraints
#' @param constraints.lambda list of constraints of length one
#' @param Ntrials list of number of trials to be used for the binomial. Defaults to NULL, will be one.
#' @param na.action eliminate NAs or keep them
#' @param weights  weights which need to be constant across multiple measurements. Negative weights are not allowed.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
#' @param control  list of parameters for controlling the fitting process. See \code{\link{mixoglmm.control}} for details.
#' @param ...  additional arguments.
#' @return an object of class mixoglmm
#' @examples
#' sum(1:10)
#'
#' \dontrun{
#' sum("a")
#' }
#' @export
#'
mixoglmm <- function(formula, data,
                   families = NULL,
                   cor_struct_gauss = NULL,
                   constraints.beta = NULL,
                   constraints.lambda = NULL,
                   Ntrials = NULL,
                   na.action, weights,
                   offset = NULL,
                   contrasts = NULL,
                   control = mixoglmm.control(),
                   ...) {
  ## call
  call <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, model response
  mt <- attr(mf, "terms")
 # mtX <- terms(formula, data = data, rhs = 1L)
  Y <- Formula::model.part(formula, data = mf, lhs = 1L)
  x <- model.matrix(formula, data = mf, rhs = 1L, contrasts = contrasts)
  # x <- if (!is.empty.model(mt))
  #   model.matrix(mt, mf, contrasts)
  # else matrix(, NROW(y), 0)
  # #attr(x, "assign") <- attrassigndefault(x, mt)
  ## offset
  offset_form <- model.offset(mf)
  offset_arg  <- offset
  if (!is.null(offset_form) & !is.null(offset_arg)) stop("Offsets should be specified either in the formula or in the argument.")
  if (!is.null(offset_form)) {
    offset <- rep(offset_form, NCOL(Y))
  }
  if (is.null(offset)) offset <- 0
  if (length(offset) == 1) offset <- rep(rep.int(offset, NROW(Y)), NCOL(Y))
  if (!is.null(offset_arg)) {
  if (length(offset_arg) != NCOL(Y)) stop(sprintf("Offsets should be a list of length equal to the number of responses. Each element of the list must be a vector of length %i", NROW(Y)))
  if (any(sapply(offset, function(x) length(x) > 0 & length(x) != NROW(Y)))) stop(sprintf("Offsets should be a list of length equal to the number of responses. Each element of the list must be a vector of length %i", NROW(Y)))
   offset <- c(sapply(offset_arg, function(x) if (length(x) == 0) rep.int(0, NROW(Y)) else x))
  }

  ## Ntrials # todo check if Ntrials specified for other than binomial
  if (is.null(Ntrials)) Ntrials <- 1
  if(length(Ntrials) == 1) Ntrials <- rep(list(rep.int(Ntrials, NROW(Y))), NCOL(Y))
  if (length(Ntrials) != NCOL(Y)) stop("The argument Ntrials must be a list of length equal to the number of responses.")
  if (any(sapply(Ntrials, function(x) length(x) > 0 & length(x) != NROW(Y)))) stop(sprintf("Ntrials must be a list of length equal to the number of responses. Each element of the list must be a vector of length %i.", NROW(Y)))
  Ntrials <- sapply(Ntrials, function(x) if (length(x) == 0) rep.int(1, NROW(Y)) else x)


  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, NROW(Y))
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  ## families
  if (is.null(families)) {
    families <- rep(list(gaussian()), NCOL(Y))
  } else {
    families <- check_families(families, x, Y)
  }
  ## sanity checks
  Y <- check_y(families, Y)
  ## correlation structure
  if (is.null(cor_struct_gauss)){
    cor_struct_gauss <- cor_ident()
  }
  cor_struct_gauss <- check_cor_structure(cor_struct_gauss, families)
  # constraints
  if (is.null(constraints.beta)) {
    constraints.beta <- rep(list(diag(NCOL(Y))), NCOL(x))
    names(constraints.beta) <- colnames(x)
    # cat("No constraints on the regression coefficients are specified. These will be set to identity and one coefficient for each response and each covariate will be estimated.")
  } else {
    check_constraints(constraints.beta, x, Y)
    constraints.beta <- constraints.beta[colnames(x)]
  }
  if (is.null(constraints.lambda)) {
    constraints.lambda <- list(diag(NCOL(Y)))
  } else {
    if (length(constraints.lambda) > 1) stop("Constraints.lambda must be a list of length one.")
  }
  ## Fit model
  fit <- mixoglmm_fit(y = Y, x = x,
                    cor_structure = cor_struct_gauss,
                    constraints.beta = constraints.beta,
                    constraints.lambda = constraints.lambda,
                    w =  weights, Ntrials = Ntrials, offset = offset,
                    families = families, control = control)
  fit$families <- families
  fit$call <- call
  fit$mf <- mf
  fit$responses <- Y
  fit$x <- x
  fit$constraints.beta <- constraints.beta
  fit$constraints.lambda <- constraints.lambda
  fit$weights <- weights
  fit$offset <- offset
  fit$formula <- formula
  fit$control <- control
 # fit$gf <- gf
  class(fit) <- "mixoglmm"

  fit
}

#' @title Print Method for class mixoglmm.
#' @description Prints regression coefficients, variance of random effects and parameters of the error structure of class mixoglmm.
#' @param x object of class \code{'mixoglmm'}
#' @param call displays function call if \code{TRUE}
#' @param ... further arguments passed to or from other methods.
#' @method print mixoglmm
#' @export
print.mixoglmm <- function(x, call = FALSE, ...) {
  if(call){
    cat("\nCall:\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }
  cat("Formula:\n")
  print(x$formula)
  cat("\n")
  cat("Log-likelihood at convergence:", - x$objective, "\n")
  pars <- x$parameters
  ## ------------------------------
  if (x$dims$Pstar == 0) {
    cat("No regression coefficients.\n\n")
  } else {
    cat("\nRegression coefficients:\n")
    beta <- pars[seq_len(x$dims$Pstar)]
    print(beta, ...)
  }
  ## ------------------------------
  cat("\nRandom effects:\n")
  cat("Random effects coefficients:\n")
  lambda <- pars[-seq_len(x$dims$Pstar + 1 + x$dims$G + x$dims$K2)]
  print(lambda, ...)
  ## ------------------------------
  cat("Random effects standard deviation:\n")
  tau <-  pars[x$dims$Pstar + 1]
  mat <- as.data.frame(cbind(tau))
  x$gf <- "Subject"
  mat[,"Grouping factor"] <-  x$gf
  mat <- mat[,c(2,1)]
  colnames(mat)[2] <- "Std.dev"
  print(mat, row.names = FALSE, ...)
  cat("Number of groups: ",  x$dims$N, "\n")
  ## ------------------------------
  cat("\nCorrelation of the Gaussian response variables, conditional on the random effects:\n")
  if (x$dims$G == 0) {
    cat("No correlation parameters.\n")
  } else {
    gamma <- pars[x$dims$Pstar + 1 + seq_len(x$dims$G)]
    print(gamma, ...)
  }
  ## ------------------------------
  cat("\nStandard deviation of the Gaussian response variables, conditional on the random effects:\n")
  omega <- pars[x$dims$Pstar + 1 + x$dims$G + seq_len(x$dims$K2)]
  #fams <- sapply(x$families, "[[", "family")
  #idnn.col <- which(fams != "gaussian")
  #names(omega2) <- colnames(x$responses)[-idnn.col]
  print(omega, ...)
}
#' @title Summary method for mixoglmm.
#' @description Summary of thresholds, regression coefficients
#' and parameters of the error structure of class \code{'mixoglmm'}.
#' @param object object of class \code{'mixoglmm'}
#' @param call displays function call if \code{TRUE}
#' @param ... further arguments passed to or from other methods.
#' @method summary mixoglmm
#' @export
summary.mixoglmm <- function(object, call = FALSE, ...)
{
  mat <- cbind.data.frame(c("nunits", object$dims$N), c("ndim", object$dims$M),
                          c("logLik", round(-object$objective,2)))
  ## extend coefficient table
  Pstar <- object$dims$Pstar
  nlambda <- object$dims$nlambda
  G <- object$dims$G
  K2 <- object$dims$K2
  ###########################
  cf <- object$parameters
  se <- sqrt(diag(object$vcov))
  se <- c(se[seq_len(Pstar + 1 + G + K2)], 0,
          se[- seq_len(Pstar + 1 + G + K2)])
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  cf[cf[, 2] == 0, 3] <- NA
  cf[cf[, 2] == 0, 4] <- NA
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  summary.output <- list()
  summary.output$formula <- object$formula
  ## ------------------------
  cat("Formula: ")
  print(summary.output$formula, ...)
  cat("\n")
  summary.output$info <- format(mat, justify="right")
  write.table(summary.output$info, row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(call){
    summary.output$call <- object$call
    cat("\nCall: ",
        paste(deparse(object$call, width.cutoff = getOption("deparse.cutoff")), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }
  ## ------------------------
  cat("\nCoefficients:\n")
  cfbeta <- cf[seq_len(Pstar), ,drop = FALSE]
  summary.output$coefficients <- printCoefmat(cfbeta, signif.legend = FALSE)
  #---------------------------------------
  cat("\nRandom effects:\n")
  cat("\tCoefficients:\n")
  cflambda <- cf[Pstar + 1 + G + K2 + seq_len(nlambda + 1), , drop = FALSE]
  summary.output$re.coefficients <- printCoefmat(cflambda, signif.legend = FALSE)
  cat("\tStandard deviation:\n")
  summary.output$re.stddev <- printCoefmat(cf[Pstar +  1, , drop = FALSE],
                                           signif.legend = FALSE)
  ## ------------------------------
  cat("\nCorrelation of the Gaussian response variables:\n")
  if (G == 0) {
    cat("No correlation parameters.\n")
  } else {
    summary.output$gauss.corr <- printCoefmat(cf[Pstar + 1 + seq_len(G), , drop = FALSE],
                                             signif.legend = FALSE)
  }
  ## ------------------------------
  cat("\nStandard deviation of the Gaussian response variables:\n")
  summary.output$gauss.stddev <- printCoefmat(cf[Pstar + 1 + G + seq_len(K2), , drop = FALSE],
                                           signif.legend = FALSE)
  class(summary.output) <- "summary.mixoglmm"
  return(invisible(summary.output))
}
# print.summary.mixoglmm <- function(summary.output, ...){
#   cat("Formula: ")
#   print(summary.output$formula, ...)
#   cat("\n")
#   write.table(summary.output$info, row.names = FALSE, col.names = FALSE, quote = FALSE)
#   cat("\nCoefficients:\n")
#   print(summary.output$coefficients)
# }

#' @title Control functions for mixoglmm()
#' @description Control arguments are set.
# #' @param start.values.beta vector of (optional) starting values for the regression coefficients.
#' @param nGHQ integer for the number of quadrature points to be used in the Gauss-Hermite quadrature.
#' @param solver.nlminb.control a list of control arguments to be passed to \code{\link{nlminb}}.
# #' @param scale If \code{scale = TRUE}, then for each response the corresponding covariates of \code{\link{class}} \code{"numeric"} are standardized before fitting,
# #'  i.e., by substracting the mean and dividing by the standard deviation.
#' @param NR.control a list of control arguments to be passed to the Newton Raphson routine for extracting the posterior modes of the random effects.
#' @seealso \code{\link{mixoglmm}}
#' @export
mixoglmm.control <- function(#start.values.beta = NULL,
                             nGHQ = 10L,
                             solver = "newuoa",
                             solver.optimx.control = list(maxit = 10000, trace = 0, kkt = FALSE),
                             NR.control = list(maxIter = 500, gradTol = 1e-6, tolerance = 1e-8, trace = 1,
                                               innerCtrl = c("giveError", "warn"),
                                               maxLineIter = 50)
                             ){
  if (is.null(solver.optimx.control$maxit)) solver.optimx.control$maxit <- 10000
  if (is.null(solver.optimx.control$trace)) solver.optimx.control$trace <- 0
  ## Controls for the Newton-Raphson Algorithm used in extracting the estimators of the random effects
  if (is.null(NR.control$maxIter)) NR.control$maxIter <- 100
  if (is.null(NR.control$tolerance)) NR.control$tolerance <- 1e-8
  if (is.null(NR.control$trace)) NR.control$trace <- 1
  if (is.null(NR.control$maxLineIter)) NR.control$maxLineIter <- 50
  list(#start.values.beta = start.values.beta,
       nGHQ = nGHQ,
       solver = solver,
       solver.optimx.control = solver.optimx.control,
       NR.control = NR.control)
}
#' @title vcov of Multivariate Ordinal Regression Models.
#' @description
#' \code{vcov} is a generic function which extracts the Godambe information matrix from objects of class \cr
#' \code{'mixoglmm'}.
#' @param object an object of class \code{'mixoglmm'}.
#' @param ... further arguments passed to or from other methods.
#' @method vcov mixoglmm
#' @export
vcov.mixoglmm <- function(object, ...) object$vcov

#' #' @title Number of units in the GLMM with mixed responses.
#' #' @description \code{nobs} is a generic function which extracts the number of units in the GLMM with mixed responses.
#' #' @param object an object of class \code{'mixoglmm'}.
#' #' @param ... further arguments passed to or from other methods.
#' #' @method nobs mixoglmm
#' #' @export
#' nobs.mixoglmm <- function(object, ...) object$dims$N
#'

#' @title Correlation Structures in cocas
#' @description Different \code{error.structures} are available in \pkg{mvord}:
#' \itemize{
#' \item general correlation structure (default) \code{cor_general(~ 1)},
#' \item general covariance structure \code{cov_general(~ 1)},
#' }
#' For more details see vignette.
#' @param formula \code{\link{formula}} object
#' @param value specifies values of the correlation  parameters.
#' For \code{cor_equi()}  a vector of correlations (in (-1,1)) of length one
#' (same correlation for all subjects) or a vector of length equal to the number of subjects.
#' For \code{cor_general()}, it can be either a vector of the lower triangular correlation matrix elements
#' (same structure for all subjects) or a matrix with number of rows equal to the number of subjects.
#' For \code{cov_general()}, it can be either a vector of the lower triangular covariance matrix elements (including the diagonal)
#' (same structure for all subjects) or a matrix with number of rows equal to the number of subjects.
#' Default is \code{value = numeric(0)} object.
#' In this case the correlation parameters are initialized with zero (and variance parameters with 1 for \code{cov_general})
#' @param fixed logical specifying whether the parameters of the error structure should not be optimized in the procedure, but will be \cr
#' fixed to the values specified in the argument \code{value}. Defaults to \code{fixed = FALSE}.
#' @export
#' @name cor_struct
cor_general <-
  ## Constructor for the cor_general class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE) {
    obj <- list(name = "cor_general", formula = formula,
                value = value, fixed = fixed)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_general", "cor_struct")
    obj
  }
#' @rdname cor_struct
#' @export
cor_equi <-
  ## Constructor for the cor_equi class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE)
  {
    obj <- list(name = "cor_equi",
                formula = formula,
                value = value, fixed = fixed)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_equi", "cor_struct")
    obj
  }
#' @rdname cor_struct
#' @export
cor_ident <-
  ## Constructor for the cor_general class
  function(formula = ~ 1, value = numeric(0), fixed = TRUE) {
    obj <- list(name = "cor_ident", formula = formula,
                value = value, fixed = fixed)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_ident", "cor_struct")
    obj
  }

##############
build_cor <-
  ## builder the correlation matrix
  function(eobj, ...) UseMethod("build_cor")

transf_cor <-
  ## returns a vector of parameters of the correlation structure to be shown in summary
  function(eobj, ...) UseMethod("transf_cor")

dtcor.cor <- function(eobj, ...) {
  ## takes the original parameters which are to be displayed in summary
  ## computes the Jacobian for the transformed parameters entering the optimizer
  ## d tpar/d par
  UseMethod("dtcor.cor")
}

#build_error_struct_fixed <-
#  ## extractor for correlation matrix
#  function(eobj, ...) UseMethod("build_error_struct_fixed")

#start_values <-
#  ## extractor for correlation matrix
#  function(eobj, ...) UseMethod("start_values")

#finalize_fun <-
#  ## finalizes the structures
#  function(eobj, ...) UseMethod("finalize_fun")

finalize <-
 function(eobj, ...) UseMethod("finalize")

#get_covariate <-
#  ## initializes the structures
#  function(eobj, ...) UseMethod("get_covariate")

init_fun <-
  ## initializes the structures
  function(eobj, ...) UseMethod("init_fun")
#####################
formula.cor_struct <-
  ## Accessor for the formula
  function(x, ...) eval(attr(x, "formula"))

#####################
init_fun.cor_general <-
  function(eobj,  y = NULL) {
    J <- NCOL(y)
    y_names <- colnames(y)
    attr(eobj, "ynames") <-
    attr(eobj, "ndim") <- J
    attr(eobj, "npar") <- J * (J - 1)/2
    attr(eobj, "start") <- double(J * (J - 1)/2)
    attr(eobj, "par_names") <- apply(combn(J,2), 2, function(x) paste0("corr ", y_names[x[1]], y_names[x[2]]))
    eobj
  }

init_fun.cor_equi <-
  function(eobj,  y = NULL) {
    attr(eobj, "ynames") <- colnames(y)
    attr(eobj, "ndim") <- NCOL(y)
    attr(eobj, "npar") <- 1
    attr(eobj, "start") <- double(1)
    attr(eobj, "par_names") <- "corr"
    eobj
  }

init_fun.cor_ident <-
  function(eobj,  y = NULL) {
    attr(eobj, "ynames") <- colnames(y)
    attr(eobj, "ndim") <- NCOL(y)
    attr(eobj, "npar") <- 0
    attr(eobj, "start") <- double(0)
    attr(eobj, "par_names") <- NULL
    eobj
  }
################
build_cor.cor_general <-
  function(eobj, tpar) {
    ## takes the transformed parameters and builds the general correlation matrix
    ## uses the spherical parameterization
    ndim <- attr(eobj, "ndim")
    angles <- pi * exp(tpar)/(1 + exp(tpar))
    cosmat <- diag(ndim)
    cosmat[lower.tri(cosmat)] <- cos(angles)
    S1 <- matrix(0, nrow = ndim, ncol = ndim)
    S1[, 1L] <- 1
    S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = TRUE)] <- sin(angles)
    tLmat <- sapply(seq_len(ndim),
                    function(j) cosmat[j, ] * cumprod(S1[j, ]))
    sigma <- crossprod(tLmat)
    return(sigma)
  }

build_cor.cor_equi <-
  function(eobj, tpar) {
    ## takes the transformed parameters and builds the general correlation matrix
    ## uses fisher z transformation
    ndim <- attr(eobj, "ndim")
    sigma <- diag(ndim)
    sigma[lower.tri(sigma)] <- sigma[upper.tri(sigma)] <- z2r(tpar)
    return(sigma)
  }

z2r <- function(z) {
  ifelse(z > 354, 1, (exp(2 * z) - 1)/(1 + exp(2 * z)))
}

build_cor.cor_ident <-
  function(eobj, tpar) {
    ## takes the transformed parameters and builds the general correlation matrix
    ## uses fisher z transformation
    ndim <- attr(eobj, "ndim")
    sigma <- diag(ndim)
    return(sigma)
  }
################
transf_cor.cor_general <-
  function(eobj, tpar) {
    ## takes the transformed parameters and builds the general correlation matrix
    ## uses the spherical parameterization
    ndim <- attr(eobj, "ndim")
    angles <- pi * exp(tpar)/(1 + exp(tpar))
    cosmat <- diag(ndim)
    cosmat[lower.tri(cosmat)] <- cos(angles)
    S1 <- matrix(0, nrow = ndim, ncol = ndim)
    S1[, 1L] <- 1
    S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = TRUE)] <- sin(angles)
    tLmat <- sapply(seq_len(ndim),
                    function(j) cosmat[j, ] * cumprod(S1[j, ]))
    sigma <- crossprod(tLmat)
    return(sigma[lower.tri(sigma)])
  }

transf_cor.cor_equi <-
  function(eobj, tpar) {
    ## takes the transformed parameters and builds the general correlation matrix
    ## uses fisher z transformation
    return(z2r(tpar))
  }

transf_cor.cor_ident <-
  function(eobj, tpar) {
    ## takes the transformed parameters and builds the general correlation matrix
    ## uses fisher z transformation
  }

################
finalize.cor_struct <- function(eobj, tpar) {
  R <- build_cor(eobj, tpar)
  colnames(R) <- rownames(R) <- attr(eobj,"ynames")

  R
}
####
backtransf_spherical <- function(par, ndim, i = NULL){
  if (is.null(i)) i <- seq_along(par)
  R <- angmat <- matrix(1, ncol = ndim , nrow = ndim)
  R[lower.tri(R)] <- R[upper.tri(R)] <- par
  l <- t(chol(R))
  angmat[-1,1] <- acos(l[-1,1])
  for (j in 2:(ndim - 1)){
    sinprod <- apply(sin(angmat[, seq_len(j-1), drop = FALSE]), 1, prod) ## denominator in division
    angmat[-seq_len(j), j] <- acos((l/sinprod)[-seq_len(j), j])
  }
  angdivpi <- angmat[lower.tri(angmat)]/pi
  log(angdivpi/(1-angdivpi))[i]
}


dtcor.cor.cor_general <-
  function(eobj, par) {
    ndim <- attr(eobj, "ndim")
    x <- lapply(seq_along(par), function(i)
      numDeriv::grad(function(par) backtransf_spherical(par, ndim, i), par))
    do.call("rbind", x)
  }

dtcor.cor.cor_equi <-
  function(eobj, par) {
    x <- 1/(1 + par) * 1/(1 - par)
    if (length(par) == 1) x else diag(x)
  }


dtcor.cor.cor_ident <-
  function(eobj, par) {
    diag(0)
  }

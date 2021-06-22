check <- function(...){
  stopifnot(...)
}

check_constraints <- function(constraints, x, y) {
  if (length(constraints) != NCOL(x))  stop(sprintf("Length of constraint list does not coincide with the number of covariates in the model.matrix.
                                                    The model.matrix contains %i columns with column names: %s.",
                                                    NCOL(x), paste(colnames(x), collapse = ", ")))
  if (!all(sapply(constraints, NROW) == NCOL(y)))  stop("The number of rows of the constraint matrices does not equal the number of responses.")
  ## todo what if the constraint matrices are specified in different order
  ## match(names(constraints), colnames(x))
  ## TODO make a function build constraints?
  ## give only warning if the intercept is forgotten?
  if (anyNA(match(names(constraints), colnames(x))))
    stop(sprintf("The names of the constraint matrices do not fit the names of the covariates. The names of the covariates in the model.matrix are: %s.",
                 paste(colnames(x), collapse = ", ")))
}



check_families <- function(families, x, y) {
  if (length(families) != NCOL(y)) stop("Number of specified families does not coincide with the number of responses.")
  if (anyNA(match(names(families), colnames(y)))) stop("The names of the specified families do not fit the names of the response variables in formula.")
  if (!all(sapply(families, "[[", "family") %in% c("gaussian", "binomial", "poisson")))  stop("One or more of the specified families are not supported. Only binomial(), poisson() or gaussian() are supported by the mixoglmm package.")
  families <- lapply(families, function(x) {
    x$loglik <- switch(x$family,
                         "binomial" = function(y, mu, w, n) {
                           #if (!(x$validmu(mu))) mu <- pmax(0, pmin(mu, 1))
                           w * dbinom(y, n, mu, log = TRUE)
                           },
                         "poisson"  = function(y, mu, w, n = NULL) w * dpois(y, mu, log = TRUE))
    x$dev.resids <- x$aic <- x$initialize <- x$simulate <- NULL
    x
  })
  families <- families[colnames(y)]
  return(families)
}

check_y <- function(family, y) {
  for (i in seq_len(NCOL(y))) {
    if ((family[[i]]$family == "binomial") && !is.numeric(y[, i])) {
      if (!is.factor(y[, i])) {
        stop("The Binomial family only allows for numeric variables or factors with two levels.")
      } else {
        if (nlevels(y[, i]) > 2) stop("The Binomial family only allows for factors with two levels.")
        y[, i] <- as.integer(y[, i] != levels(y[, i])[1L])
      }
    }
    if ((family[[i]]$family == "poisson") && !is.numeric(y[, i])) {
      stop("The Poisson family only allows for numeric variables.")
    }
    if ((family[[i]]$family == "gaussian") && !is.numeric(y[, i])) {
      stop("The Gaussian family only allows for numeric variables.")
    }
  }
  y
}

check_cor_structure <- function(cor_structure, families) {
  K2 <- sum(sapply(families, "[[", "family") == "gaussian")
  formula <- cor_structure$formula
  if (K2 == 1) {
    cat("Only one Gaussian response, setting the correlation structure to cor_ident().")
    cor_structure <- cor_ident(formula = formula)
  }
  if (K2 == 2 && cor_structure$name == "cor_general"){
    cat("Only two Gaussian responses, setting the correlation structure to cor_equi().")
    cor_structure <- cor_equi(formula = formula)
  }
  cor_structure
}

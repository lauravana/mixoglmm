dmvnorm.log.cocas1 <-
  function(x, sigmainv) {
    # computes the log likelihood of a centered multivariate normal x
    distval <- mahalanobis2(x, sigmainv)
    logdet <- tryCatch(sum(log(1 / eigen(sigmainv, symmetric = TRUE, only.values = TRUE)$values)))
    logretval <- - (ncol(x) * log(2 * pi) + logdet +  distval)/2
    sum(logretval)
  }

mahalanobis2 <-
  function(x, sigmainv) {
    # computes the (x - center)' sigmainv (x - center) where sigmainv = Sigma^-1
    setNames(rowSums((x %*% sigmainv) * x), rownames(x))
  }


mahalanobis.l <- function(x, center, cov, inverted = FALSE) {
  x <- if (is.vector(x))
    matrix(x, ncol = length(x))
  else as.matrix(x)
  x <- sweep(x, 2, center)
  if (!inverted)
    cov <- chol2inv(chol(cov))
  setNames(rowSums((x %*% cov) * x), rownames(x))
}

dmvnorm.log.cocas <- function(x, mean, sigma, inverted) {
  distval <- mahalanobis.l(x, center = mean, cov = sigma, inverted = inverted)
  if (!inverted) {
    logdet <- tryCatch(sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values)))
  } else {
    logdet <- tryCatch(sum(log(1 / eigen(sigma, symmetric = TRUE, only.values = TRUE)$values)))
  }

  logretval <- - (ncol(x) * log(2 * pi) + logdet +  distval)/2
  sum(logretval)
}
## links
## TODO include random effect in formula, otherwise it looks weird?

make.dmu.deta <- function(linkstr) {
  ## needed for the second derivative wrt to the random effects
  switch(linkstr,
         "logit"    = {
           logit_link <- make.link("logit")
           function(eta) logit_link$mu.eta(eta) * (1 - 2 * logit_link$linkinv(eta))
         },
         "probit"   = function(eta) -eta * pmax(dnorm(eta), .Machine$double.eps),
         "cauchit"  = function(eta) -2 *  eta/(pi * (1 + eta^2)^2), #-2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps),
         "cloglog"  = function(eta) pmax((1 - exp(eta)) * exp(eta - exp(eta)), .Machine$double.eps),
         # not implemeted "loglog"   = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
         "identity" = function(eta) rep.int(0, length(eta)),
         "log"      = function(eta) pmax(exp(eta), .Machine$double.eps),
         "sqrt"     = function(eta) rep.int(2, length(eta)),
         "1/mu^2"   = function(eta) 3/(4 * eta^2.5),
         "inverse"  = function(eta) 2/(eta^3))
}


make_coef_names <- function(x_names = NULL, y_names, constraints) {
  unlist(lapply(seq_along(constraints), function(i) {
    cm <- constraints[[i]]
    if (all(cm %in% c(0,1))) {
      nam <- apply(cm, 2, function(x) paste0(y_names[as.logical(x)], collapse = "."))
    } else {
      nam <- seq_len(NCOL(cm))
    }
    paste(x_names[i], nam)
  }))
}

update_families <- function(families) {
  lapply(families, function(x) {
    ## first derivative of the log likeligood wrt to mu
    x$loglik.mu <- switch(x$family,
                          "binomial" = function(y, mu, w, n) w * ((y/mu) - (n - y)/(1 - mu)),
                          "poisson"  = function(y, mu, w, n) w * (y/mu + (dpois(y, mu, log = TRUE) - y * log(mu))/mu))
    ## second derivative of the log likeligood wrt to mu
    x$dloglik.dmu <- switch(x$family,
                            "binomial" = function(y, mu, w, n) w * (-(y/mu^2) - (n - y)/(1 - mu)^2),
                            "poisson"  = function(y, mu, w, n) w * (- y/mu^2))
    ## second derivative of mu wrt to eta
    x$dmu.deta <-  make.dmu.deta(x$link)
    x
  })
}
transform_parameters2 <- function(tpar, dims, y2, x_constr, constraints.lambda,
                                  offset, cor_structure,
                                  idnn.row, idnn.col, ind.y2) {
  # print(tpar)
  beta   <-     tpar[seq_len(dims$Pstar)]
  tau2   <- exp(tpar[dims$Pstar + 1])
  gamma  <-     tpar[dims$Pstar + 1 + seq_len(dims$G)]
  omega  <- exp(tpar[dims$Pstar + 1 + dims$G + seq_len(dims$K2)])
  lambda <- tpar[dims$Pstar + 1 + dims$G + dims$K2 + seq_len(NCOL(constraints.lambda) - 1L)]
  lambda <- drop(constraints.lambda %*% c(1, lambda))
  lambda1 <- lambda[ idnn.col] # rep.int(1, dims$K1)
  lambda2 <- lambda[-idnn.col] # rep.int(1, dims$K2)

  R <- build_cor(cor_structure, gamma)
  Omega <- tcrossprod(omega) * R
  Sigma <- Omega + tau2 * tcrossprod(lambda2)
  Sigmainv <- tryCatch(chol2inv(chol(Sigma)))
  ## Xbeta
  xbeta  <- x_constr %*% beta + offset
  ## Xbeta non-normal responses
  xbeta1 <- matrix(xbeta[ idnn.row, ], ncol = dims$K1)
  ## Xbeta normal responses
  xbeta2 <- matrix(xbeta[-idnn.row, ], ncol = dims$K2)

  ## Errors of normal responses
  eps <- y2 - xbeta2
  eps[is.na(eps)] <- 0
  ## Parameters of the (normal) posterior distribution of the random effects
  lambda2TSigmainv <- crossprod(lambda2, Sigmainv)
  kappa2 <- tau2 - tau2^2 * drop(lambda2TSigmainv %*% lambda2)
  nu     <- tau2 * drop(tcrossprod(eps, lambda2TSigmainv))
  #nu <- tau2 * sapply(seq_len(dims$N), function(i) sum(eps[i, ] * lambda2TSigmainv, na.rm = T))
  output <- list(lambda1 = lambda1,
                 lambda2 = lambda2,
                 Sigma = Sigma,
                 Sigmainv = Sigmainv,
                 xbeta1 = xbeta1,
                 xbeta2 = xbeta2,
                 y2errors = eps,
                 nu = nu, kappa2 = kappa2)
}

dttau2.tau <- function(x) if (length(x) == 1) 2/x^2 else diag(2/x^2)

dtomega.omega <- function(x) if (length(x) == 1) 1/x else diag(1/x)

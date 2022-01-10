#' @title Extract random effects for mixoglmm.
#' @description extract random effects .
#' @param object object of class \code{'mixoglmm'}
#' @param method conditional modes or conditional means of the random effects
#' @export
extract_ranef <- function(object, method = c("conditional means", "conditional modes")) {
  ######################
  par <- object$res$par
  w <- object$weights
  Ntrials <- object$Ntrials
  control <- object$control
  dims <- object$dims
  fams <- sapply(object$families, "[[", "family")
  family_nn <- object$families[fams != "gaussian"]

  ## Functions connected to the family
  family_nn <- update_families(family_nn)
  y1   <- object$responses[, fams != "gaussian", drop = FALSE]
  y2   <- as.matrix(object$responses[, fams == "gaussian", drop = FALSE])

  x_constr <- do.call("cbind", lapply(seq_along(object$constraints.beta),
                                      function(p) kronecker(object$constraints.beta[[p]], object$x[, p])))

  idnn.col <- which(fams != "gaussian")
  idnn.row <- as.vector(sapply(which(fams != "gaussian"), function(i) (i-1) * object$dims$N + seq_len(object$dims$N)))
  tmp <- transform_parameters2(par, object$dims, y2, x_constr, object$constraints.lambda[[1]],
                               object$offset,
                               object$cor_structure, idnn.row, idnn.col)
  gq <- statmod::gauss.quad(control$nGHQ, kind = "hermite")
  ranef <- switch(method,
                  "conditional means" = extract_ranef_means(y1, tmp, family_nn, w, Ntrials, dims, gq),
                  "conditional modes" = extract_ranef_modes(y1, tmp, family_nn, w, Ntrials, dims,
                                                            control$NR.control, gq))
  return(ranef)
}



extract_ranef_modes <- function(y1, tmp, family_nn, w, Ntrials, dims, ctrl, gq) {
  ## gradient of log posterior density of u wrt to u
  gprime <- function(u) {
    rowSums(sapply(seq_len(dims$K1), function(k) {
      eta <- tmp$xbeta1[, k] + tmp$lambda1[k] * u
      mu <- family_nn[[k]]$linkinv(eta)
      tmp$lambda1[k] * family_nn[[k]]$loglik.mu(y1[, k], mu, w, Ntrials[, k]) *
        family_nn[[k]]$mu.eta(eta)
    }), na.rm = TRUE) - (u - tmp$nu)/tmp$kappa2
  }
  ## hessian of log posterior density of u wrt to u
  gprimeprime <- function(u) {
    rowSums(sapply(seq_len(dims$K1), function(k) {
      eta <- tmp$xbeta1[, k] + tmp$lambda1[k] * u
      mu <- family_nn[[k]]$linkinv(eta)
      tmp$lambda1[k]^2 *
        (family_nn[[k]]$dloglik.dmu(y1[, k], mu, w, Ntrials[,k]) * family_nn[[k]]$mu.eta(eta)^2 +
           family_nn[[k]]$loglik.mu(y1[, k], mu, w, Ntrials[,k]) * family_nn[[k]]$dmu.deta(eta))
    }), na.rm = TRUE) - 1/tmp$kappa2
  }
  stepFactor <- 1
  innerIter <- 0
  ctrl$maxLineIter <- 10
  u <- 0
  gradient <- gprime(u)
  maxGrad <- max(abs(gradient))
  conv <- -1  ## Convergence flag
  message <- "Iteration limit reached when updating the random effects"
  logPostDens <- log_dens_u_cond_y(u, y1, tmp, family_nn, w, Ntrials, dims, gq)
  ## Newton-Raphson algorithm:
  ctrl$innerCtrl <- "giveError"
  for(i in 1:(ctrl$maxIter + 1L)) {
    if(maxGrad < ctrl$gradTol) {
      message <- "max|gradient| < tol, so current iterate is probably solution"
      if(ctrl$trace > 0)
        cat("\nOptimizer converged after", i, "iterations!", "max|grad|:",
            maxGrad, message, fill = TRUE)
      conv <- 0
      break
    }
    D  <- gprimeprime(u)
    step <- gradient / D
    u <- u - stepFactor * step

    logPostDensTry <- log_dens_u_cond_y(u, y1, tmp, family_nn, w, Ntrials, dims, gq)
    lineIter <- 0
    ## simple line search, i.e. step halfing:
    while(logPostDensTry < logPostDens) {
      stepFactor <- stepFactor/2
      u <- u + stepFactor * step
      logPostDensTry <- log_dens_u_cond_y(u, y1, tmp, family_nn, w, Ntrials, dims, gq)
      lineIter <- lineIter + 1
      if(ctrl$trace > 0) cat(sprintf("iter: %i \t step factor: %0.4f\t value: %0.6f \t max|grad|:%0.3f\n",
                                     i+innerIter, stepFactor, logPostDensTry, maxGrad))

      if(lineIter > ctrl$maxLineIter){
        message <- "Step factor reduced below minimum when computing the random effects."
        conv <- 1
        break
      }
      innerIter <- innerIter + 1
    }
    logPostDens <- logPostDensTry
    gradient <- gprime(u)
    maxGrad <- max(abs(gradient))
    stepFactor <- min(1, 2 * stepFactor)
  }
  if(conv != 0)
    stop(message, "\n  at iteration ", i)
  return(u)
}

extract_ranef_means <- function(y1, tmp, family_nn, w, Ntrials, dims, gq) {
  ## gauss hermite for the posterior means
  Z <- rep.int(1, dims$N)
  unodes <- tcrossprod(Z, sqrt(2 * tmp$kappa2) * gq$nodes) + tmp$nu
  const <- exp(rowSums(sapply(seq_len(dims$K1),  function(k) {
    eta <-  tmp$xbeta1[, k] + tmp$lambda1[k] * unodes
    phat   <- family_nn[[k]]$linkinv(eta)
    delta <- (phat %*% gq$weights)/ sqrt(pi)
    family_nn[[k]]$loglik(y1[, k], delta, w, Ntrials[, k])
  }), na.rm  = TRUE))

  ll <- lapply(seq_len(dims$K1),  function(k) {
    eta <- tmp$xbeta1[, k] + tmp$lambda1[k] * unodes
    mu  <- family_nn[[k]]$linkinv(eta)
    ll  <- family_nn[[k]]$loglik(y1[, k], mu, w, Ntrials[, k])
    ll[is.na(ll)] <- 0
    ll
  })
  huj <- unodes * exp(Reduce("+", ll))
  hu <- drop((huj %*% gq$weights)/ sqrt(pi))
  postmeanu <- hu/const
  return(postmeanu)
}

## log posterior density
log_dens_u_cond_y <- function(u, y1, tmp, family_nn, w, Ntrials, dims, gq) {
  logf_mu <- sum(sapply(seq_len(dims$K1),  function(k) {
    eta <- tmp$xbeta1[, k] + tmp$lambda1[k] * u
    mu  <- family_nn[[k]]$linkinv(eta)
    family_nn[[k]]$loglik(y = y1[, k], mu = mu, w = w,
                          n = Ntrials[, k])
  }), na.rm = TRUE)
  Z <- rep(1, dims$N)
  unodes <- tcrossprod(Z, sqrt(2 * tmp$kappa2) * gq$nodes) + tmp$nu
  logf_delta <- sum(sapply(seq_len(dims$K1),  function(k) {
    eta <- tmp$xbeta1[, k] + tmp$lambda1[k] * unodes
    phat   <- family_nn[[k]]$linkinv(eta)
    delta <- (phat %*% gq$weights)/ sqrt(pi)
    family_nn[[k]]$loglik(y = y1[, k], mu = delta, w = w,
                          n = Ntrials[, k])
  }), na.rm = TRUE)
  logf_mu + sum(dnorm(u, tmp$nu, sqrt(tmp$kappa2), log = TRUE)) - logf_delta
}

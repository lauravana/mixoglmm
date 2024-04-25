mixoglmm_fit <- function(y, x, cor_structure,
                         constraints.beta,
                         constraints.lambda, w, Ntrials,
                         offset, families,  control) {
  fams <- sapply(families, "[[", "family")
  idnn.col <- which(fams != "gaussian")
  idn.col  <- which(fams == "gaussian")
  N <- NROW(y)
  y1   <- y[,  idnn.col, drop = FALSE]
  y2   <- as.matrix(y[, -idnn.col, drop = FALSE])
  K1 <- NCOL(y1)
  K2 <- NCOL(y2)
  K  <- K1 + K2

  x_constr <- do.call("cbind", lapply(seq_along(constraints.beta),
                                      function(p) kronecker(constraints.beta[[p]], x[, p])))
  obj <- list()

  ## correlation structure for gaussian responses
  obj$cor_structure <- init_fun(cor_structure, y = y2)

  ## Constraints lambda
  if (!is.null(constraints.lambda) & is.list(constraints.lambda)) {
    if (length(constraints.lambda) > 1) stop("Constraints.lambda must be a list of length one.")
    constraints.lambda <- constraints.lambda[[1]]
  }
  if (is.null(constraints.lambda)) {
    constraints.lambda <- diag(NCOL(y))
    constraints.lambda[idn.col[1], 1] <- 1
    constraints.lambda <- constraints.lambda[, - idn.col[1], drop = FALSE]
  }
  ## Identifiability check:
  n_error_param_allowed <- (K * (K + 1)/2 - 1)
  n_error_param_to_estimate <-
    ## correlation params in Omega
    attr(obj$cor_structure, "npar") +
    ## Lambda
    NCOL(constraints.lambda) +
    ## Tau2
    0 +
    ## sigmas
    K2
  if (n_error_param_allowed < n_error_param_to_estimate)
    stop(sprintf("No of allowed parameters is %i, to estimate we have %i. Consider imposing further constraints on lambda or choosing a more parsimonious correlation structure.",
                n_error_param_allowed, n_error_param_to_estimate))
  
  obj$dims <- list(N = N, K = K, K1 = K1, K2 = K2,
                   Pstar =  NCOL(x_constr),
                   G = attr(obj$cor_structure, "npar"),
                   nlambda = NCOL(constraints.lambda))

  idnn.row <- as.vector(sapply(idnn.col, function(i) (i-1) * N + seq_len(N)))

  x2_constr <- x_constr[-idnn.row, ]
  family_nn <- families[idnn.col]



  ## starting values - TODO
  ty <- as.matrix(y)
  ## for binomial we need to transform
  ty <- ty/Ntrials
  ## for poisson we need to transform
  ty[, fams == "poisson"] <- (ty[, fams == "poisson"])^(1/3)
  start_beta <- glm(c(ty) ~ 0 + x_constr)$coefficients
  if (is.null(control$start_values)) {
    start_values <- c(double(obj$dims$Pstar), #start_beta,#
                      double(1), # tau
                      attr(obj$cor_structure , "start"), # rho
                      double(K2),#log(apply(y2, 2, sd, na.rm = TRUE)),
                      double(obj$dims$nlambda) + 1L)
  } else {
    if (length(control$start_values) !=
        (obj$dims$Pstar + 1 + obj$dims$G + K2 + obj$dims$nlambda)) stop("Starting values in control are not of correct length: Specify starting values for beta, log(tau), correlation structure, diagonal of Omega and lambdas")
    start_values <- control$start_values
  }


  Z <- rep.int(1, N)


  ## non - normal data log lik
  Ntrials <- Ntrials[,idnn.col, drop = FALSE]

  ind_y2 <-
    apply(unique(!is.na(y2)), 1, function(x)
    list(ind.row = which(colSums(t(!is.na(y2)) == x) == NCOL(y2) ),
         ind.col = x,
         normfun = ifelse(sum(x) == 1, dnorm, dmvnorm)))

  ## Gauss-Hermite Quadrature
  gq <- statmod::gauss.quad(control$nGHQ, kind = "hermite")

  ## Optimize negative log likelihood
  obj$res <- suppressWarnings(optimx(start_values, function(par) negloglik(par,
                                                                         y1, y2, x_constr, #x2_constr,
                                                                         Z, ind_y2,
                                                                         constraints.lambda,
                                                                         obj$cor_structure,
                                                                         w, Ntrials, offset,idnn.row, idnn.col,
                                                                         family_nn,#family_nn_ll,
                                                                         obj$dims, gq),
                                        method = control$solver,
                                        hessian = FALSE,
                                        control =  control$solver.optimx.control))
  # obj$res <- nlminb(start_values, function(par) negloglik(par,
  #                                                      y1, y2, x_constr, #x2_constr,
  #                                                      Z, ind_y2,
  #                                                      constraints.lambda,
  #                                                      obj$cor_structure,
  #                                                      w, Ntrials, offset,idnn.row, idnn.col,
  #                                                      family_nn,#family_nn_ll,
  #                                                      obj$dims, gq),
  #                   control = control$solver.nlminb.control)
  obj$par <- unlist(obj$res[1:length(start_values)])
  obj$objective <- unlist(obj$res["value"])
  if (obj$res$convcode != 0){
    print(obj$res)
    warning("NO/FALSE CONVERGENCE - choose a different optimizer, increase iterations or use different starting values.")
  }
  ## Compute Hessian numerically
  tparHess <- numDeriv::hessian(function(par) negloglik(par,
                                                 y1, y2, x_constr,
                                                 Z,ind_y2,
                                                 constraints.lambda,
                                                 obj$cor_structure,
                                                 w, Ntrials, offset,idnn.row, idnn.col,
                                                 family_nn,
                                                 obj$dims, gq), obj$par)

  ##---------------------
  tpar   <- obj$par
  beta   <- tpar[seq_len(obj$dims$Pstar)]
  ttau2 <-  tpar[obj$dims$Pstar + 1]
  tau    <- exp(ttau2/2)
  gamma  <- tpar[obj$dims$Pstar + 1 + seq_len(obj$dims$G)]
  gamma  <- transf_cor(obj$cor_structure, gamma)
  omega  <- exp(tpar[obj$dims$Pstar + 1 + obj$dims$G + seq_len(obj$dims$K2)])

  lambdas <- tpar[obj$dims$Pstar + 1 + obj$dims$G + obj$dims$K2 + seq_len(obj$dims$nlambda)]
  lambda  <- drop(constraints.lambda %*% c(1, lambdas))
  names.beta   <- make_coef_names(x_names = colnames(x),
                                  y_names = colnames(y),
                                  constraints = constraints.beta)
  names.lambda <- paste0("lambda", colnames(y))#make_lambda_names(y_names = colnames(y),
                  #                      constraints = constraints.lambda)

  names.tau <- "tau"
  names.gamma <- attr(obj$cor_structure, "par_names")
  names.omega <- colnames(y)[-idnn.col]

  obj$parameters <- c(beta, tau, gamma, omega, lambda)
  names(obj$parameters) <- c(names.beta, names.tau, names.gamma, names.omega, names.lambda)
  # --------------------
  if (TRUE) {
  J <- as.matrix(Matrix::bdiag(list(diag(obj$dim$Pstar), ## d tbeta/d beta
                              dttau2.tau(tau), ## d ttau2/d tau
                              dtcor.cor(obj$cor_structure, gamma),
                              dtomega.omega(omega),# diag(dtomega.omega(omega))),## d tomega/d omega
                              diag(obj$dims$nlambda) # d tlambda/d lambda
                              )))
  H <- crossprod(J, tparHess) %*% J
  obj$vcov <- tryCatch(chol2inv(chol(H)),
                       error=function(e) {
                         warning("\nCondition number close to zero! Hessian is approximated by nearest positive semidefinite matrix.\n")
                         chol2inv(chol(Matrix::nearPD(H)$mat))
                       }
                       )
  }
  ###########################
  obj$Ntrials <- Ntrials
  obj$constraints.lambda <- constraints.lambda
  obj
}

negloglik <- function(par, y1, y2,  x_constr, #x2_constr,
                   Z, ind_y2, constraints.lambda,
                   cor_structure,
                   w, Ntrials,offset,
                   idnn.row, idnn.col,
                   family_nn,
                   dims, gq) {
  ## assign parameters

  tmp <- transform_parameters2(par, dims, y2, x_constr, constraints.lambda, offset,
                               cor_structure, idnn.row, idnn.col, ind_y2)
  unodes <- tcrossprod(Z, sqrt(2 * tmp$kappa2) * gq$nodes) + tmp$nu
  nll1 <- - sum(sapply(seq_len(dims$K1),  function(k) {
    #ynodes <- tcrossprod(Z, sqrt(2 * kappa2) * gq$nodes) + xbeta1[, i] + lambda1[i] * nu
    eta <- tmp$xbeta1[, k] + tmp$lambda1[k] * unodes
    phat   <- family_nn[[k]]$linkinv(eta)
    delta <- (phat %*% gq$weights)/ sqrt(pi)
    family_nn[[k]]$loglik(y = y1[, k], mu = delta, w = w,
                          n = Ntrials[, k])
  }), na.rm  = TRUE)

  nll2 <- - sum(sapply(ind_y2, function(id) {
    ifelse(any(id$ind.col),
           sum(w[id$ind.row] * id$normfun(tmp$y2errors[id$ind.row, id$ind.col],
                                          mean = rep.int(0, sum(id$ind.col)),
                                          tmp$Sigma[id$ind.col, id$ind.col], log = TRUE)),
           0)
  }))

  nll1 + nll2
}

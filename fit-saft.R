# load libraries
library(splines)
library(survival)
library(SimSurvNMarker)

# source  we need
library(Rcpp)
sourceCpp("splines.cpp", embeddedR = FALSE)

# estimates a spline-based AFT.
#
# Args:
#   y: observed event time or censoring time.
#   X: design matrix for the fixed effects.
#   event: event indicator.
#   n_knots: number of knots for B-splines.
#   gl_dat: Gauss–Legendre Quadrature node and weights.
#   basis_type: type of basis to use.
#   use_integrate: logical for whether to use R's integrate function or
#                  Gauss–Legendre quadrature.
#   maxit: maximum number of iterations to pass to optim.
#   check_grads: TRUE if gradients should be checked.
#   fix_boundary_knots: TRUE if the boundary knots should be fixed.
#
# Returns:
#   coefs: coefficient estimates.
#   logLik: maximum log likelihood.
#   optim: output from optim.
#   eval_basis: function evaluate the basis functions.
#   beta: estimated AFT parameters.
#   gamma: estimated splines coefficients.
saft_fit <- function(y, X, event, n_knots, gl_dat, basis_type = c("bs", "ns"),
                     use_integrate = FALSE, maxit = 1000L, check_grads = FALSE,
                     fix_boundary_knots = TRUE){
  # get the knot placement
  knots <- quantile(y[event > 0], 1:n_knots / (n_knots + 1))
  b_knots <- range(y, 0)

  # assign function to evaluate the spline basis
  basis_type <- basis_type[1]
  switch(basis_type,
    bs = {
      # assign basis matrix we will use for memory while doing the log likelihood
      # evaluations
      Bmat <- matrix(0., NROW(X), n_knots + 3L)

      # get point and assign function to evaluate the basis functions
      bs_ptr <- NULL
      set_ptr <- function(b_knots_arg = b_knots)
        bs_ptr <<- get_bs_ptr(knots = knots, boundary_knots = b_knots_arg,
                             intercept = FALSE)
      set_ptr()

      eval_basis <- function(x, gamma, ders = 0L){
        eval_spline_basis_fill(x, bs_ptr, Bmat, ders = ders)
        # account for the intercept
        gamma[1] * (ders == 0L) + drop(Bmat %*% gamma[-1])
      }

      # function to return
      eval_basis_out <- function(x)
        cbind(1, eval_spline_basis(x, bs_ptr))
    },
    ns = {
      # assign basis matrix we will use for memory while doing the log likelihood
      # evaluations
      Bmat <- matrix(0., NROW(X), n_knots + 1L)

      # get point and assign function to evaluate the basis functions
      ns_ptr <- NULL
      set_ptr <- function(b_knots_arg = b_knots)
        ns_ptr <<- get_ns_ptr(knots = knots, boundary_knots = b_knots_arg,
                              intercept = FALSE)
      set_ptr()

      eval_basis <- function(x, gamma, ders = 0L){
        eval_spline_basis_fill(x, ns_ptr, Bmat, ders = ders)
        # account for the intercept
        gamma[1] * (ders == 0L) + drop(Bmat %*% gamma[-1])
      }

      # function to return
      eval_basis_out <- function(x)
        cbind(1, eval_spline_basis(x, ns_ptr))
    },
    stop(sprintf("basis_type '%s' is not implemented", basis_type)))

  # get starting values
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop=FALSE]
  beta <- -coef(survreg(Surv(y, event) ~ X, dist = "exponential"))
  gamma <- switch(
    basis_type,
    bs = numeric(n_knots + 4),
    ns = numeric(n_knots + 2),
    stop(sprintf("basis_type '%s' is not implemented", basis_type)))
  gamma[1] <- beta["(Intercept)"]
  beta <- beta[-1]
  n_beta <- length(beta)
  n_gamma <- length(gamma)

  # move nodes and recale the nodes
  gl_node <- (gl_dat$node + 1) * .5

  # assign negative log likelihood function
  ll_func <- function(par){
    beta <- head(par, n_beta)
    gamma <- tail(par, n_gamma)

    # handle terms from the hazard
    eta <- drop(X %*% beta)
    w <- exp(eta)
    w_y <- w * y
    if(!fix_boundary_knots)
      set_ptr(range(0, w_y))

    l1 <- -sum(event * (eta + eval_basis(w_y, gamma)))

    # handle terms from the survival function
    if(!use_integrate){
      # use Gauss–Legendre quadrature
      l2 <- mapply(function(x, wei)
        sum(w * exp(eval_basis(w * y * x, gamma)) * wei * .5 * y),
        x = gl_node, w = gl_dat$weight)
      l2 <- sum(l2)
    } else {
      # use R's integrate function
      l2 <- mapply(function(w, y){
        out <- try(
          integrate(function(z) exp(drop(eval_basis_out(w * z) %*% gamma)),
                    lower = 0, upper = y,
                    rel.tol = sqrt(.Machine$double.eps))$value, silent = TRUE)
        if(inherits(out, "try-error"))
          out <- Inf
        out
      }, w = w, y = y)
      l2 <- sum(l2 * w)
    }

    # cat(sprintf("l1: %.2f l2: %.2f sum: %.2f\n", l1, l2, l1 + l2))
    l1 + l2
  }

  # assign gradient of the negative log likelihood function
  d_ll_func <- function(par){
    beta <- head(par, n_beta)
    gamma <- tail(par, n_gamma)

    if(!fix_boundary_knots)
      stop("gradient should not be called with 'fix_boundary_knots' equal to FALSE")

    # handle terms from the hazard
    eta <- drop(X %*% beta)
    w <- exp(eta)
    w_y <- w * y

    fac_beta <- event * (1 + eval_basis(w_y, gamma, ders = 1L) * w_y)
    d_beta <- -drop(fac_beta %*% X)

    eval_basis(w_y, gamma, ders = 0L)
    d_gamma <- -c(sum(event), drop(event %*% Bmat))

    # handle terms from the survival function
    survival_terms <- mapply(function(x, wei){
      w_y_x <- w_y * x
      fac <- w * exp(eval_basis(w_y_x, gamma)) * wei * .5 * y
      d_gamma <- c(sum(fac), fac %*% Bmat)

      d_log_lambda_0 <- eval_basis(w_y_x, gamma, ders = 1L)
      fac_beta <- fac * (1 + d_log_lambda_0 * w_y_x)
      d_beta_terms <- fac_beta %*% X
      c(d_beta_terms, d_gamma)
    }, x = gl_node, w = gl_dat$weight)

    survival_terms <- rowSums(survival_terms)

    c(d_beta + head(survival_terms, n_beta),
      d_gamma + tail(survival_terms, n_gamma))
  }

  # check gradients if requested
  if(check_grads){
    truth <- numDeriv::grad(ll_func, c(beta, gamma))
    func_grads <- d_ll_func(c(beta, gamma))
    stopifnot(isTRUE(all.equal(truth, func_grads, check.attributes = FALSE,
                               tolerance = 1e-4)))
  }

  # optimize and return. First find better values for gamma
  method <- if(fix_boundary_knots) "BFGS" else "Nelder-Mead"
  init <- optim(gamma, function(x) ll_func(c(beta, x)),
                function(x) tail(d_ll_func(c(beta, x)), n_gamma),
                method = method)
  gamma <- init$par

  # then do joint optimization
  res <- optim(c(beta, gamma), ll_func, d_ll_func,
               control = list(maxit = maxit), method = method)
  if(res$convergence != 0)
    warning(sprintf("convergence code is %d", res$convergence))

  if(check_grads){
    truth <- numDeriv::grad(ll_func, res$par)
    func_grads <- d_ll_func(res$par)
    stopifnot(isTRUE(all.equal(truth, func_grads, check.attributes = FALSE,
                               tolerance = 1e-4)))
  }

  if(!fix_boundary_knots)
    # make sure the boundary knots are set
    ll_func(res$par)

  coefs <- setNames(res$par, c(colnames(X), paste0("gamma", seq_len(n_gamma))))
  hess <- optimHess(coefs, ll_func,
                    if(fix_boundary_knots) d_ll_func else NULL)
  vcov <- try(solve(hess))

  structure(
    list(
      coefs = coefs, logLik = -res$value, optim = res,
      eval_basis = eval_basis_out, beta = head(coefs, n_beta),
      gamma = tail(coefs, n_gamma), vcov = vcov,
      ll = ll_func, grad = d_ll_func),
    class="saft_fit")
}

# estimates a spline-based AFT with C++.
#
# The interface is like saft_fit but the backend uses C++.
saft_fit_cpp <- function(
  y, X, event, n_knots, gl_dat, basis_type = c("bs", "ns"),
  use_integrate = FALSE, maxit = 1000L, check_grads = FALSE){
  # get the knot placement
  knots <- quantile(y[event > 0], 1:n_knots / (n_knots + 1))
  b_knots <- range(y, 0)

  basis_type <- basis_type[1]

  # get starting values
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop=FALSE]
  beta <- -coef(survreg(Surv(y, event) ~ X, dist = "exponential"))
  gamma <- switch(
    basis_type,
    bs = numeric(n_knots + 4),
    ns = numeric(n_knots + 2),
    stop(sprintf("basis_type '%s' is no implemented", basis_type)))
  gamma[1] <- beta["(Intercept)"]
  beta <- beta[-1]
  n_beta <- length(beta)
  n_gamma <- length(gamma)

  # get a pointer to the C++ object
  cpp_ptr <- get_saft_fitter_ptr(
    X = X, y = y, event = event, knots = knots, boundary_knots = b_knots,
    which_base = basis_type)

  # assign negative log likelihood function
  ll_func <- function(par){
    beta <- head(par, n_beta)
    gamma <- tail(par, n_gamma)

    eval_log_likelihood_cpp(beta = beta, gamma = gamma, nodes = gl_dat$node,
                            weights = gl_dat$weight, ptr = cpp_ptr)
  }

  # assign gradient of the negative log likelihood function
  d_ll_func <- function(par){
    beta <- head(par, n_beta)
    gamma <- tail(par, n_gamma)

    drop(eval_log_likelihood_grad_cpp(
      beta = beta, gamma = gamma, nodes = gl_dat$node, weights = gl_dat$weight,
      ptr = cpp_ptr))
  }

  # check gradients if requested
  if(check_grads){
    truth <- numDeriv::grad(ll_func, c(beta, gamma))
    func_grads <- d_ll_func(c(beta, gamma))
    stopifnot(isTRUE(all.equal(truth, func_grads, check.attributes = FALSE,
                               tolerance = 1e-4)))
  }

  # optimize and return. First find better values for gamma
  init <- optim(gamma, function(x) ll_func(c(beta, x)),
                function(x) tail(d_ll_func(c(beta, x)), n_gamma),
                method = "BFGS")
  gamma <- init$par

  # then do joint optimization
  res <- optim(c(beta, gamma), ll_func, d_ll_func,
               control = list(maxit = maxit), method = "BFGS")
  if(res$convergence != 0)
    warning(sprintf("convergence code is %d", res$convergence))

  if(check_grads){
    truth <- numDeriv::grad(ll_func, res$par)
    func_grads <- d_ll_func(res$par)
    stopifnot(isTRUE(all.equal(truth, func_grads, check.attributes = FALSE,
                               tolerance = 1e-4)))
  }

  coefs <- setNames(res$par, c(colnames(X), paste0("gamma", seq_len(n_gamma))))
  hess <- optimHess(coefs, ll_func, d_ll_func)
  vcov <- try(solve(hess))

  structure(
    list(coefs = coefs, logLik = -res$value, optim = res,
         beta = head(coefs, n_beta), gamma = tail(coefs, n_gamma), vcov = vcov,
         ll = ll_func, grad = d_ll_func),
    class="saft_fit")

}

vcov.saft_fit <- function(object, ...) object$vcov
coef.saft_fit <- function(object, ...) object$coefs
summary.saft_fit <- function(object, ...) {
    est <- coef(object)
    se <- sqrt(diag(vcov(object)))
    zval <- est/se
    ans <- c(object,
             list(table = cbind(Estimate = est,
                                       `Std. Error` = se,
                                       `z value` = zval,
                                       `Pr(>|z|)` = 2 * pnorm(abs(zval), lower.tail = FALSE))))
    class(ans) <- "summary.saft_fit"
    ans
}
print.summary.saft_fit <- function(object, ...) {
    stats::printCoefmat(object$table)
    cat("log(likelihood) =", object$logLik, "\n")
}
coef.summary.saft_fit <- function(object) object$table


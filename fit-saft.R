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
#
# Returns:
#   coefs: coefficient estimates.
#   logLik: maximum log likelihood.
#   optim: output from optim.
#   eval_basis: function evaluate the basis functions.
#   beta: estimated AFT parameters.
#   gamma: estimated splines coefficients.
#   maxit: maximum number of iterations to pass to optim.
saft_fit <- function(y, X, event, n_knots, gl_dat, basis_type = c("bs", "ns"),
                     use_integrate = FALSE, maxit = 1000L){
  # get the knot placement
  knots <- quantile(y[event > 0], 1:n_knots / (n_knots + 1))
  b_knots <- range(y, 0)

  # assign function to evaluate the spline basis
  basis_type <- basis_type[1]
  switch(basis_type,
    bs = {
      bs_ptr <- get_bs_ptr(knots = knots, boundary_knots = b_knots,
                           intercept = FALSE)
      eval_basis <- function(x){
        out <- eval_spline_basis(x, bs_ptr)
        cbind(1, out) # add intercept
      }
    },
    ns = {
      ns_ptr <- bs_ptr <- get_ns_ptr(knots = knots, boundary_knots = b_knots,
                                     intercept = FALSE)
      eval_basis <- function(x){
        out <- eval_spline_basis(x, bs_ptr)
        cbind(1, out) # add intercept
      }
    },
    stop(sprintf("basis_type '%s' is no implemented", basis_type)))

  # get starting values
  X <- X[, setdiff(colnames(X), "(Intercept)")]
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

  # move nodes and recale the nodes
  gl_node <- (gl_dat$node + 1) * .5

  # assign negative log likelihood function
  ll_func <- function(par){
    beta <- head(par, n_beta)
    gamma <- tail(par, n_gamma)

    # handle terms from the hazard
    eta <- drop(X %*% beta)
    w <- exp(eta)
    l1 <- -sum(event * (eta + drop(eval_basis(w * y) %*% gamma)))

    # handle terms from the survival function
    if(!use_integrate){
      # use Gauss–Legendre quadrature
      l2 <- mapply(function(x, wei)
        sum(w * exp(drop(eval_basis(w * y * x) %*% gamma)) * wei * .5 * y),
        x = gl_node, w = gl_dat$weight)
      l2 <- sum(l2)
    } else {
      # use R's integrate function
      l2 <- mapply(function(w, y){
        out <- try(
          integrate(function(z) exp(drop(eval_basis(w * z) %*% gamma)),
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

  # optimize and return. First find better values for gamma
  init <- optim(gamma, function(x) ll_func(c(beta, x)))
  gamma <- init$par

  # then do joint optimization
  res <- optim(c(beta, gamma), ll_func, control = list(maxit = maxit))
  if(res$convergence != 0)
    warning(sprintf("convergence code is %d", res$convergence))

  coefs <- setNames(res$par, c(colnames(X), paste0("gamma", seq_len(n_gamma))))
  list(coefs = coefs, logLik = -res$value, optim = res, eval_basis = eval_basis,
       beta = head(coefs, n_beta), gamma = tail(coefs, n_gamma))
}

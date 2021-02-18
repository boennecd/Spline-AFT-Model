# load libraries
library(splines)
library(survival)
library(SimSurvNMarker)

# prepare colon data set
colon_use <- within(
  colon[complete.cases(colon) & colon$etype==2, ], {
    time <- time/365.25
    age_c <- scale(age)
})

mf <- model.frame(time ~ rx + sex + age_c + obstruct + perfor + adhere +
                    factor(differ) + factor(extent) + surg + node4,
                  colon_use)

X <- model.matrix(terms(mf), mf)

y <- model.response(mf)
event <- colon_use$status

# get Gauss–Legendre quadrature node and weights
gl_dat <- get_gl_rule(50)

# estimates a spline-based AFT.
#
# Args:
#   y: observed event time or censoring time.
#   X: design matrix for the fixed effects.
#   event: event indicator.
#   n_knots: number of knots for B-splines.
#   gl_dat: Gauss–Legendre Quadrature node and weights
#
# Returns:
#   coefs: coefficient estimates.
#   logLik: maximum log likelihood.
#   optim: output from optim.
#   eval_basis: function evaluate the basis functions.
#   beta: estimated AFT parameters.
#   gamma: estimated splines coefficients.
saft_fit <- function(y, X, event, n_knots, gl_dat){
  # get the knot placement
  knots <- quantile(y[event > 0], 1:n_knots / (n_knots + 1))
  b_knots <- range(y, 0)

  # evaluates the spline basis
  ns_obj <- get_ns_spline(sort(c(knots, b_knots)), intercept = FALSE,
                          do_log = FALSE)
  eval_basis <- function(x){
    out <- ns_obj(x)
    cbind(1, out) # add intercept
  }

  # get starting values
  X <- X[, setdiff(colnames(X), "(Intercept)")]
  beta <- -coef(survreg(Surv(y, event) ~ X, dist = "exponential"))
  gamma <- numeric(n_knots + 2)
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

    eta <- drop(X %*% beta)
    w <- exp(eta)
    l1 <- -sum(event * (eta + drop(eval_basis(w * y) %*% gamma)))

    # use Gauss–Legendre Quadrature
    l2 <- mapply(function(x, wei)
      sum(w * exp(drop(eval_basis(w * y * x) %*% gamma)) * wei * .5 * y),
      x = gl_node, w = gl_dat$weight)
    l2 <- sum(l2)

    # use R's integrate function
    # l2 <- mapply(function(w, y){
    #   out <- try(
    #     integrate(function(z) exp(drop(eval_basis(w * z) %*% gamma)),
    #               lower = 0, upper = y)$value, silent = TRUE)
    #   if(inherits(out, "try-error"))
    #     out <- Inf
    #   out
    # }, w = w, y = y)
    # l2 <- sum(l2 * w)

    # cat(sprintf("l1: %.2f l2: %.2f sum: %.2f\n", l1, l2, l1 + l2))
    l1 + l2
  }

  # optimize and return. First find better values for gamma
  init <- optim(gamma, function(x) ll_func(c(beta, x)))
  gamma <- init$par

  # then do joint optimization
  res <- optim(c(beta, gamma), ll_func, control = list(maxit = 1000L))

  coefs <- setNames(res$par, c(colnames(X), paste0("gamma", seq_len(n_gamma))))
  list(coefs = coefs, logLik = -res$value, optim = res, eval_basis = eval_basis,
       beta = head(coefs, n_beta), gamma = tail(coefs, n_gamma))
}

# use the function
system.time(
  res <- saft_fit(y = y, X = X, event = event, n_knots = 4, gl_dat = gl_dat))
res$coefs
res$logLik

# plot the spline
jpeg("saft_fit-spline-ns.jpg", width = 600, height = 400)
ys <- seq(0, 2 * max(y), length.out = 1000)
plot(ys, res$eval_basis(ys) %*% res$gamma, type = "l", ylim = c(-5.5, -2.75),
     main = "saft_fit (ns)", xlab = "Time", ylab = "Spline", bty = "l")
grid()
dev.off()

# compare with Weibull model
survreg_res <- survreg(Surv(y, event) ~ X - 1)
coef(survreg_res)
logLik(survreg_res)
